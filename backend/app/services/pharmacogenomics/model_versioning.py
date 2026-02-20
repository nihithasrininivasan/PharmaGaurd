"""
Model Versioning Module - Track algorithm changes, enable rollback, A/B testing.

Features:
- Semantic versioning (major.minor.patch)
- Model snapshots with all parameters
- Performance metrics tracking per version
- Rollback capability
- A/B testing support
- Deterministic reproduction
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict
from datetime import datetime


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class ModelVersion:
    """Complete snapshot of risk engine configuration."""
    version: str                                # Semantic version (e.g., "2.1.3")
    timestamp: str                              # ISO8601 timestamp
    notes: str                                  # Change description

    # Algorithm parameters
    penalty_weights: Dict[str, float]
    severity_base_scores: Dict[str, float]
    confidence_component_weights: Dict[str, float]
    phenotype_modifiers: Dict[str, float]

    # Learning state
    learning_priors: Dict
    calibration_factors: Dict

    # Performance metrics
    total_predictions: int
    total_feedback: int
    accuracy_by_confidence: Dict[str, float]
    false_positive_rate: float
    false_negative_rate: float
    mean_calibration_error: float

    # Metadata
    deployed: bool = False
    baseline_version: bool = False


@dataclass
class VersionHistory:
    """History of all model versions."""
    versions: List[ModelVersion]
    current_version: str
    baseline_version: str


# ============================================================================
# Model Version Manager
# ============================================================================

class ModelVersionManager:
    """
    Manage model versions with rollback and A/B testing support.

    Provides:
    - Version snapshots
    - Version comparison
    - Rollback to previous versions
    - A/B testing between versions
    """

    def __init__(self, versions_dir: Path = Path("data/model_versions")):
        self.versions_dir = versions_dir
        self.versions_dir.mkdir(parents=True, exist_ok=True)

    def save_version(self, version: ModelVersion):
        """
        Save a model version snapshot.

        Args:
            version: ModelVersion to save
        """
        file_path = self.versions_dir / f"{version.version}.json"

        with open(file_path, 'w') as f:
            json.dump(asdict(version), f, indent=2)

        print(f"Saved model version {version.version} to {file_path}")

    def load_version(self, version_id: str) -> Optional[ModelVersion]:
        """
        Load a specific model version.

        Args:
            version_id: Version identifier (e.g., "2.1.3")

        Returns:
            ModelVersion or None if not found
        """
        file_path = self.versions_dir / f"{version_id}.json"

        if not file_path.exists():
            print(f"Version {version_id} not found")
            return None

        with open(file_path, 'r') as f:
            data = json.load(f)

        return ModelVersion(**data)

    def list_versions(self) -> List[str]:
        """
        List all available model versions.

        Returns:
            List of version identifiers sorted by timestamp
        """
        version_files = sorted(self.versions_dir.glob("*.json"))
        versions = []

        for file_path in version_files:
            version_id = file_path.stem
            versions.append(version_id)

        return versions

    def get_version_history(self) -> VersionHistory:
        """
        Get complete version history.

        Returns:
            VersionHistory with all versions
        """
        version_ids = self.list_versions()
        versions = []

        for vid in version_ids:
            version = self.load_version(vid)
            if version:
                versions.append(version)

        # Sort by timestamp
        versions.sort(key=lambda v: v.timestamp, reverse=True)

        # Determine current and baseline
        current = versions[0].version if versions else "0.0.0"
        baseline = next(
            (v.version for v in versions if v.baseline_version),
            current
        )

        return VersionHistory(
            versions=versions,
            current_version=current,
            baseline_version=baseline,
        )

    def compare_versions(
        self,
        version_a: str,
        version_b: str
    ) -> Dict:
        """
        Compare two model versions.

        Args:
            version_a: First version ID
            version_b: Second version ID

        Returns:
            Comparison dictionary with differences
        """
        v_a = self.load_version(version_a)
        v_b = self.load_version(version_b)

        if not v_a or not v_b:
            return {"error": "One or both versions not found"}

        comparison = {
            "version_a": version_a,
            "version_b": version_b,
            "performance_diff": {
                "accuracy_delta": (
                    v_b.mean_calibration_error - v_a.mean_calibration_error
                ),
                "fp_rate_delta": (
                    v_b.false_positive_rate - v_a.false_positive_rate
                ),
                "fn_rate_delta": (
                    v_b.false_negative_rate - v_a.false_negative_rate
                ),
            },
            "parameter_changes": self._diff_parameters(v_a, v_b),
        }

        return comparison

    def rollback(self, target_version: str) -> ModelVersion:
        """
        Rollback to a previous model version.

        Args:
            target_version: Version to rollback to

        Returns:
            Loaded ModelVersion

        Raises:
            ValueError: If target version not found
        """
        version = self.load_version(target_version)

        if not version:
            raise ValueError(f"Version {target_version} not found")

        print(f"Rolling back to version {target_version}")
        print(f"Notes: {version.notes}")

        # Mark as deployed
        version.deployed = True
        self.save_version(version)

        return version

    def create_new_version(
        self,
        from_version: Optional[str],
        version_increment: str,
        notes: str,
        **kwargs
    ) -> ModelVersion:
        """
        Create a new model version from an existing one.

        Args:
            from_version: Base version to copy from (None = latest)
            version_increment: "major", "minor", or "patch"
            notes: Description of changes
            **kwargs: Parameters to update

        Returns:
            New ModelVersion
        """
        # Load base version
        if from_version:
            base = self.load_version(from_version)
        else:
            history = self.get_version_history()
            if history.versions:
                base = history.versions[0]
            else:
                base = self._create_default_version()

        # Increment version number
        new_version_id = self._increment_version(base.version, version_increment)

        # Copy all parameters from base
        new_version_data = asdict(base)

        # Update with new values
        new_version_data.update(kwargs)
        new_version_data["version"] = new_version_id
        new_version_data["timestamp"] = datetime.now().isoformat()
        new_version_data["notes"] = notes
        new_version_data["deployed"] = False

        new_version = ModelVersion(**new_version_data)
        self.save_version(new_version)

        return new_version

    def _increment_version(self, current: str, increment: str) -> str:
        """Increment semantic version."""
        parts = current.split(".")
        major, minor, patch = int(parts[0]), int(parts[1]), int(parts[2])

        if increment == "major":
            major += 1
            minor = 0
            patch = 0
        elif increment == "minor":
            minor += 1
            patch = 0
        elif increment == "patch":
            patch += 1
        else:
            raise ValueError(f"Invalid increment: {increment}")

        return f"{major}.{minor}.{patch}"

    def _diff_parameters(self, v_a: ModelVersion, v_b: ModelVersion) -> Dict:
        """Compare parameters between two versions."""
        changes = {}

        # Compare penalty weights
        if v_a.penalty_weights != v_b.penalty_weights:
            changes["penalty_weights"] = {
                k: (v_a.penalty_weights.get(k), v_b.penalty_weights.get(k))
                for k in set(v_a.penalty_weights) | set(v_b.penalty_weights)
                if v_a.penalty_weights.get(k) != v_b.penalty_weights.get(k)
            }

        # Compare severity scores
        if v_a.severity_base_scores != v_b.severity_base_scores:
            changes["severity_base_scores"] = {
                k: (v_a.severity_base_scores.get(k), v_b.severity_base_scores.get(k))
                for k in set(v_a.severity_base_scores) | set(v_b.severity_base_scores)
                if v_a.severity_base_scores.get(k) != v_b.severity_base_scores.get(k)
            }

        return changes

    def _create_default_version(self) -> ModelVersion:
        """Create a default version 1.0.0."""
        return ModelVersion(
            version="1.0.0",
            timestamp=datetime.now().isoformat(),
            notes="Initial baseline version",
            penalty_weights={},
            severity_base_scores={},
            confidence_component_weights={},
            phenotype_modifiers={},
            learning_priors={},
            calibration_factors={},
            total_predictions=0,
            total_feedback=0,
            accuracy_by_confidence={},
            false_positive_rate=0.0,
            false_negative_rate=0.0,
            mean_calibration_error=0.0,
            baseline_version=True,
        )


# ============================================================================
# A/B Testing Support
# ============================================================================

class ABTestManager:
    """
    Manage A/B testing between two model versions.

    Allows:
    - Split traffic between versions
    - Track performance per version
    - Determine winner
    """

    def __init__(self, version_manager: ModelVersionManager):
        self.version_manager = version_manager
        self.active_tests: Dict[str, Dict] = {}

    def start_ab_test(
        self,
        test_id: str,
        version_a: str,
        version_b: str,
        traffic_split: float = 0.5,
    ):
        """
        Start an A/B test between two versions.

        Args:
            test_id: Unique test identifier
            version_a: Control version
            version_b: Treatment version
            traffic_split: Fraction to version_b (default 0.5 = 50/50)
        """
        self.active_tests[test_id] = {
            "version_a": version_a,
            "version_b": version_b,
            "traffic_split": traffic_split,
            "start_time": datetime.now().isoformat(),
            "results_a": [],
            "results_b": [],
        }

        print(f"Started A/B test '{test_id}': {version_a} vs {version_b}")

    def record_result(
        self,
        test_id: str,
        version: str,
        was_correct: bool,
        confidence: float,
    ):
        """Record a prediction result for A/B test."""
        if test_id not in self.active_tests:
            return

        test = self.active_tests[test_id]
        result = {"correct": was_correct, "confidence": confidence}

        if version == test["version_a"]:
            test["results_a"].append(result)
        elif version == test["version_b"]:
            test["results_b"].append(result)

    def get_test_results(self, test_id: str) -> Dict:
        """Get A/B test results and statistics."""
        if test_id not in self.active_tests:
            return {"error": "Test not found"}

        test = self.active_tests[test_id]

        # Calculate accuracy for each version
        results_a = test["results_a"]
        results_b = test["results_b"]

        accuracy_a = (
            sum(1 for r in results_a if r["correct"]) / len(results_a)
            if results_a else 0.0
        )
        accuracy_b = (
            sum(1 for r in results_b if r["correct"]) / len(results_b)
            if results_b else 0.0
        )

        # Determine winner
        if accuracy_b > accuracy_a:
            winner = test["version_b"]
            improvement = accuracy_b - accuracy_a
        else:
            winner = test["version_a"]
            improvement = accuracy_a - accuracy_b

        return {
            "test_id": test_id,
            "version_a": test["version_a"],
            "version_b": test["version_b"],
            "samples_a": len(results_a),
            "samples_b": len(results_b),
            "accuracy_a": accuracy_a,
            "accuracy_b": accuracy_b,
            "winner": winner,
            "improvement": improvement,
        }


# ============================================================================
# Convenience Functions
# ============================================================================

def create_version_manager(
    versions_dir: Path = Path("data/model_versions")
) -> ModelVersionManager:
    """Create a new model version manager."""
    return ModelVersionManager(versions_dir)


def create_ab_test_manager(
    version_manager: ModelVersionManager
) -> ABTestManager:
    """Create an A/B test manager."""
    return ABTestManager(version_manager)
