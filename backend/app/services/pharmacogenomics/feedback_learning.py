"""
Feedback Learning Module - Bayesian adaptive learning from clinical corrections.

Implements:
- Bayesian posterior updates with bounded learning
- Time-based decay of old feedback
- Confidence-weighted feedback integration
- Batch calibration from historical data
- Overfitting prevention (max delta, bounds, regularization)
"""

from __future__ import annotations

import math
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta
from collections import defaultdict


# ============================================================================
# Configuration
# ============================================================================

DEFAULT_LEARNING_RATE = 0.1         # α: How much to trust new feedback
DEFAULT_DECAY_RATE = 0.95           # β: Monthly decay factor
MIN_PRIOR = 0.80                    # Floor for prior values
MAX_PRIOR = 1.50                    # Ceiling for prior values
MAX_DELTA_PER_UPDATE = 0.10         # Maximum change per single feedback
MIN_FEEDBACK_FOR_UPDATE = 1         # Minimum feedback events before applying


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class FeedbackEvent:
    """Single clinician feedback event."""
    gene: str
    reported_diplotype: str
    correct_diplotype: str
    timestamp: datetime
    feedback_quality: float = 1.0    # Clinician confidence [0, 1]
    comments: Optional[str] = None

    def months_since(self, reference_time: datetime) -> float:
        """Calculate months elapsed since this feedback."""
        delta = reference_time - self.timestamp
        return delta.days / 30.0


@dataclass
class LearningPriors:
    """Learning priors for all genes."""
    genes: Dict[str, Dict[str, Dict[str, float]]]
    metadata: Dict
    model_version: str = "2.0.0"
    last_calibration: Optional[str] = None

    def get_diplotype_prior(self, gene: str, diplotype: str) -> float:
        """Get prior for specific gene-diplotype combination."""
        return (
            self.genes
            .get(gene, {})
            .get("diplotype_priors", {})
            .get(diplotype, 1.0)
        )

    def set_diplotype_prior(self, gene: str, diplotype: str, value: float):
        """Set prior for specific gene-diplotype combination."""
        if gene not in self.genes:
            self.genes[gene] = {
                "diplotype_priors": {},
                "variant_weights": {}
            }
        self.genes[gene]["diplotype_priors"][diplotype] = value


# ============================================================================
# Bayesian Feedback Learner
# ============================================================================

class BayesianFeedbackLearner:
    """
    Adaptive learning from clinical feedback using Bayesian updates.

    Features:
    - Bounded updates to prevent overfitting
    - Time-based decay of old feedback
    - Confidence-weighted feedback integration
    - Regularization to prevent single-sample overcorrection
    """

    def __init__(
        self,
        learning_rate: float = DEFAULT_LEARNING_RATE,
        decay_rate: float = DEFAULT_DECAY_RATE,
        min_prior: float = MIN_PRIOR,
        max_prior: float = MAX_PRIOR,
        max_delta: float = MAX_DELTA_PER_UPDATE,
    ):
        self.learning_rate = learning_rate
        self.decay_rate = decay_rate
        self.min_prior = min_prior
        self.max_prior = max_prior
        self.max_delta = max_delta

    def update_prior(
        self,
        current_prior: float,
        feedback_quality: float,
        months_since_last_update: float,
    ) -> Tuple[float, str]:
        """
        Update a diplotype prior from a single feedback event.

        Args:
            current_prior: Current prior value (default 1.0)
            feedback_quality: Clinician confidence in feedback [0, 1]
            months_since_last_update: Time since last calibration

        Returns:
            (updated_prior, explanation)
        """
        # 1. Apply time decay to current prior
        # Drift back toward neutral (1.0) over time
        decayed_prior = 1.0 + (current_prior - 1.0) * (
            self.decay_rate ** months_since_last_update
        )

        # 2. Feedback signal (positive reinforcement)
        # Higher quality → stronger signal
        # Signal range: [1.0, 1.1] for quality [0, 1]
        signal = 1.0 + (feedback_quality * 0.1)

        # 3. Bayesian weighted update
        # learning_rate controls trust in new feedback vs. old prior
        updated = (
            self.learning_rate * signal +
            (1 - self.learning_rate) * decayed_prior
        )

        # 4. Limit delta to prevent single-sample overcorrection
        delta = updated - current_prior
        if abs(delta) > self.max_delta:
            delta = math.copysign(self.max_delta, delta)
            updated = current_prior + delta
            limited = True
        else:
            limited = False

        # 5. Bound to safe range
        bounded_prior = max(self.min_prior, min(self.max_prior, updated))

        # Explanation
        explanation = (
            f"Updated from {current_prior:.3f} → {bounded_prior:.3f} "
            f"(decayed={decayed_prior:.3f}, signal={signal:.3f}, "
            f"delta={delta:+.3f}"
        )
        if limited:
            explanation += ", delta_limited"
        explanation += ")"

        return bounded_prior, explanation

    def batch_calibration(
        self,
        feedback_history: List[FeedbackEvent],
        reference_time: Optional[datetime] = None,
    ) -> LearningPriors:
        """
        Recalibrate all priors from accumulated feedback history.

        Args:
            feedback_history: All feedback events
            reference_time: Reference time for decay (default: now)

        Returns:
            Updated LearningPriors object
        """
        if reference_time is None:
            reference_time = datetime.now()

        priors = LearningPriors(
            genes={},
            metadata={
                "total_feedback_events": len(feedback_history),
                "last_updated": reference_time.isoformat(),
            },
            last_calibration=reference_time.isoformat(),
        )

        # Group feedback by (gene, correct_diplotype)
        grouped: Dict[Tuple[str, str], List[FeedbackEvent]] = defaultdict(list)
        for event in feedback_history:
            grouped[(event.gene, event.correct_diplotype)].append(event)

        # Update each gene-diplotype combination
        for (gene, diplotype), events in grouped.items():
            if len(events) < MIN_FEEDBACK_FOR_UPDATE:
                continue  # Skip if insufficient feedback

            # Start from neutral prior
            current = 1.0

            # Apply all feedback events sequentially (oldest to newest)
            events_sorted = sorted(events, key=lambda e: e.timestamp)

            for event in events_sorted:
                months_ago = event.months_since(reference_time)
                current, _ = self.update_prior(
                    current_prior=current,
                    feedback_quality=event.feedback_quality,
                    months_since_last_update=months_ago,
                )

            # Store final prior
            priors.set_diplotype_prior(gene, diplotype, current)

        return priors

    def incremental_update(
        self,
        current_priors: LearningPriors,
        new_feedback: FeedbackEvent,
    ) -> LearningPriors:
        """
        Update priors incrementally with a single new feedback event.

        Args:
            current_priors: Current learning priors
            new_feedback: New feedback event

        Returns:
            Updated LearningPriors
        """
        gene = new_feedback.gene
        diplotype = new_feedback.correct_diplotype

        # Get current prior (default 1.0)
        current = current_priors.get_diplotype_prior(gene, diplotype)

        # Calculate months since last calibration
        if current_priors.last_calibration:
            last_cal = datetime.fromisoformat(current_priors.last_calibration)
            months_since = new_feedback.months_since(last_cal)
        else:
            months_since = 0.0

        # Update prior
        updated, explanation = self.update_prior(
            current_prior=current,
            feedback_quality=new_feedback.feedback_quality,
            months_since_last_update=months_since,
        )

        # Create updated priors object
        updated_priors = LearningPriors(
            genes=current_priors.genes.copy(),
            metadata=current_priors.metadata.copy(),
            model_version=current_priors.model_version,
            last_calibration=current_priors.last_calibration,
        )

        updated_priors.set_diplotype_prior(gene, diplotype, updated)

        # Update metadata
        updated_priors.metadata["total_feedback_events"] = (
            current_priors.metadata.get("total_feedback_events", 0) + 1
        )
        updated_priors.metadata["last_updated"] = datetime.now().isoformat()
        updated_priors.metadata["last_update_explanation"] = explanation

        return updated_priors


# ============================================================================
# Persistence Layer
# ============================================================================

class LearningPriorsManager:
    """Manage loading/saving of learning priors."""

    def __init__(self, file_path: Path = Path("data/learning_priors.json")):
        self.file_path = file_path

    def load(self) -> LearningPriors:
        """Load learning priors from disk."""
        if not self.file_path.exists():
            return self._create_default_priors()

        try:
            with open(self.file_path, 'r') as f:
                data = json.load(f)

            return LearningPriors(
                genes=data.get("genes", {}),
                metadata=data.get("metadata", {}),
                model_version=data.get("model_version", "1.0.0"),
                last_calibration=data.get("last_calibration"),
            )
        except Exception as e:
            print(f"Error loading learning priors: {e}")
            return self._create_default_priors()

    def save(self, priors: LearningPriors):
        """Save learning priors to disk."""
        self.file_path.parent.mkdir(parents=True, exist_ok=True)

        data = {
            "genes": priors.genes,
            "metadata": priors.metadata,
            "model_version": priors.model_version,
            "last_calibration": priors.last_calibration,
        }

        with open(self.file_path, 'w') as f:
            json.dump(data, f, indent=2)

    def _create_default_priors(self) -> LearningPriors:
        """Create default empty priors."""
        return LearningPriors(
            genes={},
            metadata={
                "total_feedback_events": 0,
                "last_updated": None,
            },
        )


# ============================================================================
# Convenience Functions
# ============================================================================

def load_learning_priors(
    file_path: Path = Path("data/learning_priors.json")
) -> LearningPriors:
    """Load learning priors from disk."""
    manager = LearningPriorsManager(file_path)
    return manager.load()


def save_learning_priors(
    priors: LearningPriors,
    file_path: Path = Path("data/learning_priors.json")
):
    """Save learning priors to disk."""
    manager = LearningPriorsManager(file_path)
    manager.save(priors)


def get_diplotype_boost(
    priors: LearningPriors,
    gene: str,
    diplotype: str
) -> float:
    """Get feedback learning boost for a gene-diplotype combination."""
    return priors.get_diplotype_prior(gene, diplotype)
