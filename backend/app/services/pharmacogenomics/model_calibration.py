"""
Model Calibration Module - Confidence calibration, drift detection, and performance tracking.

Features:
- Confidence calibration (ensure confidence scores match real accuracy)
- Drift detection (monitor model performance degradation)
- False positive/negative tracking
- Penalty weight recalibration
- Distribution shift detection
"""

from __future__ import annotations

import math
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, asdict
from datetime import datetime
from collections import defaultdict


# ============================================================================
# Configuration
# ============================================================================

MIN_SAMPLES_FOR_CALIBRATION = 10    # Minimum samples per bin
DRIFT_ALERT_THRESHOLD = 2.0         # Standard deviations
CALIBRATION_BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class PredictionOutcome:
    """Record of a single prediction and its outcome."""
    prediction_id: str
    timestamp: datetime
    gene: str
    diplotype_predicted: str
    diplotype_actual: Optional[str]
    confidence: float
    risk_score: float
    risk_level: str
    was_correct: Optional[bool]


@dataclass
class CalibrationStats:
    """Calibration statistics for a confidence bin."""
    bin_label: str
    confidence_range: Tuple[float, float]
    total_predictions: int
    correct_predictions: int
    empirical_accuracy: float
    calibration_error: float


@dataclass
class ModelPerformanceMetrics:
    """Overall model performance metrics."""
    total_predictions: int
    total_feedback: int
    overall_accuracy: float
    accuracy_by_confidence: Dict[str, float]
    false_positive_rate: float
    false_negative_rate: float
    precision_by_severity: Dict[str, float]
    mean_calibration_error: float
    timestamp: str


# ============================================================================
# Confidence Calibrator
# ============================================================================

class ConfidenceCalibrator:
    """
    Calibrate confidence scores to match real-world accuracy.

    Tracks:
    - Prediction vs. outcome for confidence bins
    - Calibration curves
    - Calibration correction factors
    """

    def __init__(self, bins: Optional[List[float]] = None):
        self.bins = bins or CALIBRATION_BINS
        self.outcomes_by_bin: Dict[int, List[PredictionOutcome]] = defaultdict(list)

    def record_outcome(self, outcome: PredictionOutcome):
        """Record a prediction outcome for calibration."""
        bin_idx = self._get_bin_index(outcome.confidence)
        self.outcomes_by_bin[bin_idx].append(outcome)

    def get_calibration_stats(self) -> List[CalibrationStats]:
        """Calculate calibration statistics for all bins."""
        stats = []

        for bin_idx in range(len(self.bins) - 1):
            bin_lower = self.bins[bin_idx]
            bin_upper = self.bins[bin_idx + 1]

            outcomes = self.outcomes_by_bin.get(bin_idx, [])
            outcomes_with_result = [o for o in outcomes if o.was_correct is not None]

            if len(outcomes_with_result) < MIN_SAMPLES_FOR_CALIBRATION:
                continue  # Skip bins with insufficient data

            total = len(outcomes_with_result)
            correct = sum(1 for o in outcomes_with_result if o.was_correct)
            empirical_accuracy = correct / total

            # Expected accuracy = midpoint of bin
            expected_accuracy = (bin_lower + bin_upper) / 2.0
            calibration_error = abs(empirical_accuracy - expected_accuracy)

            stats.append(CalibrationStats(
                bin_label=f"{bin_lower:.1f}-{bin_upper:.1f}",
                confidence_range=(bin_lower, bin_upper),
                total_predictions=total,
                correct_predictions=correct,
                empirical_accuracy=empirical_accuracy,
                calibration_error=calibration_error,
            ))

        return stats

    def get_calibration_factor(self, confidence: float) -> float:
        """
        Get calibration correction factor for a confidence score.

        Returns:
            Correction factor (multiply confidence by this to get calibrated score)

        Example:
            If confidence=0.90 but empirical accuracy=0.85,
            return 0.85/0.90 = 0.944
        """
        bin_idx = self._get_bin_index(confidence)
        outcomes = self.outcomes_by_bin.get(bin_idx, [])
        outcomes_with_result = [o for o in outcomes if o.was_correct is not None]

        if len(outcomes_with_result) < MIN_SAMPLES_FOR_CALIBRATION:
            return 1.0  # No calibration data

        total = len(outcomes_with_result)
        correct = sum(1 for o in outcomes_with_result if o.was_correct)
        empirical_accuracy = correct / total

        # Correction factor
        # If overconfident (accuracy < confidence), factor < 1.0
        # If underconfident (accuracy > confidence), factor > 1.0
        factor = empirical_accuracy / max(0.01, confidence)

        # Bound correction to avoid extreme adjustments
        return max(0.80, min(1.20, factor))

    def mean_calibration_error(self) -> float:
        """Calculate mean absolute calibration error across all bins."""
        stats = self.get_calibration_stats()
        if not stats:
            return 0.0

        total_error = sum(s.calibration_error for s in stats)
        return total_error / len(stats)

    def recalibrate_penalty_weights(
        self,
        current_penalties: Dict[str, float]
    ) -> Dict[str, float]:
        """
        Adjust penalty weights based on calibration curve.

        If model is consistently overconfident, increase penalties.
        If model is consistently underconfident, decrease penalties.

        Returns:
            Updated penalty weights
        """
        avg_overconfidence = self._compute_average_overconfidence()

        adjusted_penalties = current_penalties.copy()

        if avg_overconfidence > 0.05:  # >5% overconfident
            # Increase all penalties by 10%
            print(f"Model overconfident by {avg_overconfidence:.2%} — increasing penalties by 10%")
            for key in adjusted_penalties:
                adjusted_penalties[key] = min(0.90, adjusted_penalties[key] * 1.10)

        elif avg_overconfidence < -0.05:  # >5% underconfident
            # Decrease all penalties by 10%
            print(f"Model underconfident by {abs(avg_overconfidence):.2%} — decreasing penalties by 10%")
            for key in adjusted_penalties:
                adjusted_penalties[key] = max(0.05, adjusted_penalties[key] * 0.90)

        return adjusted_penalties

    def _get_bin_index(self, confidence: float) -> int:
        """Get bin index for a confidence score."""
        for i in range(len(self.bins) - 1):
            if self.bins[i] <= confidence < self.bins[i + 1]:
                return i
        return len(self.bins) - 2  # Last bin

    def _compute_average_overconfidence(self) -> float:
        """
        Compute average overconfidence across all bins.

        Positive = overconfident (confidence > accuracy)
        Negative = underconfident (confidence < accuracy)
        """
        stats = self.get_calibration_stats()
        if not stats:
            return 0.0

        total_error = 0.0
        for s in stats:
            expected = (s.confidence_range[0] + s.confidence_range[1]) / 2.0
            # Positive if overconfident
            total_error += (expected - s.empirical_accuracy)

        return total_error / len(stats)


# ============================================================================
# Drift Detector
# ============================================================================

class DriftDetector:
    """
    Monitor model performance for significant drift.

    Uses statistical process control to detect:
    - Accuracy degradation
    - Distribution shifts
    - Performance anomalies
    """

    def __init__(
        self,
        baseline_accuracy: float = 0.85,
        baseline_std: float = 0.05,
        alert_threshold: float = DRIFT_ALERT_THRESHOLD,
    ):
        self.baseline_accuracy = baseline_accuracy
        self.baseline_std = baseline_std
        self.alert_threshold = alert_threshold

    def check_accuracy_drift(self, recent_accuracy: float) -> Optional[str]:
        """
        Check if current accuracy has drifted from baseline.

        Args:
            recent_accuracy: Recent accuracy (e.g., last 100 predictions)

        Returns:
            Alert message if drift detected, None otherwise
        """
        z_score = abs(recent_accuracy - self.baseline_accuracy) / self.baseline_std

        if z_score > self.alert_threshold:
            direction = "decreased" if recent_accuracy < self.baseline_accuracy else "increased"
            return (
                f"DRIFT ALERT: Accuracy {direction} to {recent_accuracy:.2%} "
                f"(baseline: {self.baseline_accuracy:.2%}, z-score={z_score:.2f})"
            )

        return None

    def check_distribution_shift(
        self,
        baseline_dist: Dict[str, int],
        current_dist: Dict[str, int]
    ) -> Tuple[float, Optional[str]]:
        """
        Detect input distribution shift using KL divergence.

        Args:
            baseline_dist: Historical phenotype/severity frequencies
            current_dist: Recent phenotype/severity frequencies

        Returns:
            (kl_divergence, alert_message)
        """
        # Normalize to probabilities
        all_keys = set(baseline_dist.keys()) | set(current_dist.keys())

        baseline_probs = []
        current_probs = []

        for key in sorted(all_keys):
            baseline_probs.append(baseline_dist.get(key, 0))
            current_probs.append(current_dist.get(key, 0))

        baseline_sum = sum(baseline_probs)
        current_sum = sum(current_probs)

        if baseline_sum == 0 or current_sum == 0:
            return 0.0, None

        baseline_probs = [p / baseline_sum for p in baseline_probs]
        current_probs = [p / current_sum for p in current_probs]

        # Calculate KL divergence
        kl_div = 0.0
        for p, q in zip(current_probs, baseline_probs):
            if p > 0 and q > 0:
                kl_div += p * math.log(p / q)

        if kl_div > 0.1:  # Threshold for significant shift
            return kl_div, (
                f"DRIFT ALERT: Input distribution shifted significantly "
                f"(KL divergence={kl_div:.4f})"
            )

        return kl_div, None


# ============================================================================
# Performance Tracker
# ============================================================================

class PerformanceTracker:
    """
    Track model performance metrics over time.

    Metrics:
    - Overall accuracy
    - False positive/negative rates
    - Precision by severity level
    - Calibration quality
    """

    def __init__(self):
        self.outcomes: List[PredictionOutcome] = []

    def record_outcome(self, outcome: PredictionOutcome):
        """Record a prediction outcome."""
        self.outcomes.append(outcome)

    def compute_metrics(self) -> ModelPerformanceMetrics:
        """Compute comprehensive performance metrics."""
        outcomes_with_result = [o for o in self.outcomes if o.was_correct is not None]

        if not outcomes_with_result:
            return ModelPerformanceMetrics(
                total_predictions=0,
                total_feedback=0,
                overall_accuracy=0.0,
                accuracy_by_confidence={},
                false_positive_rate=0.0,
                false_negative_rate=0.0,
                precision_by_severity={},
                mean_calibration_error=0.0,
                timestamp=datetime.now().isoformat(),
            )

        # Overall accuracy
        total = len(outcomes_with_result)
        correct = sum(1 for o in outcomes_with_result if o.was_correct)
        overall_accuracy = correct / total

        # Accuracy by confidence bins
        calibrator = ConfidenceCalibrator()
        for outcome in outcomes_with_result:
            calibrator.record_outcome(outcome)

        calibration_stats = calibrator.get_calibration_stats()
        accuracy_by_confidence = {
            s.bin_label: s.empirical_accuracy
            for s in calibration_stats
        }

        # False positive/negative rates
        # (Requires ground truth labels — placeholder logic)
        fp_rate = self._compute_false_positive_rate(outcomes_with_result)
        fn_rate = self._compute_false_negative_rate(outcomes_with_result)

        # Precision by severity
        precision_by_severity = self._compute_precision_by_severity(outcomes_with_result)

        # Mean calibration error
        mce = calibrator.mean_calibration_error()

        return ModelPerformanceMetrics(
            total_predictions=len(self.outcomes),
            total_feedback=len(outcomes_with_result),
            overall_accuracy=overall_accuracy,
            accuracy_by_confidence=accuracy_by_confidence,
            false_positive_rate=fp_rate,
            false_negative_rate=fn_rate,
            precision_by_severity=precision_by_severity,
            mean_calibration_error=mce,
            timestamp=datetime.now().isoformat(),
        )

    def _compute_false_positive_rate(self, outcomes: List[PredictionOutcome]) -> float:
        """Compute false positive rate (high-risk predictions that were wrong)."""
        high_risk_outcomes = [
            o for o in outcomes
            if o.risk_level in ("high", "critical")
        ]

        if not high_risk_outcomes:
            return 0.0

        false_positives = sum(1 for o in high_risk_outcomes if not o.was_correct)
        return false_positives / len(high_risk_outcomes)

    def _compute_false_negative_rate(self, outcomes: List[PredictionOutcome]) -> float:
        """Compute false negative rate (low-risk predictions that were wrong)."""
        low_risk_outcomes = [
            o for o in outcomes
            if o.risk_level in ("none", "low")
        ]

        if not low_risk_outcomes:
            return 0.0

        false_negatives = sum(1 for o in low_risk_outcomes if not o.was_correct)
        return false_negatives / len(low_risk_outcomes)

    def _compute_precision_by_severity(
        self,
        outcomes: List[PredictionOutcome]
    ) -> Dict[str, float]:
        """Compute precision for each severity level."""
        by_severity = defaultdict(list)
        for o in outcomes:
            by_severity[o.risk_level].append(o)

        precision = {}
        for severity, sevgroup in by_severity.items():
            if not sevgroup:
                continue
            correct = sum(1 for o in sevgroup if o.was_correct)
            precision[severity] = correct / len(sevgroup)

        return precision


# ============================================================================
# Persistence
# ============================================================================

class CalibrationDataManager:
    """Manage persistence of calibration data."""

    def __init__(self, file_path: Path = Path("data/calibration_data.json")):
        self.file_path = file_path

    def save_metrics(self, metrics: ModelPerformanceMetrics):
        """Save performance metrics to disk."""
        self.file_path.parent.mkdir(parents=True, exist_ok=True)

        with open(self.file_path, 'w') as f:
            json.dump(asdict(metrics), f, indent=2)

    def load_metrics(self) -> Optional[ModelPerformanceMetrics]:
        """Load performance metrics from disk."""
        if not self.file_path.exists():
            return None

        try:
            with open(self.file_path, 'r') as f:
                data = json.load(f)
            return ModelPerformanceMetrics(**data)
        except Exception as e:
            print(f"Error loading calibration data: {e}")
            return None


# ============================================================================
# Convenience Functions
# ============================================================================

def create_confidence_calibrator() -> ConfidenceCalibrator:
    """Create a new confidence calibrator."""
    return ConfidenceCalibrator()


def create_drift_detector(
    baseline_accuracy: float = 0.85,
    baseline_std: float = 0.05
) -> DriftDetector:
    """Create a new drift detector."""
    return DriftDetector(
        baseline_accuracy=baseline_accuracy,
        baseline_std=baseline_std,
    )
