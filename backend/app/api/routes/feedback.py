from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel, Field
from typing import Dict, List, Optional
import json
from pathlib import Path
from datetime import datetime
import sys
sys.path.append(str(Path(__file__).parent.parent.parent))

from app.services.pharmacogenomics.feedback_learning import (
    BayesianFeedbackLearner,
    LearningPriorsManager,
    FeedbackEvent,
    load_learning_priors,
    save_learning_priors,
)

router = APIRouter()

PRIORS_FILE = Path("data/learning_priors.json")

class FeedbackRequest(BaseModel):
    gene: str = Field(..., description="Gene symbol")
    drug: str = Field(..., description="Drug name")
    reported_diplotype: str = Field(..., description="The diplotype originally reported by the system")
    correct_diplotype: str = Field(..., description="The correct diplotype as verified by clinician")
    comments: Optional[str] = Field(None, description="Optional comments or rationale")
    feedback_quality: float = Field(1.0, ge=0.0, le=1.0, description="Clinician confidence in feedback (0-1)")

class FeedbackResponse(BaseModel):
    message: str
    status: str
    prior_update: Optional[Dict] = Field(None, description="Details of prior update")

@router.post("/", response_model=FeedbackResponse)
async def submit_feedback(feedback: FeedbackRequest):
    """
    Submit clinical feedback to improve the engine using Bayesian learning.

    Features:
    - Bayesian posterior updates with bounded learning
    - Time-based decay of old feedback
    - Confidence-weighted feedback integration
    - Overfitting prevention
    """
    # Load current learning priors
    priors_manager = LearningPriorsManager(PRIORS_FILE)
    current_priors = priors_manager.load()

    # Create feedback event
    feedback_event = FeedbackEvent(
        gene=feedback.gene,
        reported_diplotype=feedback.reported_diplotype,
        correct_diplotype=feedback.correct_diplotype,
        timestamp=datetime.now(),
        feedback_quality=feedback.feedback_quality,
        comments=feedback.comments,
    )

    # Initialize Bayesian learner
    learner = BayesianFeedbackLearner(
        learning_rate=0.1,
        decay_rate=0.95,
        min_prior=0.80,
        max_prior=1.50,
        max_delta=0.10,
    )

    # Get current prior for this diplotype
    current_prior = current_priors.get_diplotype_prior(
        feedback.gene,
        feedback.correct_diplotype
    )

    # Calculate months since last calibration
    if current_priors.last_calibration:
        last_cal = datetime.fromisoformat(current_priors.last_calibration)
        months_since = feedback_event.months_since(last_cal)
    else:
        months_since = 0.0

    # Update prior using Bayesian learning
    new_prior, explanation = learner.update_prior(
        current_prior=current_prior,
        feedback_quality=feedback.feedback_quality,
        months_since_last_update=months_since,
    )

    # Apply incremental update
    updated_priors = learner.incremental_update(
        current_priors=current_priors,
        new_feedback=feedback_event,
    )

    # Save updated priors
    priors_manager.save(updated_priors)

    return FeedbackResponse(
        status="success",
        message=f"Feedback recorded and integrated using Bayesian learning for {feedback.gene}:{feedback.correct_diplotype}",
        prior_update={
            "gene": feedback.gene,
            "diplotype": feedback.correct_diplotype,
            "previous_prior": round(current_prior, 4),
            "new_prior": round(new_prior, 4),
            "explanation": explanation,
            "total_feedback_events": updated_priors.metadata.get("total_feedback_events", 0),
        }
    )
