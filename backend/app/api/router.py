from fastapi import APIRouter
from app.api.routes import upload, pharmacogenomics, feedback, polypharmacy

api_router = APIRouter()

api_router.include_router(upload.router, prefix="/upload", tags=["Upload"])
api_router.include_router(pharmacogenomics.router, prefix="/pharmacogenomics", tags=["Pharmacogenomics"])
api_router.include_router(feedback.router, prefix="/feedback", tags=["Feedback"])
api_router.include_router(polypharmacy.router, prefix="/polypharmacy", tags=["Polypharmacy"])
