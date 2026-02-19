from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.routes import analysis
from app.core import logging  # Initialize logging

app = FastAPI(
    title="PharmaGuard Backend",
    description="Pharmacogenomics Analysis Service",
    version="1.0.0"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(analysis.router, prefix="/api/v1", tags=["Analysis"])

@app.get("/health")
async def health_check():
    return {"status": "ok"}
