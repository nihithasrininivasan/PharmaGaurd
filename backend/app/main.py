from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.router import api_router
from app.services.pharmacogenomics.cpic_loader import get_cpic_loader

app = FastAPI(
    title="PharmaGuard API",
    description="Pharmacogenomic Decision Engine with Adaptive Learning & Multi-Drug Risk Integration",
    version="2.0.0"
)

# CORS Configuration
origins = [
    "http://localhost",
    "http://localhost:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include API Router
app.include_router(api_router, prefix="/api/v1")

@app.on_event("startup")
async def startup_event():
    # Preload CPIC data
    print("Preloading CPIC data...")
    get_cpic_loader()

@app.get("/health")
async def health_check():
    return {"status": "ok", "service": "PharmaGuard"}
