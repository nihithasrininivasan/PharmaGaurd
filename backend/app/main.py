from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.router import api_router
from app.api.routes import analysis
from app.core import logging  # Initialize logging
from app.services.pharmacogenomics.cpic_loader import get_cpic_loader
from app.services.llm.groq_client import GroqClient

app = FastAPI(
    title="PharmaGuard API",
    description="Pharmacogenomic Decision Engine with Adaptive Learning & Multi-Drug Risk Integration",
    version="2.0.0"
)

# CORS Configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include API Routers
app.include_router(api_router, prefix="/api/v1")
app.include_router(analysis.router, prefix="/api/v1", tags=["Analysis"])

@app.on_event("startup")
async def startup_event():
    # Preload CPIC data
    print("Preloading CPIC data...")
    get_cpic_loader()

    # Verify Groq connectivity (non-blocking)
    try:
        client = GroqClient()
        result = await client.generate_text("Respond with OK.")
        if result:
            print("🔥 Groq LLM connected successfully.")
        else:
            print("⚠️ Groq returned empty response — check GROQ_API_KEY.")
    except Exception as e:
        print(f"⚠️ Groq not reachable ({e}). LLM explanations will use fallback. Server starting anyway.")

@app.get("/health")
async def health_check():
    return {"status": "ok", "service": "PharmaGuard"}
