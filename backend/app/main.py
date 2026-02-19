from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.routes import analysis
from app.core import logging  # Initialize logging
from app.services.llm.ollama_client import OllamaClient

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

@app.on_event("startup")
async def warmup_llm():
    """Preload the LLM model into memory on server startup."""
    try:
        client = OllamaClient()
        await client.generate_text("Warmup request. Respond with OK.")
        print("🔥 Ollama model warmed up.")
    except Exception:
        print("⚠️ Ollama warmup failed. First request may be slow.")

@app.get("/health")
async def health_check():
    return {"status": "ok"}
