import logging
import json
import time
import random
import urllib.request
import urllib.error
from typing import Optional, Dict, Any

# Configure structured logging
logger = logging.getLogger(__name__)

class OllamaClient:
    """
    Client for interacting with the local Ollama instance running Llama3.
    Uses standard library urllib for zero-dependency operation.
    TURBO MODE: Optimized for 2-3 second clinical responses.
    """
    def __init__(self, base_url: str = "http://127.0.0.1:11434", model: str = "llama3"):
        self.base_url = base_url.rstrip("/")
        self.model = model
        self.generate_endpoint = f"{self.base_url}/api/generate"

    async def generate_text(self, prompt: str) -> Optional[str]:
        """
        Generates deterministic clinical explanations.
        TURBO: num_predict=60, num_ctx=1024 for fast inference.
        """
        logger.info("Sending request to Ollama", extra={"model": self.model})
        
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "keep_alive": "15m",
            "options": {
                "num_predict": 40,
                "temperature": 0.1,
                "top_p": 0.85,
                "repeat_penalty": 1.1,
                "num_ctx": 1024
            }
        }
        
        return self._make_request(payload)

    async def generate_chat_text(self, prompt: str) -> Optional[str]:
        """
        Generates chatbot responses with slight variance.
        TURBO: Same speed params + random seed for unique answers.
        """
        logger.info("Sending CHAT request to Ollama", extra={"model": self.model})

        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "keep_alive": "15m",
            "options": {
                "num_predict": 60,
                "temperature": 0.15,
                "top_p": 0.85,
                "repeat_penalty": 1.1,
                "num_ctx": 1024,
                "seed": random.randint(1, 999999)
            }
        }

        return self._make_request(payload)

    def _make_request(self, payload: Dict[str, Any], retries: int = 2) -> Optional[str]:
        """Helper to make HTTP request with simple retry logic."""
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            self.generate_endpoint,
            data=data,
            headers={"Content-Type": "application/json"},
            method="POST"
        )

        for attempt in range(retries + 1):
            try:
                with urllib.request.urlopen(req, timeout=30) as response:
                    if response.status == 200:
                        resp_data = json.loads(response.read().decode("utf-8"))
                        generated_text = resp_data.get("response", "")
                        logger.info("Ollama request successful", extra={"response_length": len(generated_text)})
                        return generated_text
                    else:
                        logger.error(f"Ollama returned status {response.status}")
            except urllib.error.URLError as e:
                logger.error(f"Error communicating with Ollama: {e}")
                if attempt < retries:
                    time.sleep(1 * (attempt + 1))  # Simple backoff
                else:
                    return None
            except Exception as e:
                logger.error(f"Unexpected error in Ollama client: {e}")
                return None
        return None
