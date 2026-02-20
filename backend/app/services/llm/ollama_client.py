import logging
import httpx
import backoff
import random
from typing import Optional

# Configure structured logging
logger = logging.getLogger(__name__)

# Global shared HTTP client — connection reuse across all requests (TURBO: 30s timeout)
_shared_client = httpx.AsyncClient(timeout=30.0)


class OllamaClient:
    """
    Client for interacting with the local Ollama instance running Llama3.
    Uses a global shared httpx.AsyncClient for connection reuse.
    TURBO MODE: Optimized for 2-3 second clinical responses.
    """
    def __init__(self, base_url: str = "http://127.0.0.1:11434", model: str = "llama3"):
        self.base_url = base_url
        self.model = model
        self.generate_endpoint = f"{self.base_url}/api/generate"

    @backoff.on_exception(
        backoff.expo,
        (httpx.RequestError, httpx.HTTPStatusError),
        max_tries=2,
        giveup=lambda e: isinstance(e, httpx.HTTPStatusError) and e.response.status_code < 500
    )
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
        
        try:
            response = await _shared_client.post(self.generate_endpoint, json=payload)
            response.raise_for_status()
            
            data = response.json()
            generated_text = data.get("response", "")
            
            logger.info("Ollama request successful", extra={"response_length": len(generated_text)})
            return generated_text
                
        except (httpx.HTTPStatusError, httpx.RequestError) as e:
            logger.error(f"Error communicating with Ollama: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in Ollama client: {str(e)}")
            return None

    @backoff.on_exception(
        backoff.expo,
        (httpx.RequestError, httpx.HTTPStatusError),
        max_tries=2,
        giveup=lambda e: isinstance(e, httpx.HTTPStatusError) and e.response.status_code < 500
    )
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

        try:
            response = await _shared_client.post(self.generate_endpoint, json=payload)
            response.raise_for_status()

            data = response.json()
            generated_text = data.get("response", "")

            logger.info("Chat Ollama request successful", extra={"response_length": len(generated_text)})
            return generated_text

        except (httpx.HTTPStatusError, httpx.RequestError) as e:
            logger.error(f"Error communicating with Ollama (chat): {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in Ollama chat client: {str(e)}")
            return None



