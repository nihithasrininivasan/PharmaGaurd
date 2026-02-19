import logging
import httpx
import backoff
from typing import Optional

# Configure structured logging
logger = logging.getLogger(__name__)

class OllamaClient:
    """
    Client for interacting with the local Ollama instance running Llama3.
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
        Generates text using the Ollama API with retries and timeout.
        
        Args:
            prompt: The text prompt to send to the model.
            
        Returns:
            The generated text string, or None if request fails after retries.
        """
        logger.info("Sending request to Ollama", extra={"model": self.model})
        
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "keep_alive": "5m",  # Optimization: Keep model warm
            "options": {
                "num_predict": 120, # Optimization: Limit output tokens
                "temperature": 0.2, # Optimization: Deterministic
                "top_p": 0.9        # Optimization: Focus on probable tokens
            }
        }
        
        try:
             async with httpx.AsyncClient(timeout=60.0) as client:
                response = await client.post(self.generate_endpoint, json=payload)
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
