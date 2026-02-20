import logging
import httpx
import backoff
import random
import os
from typing import Optional
from dotenv import load_dotenv, find_dotenv

# Load .env file (walks up directories to find it)
load_dotenv(find_dotenv())

logger = logging.getLogger(__name__)

# Groq API config
GROQ_API_KEY = os.environ.get("GROQ_API_KEY", "")
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
GROQ_MODEL = os.environ.get("GROQ_MODEL", "llama-3.1-8b-instant")

# Shared HTTP client for connection reuse
_shared_client = httpx.AsyncClient(timeout=30.0)


class GroqClient:
    """
    Client for Groq's hosted Llama 3.1 API (OpenAI-compatible).
    Drop-in replacement for OllamaClient.
    """

    def __init__(self, api_key: str = None, model: str = None):
        self.api_key = api_key or GROQ_API_KEY
        self.model = model or GROQ_MODEL

    @backoff.on_exception(
        backoff.expo,
        (httpx.RequestError, httpx.HTTPStatusError),
        max_tries=2,
        giveup=lambda e: isinstance(e, httpx.HTTPStatusError) and e.response.status_code < 500
    )
    async def generate_text(self, prompt: str) -> Optional[str]:
        """
        Generates deterministic clinical explanations via Groq.
        Low temperature for consistent, factual responses.
        """
        logger.info("Sending request to Groq", extra={"model": self.model})

        payload = {
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 60,
            "temperature": 0.1,
            "top_p": 0.85,
        }

        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

        try:
            response = await _shared_client.post(GROQ_API_URL, json=payload, headers=headers)
            response.raise_for_status()

            data = response.json()
            generated_text = data["choices"][0]["message"]["content"]

            logger.info("Groq request successful", extra={"response_length": len(generated_text)})
            return generated_text

        except (httpx.HTTPStatusError, httpx.RequestError) as e:
            logger.error(f"Error communicating with Groq: {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in Groq client: {str(e)}")
            return None

    @backoff.on_exception(
        backoff.expo,
        (httpx.RequestError, httpx.HTTPStatusError),
        max_tries=2,
        giveup=lambda e: isinstance(e, httpx.HTTPStatusError) and e.response.status_code < 500
    )
    async def generate_chat_text(self, prompt: str) -> Optional[str]:
        """
        Generates chatbot responses with slight variance via Groq.
        """
        logger.info("Sending CHAT request to Groq", extra={"model": self.model})

        payload = {
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 80,
            "temperature": 0.15,
            "top_p": 0.85,
            "seed": random.randint(1, 999999),
        }

        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

        try:
            response = await _shared_client.post(GROQ_API_URL, json=payload, headers=headers)
            response.raise_for_status()

            data = response.json()
            generated_text = data["choices"][0]["message"]["content"]

            logger.info("Chat Groq request successful", extra={"response_length": len(generated_text)})
            return generated_text

        except (httpx.HTTPStatusError, httpx.RequestError) as e:
            logger.error(f"Error communicating with Groq (chat): {str(e)}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in Groq chat client: {str(e)}")
            return None
