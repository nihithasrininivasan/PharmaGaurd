FROM python:3.11-slim

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code and data
COPY app/ ./app/
COPY data/ ./data/

# Expose port
EXPOSE 10000

# Start server (Render uses port 10000 by default)
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "10000"]
