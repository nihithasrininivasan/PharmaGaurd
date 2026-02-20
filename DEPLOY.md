# Deploying PharmaGaurd Backend

## Docker Deployment (Recommended)

1. **Build the Image**
   ```bash
   cd backend
   docker build -t pharmaguard-backend .
   ```

2. **Run the Container**
   ```bash
   docker run -p 8000:8000 pharmaguard-backend
   ```
   The API will be available at `http://localhost:8000`.

## Cloud Hosting (Render, Railway, Heroku)

This project is ready for cloud deployment.
- **Build Command**: `pip install -r requirements.txt`
- **Start Command**: `uvicorn app.main:app --host 0.0.0.0 --port $PORT`

## VCF Validation
Strict VCF validation is enabled. Only files starting with `##fileformat=VCF` are accepted, ensuring ISGCR compliance.
