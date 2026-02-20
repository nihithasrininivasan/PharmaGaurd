from fastapi import APIRouter, File, UploadFile, HTTPException, Query
from typing import List, Optional

from app.services.vcf.pharmaguard_adapter import (
    analyze_vcf_for_drugs, 
    SUPPORTED_DRUGS_TO_GENE
)
from app.services.vcf.parser import VcfParseError

router = APIRouter()

@router.post("/", response_model=List[dict])
async def upload_vcf(
    file: UploadFile = File(...),
    drugs: Optional[List[str]] = Query(None, description="List of drugs to analyze. Defaults to all supported drugs.")
):
    """
    Upload a VCF file (uncompressed .vcf or gzip .vcf.gz) to generate a Pharmacogenomic Report.
    
    - **file**: The VCF file containing variant data.
    - **drugs**: Optional list of drug names to filter the report. If omitted, all supported drugs are analyzed.
    """
    if not (file.filename.endswith(".vcf") or file.filename.endswith(".vcf.gz")):
        raise HTTPException(status_code=400, detail="Invalid file format. Please upload .vcf or .vcf.gz")

    try:
        content = await file.read()
        
        # Determine target drugs
        target_drugs = drugs if drugs else list(SUPPORTED_DRUGS_TO_GENE.keys())
        
        # Analyze using the VCF adapter (which now uses the real engine)
        results = analyze_vcf_for_drugs(content, target_drugs)
        
        return results
        
    except VcfParseError as e:
        raise HTTPException(status_code=400, detail=f"VCF Parsing Error: {str(e)}")
    except Exception as e:
        print(f"Error processing VCF upload: {e}")
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")

