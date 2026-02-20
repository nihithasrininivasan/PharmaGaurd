import sys
import os
from pathlib import Path

# Add backend to path so we can import app
sys.path.append(os.path.join(os.getcwd(), "backend"))

from fastapi.testclient import TestClient
from app.main import app

client = TestClient(app)

def test_vcf_upload():
    vcf_path = Path("test_data/manual_test.vcf")
    if not vcf_path.exists():
        print("VCF not found. Run generate_test_vcf.py first.")
        return

    print(f"Testing upload with {vcf_path}...")
    with open(vcf_path, "rb") as f:
        # manual_test.vcf has variants for CYP2C19, CYP2D6, CYP2C9, DPYD
        # We filter for clopidogrel (CYP2C19)
        response = client.post(
            "/api/v1/upload/", 
            files={"file": ("manual.vcf", f, "text/plain")},
            params={"drugs": ["clopidogrel"]} 
        )
    
    if response.status_code == 200:
        print("Success! Reviewing response...")
        data = response.json()
        import json
        print(json.dumps(data, indent=2))
        
        # Validation
        if len(data) == 1:
            print("Verified: Exactly 1 report returned.")
        else:
            print(f"Error: Expected 1 report, got {len(data)}")
            
        report = data[0]
        if report["drug"] == "CLOPIDOGREL":
             print("Verified: Drug is CLOPIDOGREL.")
        else:
             print(f"Error: Drug mismatch {report['drug']}")
             
        if report["pharmacogenomic_profile"]["primary_gene"] == "CYP2C19":
             print("Verified: Gene is CYP2C19.")
             
        variants = report["pharmacogenomic_profile"]["detected_variants"]
        print(f"Found {len(variants)} variants for CYP2C19.")
        if len(variants) > 0:
            print("Verified: Variants detected.")
        else:
            print("Warning: No variants detected? (Check manual_test.vcf content)")

    else:
        print(f"Failed: {response.status_code}")
        print(response.text)

if __name__ == "__main__":
    test_vcf_upload()
