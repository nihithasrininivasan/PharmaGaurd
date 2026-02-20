import requests
import json
import time
from typing import List, Dict

API_URL = "http://localhost:8000/api/v1/pharmacogenomics/report"

def run_test_case(name: str, drug: str, gene: str, variants: List[Dict]):
    print(f"\n{'='*60}")
    print(f"TEST CASE: {name}")
    print(f"Drug: {drug} | Gene: {gene}")
    print(f"Variants: {len(variants)} provided")
    print(f"{'-'*60}")

    payload = {
        "drug": drug,
        "gene": gene, 
        "patient_id": "TEST_USER_001",
        # Note: The API currently accepts variants in the body but query params for drug/gene
        # We need to construct the URL with params and body
    }
    
    # Construct query params
    url = f"{API_URL}?drug={drug}&gene={gene}&patient_id=TEST_PATIENT_{name.replace(' ', '_')}"
    
    try:
        start_time = time.time()
        response = requests.post(url, json=variants)
        duration = time.time() - start_time
        
        if response.status_code == 200:
            data = response.json()
            
            # Print key results
            profile = data.get("pharmacogenomic_profile", {})
            risk = data.get("risk_assessment", {})
            rec = data.get("clinical_recommendation", {})
            
            print(f"✅ SUCCESS ({duration:.2f}s)")
            print(f"  Diplotype: {profile.get('diplotype')}")
            print(f"  Phenotype: {profile.get('phenotype')}")
            print(f"  Conf:      {risk.get('confidence_score')}")
            print(f"  Risk:      {risk.get('risk_label')}")
            print(f"  Severity:  {risk.get('severity')}")
            print(f"  Rec:       {rec.get('text')}")
            
            # Full JSON dump (optional, condensed)
            # print("\nFull Response:")
            # print(json.dumps(data, indent=2))
        else:
            print(f"❌ FAILED (Status {response.status_code})")
            print(response.text)
            
    except Exception as e:
        print(f"❌ ERROR: {e}")

def main():
    print("Running Pharmacogenomics API Tests...")
    print("Ensure server is running: uvicorn app.main:app --host 0.0.0.0 --port 8000")
    
    # Case 1: CYP2C19 *1/*1 (Normal Metabolizer) - No variants
    run_test_case(
        "CYP2C19 Normal Metabolizer",
        "clopidogrel", 
        "CYP2C19",
        [] # No variants = wildtype *1/*1
    )

    # Case 2: CYP2C19 *1/*2 (Intermediate Metabolizer)
    # *2 defined by 4 variants. Providing all to ensure robust call.
    run_test_case(
        "CYP2C19 Intermediate (*1/*2)",
        "clopidogrel",
        "CYP2C19",
        [
            {"chrom": "chr10", "pos": 94775367, "ref": "A", "alt": "G", "zygosity": "HET", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94775507, "ref": "G", "alt": "A", "zygosity": "HET", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A", "zygosity": "HET", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94842866, "ref": "A", "alt": "G", "zygosity": "HET", "quality": 100.0, "filter": "PASS"}
        ]
    )

    # Case 3: CYP2C19 *2/*2 (Poor Metabolizer)
    run_test_case(
        "CYP2C19 Poor Metabolizer (*2/*2)",
        "clopidogrel",
        "CYP2C19",
        [
            {"chrom": "chr10", "pos": 94775367, "ref": "A", "alt": "G", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94775507, "ref": "G", "alt": "A", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94842866, "ref": "A", "alt": "G", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"}
        ]
    )

    # Case 4: CYP2D6 *4/*4 (Poor Metabolizer)
    # *4 defined by 42126611 C>T (among others, this is a key SNP)
    # Note: Adjust position to match what's in your specific cache if needed. 
    # Based on standard GRCh38, but let's check one from the cache logic if possible.
    # Assuming the ETL parsed standard positions.
    # We'll use a likely position for *4 (splicing defect 1846G>A is common, or 100C>T)
    # Let's try to simulate *4 with the 1846G>A (rs3892097) equivalent if mapped, 
    # but exact position depends on the reference used in the Excel file.
    # For now, let's test CYP2C19 thoroughly as we visualized it.
    
    # Case 5: Warfarin / CYP2C9 *3 (Decreased Function)
    # Need to check if CYP2C9 *3 is in our cache. Usually 42530188 or similar in GRCh38.
    
    print("\n" + "="*60)
    print("TESTTING FEEDBACK LOOP")
    print("="*60)
    
    # Simulate correct call but low confidence initially
    # *2/*2 (Poor Metabolizer)
    variants_pm = [
            {"chrom": "chr10", "pos": 94775367, "ref": "A", "alt": "G", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94775507, "ref": "G", "alt": "A", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"},
            {"chrom": "chr10", "pos": 94842866, "ref": "A", "alt": "G", "zygosity": "HOM_ALT", "quality": 100.0, "filter": "PASS"}
    ]
    
    print("1. Initial Call (Should have honest confidence)")
    # run_test_case will print the confidence
    run_test_case("CYP2C19 Poor Metabolizer (Initial)", "clopidogrel", "CYP2C19", variants_pm)
    
    print("\n2. Submitting Feedback (Correcting/Verifying *2/*2)")
    feedback_url = "http://localhost:8000/api/v1/feedback/"
    feedback_data = {
        "gene": "CYP2C19",
        "drug": "clopidogrel",
        "reported_diplotype": "*2/*2",
        "correct_diplotype": "*2/*2", # Confirming it's correct to boost confidence
        "comments": "Clinically verified"
    }
    try:
        resp = requests.post(feedback_url, json=feedback_data)
        print(f"Feedback Status: {resp.status_code}")
        print(resp.json())
    except Exception as e:
        print(f"Feedback Error: {e}")
        
    print("\n3. Second Call (Should have BOOSTED confidence)")
    run_test_case("CYP2C19 Poor Metabolizer (Post-Feedback)", "clopidogrel", "CYP2C19", variants_pm)

    print("\nDone.")

if __name__ == "__main__":
    main()
