
import requests
import time
import os
import sys

# Configuration
API_URL = "http://localhost:8000/api/v1/upload/"
TEST_DATA_DIR = "/Users/praneshjs/Education/rift/PharmaGaurd/test_data"

# Test Scenarios
# Format: (Filename, Drug, Expected Diplotype, Expected Phenotype Pattern, Expected Risk Keywords)
SCENARIOS = [
    (
        "single_cyp2d6_hom_alt.vcf", "codeine",
        "*4/*4", "Poor Metabolizer",
        ["Avoid", "Toxic", "Contraindicated"]
    ),
    (
        "single_cyp2c19_hom_alt.vcf", "clopidogrel",
        "*2/*2", "Poor Metabolizer",
        ["Avoid", "Alternative", "Ineffective"]
    ),
    (
        "single_tpmt_het.vcf", "thioguanine",
        "*1/*3C", "Intermediate",
        ["Adjust Dosage", "Reduce Dose", "Lower Dose"]
    ),
    (
        "single_slco1b1_hom_alt.vcf", "simvastatin",
        "*5/*5", "Poor Function",
        ["High Myopathy Risk", "Alternative", "Lower Dose", "Adjust Dosage"]
    ),
     (
        "single_dpyd_het.vcf", "fluorouracil",
        "*1/*2A", "Intermediate",
        ["Adjust Dosage", "Toxicity Risk", "Reduce Dose", "Toxic"]
    )
]

def run_stress_test():
    print("=== Starting Pharmacogenomic Engine Stress Test ===\n")
    
    passed = 0
    failed = 0
    
    for filename, drug, exp_dip, exp_pheno_pattern, exp_risk_keywords in SCENARIOS:
        filepath = os.path.join(TEST_DATA_DIR, filename)
        if not os.path.exists(filepath):
            print(f"[SKIP] {filename} not found.")
            continue
            
        print(f"Testing {filename} with {drug}...")
        
        try:
            with open(filepath, 'rb') as f:
                files = {'file': (filename, f, 'text/plain')}
                params = {'drugs': [drug]}
                r = requests.post(API_URL, files=files, params=params)
            
            if r.status_code != 200:
                print(f"  [FAIL] Server returned {r.status_code}: {r.text}")
                failed += 1
                continue
                
            data = r.json()
            if not data:
                print("  [FAIL] Empty response")
                failed += 1
                continue
                
            # Assume single drug response
            res = data[0]
            pp = res['pharmacogenomic_profile']
            ra = res['risk_assessment']
            
            # Checks
            dip_check = pp.get('diplotype') == exp_dip
            pheno_check = exp_pheno_pattern.lower() in pp.get('phenotype', '').lower()
            risk_check = any(k.lower() in ra.get('risk_label', '').lower() for k in exp_risk_keywords)
            
            if dip_check and pheno_check and risk_check:
                print(f"  [PASS] {exp_dip} | {pp['phenotype']} | {ra['risk_label']}")
                passed += 1
            else:
                print(f"  [FAIL] Mismatch:")
                if not dip_check: print(f"    Diplotype: Expected {exp_dip}, Got {pp.get('diplotype')}")
                if not pheno_check: print(f"    Phenotype: Expected pattern '{exp_pheno_pattern}', Got '{pp.get('phenotype')}'")
                if not risk_check: print(f"    Risk: Expected one of {exp_risk_keywords}, Got '{ra.get('risk_label')}'")
                failed += 1

        except Exception as e:
            print(f"  [Error] {e}")
            failed += 1
            
        print("-" * 40)

    print(f"\nCompleted. Passed: {passed}, Failed: {failed}")
    if failed > 0:
        sys.exit(1)

if __name__ == "__main__":
    # Ensure server is up
    try:
        requests.get("http://localhost:8000/health")
    except:
        print("Server not running using existing process? Waiting 5s...")
        time.sleep(5)
        
    run_stress_test()
