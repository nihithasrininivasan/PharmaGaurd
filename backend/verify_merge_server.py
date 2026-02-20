import subprocess
import time
import sys
import json
import urllib.request
import urllib.parse
import os
import signal

# Minimal VCF content for testing (CYP2D6 *4/*4 dummy)
# rs1065852 (C>T) is *4 defining variant
VCF_CONTENT = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr22>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr22\t42128945\trs16947\tC\tT\t.\tPASS\t.\tGT\t1/1
chr22\t42130692\trs1065852\tG\tA\t.\tPASS\t.\tGT\t1/1
"""

def main():
    print("Starting Uvicorn server...")
    # diverse start methods to be robust
    server_process = subprocess.Popen(
        [sys.executable, "-m", "uvicorn", "app.main:app", "--port", "8001", "--no-access-log"],
        cwd=os.getcwd(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    try:
        # Wait for server to start
        print("Waiting for server to be ready...")
        for _ in range(20):
            try:
                with urllib.request.urlopen("http://127.0.0.1:8001/health", timeout=1) as response:
                    if response.status == 200:
                        print("Server is ready!")
                        break
            except Exception:
                time.sleep(1)
        else:
            print("Server failed to start.")
            stdout, stderr = server_process.communicate(timeout=1)
            print("STDOUT:", stdout.decode())
            print("STDERR:", stderr.decode())
            return

        # Prepare multipart request using standard library (since requests/httpx missing)
        # This is painful with urllib, so we'll mock the internal call instead?
        # No, better test end-to-end. I'll construct a simple multipart body manually.
        
        boundary = "---------------------------12345678901234567890123456"
        body = []
        
        # Drug field
        body.append(f"--{boundary}")
        body.append('Content-Disposition: form-data; name="drug"')
        body.append('')
        body.append('CODEINE')
        
        # Patient ID field
        body.append(f"--{boundary}")
        body.append('Content-Disposition: form-data; name="patient_id"')
        body.append('')
        body.append('TEST_PATIENT')

        # VCF File field
        body.append(f"--{boundary}")
        body.append('Content-Disposition: form-data; name="vcf"; filename="test.vcf"')
        body.append('Content-Type: text/plain')
        body.append('')
        body.append(VCF_CONTENT)
        
        body.append(f"--{boundary}--")
        body.append('')
        
        body_str = "\r\n".join(body).encode('utf-8')
        
        req = urllib.request.Request(
            "http://127.0.0.1:8001/api/v1/analyze",
            data=body_str,
            headers={
                "Content-Type": f"multipart/form-data; boundary={boundary}",
                "Content-Length": len(body_str)
            },
            method="POST"
        )
        
        print("\nSending POST /api/v1/analyze...")
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                resp_json = json.loads(response.read().decode('utf-8'))
                print("\nResponse Received:")
                print(json.dumps(resp_json, indent=2))
                
                # Validation
                ra = resp_json.get("risk_assessment", {})
                print("\nValidation:")
                print(f"Risk Score: {ra.get('risk_score')} (Expected: not None)")
                print(f"Risk Level: {ra.get('risk_level')} (Expected: not None)")
                print(f"Confidence: {ra.get('confidence_score')}")
                
                if ra.get('risk_score') is not None:
                    print("✅ SUCCESS: Risk fields are present.")
                else:
                    print("❌ FAILURE: Risk fields missing.")
                    
        except urllib.error.HTTPError as e:
            print(f"HTTP Error: {e.code}")
            print(e.read().decode())
        except Exception as e:
            print(f"Request Error: {e}")

    finally:
        print("\nStopping server...")
        server_process.terminate()
        try:
            stdout, stderr = server_process.communicate(timeout=5)
            print("\n----- SERVER STDOUT -----")
            print(stdout.decode(errors='replace'))
            print("\n----- SERVER STDERR (Errors/Tracebacks) -----")
            print(stderr.decode(errors='replace'))
        except Exception as e:
            print(f"Error capturing logs: {e}")
            try:
                server_process.kill()
            except:
                pass

if __name__ == "__main__":
    main()
