import requests
import sys

def test_api():
    url = "http://127.0.0.1:8002/api/v1/analyze"
    files = {
        'vcf': ('test.vcf', open('C:/Users/sarve/Downloads/serve/PharmaGaurd_merge/test.vcf', 'rb'), 'text/plain')
    }
    data = {
        'drug': 'CODEINE',
        'patient_id': 'TEST_PATIENT'
    }
    
    try:
        print(f"Sending request to {url}...")
        response = requests.post(url, files=files, data=data)
        print(f"Status Code: {response.status_code}")
        if response.status_code == 200:
            print("Response JSON:")
            print(response.json())
            print("SUCCESS: API reachable and working.")
        else:
            print("Error Response:")
            print(response.text)
            print("FAILURE: API returned error.")
    except Exception as e:
        print(f"EXCEPTION: {e}")
        print("FAILURE: Could not connect to API.")

if __name__ == "__main__":
    test_api()
