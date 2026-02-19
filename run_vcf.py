import subprocess, pathlib, sys

result = subprocess.run(
    [sys.executable, "-m", "services.vcf", "test_data/manual_test.vcf"],
    capture_output=True, text=True, encoding="utf-8"
)
output = result.stdout
pathlib.Path("run_output.txt").write_text(output, encoding="utf-8")
print(f"Exit code: {result.returncode}")
print(f"Output length: {len(output)} chars, {output.count(chr(10))} lines")
if result.stderr:
    print("STDERR:", result.stderr[:300])
