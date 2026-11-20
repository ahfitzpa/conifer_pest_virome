from pathlib import Path
import urllib.request

# Repo root derived from this file, not a hardcoded developer path: these
# scripts run on the cluster as well as locally.
ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "references"
REF.mkdir(parents=True, exist_ok=True)
url = "https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260729.xlsx"
out = REF / "VMR_MSL41.v1.20260729.xlsx"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=120) as r, open(out, "wb") as f:
    f.write(r.read())
print("saved", out, out.stat().st_size, "bytes")
