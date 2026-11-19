from pathlib import Path
import urllib.request

ROOT = Path(r"C:\Users\Amy\OneDrive - University College Dublin\Documents\conifer_pest_virome")
REF = ROOT / "reference_data"
REF.mkdir(exist_ok=True)
url = "https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260729.xlsx"
out = REF / "VMR_MSL41.v1.20260729.xlsx"
req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=120) as r, open(out, "wb") as f:
    f.write(r.read())
print("saved", out, out.stat().st_size, "bytes")
