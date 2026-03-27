import re
from pathlib import Path

p = Path(__file__).resolve().parents[1] / "Duran-FerrerM-evoflux/docs/Data_source_Fig.5.html"
t = p.read_text(encoding="utf-8", errors="replace")
chunks = re.findall(r'<pre class="r"><code>(.*?)</code></pre>', t, re.DOTALL)
out = Path(__file__).resolve().parent / "_fig5_r_extracted.txt"
out.write_text("\n\n# ===== CHUNK BREAK =====\n\n".join(chunks), encoding="utf-8")
print("n_chunks", len(chunks), "bytes", out.stat().st_size)
