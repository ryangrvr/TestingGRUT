"""Build a portable HTML reader with local SVG/CSS embedded and native MathML."""
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
text=(ROOT/"GRUT_WORKING_THEORY.md").read_text()
text=text.replace("## 1. What is being completed",
    "![A schematic map of GRUT's dependencies and open links](EMERGENCE_MAP.svg)\n\n## 1. What is being completed",1)
for fn in ["CLAIMS.md","SOURCES.md","AUDIT.md"]:
    text+="\n\n"+(ROOT/fn).read_text()
text=text.replace("(CLAIMS.md)","(#claim-ledger)").replace("(SOURCES.md)","(#sources-and-provenance)")
text=text.replace("(AUDIT.md)","(#audit-and-limitations)")
with tempfile.TemporaryDirectory() as temp:
    p=Path(temp)/"reader.md";p.write_text(text)
    run=subprocess.run(["pandoc",str(p),"-f","markdown+tex_math_single_backslash",
                    "--standalone","--mathml","--toc","--toc-depth=2",
                    "--embed-resources","--resource-path",str(ROOT),
                    "--css",str(ROOT/"reader.css"),"-o",str(ROOT/"GRUT_WORKING_THEORY.html")],
                   check=True,cwd=ROOT,capture_output=True,text=True)
    if run.stderr.strip():
        raise RuntimeError("Reader emitted a warning: "+run.stderr)
print("Built GRUT_WORKING_THEORY.html")
