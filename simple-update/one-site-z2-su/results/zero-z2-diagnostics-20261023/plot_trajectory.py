"""Plot existing diagnostics only; does not run Julia or optimize a state."""
import csv
import json
import math
from pathlib import Path

root = Path(__file__).resolve().parent.parent
rows = {}
for segment in (1, 2, 3):
    source = root / f"zero-z2-seed5-segment{segment}" / "trajectory.jsonl"
    for line in source.read_text().splitlines():
        data = json.loads(line)
        if data["event"] == "residual":
            rows[data["step"]] = data["r"]
out = Path(__file__).resolve().parent
target = out / "local-residuals.svg"
if target.exists():
    raise FileExistsError(target)
with (out / "local-residuals.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["step", "star", "plaquette"])
    writer.writerows((s, r["star"], r["plaquette"]) for s, r in sorted(rows.items()))
def x(step):
    return 75 + 640 * step / 700
def y(residual):
    return 65 + (0.3 - math.log10(max(residual, 1e-16))) / 16.3 * 300
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="450">',
       '<rect width="800" height="450" fill="white"/>',
       '<text x="30" y="25" font-size="17">Seed 5: 700 updates; ground-state acceptance failed</text>',
       '<text x="30" y="46" font-size="12">H=-sum(A)-sum(B), d=D=2, dt=0.05, ALS=12; same state across all segments</text>',
       '<path d="M75 60V365H715" stroke="black" fill="none"/>']
for value in (1., 1e-3, 1e-6, 1e-9, 1e-12, 1e-15):
    yy = y(value)
    svg.append(f'<text x="10" y="{yy+4}" font-size="12">{value:.0e}</text>')
for step in range(0, 701, 100):
    svg.append(f'<text x="{x(step)-10}" y="388" font-size="12">{step}</text>')
for key, color in (("star", "#2563eb"), ("plaquette", "#dc2626")):
    points = " ".join(f"{x(s):.3f},{y(r[key]):.3f}" for s, r in sorted(rows.items()))
    svg.append(f'<polyline points="{points}" stroke="{color}" stroke-width="2" fill="none"/>')
svg.extend([f'<path d="M75 {y(1e-6)}H715" stroke="gray" stroke-dasharray="5 4"/>',
            '<text x="245" y="415" font-size="13">Cumulative update count</text>',
            '<text x="75" y="440" font-size="12">Blue: star residual; red: plaquette residual. Dashed: required 1e-6 acceptance threshold.</text>',
            '</svg>'])
target.write_text("\n".join(svg), encoding="utf-8")
print(target)
