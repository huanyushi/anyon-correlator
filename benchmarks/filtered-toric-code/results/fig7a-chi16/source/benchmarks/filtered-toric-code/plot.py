"""Draw the saved Fig. 7(a) benchmark; no numerical computation is launched."""
import argparse
import csv
import json
import math
from pathlib import Path


def plot(run, output):
    settings = json.loads((run / "settings.json").read_text())
    with (run / "spectrum.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    x = lambda k: 82 + (k / math.pi + 1) * 330
    y = lambda e: 472 - e * 64
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="565" viewBox="0 0 800 565">',
           '<rect width="800" height="565" fill="white"/>',
           '<g font-family="sans-serif">',
           '<text x="82" y="30" font-size="21">Filtered toric code — Fig. 7(a)</text>',
           f'<text x="82" y="55" font-size="14">βz=0.55, βx=0.45; χ={settings["chi"]}; {settings["nk"]} momenta; {settings["num"]} modes/sector</text>',
           '<text x="82" y="77" font-size="13" fill="#2563eb">+ 1    ○ e</text>',
           '<text x="225" y="77" font-size="13" fill="#dc2626">+ m    ○ ε</text>']
    for t in range(7):
        yy = y(t)
        svg += [f'<path d="M82 {yy}H742" stroke="#e5e7eb"/>',
                f'<text x="66" y="{yy+5}" text-anchor="end" font-size="13">{t}</text>']
    for t in (-1, -.5, 0, .5, 1):
        xx = x(t * math.pi)
        svg += [f'<path d="M{xx} 472v6" stroke="black"/>',
                f'<text x="{xx}" y="501" text-anchor="middle" font-size="13">{t:g}</text>']
    svg += ['<path d="M82 88V472H742V88Z" fill="none" stroke="black"/>',
            '<text x="410" y="532" text-anchor="middle" font-size="16">Momentum k/π</text>',
            '<text transform="translate(24,280) rotate(-90)" text-anchor="middle" font-size="16">−log|λ/λ₀|</text>']
    visible = 0
    for r in rows:
        e = float(r["quasi_energy"])
        if not math.isfinite(e) or not 0 <= e <= 6:
            continue
        xx, yy = x(float(r["k"])), y(e)
        color = "#2563eb" if r["g"] == "0" else "#dc2626"
        svg.append(f'<g fill="none" stroke="{color}" stroke-width="0.9">')
        if r["q"] == "1":
            svg.append(f'<circle cx="{xx:.3f}" cy="{yy:.3f}" r="2.4"/>')
        else:
            svg.append(f'<path d="M{xx-2.4:.3f} {yy:.3f}h4.8 M{xx:.3f} {yy-2.4:.3f}v4.8"/>')
        svg.append(f'<title>g={r["g"]}, q={r["q"]}, branch={r["branch"]}, ε={e:.6g}</title></g>')
        visible += 1
    svg.append(f'<path d="M{x(0)-3} 472h6 M{x(0)} 469v6" stroke="#2563eb"><title>Fixed point, λ/λ₀=1</title></path>')
    svg += [f'<text x="82" y="556" font-size="11" fill="#555">Qualitative benchmark; {visible}/{len(rows)} modes in displayed range. Full data are retained in spectrum.csv.</text>',
            '</g></svg>']
    if output.exists():
        raise FileExistsError(output)
    output.write_text("\n".join(svg), encoding="utf-8")
    print(f"Saved {output}: {visible}/{len(rows)} modes within ε≤6")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    plot(args.run, args.out)
