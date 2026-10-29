"""Render eight saved spectra and the paper's 4×2 layout; standard library only."""
import argparse
import csv
import json
import math
from pathlib import Path

COLORS = {"blue": "#164bd4", "red": "#d62728", "green": "#188228"}
ORDER = "aebfcgdh"


def load(run):
    settings = json.loads((run / "settings.json").read_text())
    with (run / "spectrum.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if (run / "sectors.json").exists():
        styles = {s["name"]: s for s in json.loads((run / "sectors.json").read_text())}
    else:  # Read the preserved Fig.7(a) data directly.
        styles = {f"g{g}q{q}": {"color": "blue" if g == 0 else "red",
                  "marker": "plus" if q == 0 else "circle"} for g in (0, 1) for q in (0, 1)}
        for row in rows:
            row["sector"] = f'g{row["g"]}q{row["q"]}'
    return settings, rows, styles


def marker(x, y, color, kind, size=2.3):
    if kind == "circle":
        shape = f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{size}"/>'
    elif kind == "square":
        shape = f'<rect x="{x-size:.3f}" y="{y-size:.3f}" width="{2*size}" height="{2*size}"/>'
    elif kind == "dot":
        shape = f'<circle cx="{x:.3f}" cy="{y:.3f}" r="1.9" fill="{color}"/>'
    else:
        shape = f'<path d="M{x-size:.3f} {y:.3f}h{2*size} M{x:.3f} {y-size:.3f}v{2*size}"/>'
    return f'<g stroke="{color}" stroke-width="0.8" fill="none">{shape}</g>'


def panel_svg(label, run):
    settings, rows, styles = load(run)
    x = lambda k: 65 + (k / math.pi + 1) * 250
    y = lambda e: 365 - e * 48
    svg = [f'<text x="18" y="30" font-size="25" font-weight="bold">({label})</text>',
           f'<text x="85" y="28" font-size="15">βz={settings["betaz"]:g}, βx={settings["betax"]:g}; χ={settings["chi"]}</text>']
    phase = settings.get("phase", "topological")
    legend = {"topological": "blue: + 1, ○ e     red: + m, ○ ε",
              "charge_condensed": "blue: local    red: ket string    green: ket+bra string",
              "flux_condensed": "blue: + (++), ○ (+−)/(−+), □ (−−)"}[phase]
    svg.append(f'<text x="65" y="53" font-size="12">{legend}</text>')
    for e in range(7):
        svg += [f'<path d="M65 {y(e)}H565" stroke="#e5e7eb"/>',
                f'<text x="54" y="{y(e)+4}" text-anchor="end" font-size="13">{e}</text>']
    for k in (-1, -.5, 0, .5, 1):
        svg.append(f'<text x="{x(k*math.pi)}" y="389" text-anchor="middle" font-size="13">{k:g}</text>')
    svg += ['<path d="M65 77V365H565V77Z" stroke="black" fill="none"/>',
            '<text x="315" y="415" text-anchor="middle" font-size="15">Momentum k/π</text>',
            '<text transform="translate(19,220) rotate(-90)" text-anchor="middle" font-size="15">−log|λ/λ₀|</text>']
    visible = 0
    for row in rows:
        e = float(row["quasi_energy"])
        if not math.isfinite(e) or not 0 <= e <= 6:
            continue
        style = styles[row["sector"]]
        svg.append(marker(x(float(row["k"])), y(e), COLORS[style["color"]], style["marker"]))
        visible += 1
    svg.append(marker(x(0), y(0), COLORS["blue"], "dot" if phase == "charge_condensed" else "plus"))
    svg.append(f'<text x="65" y="440" font-size="11" fill="#555">{settings["nk"]} momenta; {visible} modes in 0 ≤ ε ≤ 6</text>')
    return "\n".join(svg)


def document(body, width, height):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
            f'<rect width="100%" height="100%" fill="white"/><g font-family="sans-serif">{body}</g></svg>')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, required=True, help="JSON mapping panel labels to run paths, relative to mapping file")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    runs = {k: args.runs.parent / v for k, v in json.loads(args.runs.read_text()).items()}
    args.out.mkdir(parents=True, exist_ok=True)
    combined = []
    for i, label in enumerate(ORDER):
        body = panel_svg(label, runs[label])
        (args.out / f"figure7{label}.svg").write_text(document(body, 600, 455))
        combined.append(f'<g transform="translate({600*(i%2)},{455*(i//2)})">{body}</g>')
    (args.out / "figure7.svg").write_text(document("".join(combined), 1200, 1820))
    print(f"Saved eight panels and figure7.svg in {args.out}")


if __name__ == "__main__":
    main()
