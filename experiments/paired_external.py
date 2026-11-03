"""Paired comparisons and deployment figures on EU-ADR (PA gold) and RENET2 (associated gold), one-to-one scorer.

Usage: python experiments/paired_external.py

For a - b in (minimal - PubTator3, minimal - card, card - PubTator3): P/R/F1 differences on the documents both have,
2,000 bootstrap draws over documents (seed 20261028). Deployment: the minimal prompt with every failed document
counted as an empty answer. Writes runs/paired_external.json.
"""
import importlib.util
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "single_shot" / "runs"
_spec = importlib.util.spec_from_file_location("biored_evaluate", HERE / "biored" / "evaluate.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)
SYSTEMS = {"minimal": "v6_minimal_lean_medium", "card": "v1_card_lean_medium", "pubtator3": "pubtator3"}
DATA = {"euadr": (HERE / "euadr" / "gold_euadr.json", {"PA"}), "renet2": (HERE / "renet2" / "gold_renet2.json", {"associated"})}


def counts(run_set, variant, gold):
    out = {}
    for p, g in gold.items():
        o = RUNS / run_set / variant / p / "outputs"
        if (o / "validation.json").exists() or (o / "source.json").exists():
            out[p] = B.score_paper(p, g, o)[0]
    return out


def main():
    report = {}
    for name, (path, labels) in DATA.items():
        raw = json.loads(path.read_text(encoding="utf-8"))
        gold = {p: {**g, "relations": [r for r in g["relations"] if r["type"] in labels]} for p, g in raw.items()}
        c = {s: counts(name, v, gold) for s, v in SYSTEMS.items()}
        entry = {"documents": len(gold), "valid": {s: len(c[s]) for s in c}}
        for a, b in (("minimal", "pubtator3"), ("minimal", "card"), ("card", "pubtator3")):
            common = sorted(set(c[a]) & set(c[b]))
            point = [x - y for x, y in zip(B.prf(B.total([c[a][p] for p in common])), B.prf(B.total([c[b][p] for p in common])))]
            rng, draws = random.Random(20261028), []
            for _ in range(2000):
                d = [rng.choice(common) for _ in common]
                draws.append([x - y for x, y in zip(B.prf(B.total([c[a][p] for p in d])), B.prf(B.total([c[b][p] for p in d])))])
            ci = lambda i: [round(sorted(x[i] for x in draws)[49], 3), round(sorted(x[i] for x in draws)[1949], 3)]
            entry[f"{a} - {b}"] = {"documents": len(common), **{k: {"diff": round(point[i], 3), "ci95": ci(i)}
                                                               for i, k in enumerate(("precision", "recall", "f1"))}}
        t = B.total(list(c["minimal"].values()))
        missed = sum(len(gold[p]["relations"]) for p in gold if p not in c["minimal"])
        deployed = dict(t, fn=t["fn"] + missed)
        entry["minimal_deployment"] = dict(zip(("precision", "recall", "f1"), (round(x, 3) for x in B.prf(deployed))))
        report[name] = entry
    (RUNS / "paired_external.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
