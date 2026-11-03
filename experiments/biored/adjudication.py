"""Blinded two-round adjudication package for the disagreements between two variants on BioRED Dev.

Usage: python experiments/biored/adjudication.py [variant_a] [variant_b] [--split Dev] [--seed 20261103]
  defaults: v1_card_lean_medium v6_minimal_lean_medium

Units are gene-disease pairs in abstracts both variants finished. A pair of A and a pair of B are the same item when
their gene entities share an overlapping mention and so do their disease entities. Four strata, sampled at random:
  a_only      supported by A only                          (20)
  b_only      supported by B only                          (40)
  both        supported by both                            (20)
  gold_missed BioRED gold relation supported by neither    (20)
Round 1 (relation): annotators see only the abstract with the two entities marked and judge what the paper states.
Round 2 (evidence): annotators see anonymous evidence packets (the sentences one system cited for an item, random
IDs; one packet per distinct citation set) and judge whether the packet suffices for the round-1 claim.
Neither round shows system names, statuses or gold labels. key_coordinator.csv maps items and packets back to
systems, gold and strata, with each stratum's population and sampling probability; keep it away from annotators.
This is a development-set diagnostic of the disagreement, not an estimate of overall precision or recall.
Writes experiments/biored/adjudication/<split>_<a>_vs_<b>/.
"""
import argparse
import csv
import html
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT))
from biomed.jats import parse_article  # noqa: E402
from evaluate import RUNS, overlap, read  # noqa: E402

SIZES = {"a_only": 20, "b_only": 40, "both": 20, "gold_missed": 20}
ROUND1 = ["entities_correct (yes / gene wrong / disease wrong / both wrong)",
          "paper_states (relation / hedged relation / explicitly no relation / no statement)",
          "source (this study / background / both / unclear)", "notes"]
ROUND2 = ["sufficient (sufficient / insufficient / contradicts)", "notes"]


def system_pairs(out, statuses=("supported",)):
    """Pairs with these statuses, with the mention spans of both ends and the cited sentences."""
    spans = {}
    for m in read(out / "mentions.csv"):
        spans.setdefault(m["entity_id"], []).append((m["block_id"], int(m["start"]), int(m["end"])))
    evidence = {e["evidence_id"]: e for e in read(out / "evidence.csv")}
    pairs = []
    for p in read(out / "pairs.csv"):
        if p["status"] in statuses:
            cited = [evidence[v] for v in json.loads(p["evidence_ids"] or "[]") if v in evidence]
            pairs.append({"pair_id": p["pair_id"], "status": p["status"], "gene": p["gene_name"], "disease": p["disease_name"],
                          "gene_spans": spans.get(p["gene_id"], []), "disease_spans": spans.get(p["disease_id"], []),
                          "sentences": sorted({(e["block_id"], int(e["start"]), int(e["end"])) for e in cited}),
                          "evidence": cited})
    return pairs


def same(p, q):
    return overlap(p["gene_spans"], q["gene_spans"]) and overlap(p["disease_spans"], q["disease_spans"])


def marked_text(blocks, gene_spans, disease_spans):
    """Abstract as HTML with gene mentions in <b class=g> and disease mentions in <b class=d>."""
    parts = []
    for b in blocks:
        text, marks = b["text"], sorted([(s, e, "g") for bb, s, e in gene_spans if bb == b["block_id"]] +
                                        [(s, e, "d") for bb, s, e in disease_spans if bb == b["block_id"]])
        out, pos = [], 0
        for s, e, cls in marks:
            if s < pos:
                continue  # overlapping marks: keep the first
            out += [html.escape(text[pos:s]), f'<b class="{cls}">{html.escape(text[s:e])}</b>']
            pos = e
        parts.append("<p>" + "".join(out) + html.escape(text[pos:]) + "</p>")
    return "".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("a", nargs="?", default="v1_card_lean_medium")
    parser.add_argument("b", nargs="?", default="v6_minimal_lean_medium")
    parser.add_argument("--split", default="Dev")
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    gold_all = json.loads((HERE / "gold_all.json").read_text(encoding="utf-8"))
    split = {p["paper_id"]: p["split"] for p in json.loads((HERE / "sample_all.json").read_text(encoding="utf-8"))["papers"]}
    base = RUNS / "biored_all"
    done = lambda v, p: (base / v / p / "outputs" / "validation.json").exists()
    papers = sorted(p for p in gold_all if split[p] == args.split and done(args.a, p) and done(args.b, p))
    rows_of = {v: [r for r in read(base / v / "biored_eval_rows.csv") if r["paper_id"] in papers] for v in (args.a, args.b)}
    tp_pairs = {v: {(r["paper_id"], r["our_pair"]) for r in rows_of[v] if r["outcome"] == "TP"} for v in rows_of}

    population = {k: [] for k in SIZES}
    for paper in papers:
        pa, pb = system_pairs(base / args.a / paper / "outputs"), system_pairs(base / args.b / paper / "outputs")
        matched_b = set()
        for p in pa:
            q = next((q for q in pb if same(p, q)), None)
            if q:
                matched_b.add(q["pair_id"])
            population["both" if q else "a_only"].append({"paper": paper, args.a: p, args.b: q})
        population["b_only"] += [{"paper": paper, args.a: None, args.b: q} for q in pb if q["pair_id"] not in matched_b]
        # Gold relations neither variant recovered: both scorers list them as FN.
        missed = {(r["gold_gene"], r["gold_disease"]) for r in rows_of[args.a] if r["paper_id"] == paper and r["outcome"] == "FN"}
        missed &= {(r["gold_gene"], r["gold_disease"]) for r in rows_of[args.b] if r["paper_id"] == paper and r["outcome"] == "FN"}
        mentions = gold_all[paper]["mentions"]
        for gene, disease in sorted(missed):
            concept = lambda label: label.rsplit("(", 1)[1].rstrip(")")
            spans = lambda c: [(m["block_id"], m["start"], m["end"]) for m in mentions if c in m["identifier"].split(",")]
            population["gold_missed"].append({"paper": paper, args.a: None, args.b: None, "gold": {
                "gene": gene.rsplit(" (", 1)[0], "disease": disease.rsplit(" (", 1)[0],
                "gene_spans": spans(concept(gene)), "disease_spans": spans(concept(disease))}})

    rng = random.Random(args.seed)
    items = []
    for stratum, size in SIZES.items():
        picked = rng.sample(population[stratum], min(size, len(population[stratum])))
        for item in picked:
            item.update(stratum=stratum, probability=round(len(picked) / len(population[stratum]), 4))
        items += picked
    rng.shuffle(items)
    item_ids = rng.sample(range(100, 1000), len(items))
    out_dir = HERE / "adjudication" / f"{args.split}_{args.a}_vs_{args.b}"
    out_dir.mkdir(parents=True, exist_ok=True)
    blocks_of = {p: parse_article(HERE / "papers" / f"{p}.xml")["blocks"] for p in {i["paper"] for i in items}}

    round1, packets, key = [], {}, []
    for item, number in zip(items, item_ids):
        item_id = f"R{number}"
        shown = [x for x in (item[args.a], item[args.b], item.get("gold")) if x]
        gene_spans = sorted({s for x in shown for s in x["gene_spans"]})
        disease_spans = sorted({s for x in shown for s in x["disease_spans"]})
        blocks = blocks_of[item["paper"]]
        round1.append({"item_id": item_id, "gene": shown[0]["gene"], "disease": shown[0]["disease"],
                       "abstract_html": marked_text(blocks, gene_spans, disease_spans),
                       "abstract": marked_text(blocks, gene_spans, disease_spans)})
        packet_of = {}
        for variant in (args.a, args.b):
            pair = item[variant]
            if pair:
                packet_of.setdefault(tuple(pair["sentences"]), []).append(variant)
        packet_ids = {}
        for sentences, variants in packet_of.items():
            pid = f"E{rng.randrange(10000, 100000)}"
            text = {b["block_id"]: b["text"] for b in blocks}
            packets[pid] = {"packet_id": pid, "item_id": item_id, "gene": shown[0]["gene"], "disease": shown[0]["disease"],
                            "sentences": [text[b][s:e] for b, s, e in sentences]}
            for v in variants:
                packet_ids[v] = pid
        gold_hit = {v: bool(item[v]) and (item["paper"], item[v]["pair_id"]) in tp_pairs[v] for v in (args.a, args.b)}
        key.append({"item_id": item_id, "paper_id": item["paper"], "stratum": item["stratum"],
                    "stratum_population": len(population[item["stratum"]]), "sampling_probability": item["probability"],
                    f"{args.a}_pair": item[args.a]["pair_id"] if item[args.a] else "",
                    f"{args.b}_pair": item[args.b]["pair_id"] if item[args.b] else "",
                    f"{args.a}_packet": packet_ids.get(args.a, ""), f"{args.b}_packet": packet_ids.get(args.b, ""),
                    "matches_gold_relation": "yes" if item["stratum"] == "gold_missed" or any(gold_hit.values()) else "no"})

    write_round(out_dir / "round1_relation", round1, ["item_id", "gene", "disease"], ROUND1,
                "Round 1: what does the abstract state about the marked gene and disease?")
    packet_rows = sorted(packets.values(), key=lambda p: p["packet_id"])
    for p in packet_rows:
        p["evidence_html"] = "".join(f"<p>{html.escape(s)}</p>" for s in p["sentences"])
        p["evidence"] = " | ".join(p.pop("sentences"))
    write_round(out_dir / "round2_evidence", packet_rows, ["packet_id", "item_id", "gene", "disease"], ROUND2,
                "Round 2: do these sentences alone suffice for the claim you recorded in round 1?")
    with (out_dir / "key_coordinator.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(key[0]))
        writer.writeheader()
        writer.writerows(sorted(key, key=lambda k: k["item_id"]))
    summary = {s: {"population": len(population[s]), "sampled": sum(i["stratum"] == s for i in items)} for s in SIZES}
    (out_dir / "sampling.json").write_text(json.dumps({"a": args.a, "b": args.b, "split": args.split, "papers": len(papers),
                                                        "seed": args.seed, "strata": summary, "packets": len(packets)},
                                                       indent=2), encoding="utf-8")
    print(f"{len(items)} items, {len(packets)} evidence packets from {len(papers)} {args.split} abstracts -> {out_dir}")
    print(json.dumps(summary))


def write_round(stem, rows, id_fields, answer_fields, title):
    """A CSV answer sheet (text plus empty answer columns) and an HTML page to read the items in."""
    text_field = "abstract" if "abstract" in rows[0] else "evidence"
    with stem.with_suffix(".csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=id_fields + [text_field] + answer_fields)
        writer.writeheader()
        for r in sorted(rows, key=lambda r: r[id_fields[0]]):
            plain = r[text_field]
            if text_field == "abstract":
                plain = html.unescape(plain.replace('<b class="g">', "[[").replace('<b class="d">', "{{")
                                      .replace("</b>", "]]").replace("</p><p>", " ").replace("<p>", "").replace("</p>", ""))
                plain = fix_closers(plain)
            writer.writerow({**{k: r[k] for k in id_fields}, text_field: plain, **{a: "" for a in answer_fields}})
    cards = "".join(f'<section><h2>{html.escape(r[id_fields[0]])}: {html.escape(r["gene"])} &times; {html.escape(r["disease"])}</h2>'
                    f'{r[text_field + "_html"]}</section>' for r in sorted(rows, key=lambda r: r[id_fields[0]]))
    stem.with_suffix(".html").write_text(
        f'<!doctype html><meta charset="utf-8"><title>{html.escape(title)}</title><style>body{{font:15px/1.5 system-ui;'
        f'max-width:860px;margin:2em auto;padding:0 16px}}section{{border-top:1px solid #ccc;padding:.5em 0}}'
        f'b.g{{background:#cfe3ff}}b.d{{background:#ffd9cc}}</style><h1>{html.escape(title)}</h1>'
        f'<p>Gene mentions are blue, disease mentions orange. Record answers in the CSV with the same ID.</p>{cards}',
        encoding="utf-8")


def fix_closers(text):
    """Turn the generic ']]' closer after a '{{' opener into '}}'."""
    out, stack, i = [], [], 0
    while i < len(text):
        if text.startswith("[[", i) or text.startswith("{{", i):
            stack.append(text[i:i + 2])
            out.append(text[i:i + 2])
            i += 2
        elif text.startswith("]]", i) and stack:
            out.append("}}" if stack.pop() == "{{" else "]]")
            i += 2
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


if __name__ == "__main__":
    main()
