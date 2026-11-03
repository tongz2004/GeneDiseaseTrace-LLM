"""One annotator's blind check of V6's full-text pairs that BioTriplex links no relation for (precision is otherwise unknown).

Usage: python experiments/biotriplex/precision_check.py            (build the sheet)
       python experiments/biotriplex/precision_check.py --score    (after the sheet is filled in)

Pool: V6 (run set biotriplex, full text) supported pairs whose gene and disease both overlap BioTriplex entities but
whose pair BioTriplex links no relation for (the "gold has both concepts but no relation" false positives). Items:
- 20 drawn from that pool;
- 5 positive controls: BioTriplex relations;
- 5 negative controls: BioTriplex gene-disease pairs mentioned in a common paragraph, with no BioTriplex relation and
  no V6 supported pair.
All with random.Random(20261104), shuffled, renumbered P01 ... Each item shows, in document order, up to 6 paragraphs
where both entities are mentioned (or, when none has both, up to 3 with the gene and 3 with the disease), with the
mentions marked. An entity's mentions are BioTriplex's plus V6's mentions of any entity overlapping it, for every item
alike (BioTriplex alone often never puts the two in one paragraph). No system output, evidence or gold label is shown.
Writes precision_check/sheet.html (answers exported as biotriplex_precision_answers.csv) and key.csv; --score reads
the export from precision_check/ or Downloads and writes result.json. Exploratory: one annotator.
"""
import argparse
import collections
import csv
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "biored"))
from adjudication import marked_text, parse_article  # noqa: E402
from evaluate import overlap, read  # noqa: E402
from self_check import ROUND1, STATES, write_form  # noqa: E402

OUT = HERE / "precision_check"
RUN = HERE.parent / "single_shot" / "runs" / "biotriplex" / "v6_minimal_lean_medium"
SEED, POOL, CONTROLS, MAX_PARAGRAPHS = 20261104, 20, 5, 6
EXPORT = "biotriplex_precision_answers.csv"


def concept_spans(gold):
    spans = collections.defaultdict(list)
    for m in gold["mentions"]:
        spans[m["identifier"]].append((m["block_id"], m["start"], m["end"]))
    return spans


def v6_mentions(paper, run=RUN):
    """{entity_id: mention spans} of a run's output for one article (V6 by default)."""
    ours = collections.defaultdict(list)
    for m in read(run / paper / "outputs" / "mentions.csv"):
        ours[m["entity_id"]].append((m["block_id"], int(m["start"]), int(m["end"])))
    return ours


def supported_concepts(paper, spans, related, run=RUN):
    """(covered, unannotated): concept pairs any supported pair of the run overlaps, and one concept pair per
    supported pair that overlaps BioTriplex entities at both ends but no BioTriplex relation (the scorer's
    "gold has both concepts but no relation")."""
    ours = v6_mentions(paper, run)
    covered, unannotated = set(), set()
    for p in read(run / paper / "outputs" / "pairs.csv"):
        if p["status"] != "supported":
            continue
        genes = sorted(c for c in spans if c.startswith("G") and overlap(spans[c], ours[p["gene_id"]]))
        diseases = sorted(c for c in spans if c.startswith("D") and overlap(spans[c], ours[p["disease_id"]]))
        combos = [(g, d) for g in genes for d in diseases]
        covered |= set(combos)
        if combos and not set(combos) & related:
            unannotated.add(combos[0])
    return covered, unannotated


def shown_spans(paper, spans, concept, runs=(RUN,)):
    """BioTriplex mentions of a concept plus the runs' mentions of any entity overlapping it."""
    extra = [s for run in runs for spans_o in v6_mentions(paper, run).values()
             if overlap(spans_o, spans[concept]) for s in spans_o]
    return sorted(set(spans[concept]) | set(extra))


def build():
    gold_all = json.loads((HERE / "gold_biotriplex.json").read_text(encoding="utf-8"))
    papers = json.loads((HERE / "sample_biotriplex.json").read_text(encoding="utf-8"))["evaluation"]
    pool, positives, negatives = [], [], []
    for paper in papers:
        if not (RUN / paper / "outputs" / "validation.json").exists():
            continue
        gold = gold_all[paper]
        spans = concept_spans(gold)
        related = {(r["gene"], r["disease"]) for r in gold["relations"]}
        covered, unannotated = supported_concepts(paper, spans, related)
        pool += [(paper, g, d) for g, d in sorted(unannotated)]
        positives += [(paper, g, d) for g, d in sorted(related)]
        blocks_of = lambda c: {b for b, _, _ in spans[c]}
        negatives += [(paper, g, d) for g in sorted(c for c in spans if c.startswith("G"))
                      for d in sorted(c for c in spans if c.startswith("D"))
                      if (g, d) not in related and (g, d) not in covered and blocks_of(g) & blocks_of(d)]
    rng = random.Random(SEED)
    items = ([("unannotated", x) for x in rng.sample(pool, POOL)] + [("positive", x) for x in rng.sample(positives, CONTROLS)]
             + [("negative", x) for x in rng.sample(negatives, CONTROLS)])
    rng.shuffle(items)
    rows, key = [], []
    for n, (kind, (paper, g, d)) in enumerate(items, 1):
        gold, spans = gold_all[paper], concept_spans(gold_all[paper])
        article = parse_article(HERE / "papers" / f"{paper}.xml")
        gs, ds = shown_spans(paper, spans, g), shown_spans(paper, spans, d)
        gb, db = {b for b, _, _ in gs}, {b for b, _, _ in ds}
        both = [b for b in article["blocks"] if b["block_id"] in gb & db][:MAX_PARAGRAPHS]
        shown = both or ([b for b in article["blocks"] if b["block_id"] in gb][:3]
                         + [b for b in article["blocks"] if b["block_id"] in db][:3])
        shown = sorted(shown, key=lambda b: b["block_id"])
        body = "".join(f'<p style="color:#666;font-size:.85em;margin-bottom:0">{b["section"] or "Abstract"}</p>'
                       + marked_text([b], gs, ds) for b in shown)
        texts = gold["entity_texts"]
        item_id = f"P{n:02d}"
        rows.append({"item_id": item_id, "gene": " / ".join(texts[g][:3]), "disease": " / ".join(texts[d][:3]),
                     "abstract_html": body})
        key.append({"item_id": item_id, "kind": kind, "paper_id": paper, "gene": g, "disease": d,
                    "paragraphs": "both" if both else "separate", "pool_size": len(pool) if kind == "unannotated" else ""})
    OUT.mkdir(exist_ok=True)
    write_form(OUT / "sheet.html", rows, scope="这些段落", store="biotriplex_precision_v1", download=EXPORT)
    with (OUT / "key.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(key[0]))
        writer.writeheader()
        writer.writerows(key)
    print(f"{len(rows)} items -> {OUT}; pool {len(pool)} unannotated pairs, {len(positives)} relations, "
          f"{len(negatives)} co-mentioned unrelated pairs; {collections.Counter(k['kind'] for k in key)}")


def score(path):
    if not path:
        path = next(p for p in (OUT / EXPORT, Path.home() / "Downloads" / EXPORT) if p.exists())
    print(f"answers: {path}")
    key = {r["item_id"]: r for r in csv.DictReader((OUT / "key.csv").open(encoding="utf-8-sig"))}
    report = {}
    for a in csv.DictReader(Path(path).open(encoding="utf-8-sig")):
        state = a[ROUND1[1]].split(":")[0].strip().lower()
        if state not in STATES:
            print(f"{a['item_id']}: unreadable answer {a[ROUND1[1]]!r}")
            continue
        e = report.setdefault(key[a["item_id"]]["kind"], {"items": 0, "states": collections.Counter(), "entities_wrong": 0})
        e["items"] += 1
        e["states"][state] += 1
        e["entities_wrong"] += a[ROUND1[0]].strip().lower() not in ("", "yes")
    for e in report.values():
        e["relation_or_hedged"] = e["states"]["relation"] + e["states"]["hedged relation"]
        e["states"] = dict(e["states"])
    (OUT / "result.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--score", nargs="?", const="", metavar="ANSWERS_CSV")
    args = parser.parse_args()
    score(args.score) if args.score is not None else build()
