"""Blind, mixed check of pairs whose disease is not a BioTriplex disease: is it a disease at all, and is a relation stated?

Usage: python experiments/biotriplex/scope_check.py            (build the sheet)
       python experiments/biotriplex/scope_check.py --score    (after the sheet is filled in)

The scorer's "disease not a gold disease concept" false positives: supported pairs whose gene overlaps a BioTriplex
gene but whose disease overlaps no BioTriplex disease. Pools from V6 (Claude) and BioREx full text on the articles
both have a valid result for; one item per (article, gene concept, disease mention spans), merged across systems when
the two disease entities share a mention. 15 drawn per system (random.Random(20261106)), plus 5 BioTriplex relations
and 5 co-mentioned unrelated BioTriplex pairs not used in precision_check/ or precision_compare/. Shuffled, numbered
S01 ... Gene mentions as in precision_check.shown_spans; the disease mentions are those of the reporting system's
entity plus any entity of the other system sharing a mention, so both sides are marked alike. The annotator answers
whether the marked texts are a gene and a disease (a disorder or clinical phenotype) and what the paragraphs state.
Writes scope_check/sheet.html (export biotriplex_scope_answers.csv) and key.csv; --score writes result.json.
Exploratory: one annotator.
"""
import argparse
import collections
import csv
import json
import random
from pathlib import Path

import precision_check as pc
import precision_compare as cmp

HERE = Path(__file__).resolve().parent
OUT = HERE / "scope_check"
SEED, PER_SYSTEM, CONTROLS = 20261106, 15, 5
EXPORT = "biotriplex_scope_answers.csv"


def out_of_scope(paper, spans, run):
    """[(gene concept, disease mention spans, disease name)] for supported pairs with a gold gene, no gold disease."""
    ours = pc.v6_mentions(paper, run)
    names = {e["entity_id"]: e.get("canonical_name") or e.get("name") or "" for e in pc.read(run / paper / "outputs" / "entities.csv")}
    items = []
    for p in pc.read(run / paper / "outputs" / "pairs.csv"):
        if p["status"] != "supported":
            continue
        genes = sorted(c for c in spans if c.startswith("G") and pc.overlap(spans[c], ours[p["gene_id"]]))
        diseases = [c for c in spans if c.startswith("D") and pc.overlap(spans[c], ours[p["disease_id"]])]
        if genes and not diseases and ours[p["disease_id"]]:
            items.append((genes[0], tuple(sorted(ours[p["disease_id"]])), p.get("disease_name") or names.get(p["disease_id"], "")))
    return items


def build():
    gold_all = json.loads((HERE / "gold_biotriplex.json").read_text(encoding="utf-8"))
    papers = json.loads((HERE / "sample_biotriplex.json").read_text(encoding="utf-8"))["evaluation"]
    used = set()
    for key_file in (pc.OUT / "key.csv", cmp.OUT / "key.csv"):
        used |= {(r["paper_id"], r["gene"], r["disease"]) for r in csv.DictReader(key_file.open(encoding="utf-8-sig"))}
    pools, positives, negatives = {s: [] for s in cmp.SYSTEMS}, [], []
    for paper in papers:
        if not all((run / paper / "outputs" / ("validation.json" if s == "v6" else "source.json")).exists()
                   for s, run in cmp.SYSTEMS.items()):
            continue
        gold = gold_all[paper]
        spans = pc.concept_spans(gold)
        related = {(r["gene"], r["disease"]) for r in gold["relations"]}
        covered = set()
        for s, run in cmp.SYSTEMS.items():
            covered |= pc.supported_concepts(paper, spans, related, run)[0]
            seen = set()
            for g, dspans, name in out_of_scope(paper, spans, run):
                if (g, dspans) not in seen:
                    seen.add((g, dspans))
                    pools[s].append((paper, g, dspans, name))
        positives += [(paper, g, d) for g, d in sorted(related) if (paper, g, d) not in used]
        blocks_of = lambda c: {b for b, _, _ in spans[c]}
        negatives += [(paper, g, d) for g in sorted(c for c in spans if c.startswith("G"))
                      for d in sorted(c for c in spans if c.startswith("D"))
                      if (g, d) not in related and (g, d) not in covered and blocks_of(g) & blocks_of(d)
                      and (paper, g, d) not in used]
    rng = random.Random(SEED)
    drawn = []  # [(paper, gene concept, disease spans, name, systems)]
    for s in cmp.SYSTEMS:
        for paper, g, dspans, name in rng.sample(pools[s], PER_SYSTEM):
            match = next((x for x in drawn if x[0] == paper and x[1] == g and pc.overlap(list(x[2]), list(dspans))), None)
            if match:
                match[4].add(s)
            else:
                drawn.append([paper, g, dspans, name, {s}])
    items = [("+".join(sorted(x[4])), x[0], x[1], x[2], x[3]) for x in drawn]
    items += [("positive", p, g, d, None) for p, g, d in rng.sample(positives, CONTROLS)]
    items += [("negative", p, g, d, None) for p, g, d in rng.sample(negatives, CONTROLS)]
    rng.shuffle(items)
    runs = tuple(cmp.SYSTEMS.values())
    rows, key = [], []
    for n, (kind, paper, g, d, name) in enumerate(items, 1):
        gold, spans = gold_all[paper], pc.concept_spans(gold_all[paper])
        article = pc.parse_article(HERE / "papers" / f"{paper}.xml")
        gs = pc.shown_spans(paper, spans, g, runs)
        if name is None:  # control: a BioTriplex disease concept
            ds, dname = pc.shown_spans(paper, spans, d, runs), " / ".join(gold["entity_texts"][d][:3])
        else:  # the reporting system's disease mentions plus any other system's entity sharing one
            ds = sorted(set(d) | {s for run in runs for o in pc.v6_mentions(paper, run).values()
                                  if pc.overlap(o, list(d)) for s in o})
            dname = name
        gb, db = {b for b, _, _ in gs}, {b for b, _, _ in ds}
        both = [b for b in article["blocks"] if b["block_id"] in gb & db][:pc.MAX_PARAGRAPHS]
        shown = both or ([b for b in article["blocks"] if b["block_id"] in gb][:3]
                         + [b for b in article["blocks"] if b["block_id"] in db][:3])
        body = "".join(f'<p style="color:#666;font-size:.85em;margin-bottom:0">{b["section"] or "Abstract"}</p>'
                       + pc.marked_text([b], gs, ds) for b in sorted(shown, key=lambda b: b["block_id"]))
        if name is not None:
            # Systems name entities in their own style (PubTator3 uses MeSH headings), which would reveal the source;
            # show the first marked surface text instead, as for the controls.
            text = {b["block_id"]: b["text"] for b in shown}
            first = next(((blk, s, e) for blk, s, e in ds if blk in text), None)
            dname = text[first[0]][first[1]:first[2]] if first else dname
        item_id = f"S{n:02d}"
        rows.append({"item_id": item_id, "gene": " / ".join(gold["entity_texts"][g][:3]), "disease": dname,
                     "abstract_html": body})
        key.append({"item_id": item_id, "kind": kind, "paper_id": paper, "gene": g,
                    "disease": d if name is None else json.dumps(d), "paragraphs": "both" if both else "separate",
                    **{f"pool_{s}": len(pools[s]) for s in cmp.SYSTEMS}})
    OUT.mkdir(exist_ok=True)
    pc.write_form(OUT / "sheet.html", rows, scope="这些段落", store="biotriplex_scope_v1", download=EXPORT)
    with (OUT / "key.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(key[0]))
        writer.writeheader()
        writer.writerows(key)
    print(f"{len(rows)} items -> {OUT}; pools {{{', '.join(f'{s}: {len(p)}' for s, p in pools.items())}}}; "
          f"{collections.Counter(k['kind'] for k in key)}")


def score(path):
    if not path:
        path = next(p for p in (OUT / EXPORT, Path.home() / "Downloads" / EXPORT) if p.exists())
    print(f"answers: {path}")
    key = {r["item_id"]: r for r in csv.DictReader((OUT / "key.csv").open(encoding="utf-8-sig"))}
    report = collections.defaultdict(lambda: {"items": 0, "entities": collections.Counter(), "states": collections.Counter(),
                                              "disease_ok_and_relation": 0})
    for a in csv.DictReader(Path(path).open(encoding="utf-8-sig")):
        state = a[pc.ROUND1[1]].split(":")[0].strip().lower()
        entity = a[pc.ROUND1[0]].strip().lower()
        if state not in pc.STATES:
            print(f"{a['item_id']}: unreadable answer {a[pc.ROUND1[1]]!r}")
            continue
        kind = key[a["item_id"]]["kind"]
        for group in (kind.split("+") if kind not in ("positive", "negative") else [kind]):
            e = report[group]
            e["items"] += 1
            e["entities"][entity or "(blank)"] += 1
            e["states"][state] += 1
            e["disease_ok_and_relation"] += entity in ("yes", "gene wrong") and state in ("relation", "hedged relation")
    for e in report.values():
        e["entities"], e["states"] = dict(e["entities"]), dict(e["states"])
    (OUT / "result.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--score", nargs="?", const="", metavar="ANSWERS_CSV")
    args = parser.parse_args()
    score(args.score) if args.score is not None else build()
