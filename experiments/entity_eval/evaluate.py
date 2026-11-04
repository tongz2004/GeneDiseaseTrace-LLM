"""Task 1: how well are genes and diseases extracted? Stored outputs only, no model calls.

Usage: python experiments/entity_eval/evaluate.py

Rules (fixed before scoring):
- Corpora: BioRED Test (100 abstracts, exhaustive gene and disease annotation) gives precision, recall and F1.
  BioTriplex full text (locked 30 and 49 further articles) gives recall only: entities there are annotated selectively, so
  an unannotated mention is not an error.
- Unit, primary: the document-level concept. A gold concept is found when one of its mentions overlaps a mention of a
  predicted entity of the same class; a predicted entity is correct when one of its mentions overlaps a gold mention of
  that class. Identifiers are not compared. Secondary: the mention, with the same overlap rule.
- Classes: gene; disease. BioRED's class is "disease or phenotypic feature", so symptom, phenotype and clinical-event
  entities count as disease for the LLM runs. The disease_strict view counts only entities typed disease.
- Systems: V6, V1 (BioRED only) and PubTator3. Valid documents only; each table gives the document count.
- Intervals: 2,000 document-level bootstrap draws, seed 20261028. PubTator3's entity recognizer was trained on data
  including BioRED, so it is in-domain there.
Writes result.json.
"""
import collections
import csv
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
RUNS = EXPERIMENTS / "single_shot" / "runs"
GOLD_TYPE = {"GeneOrGeneProduct": "gene", "DiseaseOrPhenotypicFeature": "disease"}
PHENOTYPE = {"symptom", "phenotype", "clinical_event"}
SEED, DRAWS = 20261028, 2000


def read(path):
    if not path.exists():
        return []
    with path.open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load(name):
    return json.loads((EXPERIMENTS / name).read_text(encoding="utf-8"))


def hit(span, spans):
    return any(b == span[0] and s < span[2] and e > span[1] for b, s, e in spans)


def predicted(out, view):
    """{class: {entity_id: mention spans}} for one document."""
    kinds = {}
    for e in read(out / "entities.csv") + read(out / "auxiliary_entities.csv"):
        t = e["entity_type"]
        if t == "gene":
            kinds[e["entity_id"]] = "gene"
        elif t == "disease" or (t in PHENOTYPE and view == "disease"):
            kinds[e["entity_id"]] = "disease"
    found = collections.defaultdict(lambda: collections.defaultdict(list))
    for m in read(out / "mentions.csv"):
        if m["entity_id"] in kinds:
            found[kinds[m["entity_id"]]][m["entity_id"]].append((m["block_id"], int(m["start"]), int(m["end"])))
    return found


def gold_concepts(gold):
    concepts = collections.defaultdict(lambda: collections.defaultdict(list))
    for m in gold["mentions"]:
        kind = GOLD_TYPE.get(m["type"])
        if kind:
            for c in str(m["identifier"]).split(","):
                concepts[kind][c].append((m["block_id"], m["start"], m["end"]))
    return concepts


def counts(gold, out, view):
    """{class: [concepts found, gold concepts, entities correct, entities, mentions found, gold mentions, mentions correct, mentions]}"""
    ours, theirs, result = predicted(out, view), gold_concepts(gold), {}
    for kind in ("gene", "disease"):
        pred, ref = ours[kind], theirs[kind]
        pred_spans = [s for v in pred.values() for s in v]
        ref_spans = sorted({s for v in ref.values() for s in v})
        result[kind] = [sum(any(hit(s, pred_spans) for s in v) for v in ref.values()), len(ref),
                        sum(any(hit(s, ref_spans) for s in v) for v in pred.values()), len(pred),
                        sum(hit(s, pred_spans) for s in ref_spans), len(ref_spans),
                        sum(hit(s, ref_spans) for s in pred_spans), len(pred_spans)]
    return result


def prf(c):
    recall, precision = c[0] / max(1, c[1]), c[2] / max(1, c[3])
    return precision, recall, 2 * precision * recall / max(1e-9, precision + recall)


def summary(docs, kind, offset=0):
    """Point estimates and intervals; offset 4 selects the mention counts."""
    rows = [d[kind][offset:offset + 4] for d in docs.values()]
    total = lambda rs: [sum(x) for x in zip(*rs)]
    point = prf(total(rows))
    rng, draws = random.Random(SEED), []
    for _ in range(DRAWS):
        draws.append(prf(total(rng.choices(rows, k=len(rows)))))
    ci = lambda i: [round(sorted(d[i] for d in draws)[int(.025 * DRAWS)], 3), round(sorted(d[i] for d in draws)[int(.975 * DRAWS) - 1], 3)]
    return {"precision": round(point[0], 3), "recall": round(point[1], 3), "f1": round(point[2], 3),
            "precision_ci": ci(0), "recall_ci": ci(1), "f1_ci": ci(2), "gold": total(rows)[1]}


def corpus(papers, gold, view):
    docs = {}
    for paper, out in papers.items():
        if (out / "validation.json").exists() or (out / "source.json").exists():
            docs[paper] = counts(gold[paper], out, view)
    return docs


def main():
    split = {p["paper_id"]: p["split"] for p in load("biored/sample_all.json")["papers"]}
    gold = load("biored/gold_all.json")
    rest = load("biored/gold_test_rest.json")
    gold.update(rest)
    test = {**{p: "biored_all" for p, s in split.items() if s == "Test"}, **{p: "biored_testrest" for p in rest}}
    report = {"BioRED Test": {}}
    for system in ("v6_minimal_lean_medium", "v1_card_lean_medium", "pubtator3"):
        entry = report["BioRED Test"][system] = {}
        for view in ("disease", "disease_strict"):
            docs = corpus({p: RUNS / rs / system / p / "outputs" for p, rs in test.items()}, gold, view)
            entry["documents"] = len(docs)
            if view == "disease":
                entry["gene"] = summary(docs, "gene")
                entry["gene_mentions"] = summary(docs, "gene", 4)
                entry["disease"] = summary(docs, "disease")
                entry["disease_mentions"] = summary(docs, "disease", 4)
            else:
                entry["disease_strict"] = summary(docs, "disease")
    for name, run_set, gold_file, sample_file in (("BioTriplex locked", "biotriplex", "gold_biotriplex.json", "sample_biotriplex.json"),
                                                  ("BioTriplex further", "biotriplex2", "gold_biotriplex2.json", "sample_biotriplex2.json")):
        g, papers = load(f"biotriplex/{gold_file}"), load(f"biotriplex/{sample_file}")["evaluation"]
        report[name] = {}
        for system in ("v6_minimal_lean_medium", "pubtator3"):
            docs = corpus({p: RUNS / run_set / system / p / "outputs" for p in papers}, g, "disease")
            report[name][system] = {"documents": len(docs),
                                    **{k: {x: summary(docs, k)[x] for x in ("recall", "recall_ci", "gold")} for k in ("gene", "disease")}}
    (HERE / "result.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
