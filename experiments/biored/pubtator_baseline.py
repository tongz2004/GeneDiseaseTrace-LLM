"""PubTator3 on a benchmark sample, written in our output layout so the benchmark's scorer uses the same rules.

Usage: python experiments/biored/pubtator_baseline.py [dataset]     (default biored50)
  biored50  BioRED 50-abstract test sample  -> runs/biored50/pubtator3;  then evaluate.py pubtator3
  euadr     EU-ADR 100 abstracts            -> runs/euadr/pubtator3;     then ../euadr/evaluate.py pubtator3
  renet2    RENET2 40 full-text articles    -> runs/renet2/pubtator3;    then ../renet2/evaluate.py pubtator3
  biotriplex BioTriplex 30 full-text articles -> runs/biotriplex/pubtator3; then ../biotriplex/evaluate.py pubtator3
RENET2 and BioTriplex use the full-text export (full=true) so PubTator3's mentions in the main text can align too; its relation
predictions are the same either way (they come from the abstract, checked on BioRED).

For each sampled PMID the public PubTator3 abstract export (biocjson without full=true, cached under
data/pubtator_abstract/; title and abstract only, the same text BioRED annotates and our method reads) is aligned to the same JATS
blocks our method reads (biomed.pubtator: gene/disease mentions by surface form in the source passage). A gene-disease
pair is a positive prediction when PubTator3 lists a relation for it; everything else is no prediction.
PubTator3 has no evidence, source (novel/background) or uncertain state, so those fields stay empty.
Writes runs/biored50/pubtator3/<paper>/outputs/{entities,mentions,pairs}.csv and source.json.
"""
import csv
import json
import sys
import time
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT))
from biomed.jats import parse_article  # noqa: E402
from biomed.io import request_json  # noqa: E402
from biomed.pubtator import API, extract_entities, relation_predictions  # noqa: E402

RUNS = HERE.parent / "single_shot" / "runs"
# dataset: (sample file, folder of our JATS inputs, full-text export?)
DATASETS = {"biored50": (HERE / "sample_50.json", HERE / "papers", False),
            "biored_all": (HERE / "sample_all.json", HERE / "papers", False),
            "biored_testrest": (HERE / "sample_test_rest.json", HERE / "papers", False),
            "euadr": (HERE.parent / "euadr" / "sample_euadr.json", HERE.parent / "euadr" / "papers", False),
            "renet2": (HERE.parent / "renet2" / "sample_renet2_40.json", HERE.parent / "renet2" / "papers", True),
            "biotriplex": (HERE.parent / "biotriplex" / "sample_biotriplex.json", HERE.parent / "biotriplex" / "papers", True),
            "biotriplex2": (HERE.parent / "biotriplex" / "sample_biotriplex2.json", HERE.parent / "biotriplex" / "papers", True)}


def fetch_abstract(pmid, full=False):
    cache = HERE / "data" / ("pubtator_full" if full else "pubtator_abstract") / f"{pmid}.json"
    if not cache.exists():
        data = request_json(API + "?" + urllib.parse.urlencode({"pmids": pmid, **({"full": "true"} if full else {})}))
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(data), encoding="utf-8")
        time.sleep(0.55)
    data = json.loads(cache.read_text(encoding="utf-8"))
    documents = data.get("PubTator3", data.get("documents", [])) if isinstance(data, dict) else data
    return next(d for d in documents if str(d.get("pmid", d.get("id"))) == str(pmid))


def write(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(dataset="biored50"):
    sample_file, papers_dir, full = DATASETS[dataset]
    sample = json.loads(sample_file.read_text(encoding="utf-8"))
    for paper in sample["papers"]:
        pmid = paper.get("pmid") or paper["paper_id"].removeprefix("EUADR")
        article = parse_article(papers_dir / f"{paper['paper_id']}.xml")
        document = fetch_abstract(pmid, full)
        entities, rejected = extract_entities(article, document)
        predictions = relation_predictions(document, entities)
        out = RUNS / dataset / "pubtator3" / paper["paper_id"] / "outputs"
        out.mkdir(parents=True, exist_ok=True)
        write(out / "entities.csv", [{"entity_id": e["entity_id"], "entity_type": e["entity_type"],
                                      "canonical_name": e["canonical_name"], "external_id": e["external_id"]} for e in entities],
              ["entity_id", "entity_type", "canonical_name", "external_id"])
        write(out / "mentions.csv", [{"entity_id": e["entity_id"], "block_id": m["block_id"], "start": m["start"], "end": m["end"],
                                      "surface": m["text"]} for e in entities for m in e["mentions"]],
              ["entity_id", "block_id", "start", "end", "surface"])
        names = {e["entity_id"]: e["canonical_name"] for e in entities}
        write(out / "pairs.csv", [{"pair_id": f"{g}-{d}", "gene_id": g, "disease_id": d, "gene_name": names[g],
                                   "disease_name": names[d], "status": "supported",
                                   "relation_type": ", ".join(sorted({p["relation_type"] for p in preds})),
                                   "source_scope": "", "reason": f"PubTator3 prediction (score {max(float(p['source_score'] or 0) for p in preds):.2f})"}
                                  for (g, d), preds in sorted(predictions.items())],
              ["pair_id", "gene_id", "disease_id", "gene_name", "disease_name", "status", "relation_type", "source_scope", "reason"])
        (out / "source.json").write_text(json.dumps({"source": API, "pmid": pmid, "full_text_export": full, "written": datetime.now(timezone.utc).isoformat(),
                                                     "unaligned_annotations": len(rejected)}, indent=2), encoding="utf-8")
        print(f"{paper['paper_id']}: {sum(e['entity_type'] == 'gene' for e in entities)} genes, "
              f"{sum(e['entity_type'] == 'disease' for e in entities)} diseases, {len(predictions)} relations, "
              f"{len(rejected)} unaligned")


if __name__ == "__main__":
    main(*sys.argv[1:2])
