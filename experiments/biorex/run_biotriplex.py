"""BioREx (the supervised relation model inside PubTator3) applied paragraph by paragraph to BioTriplex full texts.

Usage: python experiments/biorex/run_biotriplex.py [--abstract]
       then python experiments/biotriplex/evaluate.py biorex_fulltext biotriplex   (or biorex_abstract ... --abstract)

Rules (fixed before the run):
- Text: the paragraphs V6 read (experiments/biotriplex/papers/BTX<pmcid>.xml, the same kept blocks, boilerplate
  sections dropped); --abstract keeps only the abstract blocks and writes the BTXA articles.
- Entities: PubTator3's full-text gene and disease annotations as aligned to those blocks by the PubTator3 baseline
  (runs/biotriplex/pubtator3/<paper>/outputs), so BioREx and the PubTator3 baseline see the same entities.
- Each paragraph is one BioREx document (BioREx was trained on abstracts; paragraphs are of similar length, and
  longer ones are truncated at 512 tokens as BioREx does). Only gene-disease instances are predicted.
- A gene-disease pair is supported when BioREx gives any label other than a "None" label in at least one paragraph.
BioREx runs in experiments/biorex/.venv (TensorFlow 2.18, CPU) through predict.py and BioREx's own converters.
Writes runs/biotriplex/biorex_{fulltext,abstract}/<paper>/outputs/{entities,mentions,pairs}.csv and source.json.
"""
import argparse
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "experiments" / "single_shot"))
from common import load  # noqa: E402

RUNS = PROJECT / "experiments" / "single_shot" / "runs" / "biotriplex"
BIOREX = HERE / "vendor" / "BioREx"
VENV = HERE / ".venv" / "Scripts" / "python.exe"
TYPE = {"gene": "GeneOrGeneProduct", "disease": "DiseaseOrPhenotypicFeature"}
NEGATIVE = {"None"} | {l for l in json.loads((BIOREX / "pretrained_model_biolinkbert" / "config.json").read_text())["id2label"].values()
                       if l.startswith("None")}


def read(path):
    with open(path, encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--abstract", action="store_true")
    args = parser.parse_args()
    sample = json.loads((PROJECT / "experiments" / "biotriplex" / "sample_biotriplex.json").read_text(encoding="utf-8"))
    name = "biorex_abstract" if args.abstract else "biorex_fulltext"
    work = HERE / "work" / name
    work.mkdir(parents=True, exist_ok=True)
    docs, doc_of = [], {}
    for paper in sample["evaluation"]:
        target = paper.replace("BTX", "BTXA") if args.abstract else paper
        article, kept, _, _ = load(paper)
        blocks = [b for b in kept if b["content_type"] != "section_heading"]
        if args.abstract:
            abstract_ids = set(json.loads((PROJECT / "experiments" / "biotriplex" / "gold_biotriplex.json")
                                          .read_text(encoding="utf-8"))[paper]["abstract_blocks"])
            blocks = [b for b in kept if b["block_id"] in abstract_ids]
        pt = RUNS / "pubtator3" / paper / "outputs"
        entity = {e["entity_id"]: e for e in read(pt / "entities.csv") if e["entity_type"] in TYPE}
        mentions = collections.defaultdict(list)
        for m in read(pt / "mentions.csv"):
            if m["entity_id"] in entity:
                mentions[m["block_id"]].append(m)
        for n, b in enumerate(blocks):
            ms = mentions.get(b["block_id"], [])
            if not any(entity[m["entity_id"]]["entity_type"] == "gene" for m in ms) or \
                    not any(entity[m["entity_id"]]["entity_type"] == "disease" for m in ms):
                continue
            doc_id = f"{paper.removeprefix('BTX')}{n:04d}"
            doc_of[doc_id] = (target, paper)
            text = b["text"].replace("|", " ").replace("\n", " ")
            lines = [f"{doc_id}|t|{text}", f"{doc_id}|a|"]
            for m in sorted(ms, key=lambda m: int(m["start"])):
                e = entity[m["entity_id"]]
                lines.append("\t".join([doc_id, m["start"], m["end"], text[int(m["start"]):int(m["end"])],
                                        TYPE[e["entity_type"]], m["entity_id"]]))
            docs.append("\n".join(lines) + "\n")
    (work / "input.pubtator").write_text("\n".join(docs) + "\n", encoding="utf-8", newline="\n")
    env = {"PYTHONIOENCODING": "utf-8", "TF_CPP_MIN_LOG_LEVEL": "2", "TF_USE_LEGACY_KERAS": "1", "CUDA_VISIBLE_DEVICES": ""}
    run = lambda cmd: subprocess.run([str(VENV)] + cmd, cwd=BIOREX, check=True, env={**__import__("os").environ, **env})
    run(["src/dataset_format_converter/convert_pubtator_2_tsv.py", "--exp_option", "biored_pred",
         "--in_pubtator_file", str(work / "input.pubtator"), "--out_tsv_file", str(work / "all.tsv")])
    # Keep gene-disease instances only (the converter pairs every two entities).
    rows = [l for l in (work / "all.tsv").read_text(encoding="utf-8").splitlines()
            if {l.split("\t")[1], l.split("\t")[2]} == set(TYPE.values())]
    (work / "processed.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(docs)} paragraphs, {len(rows)} gene-disease instances", flush=True)
    run([str(HERE / "predict.py"), str(work / "processed.tsv"), str(work / "pred")])
    labels = list(json.loads((BIOREX / "pretrained_model_biolinkbert" / "config.json").read_text())["id2label"].values())
    supported = collections.defaultdict(lambda: collections.defaultdict(list))
    for line, logits in zip(rows, (work / "pred" / "test_results.tsv").read_text().splitlines()):
        f = line.split("\t")
        values = [float(x) for x in logits.split("\t")]
        label = labels[max(range(len(values)), key=values.__getitem__)]
        if label not in NEGATIVE:
            target, paper = doc_of[f[0]]
            gene, disease = (f[3], f[4]) if f[1] == TYPE["gene"] else (f[4], f[3])
            supported[target][(gene, disease)].append(label)
    for paper in sample["evaluation"]:
        target = paper.replace("BTX", "BTXA") if args.abstract else paper
        pt = RUNS / "pubtator3" / paper / "outputs"
        out = RUNS / name / target / "outputs"
        out.mkdir(parents=True, exist_ok=True)
        entities = [e for e in read(pt / "entities.csv") if e["entity_type"] in TYPE]
        write(out / "entities.csv", entities, ["entity_id", "entity_type", "canonical_name", "external_id"])
        ids = {e["entity_id"] for e in entities}
        write(out / "mentions.csv", [m for m in read(pt / "mentions.csv") if m["entity_id"] in ids],
              ["entity_id", "block_id", "start", "end", "surface"])
        names = {e["entity_id"]: e["canonical_name"] for e in entities}
        write(out / "pairs.csv", [{"pair_id": f"{g}-{d}", "gene_id": g, "disease_id": d, "gene_name": names[g],
                                   "disease_name": names[d], "status": "supported",
                                   "relation_type": ", ".join(sorted(set(labs))), "source_scope": "",
                                   "reason": f"BioREx, {len(labs)} paragraph(s)"}
                                  for (g, d), labs in sorted(supported[target].items())],
              ["pair_id", "gene_id", "disease_id", "gene_name", "disease_name", "status", "relation_type", "source_scope", "reason"])
        (out / "source.json").write_text(json.dumps({"source": "BioREx pretrained_model_biolinkbert, per paragraph",
                                                     "input": "abstract" if args.abstract else "full text",
                                                     "entities": "PubTator3 full-text annotations aligned by pubtator_baseline.py"},
                                                    indent=2), encoding="utf-8")
    print(f"{name}: {sum(len(v) for v in supported.values())} supported pairs over {len(sample['evaluation'])} articles")


if __name__ == "__main__":
    main()
