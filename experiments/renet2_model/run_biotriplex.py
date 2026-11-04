"""RENET2's released full-text models on the locked BioTriplex set, as a supervised full-text comparator.

Usage: python experiments/renet2_model/run_biotriplex.py [--abstract]
       then python experiments/biotriplex/evaluate.py renet2_fulltext biotriplex   (renet2_abstract ... --abstract)

Rules (fixed before the run):
- Models: models/ft_models from renet2_data_models.tar.gz, all 10, default decision rule: each model's stored threshold,
  and a pair is positive when at least int(10 / 2) = 5 models vote for it (RENET2's own cutoff). Nothing is tuned.
- Text: the paragraphs V6 and BioREx read (section headings and boilerplate sections dropped), split with our sentence
  splitter (RENET2's geniass needs Ruby) and written in RENET2's format: title, abstract, body, one space between
  sentences. RENET2 reads at most 400 sentences per document; articles beyond that are truncated by the model.
- Entities: PubTator3 gene and disease annotations aligned by the PubTator3 baseline, the same as for BioREx. Entity IDs
  are replaced by indices (genes numeric, diseases D000001 ...), so RENET2 sees one ID per entity.
- --abstract keeps only the abstract paragraphs and writes the BTXA articles.
RENET2's last step (attaching entity names) fails under current pandas after the per-model files are written, so the
vote is applied here to gda_rst_XX.tsv, as RENET2 does. RENET2's models were trained on the 500 RENET2 articles, none of
which is in BioTriplex.
Writes runs/biotriplex/renet2_{fulltext,abstract}/<paper>/outputs/{entities,mentions,pairs}.csv and source.json.
"""
import argparse
import collections
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT / "experiments" / "single_shot"))
from common import load, split_sentences  # noqa: E402

RUNS = PROJECT / "experiments" / "single_shot" / "runs" / "biotriplex"
VENV = HERE / ".venv" / "Scripts" / "python.exe"
BIOTRIPLEX = PROJECT / "experiments" / "biotriplex"
MODELS, CUTOFF, TYPES = 10, 5, {"gene": "Gene", "disease": "Disease"}


def read(path):
    with open(path, encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write(path, rows, fields):
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build(paper, abstract_ids, only_abstract):
    """(title, abstract, body sentence lists, annotation lines, id map, sentence count) for one article."""
    _, kept, _, _ = load(paper)
    blocks = [b for b in kept if b["content_type"] != "section_heading"
              and (not only_abstract or b["block_id"] in abstract_ids)]
    pt = RUNS / "pubtator3" / paper / "outputs"
    entity = {e["entity_id"]: e for e in read(pt / "entities.csv") if e["entity_type"] in TYPES}
    ids, by_block = {}, collections.defaultdict(list)
    for n, e in enumerate(sorted(entity), 1):
        ids[e] = str(n) if entity[e]["entity_type"] == "gene" else f"D{n:06d}"
    for m in read(pt / "mentions.csv"):
        if m["entity_id"] in entity:
            by_block[m["block_id"]].append(m)
    parts, anns, pos, dropped = {"title": [], "abstract": [], "body": []}, [], 0, 0
    for i, b in enumerate(blocks):
        part = "title" if i == 0 else "abstract" if b["block_id"] in abstract_ids else "body"
        spans = []
        for s, e in split_sentences(b["text"]):
            segment = b["text"][s:e]
            lead = len(segment) - len(segment.lstrip())
            sentence = segment.strip()
            if sentence:
                parts[part].append(sentence)
                spans.append((s + lead, s + lead + len(sentence), pos))
                pos += len(sentence) + 1
        for m in sorted(by_block.get(b["block_id"], ()), key=lambda m: int(m["start"])):
            s, e = int(m["start"]), int(m["end"])
            hit = next((x for x in spans if x[0] <= s and e <= x[1]), None)
            if not hit:
                dropped += 1
                continue
            offset = hit[2] + s - hit[0]
            anns.append((offset, offset + e - s, b["text"][s:e], TYPES[entity[m["entity_id"]]["entity_type"]],
                         ids[m["entity_id"]], i + 1))
    return parts, anns, ids, dropped


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--abstract", action="store_true")
    args = parser.parse_args()
    sample = json.loads((BIOTRIPLEX / "sample_biotriplex.json").read_text(encoding="utf-8"))
    gold = json.loads((BIOTRIPLEX / "gold_biotriplex.json").read_text(encoding="utf-8"))
    name = "renet2_abstract" if args.abstract else "renet2_fulltext"
    work = HERE / "work" / name
    work.mkdir(parents=True, exist_ok=True)
    docs, sents, anns_out, id_of, stats = [], [], [], {}, []
    for paper in sample["evaluation"]:
        pid = paper.removeprefix("BTX")
        parts, anns, ids, dropped = build(paper, set(gold[paper]["abstract_blocks"]), args.abstract)
        all_sentences = parts["title"] + parts["abstract"] + parts["body"]
        docs += [pid, " ".join(parts["title"]), " ".join(parts["abstract"]), " ".join(parts["body"]), ""]
        sents += [pid] + all_sentences + [""]
        anns_out += ["\t".join([pid, str(a), str(b), t, ty, i, str(s)]) for a, b, t, ty, i, s in sorted(anns)] + [""]
        id_of[pid] = {v: k for k, v in ids.items()}
        stats.append((paper, len(all_sentences), len(anns), dropped))
    (work / "docs.txt").write_text("\n".join(docs) + "\n", encoding="utf-8", newline="\n")
    (work / "sentences.txt").write_text("\n".join(sents) + "\n", encoding="utf-8", newline="\n")
    (work / "anns_n.txt").write_text("\n".join(anns_out) + "\n", encoding="utf-8", newline="\n")
    (work / "anns.txt").write_text("\n".join(anns_out) + "\n", encoding="utf-8", newline="\n")
    (work / "labels.txt").write_text("pmid,geneId,diseaseId,label\n0,0,D0,0.0\n", encoding="utf-8", newline="\n")
    over = sum(s[1] > 400 for s in stats)
    print(f"{len(stats)} articles; {sum(s[1] for s in stats)} sentences; {over} over 400 sentences; "
          f"{sum(s[2] for s in stats)} mentions placed, {sum(s[3] for s in stats)} dropped (cross a sentence)", flush=True)
    env = {**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"}
    with (work / "log.txt").open("w", encoding="utf-8") as log:
        subprocess.run([str(VENV), str(HERE / "run_predict.py"), "--raw_data_dir", str(work) + "/", "--model_dir",
                        str(HERE / "weights" / "models" / "ft_models") + "/", "--gda_fn_d", str(work) + "/",
                        "--models_number", str(MODELS), "--batch_size", "4", "--no_cache_file", "--use_cuda"],
                       env=env, stdout=log, stderr=log, cwd=HERE)
    votes = collections.Counter()
    pairs = set()
    for i in range(1, MODELS + 1):
        for r in csv.DictReader((work / f"gda_rst_{i:02d}.tsv").open(encoding="utf-8"), delimiter="\t"):
            key = (r["pmid"], r["geneId"], r["diseaseId"])
            pairs.add(key)
            votes[key] += int(float(r["pred"]))
    positive = collections.defaultdict(list)
    for key in pairs:
        if votes[key] >= CUTOFF:
            positive[key[0]].append((id_of[key[0]][key[1]], id_of[key[0]][key[2]], votes[key]))
    for paper in sample["evaluation"]:
        pid = paper.removeprefix("BTX")
        target = paper.replace("BTX", "BTXA") if args.abstract else paper
        pt = RUNS / "pubtator3" / paper / "outputs"
        out = RUNS / name / target / "outputs"
        out.mkdir(parents=True, exist_ok=True)
        entities = [e for e in read(pt / "entities.csv") if e["entity_type"] in TYPES]
        write(out / "entities.csv", entities, ["entity_id", "entity_type", "canonical_name", "external_id"])
        keep = {e["entity_id"] for e in entities}
        write(out / "mentions.csv", [m for m in read(pt / "mentions.csv") if m["entity_id"] in keep],
              ["entity_id", "block_id", "start", "end", "surface"])
        names = {e["entity_id"]: e["canonical_name"] for e in entities}
        write(out / "pairs.csv", [{"pair_id": f"{g}-{d}", "gene_id": g, "disease_id": d, "gene_name": names[g],
                                   "disease_name": names[d], "status": "supported", "relation_type": "association",
                                   "source_scope": "", "reason": f"RENET2, {v} of {MODELS} models"}
                                  for g, d, v in sorted(positive[pid])],
              ["pair_id", "gene_id", "disease_id", "gene_name", "disease_name", "status", "relation_type", "source_scope", "reason"])
        (out / "source.json").write_text(json.dumps({"source": "RENET2 ft_models, 10-model vote, cutoff 5",
                                                     "input": "abstract" if args.abstract else "full text"}, indent=2), encoding="utf-8")
    print(f"{name}: {len(pairs)} candidate pairs, {sum(len(v) for v in positive.values())} positive")


if __name__ == "__main__":
    main()
