"""Second pass over a finished run: ask the model to judge the gene–target pairs it left unlisted although both are
mentioned in the same sentence.

Usage: python experiments/single_shot/second_pass.py --run-set biored_all --source v1_card_lean_medium \
           --target v1_pass2_medium --papers PMID... [--parallel 3]

Candidates are found by script, not by the model: gene and target entities (from the first answer) whose mentions
share a sentence, and whose pair the first answer left unlisted (status nf, reason 'unlisted pair').
For each paper with candidates, one extra request gets the same card and format as the source variant (cached
system prompt), the paper, the first answer and the candidate list, and returns only P lines for those pairs plus EV
lines for new evidence (IDs W1, W2, ...). The new lines are appended to the first answer, which is then expanded and
validated as usual. If the merged answer fails the checker, the paper keeps its first answer (counted as fallback).
The source run is never modified; results go to runs/<run-set>/<target>/.
"""
import argparse
import csv
import json
import shutil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import batch
from common import RUNS, load
from expand import expand, validate_run
from pack import paper_part
from run import call

INSTRUCTIONS = """# Second pass: pairs not yet judged
Your earlier answer for this paper is below. It left the gene–target pairs listed after it unjudged, although each
gene and target are mentioned in the same sentence. Judge each of these pairs now, by the same decision card.
Output only:
- exactly one P line per listed pair (status sup / neg / unc / nf; same field rules and codes as before), using the
  entity IDs given;
- EV lines for any new evidence you cite, with evidence IDs W1, W2, ... (you may also cite evidence IDs from your
  earlier answer).
No other lines."""


def read(path):
    with path.open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def candidates(paper_id, out):
    """Unlisted gene–target pairs whose mentions share a sentence: [(gene_id, target_id, gene, target, sentence IDs)]."""
    _, _, _, sentences = load(paper_id)
    spans = defaultdict(list)
    for m in read(out / "mentions.csv"):
        spans[m["entity_id"]].append((m["block_id"], int(m["start"]), int(m["end"])))
    in_sentence = lambda eid: {sid for sid, (b, s, e) in sentences.items() for bb, a, z in spans[eid] if bb == b and a < e and s < z}
    found = []
    for p in read(out / "pairs.csv"):
        if p["status"] == "not_found" and "unlisted" in p["reason"]:
            shared = sorted(in_sentence(p["gene_id"]) & in_sentence(p["disease_id"]))
            if shared:
                found.append((p["gene_id"], p["disease_id"], p["gene_name"], p["disease_name"], shared))
    return found


def run_one(args, variant_config, fixed, paper_id):
    source = RUNS / args.run_set / args.source / paper_id
    target = RUNS / args.run_set / args.target / paper_id
    if (target / "outputs" / "validation.json").exists():
        return paper_id, "skipped (done)", 0
    shutil.rmtree(target, ignore_errors=True)
    shutil.copytree(source, target)
    meta = json.loads((target / "run.json").read_text(encoding="utf-8"))
    found = candidates(paper_id, source / "outputs")
    meta["second_pass"] = {"candidates": len(found)}
    if not found:
        (target / "run.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
        return paper_id, "no candidates", 0
    first = (source / "response.txt").read_text(encoding="utf-8")
    listing = "\n".join(f"{g}|{t}|{gn} / {tn}|{', '.join(s)}" for g, t, gn, tn, s in found)
    prompt = (paper_part(paper_id)[0] + "\n\n" + INSTRUCTIONS + "\n\n# Your earlier answer\n" + first +
              "\n\n# Pairs to judge (gene_id|target_id|names|sentences where both appear)\n" + listing)
    workdir = target / "empty_workdir"
    workdir.mkdir(exist_ok=True)
    text, usage = call(batch.MODEL, variant_config["effort"], prompt, workdir, fixed)
    meta["calls"].append({**usage, "second_pass": True})
    if usage.get("exit_code") or not text.strip():
        shutil.rmtree(target, ignore_errors=True)
        return paper_id, "not started (second-pass call failed or empty)", 0
    wanted = {(g, t) for g, t, *_ in found}
    added = [line for line in text.splitlines()
             if line.startswith("EV|W") or (line.startswith("P|") and tuple(line.split("|")[1:3]) in wanted)]
    (target / "response_second_pass.txt").write_text(text, encoding="utf-8")
    (target / "response.txt").write_text(first.rstrip("\n") + "\n" + "\n".join(added) + "\n", encoding="utf-8")
    problems = expand(target)
    if problems:
        (target / "response.txt").write_text(first, encoding="utf-8")
        expand(target)
        meta["second_pass"].update(fallback=True, problems=problems[:5])
    meta["second_pass"]["lines_added"] = len(added)
    (target / "run.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    _, failed = validate_run(target)
    note = "fallback to first answer" if problems else f"{sum(l.startswith('P|') for l in added)}/{len(found)} pairs judged"
    return paper_id, ("ok, " if not failed else "validation failed, ") + note, len(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-set", required=True)
    parser.add_argument("--source", required=True, help="finished variant to start from")
    parser.add_argument("--target", required=True, help="new variant folder for the merged results")
    parser.add_argument("--papers", nargs="+", required=True)
    parser.add_argument("--parallel", type=int, default=3)
    args = parser.parse_args()
    variant_config, fixed, _ = batch.variant_files(args.source)
    with ThreadPoolExecutor(max_workers=args.parallel) as pool:
        for future in as_completed(pool.submit(run_one, args, variant_config, fixed, p) for p in args.papers):
            paper_id, result, n = future.result()
            print(f"{paper_id}: {result}", flush=True)


if __name__ == "__main__":
    main()
