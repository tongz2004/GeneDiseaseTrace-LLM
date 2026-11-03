"""Offline replay of a finished run's first answers through deterministic pre-processing (no model calls).

Usage: python experiments/single_shot/normalize.py --run-set biored_all --variant v6_minimal_lean_medium [--filter] [--papers ...]

The model's first answer (response_first.txt when a repair call was made, else response.txt; both are the model's own
saved text) is rewritten before expansion, touching P records only. Two steps of different nature, reported apart:
- canonicalize (always): representation only, no prediction is removed.
  * a non-supported pair gets relation type '-' as the data contract requires (the format has no field for a tested
    type, so the model's original code is kept in normalize.json);
  * a supported pair's main relation type and source are recomputed from its claims (display priority causal >
    correlation > association; 'both' when background and this study support that type), as expand.py checks; the
    claims are taken as authoritative because they carry the evidence IDs, and every disagreement is logged.
- filter_invalid (--filter): a prediction filter. A pair whose endpoints are not a primary gene and a primary target of
  the same unit is dropped (e.g. a gene paired with a symptom entity) and logged. This changes system output; under the
  BioRED scorer such pairs can never be true positives, so it does not change BioRED scores.
An answer with lines the parser rejects (refusals, chatter, wrong field counts) is left unchanged, so it still fails.
Prompts are unchanged, but a pipeline with this pre-processing is a new system version. The source run is never
modified. Results go to runs/<run-set>/<variant>_canon[_filter]/<paper>/ with normalize.json; papers that still fail
keep problems.txt and count as failures (no repair call is replayed). Cost and calls are the first call's only, so the
difference from the source run estimates the avoidable repair cost; it is not money saved.
"""
import argparse
import json
import shutil
from collections import Counter

import batch
from common import RUNS
from expand import PRIORITY, RTYPE, SCOPE, STATUS, WIDTH, expand, parse, split_list, validate_run

GENE_SIDE = {"gene", "ambiguous_gene"}
TARGET_SIDE = {"disease", "clinical_event", "ambiguous_disease"}
RTYPE_CODE = {v: k for k, v in RTYPE.items()}
SCOPE_CODE = {v: k for k, v in SCOPE.items()}


def main_type_and_scope(claims):
    """(relation code, scope code) implied by the supported claims with evidence, or None."""
    supported = []
    for claim in split_list(claims, ";"):
        parts = claim.split(":")
        if len(parts) == 4 and parts[0] == "sup" and parts[1] in RTYPE and parts[2] in SCOPE and parts[3]:
            supported.append((RTYPE[parts[1]], SCOPE[parts[2]]))
    if not supported:
        return None
    top = max(supported, key=lambda c: PRIORITY.get(c[0], 0))[0]
    scopes = {s for t, s in supported if t == top}
    scope = "两者兼有" if {"背景引用", "本研究结果"} <= scopes or "两者兼有" in scopes else next(iter(scopes))
    return RTYPE_CODE[top], SCOPE_CODE[scope]


def normalize(text, filter_invalid=False):
    """(normalized text, [(kind, detail)]); the text comes back unchanged when the parser rejects any line."""
    problems = []
    rec = parse(text, problems)
    if problems:
        return text, []
    units = [u[0] for u in rec["UNIT"]] or ["u1"]
    entities = {e[0]: (e[1], e[8] or units[0]) for e in rec["E"]}
    changes, pairs = [], []
    for p in rec["P"]:
        gid, tid, status, rtype, scope, claims = p[:6]
        g, t = entities.get(gid), entities.get(tid)
        if filter_invalid and (not g or not t or g[0] not in GENE_SIDE or t[0] not in TARGET_SIDE or g[1] != t[1]):
            changes.append(("filter: dropped pair", f"P {gid}-{tid} ({g[0] if g else '?'} x {t[0] if t else '?'})"))
            continue
        if status != "sup" and status in STATUS and rtype != "-":
            changes.append(("canonicalize: type set to -", f"P {gid}-{tid}: {rtype} -> - ({status})"))
            p[3] = "-"
        elif status == "sup" and (implied := main_type_and_scope(claims)) and implied != (rtype, scope):
            changes.append(("canonicalize: main type/scope recomputed", f"P {gid}-{tid}: {rtype}/{scope} -> {implied[0]}/{implied[1]}"))
            p[3], p[4] = implied
        pairs.append(p)
    rec["P"] = pairs
    lines = ["|".join([kind] + fields) for kind in WIDTH for fields in rec[kind]]
    return "\n".join(lines) + "\n", changes


def run_paper(source, target, paper_id, filter_invalid):
    src, dst = source / paper_id, target / paper_id
    first = src / ("response_first.txt" if (src / "response_first.txt").exists() else "response.txt")
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    meta = json.loads((src / "run.json").read_text(encoding="utf-8"))
    meta.update(normalized_from=str(first.relative_to(RUNS)), calls=meta["calls"][:1])
    (dst / "run.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    text, changes = normalize(first.read_text(encoding="utf-8"), filter_invalid)
    (dst / "response.txt").write_text(text, encoding="utf-8")
    (dst / "normalize.json").write_text(json.dumps(changes, ensure_ascii=False, indent=2), encoding="utf-8")
    problems = expand(dst)
    if problems:
        (dst / "problems.txt").write_text("\n".join(problems), encoding="utf-8")
        return "failed", changes, problems
    _, failed = validate_run(dst)
    return ("ok" if not failed else "validation failed"), changes, []


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-set", required=True)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--papers", nargs="+")
    parser.add_argument("--filter", action="store_true", help="also drop pairs with invalid endpoints (filter_invalid)")
    args = parser.parse_args()
    source = RUNS / args.run_set / args.variant
    target = RUNS / args.run_set / f"{args.variant}_canon{'_filter' if args.filter else ''}"
    papers = args.papers or sorted(p.name for p in source.iterdir() if (p / "run.json").exists())
    results, kinds, still = Counter(), Counter(), Counter()
    for paper_id in papers:
        result, changes, problems = run_paper(source, target, paper_id, args.filter)
        results[result] += 1
        kinds.update(kind for kind, _ in changes)
        still.update(p.split(":", 1)[-1].strip()[:70] for p in problems)
        if result != "ok":
            print(f"{paper_id}: {result}: {problems[:1]}")
    print(f"\n{target.name}: {dict(results)} of {len(papers)}; changes: {dict(kinds)}")
    if still:
        print("remaining problems:", dict(still.most_common(8)))
    batch.OUT_ROOT = RUNS / args.run_set
    batch.summarise(target.name, papers)  # summary.csv; cost and calls are the first call's only


if __name__ == "__main__":
    main()
