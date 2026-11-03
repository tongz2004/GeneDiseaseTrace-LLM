"""LLM judges on the blinded adjudication package, standing in for human annotators (a stand-in, not a gold standard).

Usage: python experiments/biored/llm_adjudicate.py [--models GPT-6.1_sol GPT-6_Luna] [--limit N] [--workers 6]
       python experiments/biored/llm_adjudicate.py --score

Each judge sees exactly what an annotator sees, one item per call, with every tool disabled:
- round 1: the abstract with the two entities marked, answered per the annotator guide (README.md);
- round 2: an evidence packet plus the judge's own round-1 conclusion for that item (not the abstract).
Answers go to <package>/llm/<model>/round{1,2}.jsonl; reruns skip items already answered.
--score joins the answers with key_coordinator.csv and writes <package>/llm/agreement.json:
judge-vs-judge Cohen's kappa, and per stratum how often each judge finds a stated relation. The strata whose
pairs match a BioRED gold relation check the judges against the gold standard before their verdict on the
disputed pairs (V6 only, no gold relation) is read.
"""
import argparse
import csv
import json
import sys
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from evaluate import HERE

sys.path.insert(0, str(HERE.parent / "single_shot"))
import run  # noqa: E402

run.MODELS["GPT-6_Luna"] = ("codex", "gpt-6-luna")
PACKAGE = HERE / "adjudication" / "Dev_v1_card_lean_medium_vs_v6_minimal_lean_medium"
OUT = PACKAGE / "llm"
STATES = ("relation", "hedged relation", "explicitly no relation", "no statement")
ENTITIES = ("yes", "gene wrong", "disease wrong", "both wrong")
SOURCES = ("this study", "background", "both", "unclear")
VERDICTS = ("sufficient", "insufficient", "contradicts")

ROUND1 = """You are annotating a biomedical abstract for a study of gene-disease relation extraction.
The gene is marked [[like this]] and the disease {{like this}} (the gene may be unmarked if its mention sits inside the
disease mention; the gene name below is authoritative). Judge only what the abstract says. Do not use outside sources
or biological background knowledge.

Gene: {gene}
Disease: {disease}

Abstract:
{abstract}

Answer four questions.
1. entities_correct: is the marked text really this gene and this disease? One of: yes / gene wrong / disease wrong / both wrong
2. paper_states: what does the abstract say about the relation between the two? One of:
   - relation: states an association, correlation or causal link (including a negative correlation, a treatment or protective effect);
   - hedged relation: only speculative wording (may, might, suggest, possible and the like);
   - explicitly no relation: states there is no association;
   - no statement: they merely co-occur, or the abstract states no relation between them.
3. source: is that statement a result of this study, background knowledge, both, or unclear? One of: this study / background / both / unclear
4. notes: anything you were unsure about, in one sentence.

Output one JSON object and nothing else:
{{"entities_correct": "...", "paper_states": "...", "source": "...", "notes": "..."}}"""

ROUND2 = """You are checking whether a set of sentences cited from a biomedical abstract supports a conclusion.
Earlier, reading the whole abstract, you concluded about the gene and disease below: "{state}".

Gene: {gene}
Disease: {disease}

Cited sentences (separated by |):
{evidence}

Using only these sentences, answer one of:
- sufficient: these sentences alone are enough to reach your earlier conclusion;
- insufficient: other parts of the abstract are needed to reach it;
- contradicts: these sentences contradict your earlier conclusion.

Output one JSON object and nothing else:
{{"sufficient": "...", "notes": "one sentence"}}"""


def rows(name):
    return list(csv.DictReader((PACKAGE / name).open(encoding="utf-8-sig")))


def answered(path):
    """Latest line per id; a saved raw answer the current parser accepts fills an answer recorded as failed."""
    out = {}
    for line in path.open(encoding="utf-8") if path.exists() else ():
        x = json.loads(line)
        if not x["answer"]:
            field, allowed = ("paper_states", STATES) if path.stem == "round1" else ("sufficient", VERDICTS)
            x["answer"] = parse(x["raw"], field, allowed)
        out[x["id"]] = x
    return out


def parse(text, field, allowed):
    try:
        data = json.loads(text[text.index("{"):text.rindex("}") + 1])
    except ValueError:
        return None
    # Judges sometimes copy the option with its gloss ("relation: states an association ..."); keep the label.
    value = str(data.get(field, "")).split(":")[0].strip().lower()
    if value not in allowed:
        return None
    data[field] = value
    return data


def ask(model, prompt, field, allowed):
    """Two attempts; returns (answer or None, raw text)."""
    for _ in range(2):
        with tempfile.TemporaryDirectory() as work:
            text, usage = run.call(model, "medium", prompt, Path(work))
        data = parse(text, field, allowed)
        if data:
            return data, text, usage
    return None, text, usage


def judge(model, limit, workers):
    out = OUT / model
    out.mkdir(parents=True, exist_ok=True)
    jobs = [("round1.jsonl", r["item_id"], ROUND1.format(gene=r["gene"], disease=r["disease"], abstract=r["abstract"]),
             "paper_states", STATES) for r in rows("round1_relation.csv")[:limit]]
    run_jobs(model, out, jobs, workers)
    first = answered(out / "round1.jsonl")
    jobs = [("round2.jsonl", r["packet_id"], ROUND2.format(state=first[r["item_id"]]["answer"]["paper_states"],
                                                         gene=r["gene"], disease=r["disease"], evidence=r["evidence"]),
             "sufficient", VERDICTS)
            for r in rows("round2_evidence.csv") if first.get(r["item_id"], {}).get("answer")]
    run_jobs(model, out, jobs[:limit], workers)


def run_jobs(model, out, jobs, workers):
    path = out / jobs[0][0] if jobs else None
    if not path:
        return
    done = {k for k, v in answered(path).items() if v.get("answer")}
    todo = [j for j in jobs if j[1] not in done]

    def one(job):
        _, key, prompt, field, allowed = job
        answer, text, usage = ask(model, prompt, field, allowed)
        return {"id": key, "answer": answer, "raw": text, "seconds": usage["seconds"], "exit_code": usage["exit_code"]}

    with ThreadPoolExecutor(workers) as pool, path.open("a", encoding="utf-8") as f:
        for result in pool.map(one, todo):
            f.write(json.dumps(result, ensure_ascii=False) + "\n")
            f.flush()
            print(f"{model} {path.stem} {result['id']}: {(result['answer'] or {}).get(jobs[0][3], 'FAILED')}")


def kappa(a, b):
    n = len(a)
    if not n:
        return None
    po = sum(x == y for x, y in zip(a, b)) / n
    ca, cb = Counter(a), Counter(b)
    pe = sum(ca[k] * cb[k] for k in set(ca) | set(cb)) / n / n
    return round((po - pe) / (1 - pe), 3) if pe < 1 else None


def score(models):
    key = {r["item_id"]: r for r in rows("key_coordinator.csv")}
    packets = {r["packet_id"]: r["item_id"] for r in rows("round2_evidence.csv")}
    r1 = {m: {k: v["answer"] for k, v in answered(OUT / m / "round1.jsonl").items() if v.get("answer")} for m in models}
    r2 = {m: {k: v["answer"] for k, v in answered(OUT / m / "round2.jsonl").items() if v.get("answer")} for m in models}
    report = {"note": "LLM judges stand in for human annotators; this is not a gold standard.", "models": models}
    common = sorted(set.intersection(*(set(r1[m]) for m in models)))
    states = {m: [r1[m][i]["paper_states"] for i in common] for m in models}
    stated = {m: ["relation" if s in ("relation", "hedged relation") else "none" for s in states[m]] for m in models}
    a, b = models
    report["round1_items_both"] = len(common)
    report["kappa_paper_states_4way"] = kappa(states[a], states[b])
    report["kappa_relation_stated_or_hedged_vs_not"] = kappa(stated[a], stated[b])
    report["raw_agreement_4way"] = round(sum(x == y for x, y in zip(states[a], states[b])) / len(common), 3)
    groups = {}
    for i in common:
        g = f"{key[i]['stratum']} / gold={key[i]['matches_gold_relation']}"
        groups.setdefault(g, []).append(i)
    report["by_stratum"] = {}
    for g, ids in sorted(groups.items()):
        entry = {"items": len(ids)}
        for m in models:
            c = Counter(r1[m][i]["paper_states"] for i in ids)
            entry[m] = {"states": dict(c), "relation_or_hedged": round((c["relation"] + c["hedged relation"]) / len(ids), 3),
                        "entities_wrong": sum(r1[m][i]["entities_correct"] != "yes" for i in ids)}
        both = [i for i in ids if all(r1[m][i]["paper_states"] in ("relation", "hedged relation") for m in models)]
        entry["both_judges_relation_or_hedged"] = len(both)
        report["by_stratum"][g] = entry
    report["disputed_items"] = [{"item": i, "stratum": key[i]["stratum"], "gold": key[i]["matches_gold_relation"],
                                 **{m: r1[m][i]["paper_states"] for m in models}}
                                for i, x, y in zip(common, states[a], states[b]) if x != y]
    report["round2"] = {}
    for m in models:
        by_system = {}
        for pid, ans in r2[m].items():
            item = key[packets[pid]]
            system = "v1" if pid == item["v1_card_lean_medium_packet"] else "v6"
            by_system.setdefault(system, Counter())[ans["sufficient"]] += 1
        report["round2"][m] = {s: dict(c) for s, c in by_system.items()}
    (OUT / "agreement.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "disputed_items"}, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--models", nargs=2, default=["GPT-6.1_sol", "GPT-6_Luna"])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--score", action="store_true")
    args = parser.parse_args()
    if args.score:
        return score(args.models)
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda m: judge(m, args.limit, args.workers), args.models))


if __name__ == "__main__":
    main()
