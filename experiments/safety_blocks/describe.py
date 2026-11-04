"""Which full texts did the provider's safety filter stop? A description of stored inputs and run records, with no model calls.

Usage: python experiments/safety_blocks/describe.py

Articles: every full text run with the minimal prompt (V6, Claude Opus 5.5) at or after November 3: BioTriplex locked (30),
BioTriplex further (49), the original 100 papers. Blocked = no valid result after one unchanged retry, with the last
attempt ending in a safety flag or a refusal citing a safety classifier. Articles that failed on format only are left out.

Term groups were fixed before any count was seen. For each group: terms per 1,000 words of the text the model read,
median in blocked and completed articles, and the AUC (chance that a blocked article has the higher density).
With several groups and no correction, this is exploratory and only describes the articles. It is meant for reporting
possible false positives to the provider and for describing the failure in the paper, not for changing text or prompts
to get past the filter.
Writes summary.json and blocked_articles.csv (IDs and titles of published open-access articles).
"""
import collections
import csv
import json
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
PROJECT = EXPERIMENTS.parent
sys.path.insert(0, str(EXPERIMENTS / "single_shot"))
from common import load  # noqa: E402

RUNS = EXPERIMENTS / "single_shot" / "runs"
VARIANT = "v6_minimal_lean_medium"
LOCKED_BLOCKED = {"BTX3320601", "BTX4302655", "BTX4550874", "BTX5444857"}  # README: blocked in both attempts
GROUPS = {
    "infection": r"infect\w*|pathogen\w*|virus\w*|viral|bacteri\w*|sepsis|septic|parasit\w*|fung\w+|microb\w*",
    "toxin_or_poison": r"toxin\w*|toxic\w*|poison\w*|venom\w*|lethal\w*",
    "vaccine_or_immunity": r"vaccin\w*|immuni[sz]\w*|antibod\w*|antigen\w*|inflamm\w*",
    "animal_experiment": r"mice|mouse|murine|rats?|rodent\w*|knockout|transgenic",
    "gene_editing": r"crispr\w*|cas9|knock-?down|sirna|shrna|plasmid\w*|vector\w*|transfect\w*",
    "dosing_or_drug": r"mg/kg|dose\w*|dosing|inject\w*|administ\w*|drug\w*",
    "cancer": r"cancer\w*|tumou?r\w*|carcinoma\w*|neoplas\w*|metasta\w*",
}


def articles():
    """[(set, paper, run output directory)]."""
    out = []
    for name, run_set, sample in (("locked", "biotriplex", "sample_biotriplex.json"), ("further", "biotriplex2", "sample_biotriplex2.json")):
        for p in json.loads((EXPERIMENTS / "biotriplex" / sample).read_text(encoding="utf-8"))["evaluation"]:
            out.append((name, p, RUNS / run_set / VARIANT / p))
    with (PROJECT / "outputs" / "inventory.csv").open(encoding="utf-8-sig") as stream:
        out += [("original100", r["paper_id"], RUNS / "full100" / VARIANT / r["paper_id"]) for r in csv.DictReader(stream)]
    return out


def outcomes():
    """{paper: last logged line} for runs whose logs we keep."""
    last = {}
    for name in ("biotriplex2_batch_first.log", "biotriplex2_batch_retry.log", "full100_v6_batch_first.log", "full100_v6_batch_rest.log"):
        for line in (RUNS / name).read_text(encoding="utf-8", errors="replace").splitlines():
            m = re.match(r"((?:BTX|PMC)\d+): (.*)", line)
            if m:
                last[m.group(1)] = m.group(2)
    return last


def status(paper, run, last):
    if (run / "outputs" / "validation.json").exists():
        return "completed"
    first = run / "response_first.txt"
    if first.exists() and "safety classifier" in first.read_text(encoding="utf-8", errors="replace"):
        return "blocked"
    if paper in LOCKED_BLOCKED or re.search(r"safeguards|safety|withheld|not going to continue|won.t continue|Understood: I won", last.get(paper, "")):
        return "blocked"
    return "format failure"


def median(values):
    values = sorted(values)
    return round(values[len(values) // 2], 2) if values else None


def auc(a, b):
    wins = sum((x > y) + .5 * (x == y) for x in a for y in b)
    return wins / (len(a) * len(b))


def main():
    last, rows = outcomes(), []
    for name, paper, run in articles():
        article, kept, _, _ = load(paper)
        text = " ".join(b["text"] for b in kept)
        words = max(1, len(text.split()))
        density = {g: 1000 * len(re.findall(rf"\b(?:{p})\b", text, flags=re.I)) / words for g, p in GROUPS.items()}
        rows.append({"set": name, "paper": paper, "status": status(paper, run, last), "words": words,
                     "title": article.get("title") or kept[0]["text"][:160], **density})
    report = {"by_set": {}, "groups": {}}
    for name in ("locked", "further", "original100"):
        report["by_set"][name] = dict(collections.Counter(r["status"] for r in rows if r["set"] == name))
    report["all"] = dict(collections.Counter(r["status"] for r in rows))
    blocked = [r for r in rows if r["status"] == "blocked"]
    done = [r for r in rows if r["status"] == "completed"]
    report["words_median"] = {"blocked": median([r["words"] for r in blocked]), "completed": median([r["words"] for r in done])}
    rng = random.Random(20261028)
    for g in list(GROUPS) + ["words"]:
        a, b = [r[g] for r in blocked], [r[g] for r in done]
        draws = sorted(auc(rng.choices(a, k=len(a)), rng.choices(b, k=len(b))) for _ in range(1000))
        report["groups"][g] = {"blocked_median": median(a), "completed_median": median(b), "auc": round(auc(a, b), 3),
                               "auc_ci95": [round(draws[24], 3), round(draws[974], 3)]}
    (HERE / "summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    with (HERE / "blocked_articles.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["set", "paper", "words", "title"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(blocked)
    print(json.dumps(report, indent=1))
    for r in blocked:
        print(r["set"], r["paper"], r["words"], r["title"][:95])


if __name__ == "__main__":
    main()
