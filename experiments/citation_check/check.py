"""Do the sentences the minimal prompt cites support the relation it reports? Automatic check plus one annotator's blind check.

Usage: python experiments/citation_check/check.py            (automatic check, build the sheets)
       python experiments/citation_check/check.py --score    (after the sheets are filled in)

Rules (fixed before any sheet was built or answered):
- Systems and text: V6 on BioRED Test (100 abstracts) and on BioTriplex full text (locked 30), valid runs only. Unit: one
  supported pair with every sentence cited by its supported claims.
- Automatic check, all supported pairs: does one cited sentence hold mentions of both entities, do the cited sentences
  hold both only jointly, or does at least one entity have no mention in them? Mentions are the run's own. This is a
  coverage count: co-mention is neither necessary (pronouns) nor sufficient (co-occurrence) for support.
- Sheet: 30 pairs per corpus, stratified by the automatic check: 15 with both entities in one cited sentence and 15
  without, drawn with random.Random(20261107), at most 3 per article and stratum. Only the cited sentences are shown.
  An article-balanced draw was built first and discarded before anyone saw it: pairs without co-mention cluster in
  articles with many pairs, so it held 4 of 30 such BioRED pairs where the population has 42%.
- False-citation controls: 6 further pairs per corpus, shown with one uncited sentence of the same document that mentions
  exactly one of the two entities. They test that "not supported" is used; they are easy because the other entity is
  absent, so they bound gross rubber-stamping only. A gold-based positive control is not possible: neither corpus
  locates relations at sentence level.
- Retest: 20 of the 60 pairs in a second sheet with new numbers and a new order, to be answered at least 3 days after
  the first. It measures one annotator's consistency and does not replace a second annotator.
- Primary figure: share of supported pairs whose citations are "sufficient", per corpus: the two strata's shares
  weighted by their sizes in the automatic check, with a bootstrap interval over articles within strata. gold_match in
  the key records whether a gold relation overlaps both entities and is used descriptively.
Writes auto.json, sheet.html, retest.html, key.csv and, with --score, result.json. Exploratory: one annotator.
"""
import argparse
import collections
import csv
import html
import json
import random
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENTS = HERE.parent
sys.path[:0] = [str(EXPERIMENTS / "biored"), str(EXPERIMENTS / "single_shot")]
from common import load  # noqa: E402
from evaluate import overlap, read  # noqa: E402
from self_check import FORM  # noqa: E402

RUNS = EXPERIMENTS / "single_shot" / "runs"
VARIANT = "v6_minimal_lean_medium"
SEED, PER_STRATUM, MAX_PER_ARTICLE, FALSE_CITATIONS, RETEST, BOOTSTRAP = 20261107, 15, 3, 6, 20, 2000
EXPORT, RETEST_EXPORT = "citation_answers.csv", "citation_retest_answers.csv"
QUESTION = "Support (sufficient/needs context/hedged/not supported)"
OPTIONS = ["sufficient", "needs context", "hedged", "not supported"]
INSTRUCTIONS = """<p>基因标为<b class="g">蓝色</b>，疾病标为<b class="d">橙色</b>。每条只显示一组句子。只根据这些句子判断它们能否支持标题里的基因与疾病有关系，不要查资料，也不要凭生物学常识。</p>
<ul style="font-size:.92em">
<li><b>sufficient</b>：这些句子本身就陈述了两者的关联、相关或因果（包括负相关、治疗或保护作用）</li>
<li><b>needs context</b>：看起来在说这件事，但要靠句子以外的内容才能确定（如代词指代、缩写、只说了一半）</li>
<li><b>hedged</b>：只用了推测性措辞（may、might、suggest、possible 等）</li>
<li><b>not supported</b>：没有陈述两者的关系，或明确说没有关系，或只是同时出现</li></ul>
<p style="font-size:.92em">有的条目可能只标出一个实体，这时也照样判断。</p>
"""


def corpora():
    """{corpus: {paper: run output directory}} for valid runs."""
    split = {p["paper_id"]: p["split"] for p in json.loads((EXPERIMENTS / "biored" / "sample_all.json").read_text(encoding="utf-8"))["papers"]}
    biored = {p: RUNS / "biored_all" / VARIANT / p / "outputs" for p, s in split.items() if s == "Test"}
    rest = json.loads((EXPERIMENTS / "biored" / "gold_test_rest.json").read_text(encoding="utf-8"))
    biored.update({p: RUNS / "biored_testrest" / VARIANT / p / "outputs" for p in rest})
    sample = json.loads((EXPERIMENTS / "biotriplex" / "sample_biotriplex.json").read_text(encoding="utf-8"))
    biotriplex = {p: RUNS / "biotriplex" / VARIANT / p / "outputs" for p in sample["evaluation"]}
    return {name: {p: out for p, out in sorted(runs.items()) if (out / "validation.json").exists()}
            for name, runs in (("BioRED", biored), ("BioTriplex", biotriplex))}


def gold_sets():
    gold = json.loads((EXPERIMENTS / "biored" / "gold_all.json").read_text(encoding="utf-8"))
    gold.update(json.loads((EXPERIMENTS / "biored" / "gold_test_rest.json").read_text(encoding="utf-8")))
    gold.update(json.loads((EXPERIMENTS / "biotriplex" / "gold_biotriplex.json").read_text(encoding="utf-8")))
    return gold


def pairs_of(paper, out, gold):
    """Supported pairs of one article with their cited sentences and whether each cited sentence mentions each entity."""
    sentences = load(paper)[3]
    cited = {e["evidence_id"]: json.loads(e["sentence_ids"]) for e in read(out / "evidence.csv")}
    mentions = collections.defaultdict(list)
    for m in read(out / "mentions.csv"):
        mentions[m["entity_id"]].append((m["block_id"], int(m["start"]), int(m["end"])))
    spans = collections.defaultdict(list)
    for m in gold["mentions"]:
        spans[m["identifier"]].append((m["block_id"], m["start"], m["end"]))
    pairs = []
    for p in read(out / "pairs.csv"):
        if p["status"] != "supported":
            continue
        ids = [s for c in json.loads(p["claims"]) if c["status"] == "supported" for v in c["evidence_ids"] for s in cited.get(v, [])]
        ids = sorted(set(ids) & set(sentences), key=lambda s: (s.split(".")[0], int(s.split(".s")[1])))
        holds = lambda entity, sid: any(b == sentences[sid][0] and s < sentences[sid][2] and e > sentences[sid][1]
                                        for b, s, e in mentions[entity])
        both = [(holds(p["gene_id"], s), holds(p["disease_id"], s)) for s in ids]
        pairs.append({"paper": paper, "pair": p, "ids": ids, "sentences": sentences, "mentions": mentions,
                      "class": "uncited" if not ids else "one sentence" if any(g and d for g, d in both)
                      else "jointly" if any(g for g, _ in both) and any(d for _, d in both) else "entity absent",
                      "gold_match": any(overlap(spans[r["gene"]], mentions[p["gene_id"]])
                                        and overlap(spans[r["disease"]], mentions[p["disease_id"]])
                                        for r in gold["relations"])})
    return pairs


def stratum(automatic_class):
    return "co-mention" if automatic_class == "one sentence" else "no co-mention"


def draw(pairs, n, rng):
    """n random pairs, at most MAX_PER_ARTICLE from one article."""
    pool, per_article, drawn = list(pairs), collections.Counter(), []
    rng.shuffle(pool)
    for p in pool:
        if len(drawn) < n and per_article[p["paper"]] < MAX_PER_ARTICLE:
            per_article[p["paper"]] += 1
            drawn.append(p)
    return drawn


def false_citation(item, rng):
    """One uncited sentence of the document that mentions exactly one of the pair's entities, or None."""
    p, sentences, mentions = item["pair"], item["sentences"], item["mentions"]
    holds = lambda entity, sid: any(b == sentences[sid][0] and s < sentences[sid][2] and e > sentences[sid][1]
                                    for b, s, e in mentions[entity])
    options = [s for s in sorted(sentences) if s not in item["ids"] and holds(p["gene_id"], s) != holds(p["disease_id"], s)]
    return rng.choice(options) if options else None


def shown(item, ids):
    """The sentences as HTML with the pair's mentions marked."""
    blocks = {b["block_id"]: b["text"] for b in load(item["paper"])[1]}
    p, parts = item["pair"], []
    for sid in ids:
        block, start, end = item["sentences"][sid]
        marks = sorted({(s, e, cls) for cls, entity in (("g", p["gene_id"]), ("d", p["disease_id"]))
                        for b, s, e in item["mentions"][entity] if b == block and start <= s and e <= end})
        text, at = "", start
        for s, e, cls in marks:
            if s >= at:
                text += html.escape(blocks[block][at:s]) + f'<b class="{cls}">{html.escape(blocks[block][s:e])}</b>'
                at = e
        parts.append(f"<p>{text}{html.escape(blocks[block][at:end])}</p>")
    return "".join(parts)


def write_sheet(path, rows, store, download):
    data = [{"item_id": r["item_id"], "gene": html.escape(r["gene"]), "disease": html.escape(r["disease"]),
             "abstract_html": r["html"]} for r in rows]
    page = re.sub(r"<p>基因标为.*?</ul>\n", lambda m: INSTRUCTIONS, FORM, flags=re.S)
    page = page.replace("__ITEMS__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/"))
    page = page.replace("__STORE__", store).replace("__DOWNLOAD__", download).replace("Self-check annotation", "Citation check")
    question = [{"column": QUESTION, "title": "这些句子能支持两者有关系吗？", "options": OPTIONS}]
    path.write_text(page.replace("__QUESTIONS__", json.dumps(question, ensure_ascii=False)), encoding="utf-8")


def build():
    gold = gold_sets()
    rng = random.Random(SEED)
    auto, rows, key = {}, [], []
    for corpus, runs in corpora().items():
        pairs = [x for paper, out in runs.items() for x in pairs_of(paper, out, gold[paper])]
        count = collections.Counter(x["class"] for x in pairs)
        auto[corpus] = {"articles": len(runs), "articles_with_pairs": len({x["paper"] for x in pairs}), "supported_pairs": len(pairs),
                        "classes": dict(count), "share_one_sentence": round(count["one sentence"] / len(pairs), 3),
                        "cited_sentences_median": sorted(len(x["ids"]) for x in pairs)[len(pairs) // 2],
                        "cited_sentences_max": max(len(x["ids"]) for x in pairs),
                        "by_gold_match": {str(m): dict(collections.Counter(x["class"] for x in pairs if x["gold_match"] == m))
                                          for m in (True, False)}}
        citable = [x for x in pairs if x["ids"]]
        real = [x for name in ("co-mention", "no co-mention")
                for x in draw([y for y in citable if stratum(y["class"]) == name], PER_STRATUM, rng)]
        spare = [x for x in citable if not any(x is y for y in real)]
        rng.shuffle(spare)
        controls = []
        for item in spare:
            sid = false_citation(item, rng)
            if sid and len(controls) < FALSE_CITATIONS:
                controls.append((item, sid))
        if len(controls) < FALSE_CITATIONS:
            raise SystemExit(f"{corpus}: only {len(controls)} false-citation controls")
        rows += [{"kind": "real", "corpus": corpus, "item": x, "ids": x["ids"]} for x in real]
        rows += [{"kind": "false citation", "corpus": corpus, "item": x, "ids": [sid]} for x, sid in controls]
    rng.shuffle(rows)
    for n, r in enumerate(rows, 1):
        r.update(item_id=f"C{n:02d}", gene=r["item"]["pair"]["gene_name"], disease=r["item"]["pair"]["disease_name"],
                 html=shown(r["item"], r["ids"]))
    retest = rng.sample([r for r in rows if r["kind"] == "real"], RETEST)
    rng.shuffle(retest)
    retest_rows = [{**r, "item_id": f"R{n:02d}", "retest_of": r["item_id"]} for n, r in enumerate(retest, 1)]
    write_sheet(HERE / "sheet.html", rows, "citation_check_v1", EXPORT)
    write_sheet(HERE / "retest.html", retest_rows, "citation_retest_v1", RETEST_EXPORT)
    with (HERE / "key.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["item_id", "kind", "corpus", "paper_id", "pair_id", "sentences", "automatic_class",
                                                    "gold_match", "retest_of"])
        writer.writeheader()
        for r in rows + retest_rows:
            writer.writerow({"item_id": r["item_id"], "kind": r["kind"], "corpus": r["corpus"], "paper_id": r["item"]["paper"],
                             "pair_id": r["item"]["pair"]["pair_id"], "sentences": " ".join(r["ids"]),
                             "automatic_class": r["item"]["class"], "gold_match": r["item"]["gold_match"],
                             "retest_of": r.get("retest_of", "")})
    (HERE / "auto.json").write_text(json.dumps(auto, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(auto, ensure_ascii=False, indent=2))
    print(f"{len(rows)} items, {len(retest_rows)} retest; shown sentences per item: max {max(len(r['ids']) for r in rows)}, "
          f"words max {max(len(re.sub('<[^>]+>', '', r['html']).split()) for r in rows)}, "
          f"total words {sum(len(re.sub('<[^>]+>', '', r['html']).split()) for r in rows)}")


def answers(name):
    path = next(p for p in (HERE / name, Path.home() / "Downloads" / name) if p.exists())
    print(f"answers: {path}")
    return {a["item_id"]: a[QUESTION].strip() for a in csv.DictReader(path.open(encoding="utf-8-sig"))}


def weighted(items, weights):
    """Share sufficient, each stratum weighted by its population share; items are (stratum, answer) pairs."""
    total = 0
    for name, weight in weights.items():
        group = [a == "sufficient" for k, a in items if k == name]
        total += weight * sum(group) / len(group)
    return total


def weighted_answers(real, weights):
    """Each answer's share of supported pairs, strata weighted by population share."""
    shares = {}
    for option in OPTIONS:
        shares[option] = 0
        for name, weight in weights.items():
            group = [a for k, a in real if stratum(k["automatic_class"]) == name]
            shares[option] += weight * group.count(option) / len(group)
    return {o: round(v, 3) for o, v in shares.items()}


def interval(items, weights, rng):
    """Bootstrap interval over articles within strata."""
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for k, a in items:
        groups[stratum(k["automatic_class"])][k["paper_id"]].append(a)
    draws = []
    for _ in range(BOOTSTRAP):
        sample = [(name, a) for name, papers in groups.items()
                  for p in rng.choices(sorted(papers), k=len(papers)) for a in papers[p]]
        draws.append(weighted(sample, weights))
    draws.sort()
    return [round(draws[int(.025 * BOOTSTRAP)], 3), round(draws[int(.975 * BOOTSTRAP) - 1], 3)]


def score():
    key = {r["item_id"]: r for r in csv.DictReader((HERE / "key.csv").open(encoding="utf-8-sig"))}
    first, rng, report = answers(EXPORT), random.Random(SEED), {}
    auto = json.loads((HERE / "auto.json").read_text(encoding="utf-8"))
    for corpus in ("BioRED", "BioTriplex"):
        real = [(k, first[i]) for i, k in key.items() if k["corpus"] == corpus and k["kind"] == "real" and first.get(i)]
        controls = [first[i] for i, k in key.items() if k["corpus"] == corpus and k["kind"] == "false citation" and first.get(i)]
        share = auto[corpus]["share_one_sentence"]
        weights = {"co-mention": share, "no co-mention": 1 - share}
        report[corpus] = {"items": len(real), "answers": dict(collections.Counter(a for _, a in real)),
                          "sufficient_weighted": round(weighted([(stratum(k["automatic_class"]), a) for k, a in real], weights), 3),
                          "ci95": interval(real, weights, rng),
                          "weighted_answers": weighted_answers(real, weights),
                          "by_automatic_class": {c: dict(collections.Counter(a for k, a in real if k["automatic_class"] == c))
                                                 for c in sorted({k["automatic_class"] for k, _ in real})},
                          "by_gold_match": {m: dict(collections.Counter(a for k, a in real if k["gold_match"] == m))
                                            for m in ("True", "False")},
                          "false_citations": dict(collections.Counter(controls))}
    if any(p.exists() for p in (HERE / RETEST_EXPORT, Path.home() / "Downloads" / RETEST_EXPORT)):
        again = answers(RETEST_EXPORT)
        both = [(first[key[i]["retest_of"]], a) for i, a in again.items() if a and first.get(key[i]["retest_of"])]
        observed = sum(a == b for a, b in both) / len(both)
        chance = sum(sum(a == o for a, _ in both) * sum(b == o for _, b in both) for o in OPTIONS) / len(both) ** 2
        report["retest"] = {"items": len(both), "same_answer": sum(a == b for a, b in both),
                            "kappa": round((observed - chance) / (1 - chance), 3) if chance < 1 else None}
    (HERE / "result.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--score", action="store_true")
    score() if parser.parse_args().score else build()
