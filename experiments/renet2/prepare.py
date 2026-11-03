"""RENET2 full-text articles as input for the frozen single-shot method, plus a pair-level gold standard.

Usage: python experiments/renet2/prepare.py

Input: data/data/ft_data/ from renet2_data_models.tar.gz (Su et al., NAR Genomics and Bioinformatics 2021; authors'
server, https://www.bio8.cs.hku.hk/RENET2/). Only the full-text data folders were kept. Data licence unstated: kept
local, not redistributed.

Text: each article's title, abstract and main text as given in docs.txt. The main text there is a truncated part of
the article (introduction and part of the results; median about 5,400 characters per article), and RENET2's labels
were made on this same text, so the method reads exactly that.
Gold: labels.txt, one row per (PMID, NCBI Gene ID, MeSH ID): 1.0 associated (predicted, then confirmed by experts),
0.5 ambiguous, 0.0 non-associated (not checked by hand: three tools agreed). Mentions: anns.txt (PubTator Central),
offsets into title + " " + abstract + " " + main text.
Sample: articles with at least one associated pair and at least 1,000 characters of main text; 40 drawn with
random.Random(20261028) from the sorted PMIDs, frozen in sample_renet2_40.json.
Output: papers/RENET<pmid>.xml (title, abstract and main text as minimal JATS) and gold_renet2.json.
"""
import csv
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT))
from biomed.jats import parse_article  # noqa: E402

_spec = importlib.util.spec_from_file_location("euadr_prepare", HERE.parent / "euadr" / "prepare.py")
_euadr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_euadr)
offset_map = _euadr.offset_map  # difflib alignment of original text to parsed blocks

DATA = HERE / "data" / "data" / "ft_data"
PAPERS = HERE / "papers"
SAMPLE = HERE / "sample_renet2_40.json"
SEED, SIZE, MIN_BODY = 20261028, 40, 1000
TYPE = {"Gene": "GeneOrGeneProduct", "Disease": "DiseaseOrPhenotypicFeature"}
LABEL = {"1.0": "associated", "0.5": "ambiguous", "0.0": "non-associated"}


def read_docs():
    lines = (DATA / "docs.txt").read_text(encoding="utf-8").split("\n")
    docs, i = {}, 0
    while i + 3 < len(lines):
        if lines[i].strip().isdigit():
            docs[lines[i].strip()] = (lines[i + 1], lines[i + 2], lines[i + 3])
            i += 5
        else:
            i += 1
    return docs


def write_jats(paper_id, pmid, title, abstract, body):
    PAPERS.mkdir(exist_ok=True)
    xml = ('<?xml version="1.0" encoding="UTF-8"?>\n<article article-type="research-article"><front><article-meta>'
           f'<article-id pub-id-type="pmid">{pmid}</article-id>'
           f'<title-group><article-title>{escape(title)}</article-title></title-group>'
           f'<abstract><p>{escape(abstract)}</p></abstract></article-meta></front>'
           f'<body><p>{escape(body)}</p></body></article>\n')
    path = PAPERS / f"{paper_id}.xml"
    path.write_text(xml, encoding="utf-8")
    return path


def main():
    docs = read_docs()
    labels = list(csv.DictReader((DATA / "labels.txt").open(encoding="utf-8")))
    positive = {r["pmid"] for r in labels if r["label"] == "1.0"}
    eligible = sorted((p for p in positive if p in docs and len(docs[p][2]) >= MIN_BODY), key=int)
    picked = sorted(random.Random(SEED).sample(eligible, SIZE), key=int)
    files = [DATA / n for n in ("docs.txt", "anns.txt", "labels.txt")]
    record = {"source": "renet2_data_models.tar.gz, data/ft_data (Su et al. 2021)",
              "data_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
              "criteria": f">= 1 associated pair and >= {MIN_BODY} characters of main text", "seed": SEED,
              "articles": len(docs), "eligible": len(eligible), "papers": []}
    anns = {}
    for line in (DATA / "anns.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 6 and parts[4] in TYPE:
            anns.setdefault(parts[0], []).append(parts)
    gold, mismatches = {}, []
    for pmid in picked:
        title, abstract, body = docs[pmid]
        paper_id = f"RENET{pmid}"
        article = parse_article(write_jats(paper_id, pmid, title, abstract, body))
        blocks = {b["block_id"]: b["text"] for b in article["blocks"]}
        parts = [("b0001", title), ("b0002", abstract), ("b0003", body)]
        if list(blocks) != [b for b, _ in parts]:
            raise SystemExit(f"{paper_id}: unexpected blocks {list(blocks)}")
        starts, pos = {}, 0
        for block_id, text in parts:
            starts[block_id] = pos
            pos += len(text) + 1
        maps = {b: offset_map(t, blocks[b]) for b, t in parts}
        mentions = []
        for a in anns.get(pmid, []):
            start, end = int(a[1]), int(a[2])
            block_id = next(b for b, t in reversed(parts) if start >= starts[b])
            lo, hi = start - starts[block_id], end - starts[block_id]
            text = dict(parts)[block_id]
            if text[lo:hi] != a[3]:
                mismatches.append(f"{paper_id}: {a[3]!r} not at its offset")
                continue
            s, e = maps[block_id][lo], maps[block_id][hi]
            for identifier in a[5].replace(",", ";").split(";"):
                if identifier and identifier != "-":
                    mentions.append({"identifier": identifier, "type": TYPE[a[4]], "text": blocks[block_id][s:e],
                                     "block_id": block_id, "start": s, "end": e})
        relations = [{"gene": r["geneId"], "disease": r["diseaseId"], "type": LABEL[r["label"]], "novel": ""}
                     for r in labels if r["pmid"] == pmid]
        gold[paper_id] = {"blocks": list(blocks), "mentions": mentions, "relations": relations}
        record["papers"].append({"paper_id": paper_id, "pmid": pmid, "characters": len(title) + len(abstract) + len(body),
                                 **{k: sum(r["type"] == k for r in relations) for k in LABEL.values()}})
    if SAMPLE.exists():
        frozen = json.loads(SAMPLE.read_text(encoding="utf-8"))
        if [p["pmid"] for p in frozen["papers"]] != picked or frozen["data_sha256"] != record["data_sha256"]:
            raise SystemExit("sample_renet2_40.json no longer matches the data and seed; refusing to redraw")
    else:
        SAMPLE.write_text(json.dumps(record, indent=2), encoding="utf-8")
    (HERE / "gold_renet2.json").write_text(json.dumps(gold, ensure_ascii=False, indent=1), encoding="utf-8")
    rel = [r for g in gold.values() for r in g["relations"]]
    chars = sorted(p["characters"] for p in record["papers"])
    print(f"RENET2: {len(docs)} articles, {len(eligible)} eligible; sample {len(picked)} (seed {SEED})")
    print(f"characters per article: median {chars[len(chars) // 2]}, range {chars[0]}-{chars[-1]}")
    print("gold pairs:", {k: sum(r['type'] == k for r in rel) for k in LABEL.values()})
    print(f"mentions placed: {sum(len(g['mentions']) for g in gold.values())}; offset mismatches: {len(mismatches)}")
    for m in mismatches[:10]:
        print("  " + m)


if __name__ == "__main__":
    main()
