"""EU-ADR gene–disease subset as input for the frozen single-shot method, plus a pair-level gold standard.

Usage: python experiments/euadr/prepare.py

Input: data/src_*.json, the EU-ADR corpus (van Mulligen et al., J Biomed Inform 2012) in the BigBio source schema,
fetched as JSON rows from the Hugging Face dataset server (bigbio/euadr, config euadr_source). Licence unstated: the
data stays local and is not redistributed.

Scope: the 100 abstracts of the Target-Disorder subset; consensus annotations only (flag True: kept by at least two
annotators). Gold relations are gene–disease pairs (target type "Genes & Molecular Sequences", disorder type
"Diseases & Disorders"); SNP–disorder relations are counted but not scored, as in the BioRED evaluation.
Labels: PA positive association, SA speculative association, NA negative association (the text states there is
none). Relations are annotated per sentence; a concept pair takes PA if any sentence gives PA, else SA, else NA.
Concepts are keyed by their first database ID (UniProt / UMLS), else by lower-cased text.

Output: papers/EUADR<pmid>.xml (title + abstract as minimal JATS, read by the unchanged pipeline through
single_shot/common.paper_path), gold_euadr.json (mentions mapped onto our blocks, relations with labels) and
sample_euadr.json (the frozen abstract list and counts).
"""
import collections
import difflib
import hashlib
import json
import sys
from pathlib import Path
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT))
from biomed.jats import parse_article  # noqa: E402

def offset_map(original, parsed):
    """Position in the parsed block for every position in the original text. Unlike BioRED, parsing can change more
    than whitespace here (some abstracts carry mis-encoded characters such as 'KrÃ¼ppel', which Unicode
    normalization rewrites), so the two strings are aligned with difflib."""
    mapping = [None] * (len(original) + 1)
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, original, parsed, autojunk=False).get_opcodes():
        for k in range(i1, i2 + 1):
            mapping[k] = j1 + (k - i1 if tag == "equal" else min(k - i1, j2 - j1))
    return mapping

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", errors="replace")

DATA, PAPERS = HERE / "data", HERE / "papers"
GENE, DISEASE, SNP = "Genes & Molecular Sequences", "Diseases & Disorders", "SNP & Sequence variations"
TYPE = {GENE: "GeneOrGeneProduct", DISEASE: "DiseaseOrPhenotypicFeature", SNP: "SequenceVariant"}
RANK = {"PA": 0, "SA": 1, "NA": 2}


def concept_key(ids_field, text):
    ids = [i for i in json.loads(ids_field.replace("'", '"')) if not i.startswith("sda/")]
    return ids[0] if ids else "text:" + text.lower().strip()


def locate(full, text, start, end):
    """Offsets as annotated, or the nearest occurrence of the annotated text when the offsets drifted."""
    if not text or full[start:end] == text:
        return start, end
    hits, i = [], full.find(text)
    while i != -1:
        hits.append(i)
        i = full.find(text, i + 1)
    if not hits:
        return None
    best = min(hits, key=lambda h: abs(h - start))
    return best, best + len(text)


def write_jats(paper_id, pmid, title, abstract):
    PAPERS.mkdir(exist_ok=True)
    xml = ('<?xml version="1.0" encoding="UTF-8"?>\n<article article-type="research-article"><front><article-meta>'
           f'<article-id pub-id-type="pmid">{pmid}</article-id>'
           f'<title-group><article-title>{escape(title)}</article-title></title-group>'
           f'<abstract><p>{escape(abstract)}</p></abstract></article-meta></front><body/></article>\n')
    path = PAPERS / f"{paper_id}.xml"
    path.write_text(xml, encoding="utf-8")
    return path


def main():
    files = sorted(DATA.glob("src_*.json"))
    docs = [r["row"] for f in files for r in json.loads(f.read_text(encoding="utf-8"))["rows"]]
    docs = sorted((d for d in docs if d["annotations"] and d["annotations"][0].split("\t")[0] == "Target-Disorder"),
                  key=lambda d: int(d["pmid"]))
    gold, stats = {}, collections.Counter()
    for doc in docs:
        pmid, title, abstract = doc["pmid"], doc["title"], doc["abstract"]
        paper_id = f"EUADR{pmid}"
        full = title + " " + abstract
        rows = [line.rstrip("\n").split("\t") for line in doc["annotations"]]
        concepts = {x[8]: x for x in rows if x[2] == "concept" and x[1] == "True"}
        all_concepts = {x[8]: x for x in rows if x[2] == "concept"}
        article = parse_article(write_jats(paper_id, pmid, title, abstract))
        blocks = {b["block_id"]: b["text"] for b in article["blocks"]}
        texts, starts = {"b0001": title, "b0002": abstract}, {"b0001": 0, "b0002": len(title) + 1}
        maps = {b: offset_map(texts[b], blocks[b]) for b in blocks}

        def to_block(start, end):
            block_id = "b0001" if start < len(title) else "b0002"
            lo, hi = start - starts[block_id], end - starts[block_id]
            if lo < 0 or hi > len(texts[block_id]):
                return None
            return block_id, maps[block_id][lo], maps[block_id][hi]

        mentions, seen = [], set()

        def add_mention(key, kind, text, start, end):
            spot = locate(full, text, start, end)
            placed = to_block(*spot) if spot else None
            if not placed:
                stats["mentions dropped (text not in abstract)"] += 1
                return
            item = (key, *placed)
            if item not in seen:
                seen.add(item)
                mentions.append({"identifier": key, "type": TYPE[kind], "text": full[spot[0]:spot[1]],
                                 "block_id": placed[0], "start": placed[1], "end": placed[2]})

        for x in concepts.values():
            if x[9] in TYPE:
                add_mention(concept_key(x[7], x[3]), x[9], x[3], int(x[4]), int(x[5]))
        labels = collections.defaultdict(set)
        for x in rows:
            if x[2] != "relation" or x[1] != "True":
                continue
            a, b = all_concepts.get(x[3]), all_concepts.get(x[4])
            if not a or not b:
                stats["relations with unknown argument"] += 1
                continue
            kinds = (a[9], b[9])
            if kinds == (SNP, DISEASE):
                stats[f"SNP–disease {x[10]} (not scored)"] += 1
                continue
            if kinds != (GENE, DISEASE):
                continue
            gene, disease = concept_key(x[5], a[3]), concept_key(x[6], b[3])
            for key, kind, span, c in ((gene, GENE, x[7], a), (disease, DISEASE, x[8], b)):
                s, e = map(int, span.split(":"))
                add_mention(key, kind, c[3] or full[s:e], s, e)
            labels[(gene, disease)].add(x[10])
            stats[f"gene–disease relations {x[10]}"] += 1
        relations = []
        for (gene, disease), found in labels.items():
            label = min(found, key=RANK.get)
            if len(found) > 1:
                stats["pairs with mixed labels"] += 1
            relations.append({"gene": gene, "disease": disease, "type": label, "novel": "", "labels": sorted(found)})
            stats[f"gene–disease pairs {label}"] += 1
        gold[paper_id] = {"blocks": list(blocks), "mentions": mentions, "relations": relations}
    (HERE / "gold_euadr.json").write_text(json.dumps(gold, ensure_ascii=False, indent=1), encoding="utf-8")
    sample = {"source": "bigbio/euadr (euadr_source), Hugging Face dataset server; EU-ADR corpus, van Mulligen et al. 2012",
              "data_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
              "criteria": "all 100 Target-Disorder abstracts; consensus annotations; gene–disease pairs scored",
              "papers": [{"paper_id": p, "gold_pairs": len(g["relations"]),
                          "labels": dict(collections.Counter(r["type"] for r in g["relations"]))} for p, g in gold.items()]}
    (HERE / "sample_euadr.json").write_text(json.dumps(sample, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{len(gold)} abstracts; {sum(1 for g in gold.values() if g['relations'])} with a gene–disease pair")
    for k, v in sorted(stats.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
