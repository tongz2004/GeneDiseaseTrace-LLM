"""BioTriplex full-text articles as input for the frozen single-shot method, plus a pair-level gold standard.

Usage: python experiments/biotriplex/prepare.py

Input: data/splits.json from https://github.com/PanagiotisFytas/BioTriplex ("Preprocessed BioTriplex"; commit in
data/SOURCE_COMMIT). Fytas et al., Bioinformatics 2026 (doi:10.1093/bioinformatics/btag037). 100 PMC open-access
articles split by article (train 70 / val 15 / test 15), annotated per subsection with gene and disease mentions,
relation phrases and gene-disease-relation links. Code is MIT; the data licence is not stated separately, so the data,
generated articles and gold stay local.

Articles: each subsection becomes a JATS <sec> (the abstract, whose text starts with a "Title:" line, becomes the
<abstract>), one <p> per non-empty line. A " Rev" suffix on a subsection name (a revised annotation that replaces the
original) is dropped. Subsection order is not stored, so numbered names are ordered by number and the rest by a fixed
rank (introduction, methods, results, discussion, conclusion, others). BTXA<pmcid> holds the abstract alone, with the
same block IDs as the abstract of BTX<pmcid>, for the abstract-versus-full-text comparison.

Gold (frozen rules, see README):
- entities have no database IDs, so mentions of one type in one article are merged when their texts match after
  lower-casing and whitespace/dash normalisation (the "..." joining a discontinuous mention is dropped), or when one directly follows the other in parentheses
  ("long form (SF)" or "SF (long form)"); each merged set gets a synthetic identifier (G1, D1, ...);
- a gene-disease pair is a relation when at least one link between their mentions carries a relation other than
  "no relation" (21 types collapsed); pairs linked only by "no relation" are kept apart as explicit negatives;
- each relation records the BioTriplex location (background / result) of its links and the subsections they are in.
Annotation offsets that do not land on their text are searched for within 5 characters (and trailing spaces
trimmed); the rest are dropped and counted.
Sample: evaluation set = val + test (30 articles), fixed before any run; tuning set = 5 train articles with relations
(random.Random(20261103)), used only to check these rules. Written to sample_biotriplex.json, with each evaluation
article's PMID (NCBI ID converter, cached in data/pmids.json) for the PubTator3 baseline.
Output: papers/BTX<pmcid>.xml, papers/BTXA<pmcid>.xml, gold_biotriplex.json (same gold under both IDs).
"""
import hashlib
import importlib.util
import json
import random
import re
import sys
from pathlib import Path
from xml.sax.saxutils import escape

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
sys.path.insert(0, str(PROJECT))
from biomed.io import request_json  # noqa: E402
from biomed.jats import parse_article  # noqa: E402

_spec = importlib.util.spec_from_file_location("euadr_prepare", HERE.parent / "euadr" / "prepare.py")
_euadr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_euadr)
offset_map = _euadr.offset_map  # difflib alignment of original text to parsed blocks

DATA = HERE / "data" / "splits.json"
PAPERS = HERE / "papers"
SAMPLE = HERE / "sample_biotriplex.json"
SEED, TUNE = 20261103, 5
TYPE = {"genes": "GeneOrGeneProduct", "diseases": "DiseaseOrPhenotypicFeature"}
RANK = [("abstract", 0), ("introduction", 1), ("background", 1), ("method", 2), ("material", 2), ("subject", 2),
        ("result", 3), ("discussion", 4), ("conclusion", 5)]


def clean_name(name):
    return re.sub(r"\s+Rev$", "", name.strip(), flags=re.I)


def order(names):
    def key(name):
        number = re.match(r"(\d+)\.", name)
        if number:
            return (1, int(number.group(1)), name)
        low = name.lower()
        return (1, next((r for w, r in RANK if w in low), 6) * 100, name)
    return sorted(names, key=key)


def lines_of(text):
    """[(start, line)] for the non-empty lines of a subsection text."""
    out, pos = [], 0
    for line in text.split("\n"):
        if line.strip():
            out.append((pos, line))
        pos += len(line) + 1
    return out


def norm(text):
    text = re.sub(r"\s*(\.\.\.|…)\s*", " ", text.lower())  # discontinuous mentions: "drug ... addiction"
    return re.sub(r"\s+", " ", re.sub(r"[‐‑–−-]", "-", text)).strip()


def locate(text, start, end, want):
    if text[start:end] == want:
        return start, end
    stripped = want.strip()
    for shift in sorted(range(-5, 6), key=abs):
        s = start + shift + (len(want) - len(want.lstrip()))
        if text[s:s + len(stripped)] == stripped:
            return s, s + len(stripped)
    return None


def write_jats(paper_id, pmcid, title, sections):
    """sections: [(name, [(start, line)], is_abstract)]. Returns the emitted paragraphs in order: (name, start, line)."""
    PAPERS.mkdir(exist_ok=True)
    emitted, abstract, body = [], "", ""
    for name, lines, is_abstract in sections:
        paragraphs = "".join(f"<p>{escape(line)}</p>" for _, line in lines)
        emitted += [(name, start, line) for start, line in lines]
        if is_abstract:
            abstract += paragraphs
        else:
            body += f"<sec><title>{escape(name)}</title>{paragraphs}</sec>"
    xml = ('<?xml version="1.0" encoding="UTF-8"?>\n<article article-type="research-article"><front><article-meta>'
           f'<article-id pub-id-type="pmc">PMC{pmcid}</article-id>'
           f'<title-group><article-title>{escape(title)}</article-title></title-group>'
           f'<abstract>{abstract}</abstract></article-meta></front><body>{body}</body></article>\n')
    (PAPERS / f"{paper_id}.xml").write_text(xml, encoding="utf-8")
    return emitted


def place(paper_id, emitted):
    """{(section, original line start): (block_id, offset map, line)} for the paragraphs of a written article."""
    blocks = [b for b in parse_article(PAPERS / f"{paper_id}.xml")["blocks"] if b["content_type"] == "paragraph"]
    if len(blocks) != len(emitted):
        raise SystemExit(f"{paper_id}: {len(emitted)} paragraphs written, {len(blocks)} parsed")
    return {(name, start): (b["block_id"], offset_map(line, b["text"]), line)
            for (name, start, line), b in zip(emitted, blocks)}


def pmid_map(pmcids):
    """{pmcid digits: PMID} from the NCBI ID converter, cached in data/pmids.json."""
    cache = HERE / "data" / "pmids.json"
    known = json.loads(cache.read_text(encoding="utf-8")) if cache.exists() else {}
    missing = [p for p in pmcids if p not in known]
    if missing:
        url = ("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?format=json&tool=biotriplex_prepare&ids="
               + ",".join("PMC" + p for p in missing))
        for r in request_json(url)["records"]:
            known[r["pmcid"].removeprefix("PMC")] = r["pmid"]
        cache.write_text(json.dumps(known, indent=1), encoding="utf-8")
    return known


class Merge:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def build(pmcid, record):
    data = record["section_data"]
    names = {clean_name(n): n for n in data}
    abstract_name = next(n for n in names if n.lower() == "abstract")
    ordered = [abstract_name] + [n for n in order(names) if n != abstract_name]
    title = ""  # the title is already the abstract's first line ("Title: ..."); repeating it would duplicate mentions
    sections = [(n, lines_of(data[names[n]]["text"]), n == abstract_name) for n in ordered]
    full = place(f"BTX{pmcid}", write_jats(f"BTX{pmcid}", pmcid, title, sections))
    place(f"BTXA{pmcid}", write_jats(f"BTXA{pmcid}", pmcid, title, sections[:1]))  # same abstract block IDs

    stats = {"repaired": 0, "dropped": 0}
    mentions, merge, by_key = [], Merge(), {}
    for name in ordered:
        x, starts = data[names[name]], [s for s, _ in lines_of(data[names[name]]["text"])]
        for kind in TYPE:
            for mid, m in x[kind].items():
                key = (name, kind, mid)
                for seg in re.split(r"[,;]", m["spans"]):
                    a, b = map(int, seg.split("~"))
                    want = m["text"] if "," not in m["spans"] else x["text"][a:b]
                    found = locate(x["text"], a, b, want)
                    if not found:
                        stats["dropped"] += 1
                        continue
                    stats["repaired"] += found != (a, b)
                    a, b = found
                    line_start = max((s for s in starts if s <= a), default=None)
                    block_id, mapping, line = full[(name, line_start)]
                    if b - line_start > len(line):
                        stats["dropped"] += 1  # crosses a line break
                        continue
                    s, e = mapping[a - line_start], mapping[b - line_start]
                    mentions.append({"key": key, "kind": kind, "type": TYPE[kind], "text": m["text"].strip(),
                                     "block_id": block_id, "start": s, "end": e, "section": name, "orig": (a, b)})
                merge.find(key)
                by_key[key] = m
    # Merge by normalised text, then "A (B)" / "B (A)" pairs of the same type in the same subsection.
    first_by_text = {}
    for m in mentions:
        t = (m["kind"], norm(m["text"]))
        merge.union(m["key"], first_by_text.setdefault(t, m["key"]))
    for name in ordered:
        text = data[names[name]]["text"]
        spans = sorted(((m["orig"], m) for m in mentions if m["section"] == name), key=lambda x: x[0])
        for (a, b), m in spans:
            for (c, d), n in spans:
                if n["kind"] == m["kind"] and n["key"] != m["key"] and 0 < c - b <= 3 \
                        and re.fullmatch(r"\s*\(\s*", text[b:c]) and text[d:d + 2].lstrip().startswith(")"):
                    merge.union(m["key"], n["key"])
    ids, counters = {}, {"genes": 0, "diseases": 0}
    for m in mentions:
        root = merge.find(m["key"])
        if root not in ids:
            counters[m["kind"]] += 1
            ids[root] = f"{'G' if m['kind'] == 'genes' else 'D'}{counters[m['kind']]}"
        m["identifier"] = ids[root]
    cluster = {k: ids[merge.find(k)] for k in by_key if merge.find(k) in ids}

    pairs = {}
    for name in ordered:
        x = data[names[name]]
        for link in x["entity_linkings"].values():
            g, d = cluster.get((name, "genes", link["geneID"])), cluster.get((name, "diseases", link["diseaseID"]))
            relation = x["relations"].get(link["relationID"])
            if not (g and d and relation):
                stats["dropped"] += 1
                continue
            p = pairs.setdefault((g, d), {"types": set(), "locations": set(), "sections": set()})
            p["types"].add(relation["relation"])
            p["locations"].add(relation.get("location", ""))
            p["sections"].add(name)
    relations, negatives = [], []
    for (g, d), p in sorted(pairs.items()):
        row = {"gene": g, "disease": d, "type": "relation", "novel": "", "relation_types": sorted(p["types"]),
               "locations": sorted(p["locations"]), "sections": sorted(p["sections"]),
               "abstract_only": p["sections"] == {abstract_name}, "body_only": abstract_name not in p["sections"]}
        (negatives if p["types"] == {"no relation"} else relations).append(row)
    blocks = sorted({b for b, _, _ in full.values()})
    gold = {"blocks": blocks, "mentions": [{k: m[k] for k in ("identifier", "type", "text", "block_id", "start", "end")}
                                           for m in mentions],
            "relations": relations, "explicit_negatives": negatives,
            "abstract_blocks": sorted(b for (n, _), (b, _, _) in full.items() if n == abstract_name)}
    texts = {}
    for m in mentions:
        texts.setdefault(m["identifier"], set()).add(m["text"])
    gold["entity_texts"] = {k: sorted(v) for k, v in sorted(texts.items())}
    return gold, stats, sum(len(x["text"]) for x in data.values())


def main():
    records = json.loads(DATA.read_text(encoding="utf-8"))
    split = {p: r["split"] for p, r in records.items()}
    evaluation = sorted(p for p in records if split[p] in ("val", "test"))
    with_relations = sorted(p for p in records if split[p] == "train"
                            and any(x["entity_linkings"] for x in records[p]["section_data"].values()))
    tune = sorted(random.Random(SEED).sample(with_relations, TUNE))
    sample = {"source": "BioTriplex, Preprocessed BioTriplex/splits.json", "commit": (HERE / "data" / "SOURCE_COMMIT").read_text().strip(),
              "data_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(), "seed": SEED,
              "evaluation": [f"BTX{p}" for p in evaluation], "evaluation_abstract": [f"BTXA{p}" for p in evaluation],
              "tuning": [f"BTX{p}" for p in tune], "split": {f"BTX{p}": split[p] for p in evaluation + tune}}
    if SAMPLE.exists():
        frozen = json.loads(SAMPLE.read_text(encoding="utf-8"))
        if {k: frozen[k] for k in ("data_sha256", "evaluation", "tuning")} != {k: sample[k] for k in ("data_sha256", "evaluation", "tuning")}:
            raise SystemExit("sample_biotriplex.json no longer matches the data and seed; refusing to redraw")
    # PMIDs for the PubTator3 baseline (experiments/biored/pubtator_baseline.py reads "papers").
    pmids = pmid_map(evaluation + tune)
    sample["papers"] = [{"paper_id": f"BTX{p}", "pmid": pmids[p]} for p in evaluation]
    SAMPLE.write_text(json.dumps(sample, indent=2), encoding="utf-8")
    gold, totals, report = {}, {"repaired": 0, "dropped": 0}, []
    for p in evaluation + tune:
        g, stats, chars = build(p, records[p])
        gold[f"BTX{p}"] = gold[f"BTXA{p}"] = g
        for k in totals:
            totals[k] += stats[k]
        report.append((f"BTX{p}", split[p], chars, len(g["relations"]), sum(r["body_only"] for r in g["relations"]),
                       len(g["explicit_negatives"])))
    (HERE / "gold_biotriplex.json").write_text(json.dumps(gold, ensure_ascii=False, indent=1), encoding="utf-8")
    for subset in ("evaluation", "tuning"):
        rows = [r for r in report if r[0] in sample[subset]]
        print(f"{subset}: {len(rows)} articles, {sum(r[3] > 0 for r in rows)} with relations; relations {sum(r[3] for r in rows)} "
              f"(only outside the abstract {sum(r[4] for r in rows)}); explicit negatives {sum(r[5] for r in rows)}; "
              f"characters median {sorted(r[2] for r in rows)[len(rows) // 2]}")
    print(f"annotation offsets repaired {totals['repaired']}, dropped (mentions or links) {totals['dropped']}")


if __name__ == "__main__":
    main()
