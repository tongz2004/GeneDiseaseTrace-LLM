"""One annotator's check of the disputed pairs: those V6 supported, V1 did not, and BioRED does not annotate.

Usage: python experiments/biored/self_check.py                       (build the sheet)
       python experiments/biored/self_check.py --score [ANSWERS_CSV]   (after the sheet is filled in)

sheet.html takes the answers by clicking, keeps them in the browser and exports sheet_answers.csv; --score reads
that export (from self_check/ or Downloads) by default, else a filled-in sheet.csv.

The 16 disputed items (adjudication stratum b_only, matches_gold_relation = no) are mixed with controls so the
annotator cannot tell them apart:
- 5 positives: package items that match a BioRED gold relation;
- 5 negatives: a gene and a disease mentioned in the same Dev abstract, neither of which is in any BioRED relation
  there, and not supported by either variant.
Items get new shuffled IDs (U01 ...). The controls check the annotator against BioRED before the disputed items
are read. Writes <package>/self_check/ with sheet.{csv,html} and key.csv (keep the key closed until scoring).
"""
import argparse
import csv
import html
import json
import random
from collections import Counter
from pathlib import Path

from adjudication import ROUND1, marked_text, parse_article, same, system_pairs, write_round
from evaluate import RUNS
from llm_adjudicate import PACKAGE, STATES, rows

HERE = PACKAGE.parents[1]
OUT = PACKAGE / "self_check"
VARIANTS = ("v1_card_lean_medium", "v6_minimal_lean_medium")
SEED = 20261104


def negatives(papers, gold_all, rng, n):
    candidates = []
    for paper in papers:
        mentions = gold_all[paper]["mentions"]
        related = {(r["gene"], r["disease"]) for r in gold_all[paper]["relations"]}
        concepts = lambda t: sorted({c for m in mentions if m["type"] == t for c in m["identifier"].split(",") if c not in ("", "-")})
        spans = lambda c: [(m["block_id"], m["start"], m["end"]) for m in mentions if c in m["identifier"].split(",")]
        name = lambda c: next(m["text"] for m in mentions if c in m["identifier"].split(","))
        supported = [p for v in VARIANTS for p in system_pairs(RUNS / "biored_all" / v / paper / "outputs")]
        for g in concepts("GeneOrGeneProduct"):
            for d in concepts("DiseaseOrPhenotypicFeature"):
                pair = {"paper": paper, "gene": name(g), "disease": name(d), "gene_spans": spans(g), "disease_spans": spans(d)}
                # Neither end in any BioRED relation of the abstract: a pair that merely lacks a relation can still
                # be a related pair under another concept ID (a gene named by its protein family, say).
                if g not in {r[0] for r in related} and d not in {r[1] for r in related}                         and not any(same(pair, s) for s in supported):
                    candidates.append(pair)
    return rng.sample(candidates, n)


def build(html_of):
    rng = random.Random(SEED)
    key = {r["item_id"]: r for r in rows("key_coordinator.csv")}
    sheet = {r["item_id"]: r for r in rows("round1_relation.csv")}
    disputed = [i for i, k in key.items() if k["stratum"] == "b_only" and k["matches_gold_relation"] == "no"]
    positives = rng.sample(sorted(i for i, k in key.items() if k["matches_gold_relation"] == "yes"), 5)
    gold_all = json.loads((HERE / "gold_all.json").read_text(encoding="utf-8"))
    papers = sorted({k["paper_id"] for k in key.values()})
    items = [{"kind": "disputed", "source": i} for i in disputed] + [{"kind": "positive", "source": i} for i in positives]
    items += [{"kind": "negative", "source": "", **p} for p in negatives(papers, gold_all, rng, 5)]
    rng.shuffle(items)
    out_rows, key_rows = [], []
    for n, item in enumerate(items, 1):
        uid = f"U{n:02d}"
        if item["kind"] == "negative":
            blocks = parse_article(HERE / "papers" / f"{item['paper']}.xml")["blocks"]
            text = marked_text(blocks, item["gene_spans"], item["disease_spans"])
            gene, disease, paper = item["gene"], item["disease"], item["paper"]
        else:
            r = sheet[item["source"]]
            gene, disease, paper = r["gene"], r["disease"], key[item["source"]]["paper_id"]
            text = html_of[item["source"]]
        out_rows.append({"item_id": uid, "gene": gene, "disease": disease, "abstract": text, "abstract_html": text})
        key_rows.append({"item_id": uid, "kind": item["kind"], "package_item": item["source"], "paper_id": paper})
    OUT.mkdir(exist_ok=True)
    if not (OUT / "sheet.csv").exists():  # never overwrite a sheet that may hold answers
        write_round(OUT / "sheet", out_rows, ["item_id", "gene", "disease"], ROUND1,
                    "Self-check: what does the abstract state about the marked gene and disease?")
    write_form(OUT / "sheet.html", out_rows)
    with (OUT / "key.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(key_rows[0]))
        writer.writeheader()
        writer.writerows(key_rows)
    print(f"{len(out_rows)} items -> {OUT}; {Counter(k['kind'] for k in key_rows)}")


FORM = """<!doctype html><meta charset="utf-8"><title>Self-check annotation</title>
<style>
body{font:15px/1.55 system-ui,sans-serif;max-width:880px;margin:0 auto;padding:0 16px 4em;background:#fff;color:#222}
header{position:sticky;top:0;background:#fff;border-bottom:1px solid #ccc;padding:.6em 0;z-index:1;display:flex;gap:1em;align-items:center;flex-wrap:wrap}
section{border-top:1px solid #ddd;padding:.8em 0}section.done h2{color:#2a7a2a}
h2{font-size:1.05em;margin:.2em 0}b.g{background:#cfe3ff}b.d{background:#ffd9cc}
fieldset{border:0;margin:.4em 0;padding:0}legend{font-weight:600;font-size:.92em;margin-bottom:.15em}
label{display:inline-block;margin:0 .9em .2em 0;cursor:pointer}input[type=text]{width:100%;box-sizing:border-box;padding:.3em}
button{font:inherit;padding:.35em .9em;cursor:pointer}.abs{background:#fafafa;padding:.2em .8em;border-radius:4px}
</style>
<header><b>自查标注</b><span id="progress"></span><button id="export">导出 CSV</button>
<span style="font-size:.85em;color:#555">答案自动保存在本浏览器；全部填完后点"导出 CSV"。</span></header>
<p>基因标为<b class="g">蓝色</b>，疾病标为<b class="d">橙色</b>。只根据__SCOPE__写了什么来判断，不要查资料，也不要凭生物学常识。</p>
<ul style="font-size:.92em">
<li><b>relation</b>：明确陈述了关联、相关或因果（包括负相关、治疗或保护作用）</li>
<li><b>hedged relation</b>：只用了推测性措辞（may、might、suggest、possible 等）</li>
<li><b>explicitly no relation</b>：明确说没有关联</li>
<li><b>no statement</b>：只是同时出现，或没有陈述两者的关系</li></ul>
<div id="items"></div>
<script>
const ITEMS = __ITEMS__;
const QUESTIONS = __QUESTIONS__;
const STORE = "__STORE__";
let answers = {};
try { answers = JSON.parse(localStorage.getItem(STORE) || "{}"); } catch (e) {}
const save = () => { try { localStorage.setItem(STORE, JSON.stringify(answers)); } catch (e) {} };
const box = document.getElementById("items");
for (const it of ITEMS) {
  const a = answers[it.item_id] = answers[it.item_id] || {};
  const s = document.createElement("section");
  s.innerHTML = `<h2>${it.item_id}：${it.gene} × ${it.disease}</h2><div class="abs">${it.abstract_html}</div>`;
  for (const q of QUESTIONS) {
    const f = document.createElement("fieldset");
    f.innerHTML = `<legend>${q.title}</legend>` + q.options.map(o =>
      `<label><input type="radio" name="${it.item_id}-${q.column}" value="${o}"${a[q.column] === o ? " checked" : ""}> ${o}</label>`).join("");
    f.addEventListener("change", e => { a[q.column] = e.target.value; save(); update(); });
    s.appendChild(f);
  }
  const n = document.createElement("input");
  n.type = "text"; n.placeholder = "notes（可留空）"; n.value = a.notes || "";
  n.addEventListener("input", () => { a.notes = n.value; save(); });
  s.appendChild(n);
  s.dataset.id = it.item_id;
  box.appendChild(s);
}
function update() {
  let done = 0;
  for (const s of box.children) {
    const ok = QUESTIONS.every(q => answers[s.dataset.id][q.column]);
    s.classList.toggle("done", ok); done += ok;
  }
  document.getElementById("progress").textContent = `已完成 ${done} / ${ITEMS.length}`;
}
update();
document.getElementById("export").onclick = () => {
  const cols = ["item_id", "gene", "disease", ...QUESTIONS.map(q => q.column), "notes"];
  const esc = v => `"${String(v ?? "").replace(/"/g, '""')}"`;
  const lines = [cols.map(esc).join(",")].concat(ITEMS.map(it => {
    const a = answers[it.item_id];
    return [it.item_id, it.gene, it.disease, ...QUESTIONS.map(q => a[q.column] || ""), a.notes || ""].map(esc).join(",");
  }));
  const blob = new Blob(["\\ufeff" + lines.join("\\r\\n")], {type: "text/csv"});
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob); link.download = "__DOWNLOAD__"; link.click();
};
</script>
"""


def write_form(path, items, scope="摘要", store="self_check_answers_v2", download="sheet_answers.csv"):
    """The sheet as a page with clickable answers, saved in the browser (under `store`) and exported as `download`.
    scope names the text shown ("摘要" here; other sheets show paragraphs)."""
    titles = ["实体标得对吗？", f"{scope}怎么描述两者的关系？", "这个陈述来自哪里？"]
    questions = [{"column": c, "title": t, "options": [o.strip() for o in c.split("(", 1)[1].rstrip(")").split("/")]}
                 for c, t in zip(ROUND1[:3], titles)]
    data = [{"item_id": r["item_id"], "gene": html.escape(r["gene"]), "disease": html.escape(r["disease"]),
             "abstract_html": r["abstract_html"]} for r in items]
    page = FORM.replace("__ITEMS__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/"))
    page = page.replace("__SCOPE__", scope).replace("__STORE__", store).replace("__DOWNLOAD__", download)
    path.write_text(page.replace("__QUESTIONS__", json.dumps(questions, ensure_ascii=False)), encoding="utf-8")


def package_html():
    """{item_id: marked abstract HTML} parsed back out of round1_relation.html."""
    page = (PACKAGE / "round1_relation.html").read_text(encoding="utf-8")
    out = {}
    for section in page.split("<section>")[1:]:
        head, body = section.split("</h2>", 1)
        out[head.removeprefix("<h2>").split(":", 1)[0]] = body.removesuffix("</section>").split("</section>")[0]
    return out


def score(path):
    """path: the filled sheet; by default the page's export (sheet_answers.csv here or in Downloads), else sheet.csv."""
    if not path:
        path = next(p for p in (OUT / "sheet_answers.csv", Path.home() / "Downloads" / "sheet_answers.csv", OUT / "sheet.csv")
                    if p.exists())
    print(f"answers: {path}")
    key = {r["item_id"]: r for r in csv.DictReader((OUT / "key.csv").open(encoding="utf-8-sig"))}
    answers = list(csv.DictReader(Path(path).open(encoding="utf-8-sig")))
    field = ROUND1[1]
    report = {}
    for a in answers:
        state = a[field].split(":")[0].strip().lower()
        if state not in STATES:
            print(f"{a['item_id']}: unreadable answer {a[field]!r}")
            continue
        entry = report.setdefault(key[a["item_id"]]["kind"], {"items": 0, "states": Counter(), "entities_wrong": 0})
        entry["items"] += 1
        entry["states"][state] += 1
        entry["entities_wrong"] += a[ROUND1[0]].strip().lower() not in ("", "yes")
    for kind, e in report.items():
        e["relation_or_hedged"] = e["states"]["relation"] + e["states"]["hedged relation"]
        e["states"] = dict(e["states"])
    (OUT / "result.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--score", nargs="?", const="", metavar="ANSWERS_CSV")
    args = parser.parse_args()
    if args.score is not None:
        score(args.score)
    else:
        build(package_html())
