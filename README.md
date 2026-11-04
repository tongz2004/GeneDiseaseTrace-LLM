<p align="right"><b>English</b> | <a href="README.zh-CN.md">简体中文</a></p>

# GeneDiseaseTrace-LLM

Gene–disease relation extraction from biomedical papers with an LLM, with sentence citations and a full evaluation
on abstracts and full texts. Each reported relation cites the sentences behind it and is labeled as the paper's own
finding or as background knowledge; pairs with weak evidence are flagged *uncertain*. A deterministic check confirms that
every citation resolves; whether a citation supports the relation still needs human review.

**Project page (results, method and interactive graphs):** https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/

## The task
For each of 100 PMC full-text papers (JATS XML): extract all genes and diseases, decide for every gene × disease
pair whether the paper gives evidence of an association, correlation or causal relationship, and build one
gene–disease bipartite graph per paper.

## Method
1. **Packing.** A script splits the paper into numbered blocks and sentences with stable IDs.
2. **One request per paper.** A 108-word task prompt and a 595-word output format go to the LLM in a single tool-less
   request. The model answers with short lines that cite sentence IDs instead of copying text. This minimal prompt is the default system.
3. **Deterministic post-processing.** A script rebuilds quotes and character offsets, completes the gene × disease
   matrix and writes entity, mention, evidence and pair tables plus GraphML. A validator checks the surface forms and
   citations; only if it finds problems is the error list sent back once for repair.

**Earlier design.** The project began with decision rules (R1–R16) and 35 boundary cases (`skill_snapshot/`), condensed
into a 1,045-word condensed rule set. The evaluation below replaced it with the minimal prompt as the default, because it lowered
recall on BioRED.

**Names.** These are parallel conditions, not successive versions. Folder names under `experiments/single_shot/variants/` keep their original IDs so that run records stay valid:

| Name | Folder |
|---|---|
| Minimal prompt | `v6_minimal_lean_medium` |
| Minimal prompt + format rules | `v6f_format_rules_medium` |
| Condensed rule set | `v1_card_lean_medium` |
| Full rule set | `v2_full_lean_medium` |

## Results
All systems use one-to-one matching, so one gold relation can credit only one prediction. Brackets are 95% document-level
bootstrap intervals. The [project page](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/) has the full tables.

| Evaluation | Result |
|---|---|
| BioRED Test, 100 abstracts | Minimal prompt: precision 0.680, recall 0.560, F1 0.614. PubTator3: F1 0.680; paired difference −0.064 [−0.151, 0.033]. PubTator3's relation model was trained on data including BioRED. |
| Condensed rule set ablation, BioRED Dev+Test | The rule set lowered recall (0.418 vs 0.604) and F1 (0.551 vs 0.665). This describes one rule set, which changed many rules and the length at once. |
| BioTriplex, locked set of 30 full texts | Claude Opus 5.5 recall rose from 0.366 on abstracts to 0.813 on full text (25 valid articles); F1 0.525 against 0.293 for BioREx and 0.361 for RENET2 on all 30. GPT-6.1 Sol: recall 0.858, F1 0.433 (27 valid). |
| BioTriplex, 49 further full texts | The prespecified comparison, with failures counted as empty, gave a recall gain of +0.322 [0.133, 0.477]. On the 35 valid articles: recall 0.829, F1 0.474; F1 exceeded RENET2's by 0.262 [0.148, 0.365]. These articles were chosen after the locked results were known, so they are further validation and not a held-out test. |
| Citation check, 30 pairs per corpus | Cited sentences alone were judged sufficient for 58.5% [45.8, 70.9] of supported pairs on BioRED and 75.4% [58.0, 89.1] on BioTriplex; 5.6% on each were judged unsupported. One annotator. |
| Task 1: gene and disease extraction, BioRED Test | Document-concept F1: genes 0.968 (PubTator3 0.917), diseases 0.874 (PubTator3 0.958). PubTator3's recognizer was trained on data including BioRED. |
| Task 2: the 100 original papers | The minimal prompt completed 81 papers with 1,180 supported relations; the condensed rules reported 881 on the same papers. The other 19 were stopped by the safety filter and keep their labeled condensed-rules result. No gold labels, so these are output counts only. |
| Prompts and models, 100 papers | Format rules cut repairs from 62 to 9 papers and cost from $0.64 to $0.43 per paper, with 6% fewer supported relations; their recall is not yet measured. Opus 4.8 had no safety stops but 42% fewer supported relations than Opus 5.5. See the [comparison](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/#compare). |
| Cost | At list prices, including repair requests: about $0.40 per BioTriplex full text, $0.64 per original paper and $0.07 per abstract. |

**Limitations.**
- Safety filters left 37 of 179 full texts without a Claude result (4 of 30, 14 of 49 and 19 of 100). We do not know what triggers them. Valid-article results apply only to articles that were not blocked. With failures counted as empty, recall on the 49 articles is 0.632, against 0.848 for BioREx, which runs locally.
- All blinded checks used one annotator who knew the hypotheses; there is no second annotator.
- BioTriplex annotation is not exhaustive, so measured precision can understate textual correctness.
- F1 comparisons depend on entity scope and training text, so they should not be read as better relation judgment.
- Entities are not mapped to database IDs (NCBI Gene / MeSH), and the model may have seen these texts in pretraining.

Graphs for the 100 papers, for both the minimal prompt and the condensed rules, are [here](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/graphs/). The [October result pages](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/results/ab/ab_summary.html) used the condensed rules and an earlier scorer.

## Repository layout
| Path | Contents |
|---|---|
| `experiments/single_shot/` | **The method.** Card, output format, packer, runner, expander, batch runner, frozen A/B variants, PubTator3 comparison and page renderers; `runs/` holds every run's outputs. Notes (Chinese): `README.md` |
| `experiments/biored/` | **The benchmark.** BioRED adapter (`prepare.py`), frozen samples and gold, scorer (`evaluate.py`), PubTator3 baseline, results page. Notes (Chinese): `README.md` |
| `experiments/biotriplex/` | Full-text benchmark: adapter, locked set, second batch protocol, scorer, paired comparison, blinded checks |
| `experiments/citation_check/` | Citation sufficiency: automatic count, blinded sheets, scorer |
| `experiments/entity_eval/` | Task 1: gene and disease extraction scores from stored outputs |
| `experiments/safety_blocks/` | Description of the full texts stopped by the safety filter |
| `experiments/biorex/`, `experiments/renet2_model/` | Local runs of the supervised comparators (third-party code and weights are downloaded, not stored) |
| `experiments/euadr/`, `experiments/renet2/` | EU-ADR and RENET2 benchmark adapters |
| `paper/` | Manuscript in preparation (LaTeX) and decision memo |
| `skill_snapshot/` | The decision rules: skill 1.3.0 with data contract and boundary cases (exported from a separate repository) |
| `toolkit/` | Contract validator and shared helpers (also the earlier agent-mode toolkit) |
| `biomed/` | JATS parser and PubTator3 client shared by everything |
| `baseline/` | The first, LLM-free baseline: PubTator3 annotations aligned to the full text (`run.ps1`, config, audit and figure scripts) |
| `site_builder/` | Builds the public project page (`build_site.py`; bilingual main page template, graph viewers) |
| `tests/` | Unit tests for the parser, pipeline and baseline |
| `archive/` | Development history: manual annotation of paper 1, blind tests, rule proposals, one-off scripts |

## Reproducing
Python 3.12; the Claude CLI must be installed and signed in for model runs. Benchmark texts and gold are not included.

```bash
# BioRED: put BIORED.zip from https://ftp.ncbi.nlm.nih.gov/pub/lu/BioRED/ in experiments/biored/data/ and unzip it
python experiments/biored/prepare.py --all
python experiments/single_shot/batch.py --variant v6_minimal_lean_medium --run-set biored_all --papers <ids from sample_all.json>
python experiments/biored/evaluate.py v6_minimal_lean_medium biored_all

# BioTriplex: put splits.json in experiments/biotriplex/data/ (see experiments/biotriplex/README.md)
python experiments/biotriplex/prepare.py
python experiments/single_shot/batch.py --variant v6_minimal_lean_medium --run-set biotriplex --papers <ids from sample_biotriplex.json>
python experiments/biotriplex/evaluate.py v6_minimal_lean_medium biotriplex
python experiments/biotriplex/compare.py v6_minimal_lean_medium v6_minimal_lean_medium:abstract

# Second batch: protocol in experiments/biotriplex/batch2_protocol.md
python experiments/biotriplex/prepare.py --batch2

# Tests and the public page
python -m unittest discover -s tests
python site_builder/build_site.py
```

The 100-paper runs use `--variant v6_minimal_lean_medium --all --run-set full100` (October: `v1_card_lean_medium`; Opus 4.8: add `--model Claude_Opus_4.8 --run-set full100_opus48`), with inputs in `.Input/100papers/`.
The PubTator3 comparison reads the baseline's aligned annotations in `outputs/`; build them with
`.\baseline\run.ps1 -Command run` (see `baseline/README.md`).

## Data and credits
BioRED: Luo et al., *Briefings in Bioinformatics* 2022 (public domain, NCBI). BioTriplex: Collins et al.,
*Bioinformatics* 2026. PubTator3: NCBI. BioREx: Lai et al. 2023. RENET2: Su et al. 2021. Input papers, benchmark texts and
full-text prompts are not redistributed here.
