<p align="right"><a href="README.md">English</a> | <b>简体中文</b></p>

# GeneDiseaseTrace-LLM

用 LLM 从生物医学论文中抽取基因–疾病关系，附句子引用，并在摘要和全文上做了完整评测。
每条报告的关系都引用支撑它的句子，并标明是论文自身的发现还是背景知识；证据不足的配对标为 *uncertain*。确定性检查可以确认每条引用都能解析，但引用是否真的支持这条关系，仍需要人工复核。

**项目主页（结果、方法和交互图）：** https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/?lang=zh

## 任务
对 100 篇 PMC 全文论文（JATS XML）逐篇处理：抽取所有基因和疾病；对每一个基因 × 疾病配对，判断论文是否提供了关联、相关或因果关系的证据；最后为每篇论文构建一张基因–疾病二部图。

## 方法
1. **打包。** 脚本把论文切分成带固定编号的段落和句子。
2. **每篇论文一次请求。** 108 词的任务提示词和 595 词的输出格式一起发给 LLM，不调用任何工具。模型用简短的行作答，引用句子编号而不复制原文。这个最小提示词（V6）是默认系统。
3. **确定性的后处理。** 脚本重建引文和字符位置，补全基因 × 疾病矩阵，并写出实体、提及、证据和配对表以及 GraphML。校验器核对实体文字和引用；只有发现问题时，才把错误清单发回模型修正一次。

**早期设计。** 项目最初制定了判定规则（R1–R16）和 35 个边界案例（`skill_snapshot/`），并压缩成 1,045 词的判定卡（V1）。下面的评测把默认系统换成了 V6，因为判定卡在 BioRED 上降低了召回率。

## 结果
所有系统都用一对一匹配，一个金标关系只能给一个预测计分。方括号是文章层 bootstrap 的 95% 区间。完整表格见[项目主页](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/?lang=zh)。

| 评测 | 结果 |
|---|---|
| BioRED Test，100 篇摘要 | 最小提示词：精确率 0.680、召回率 0.560、F1 0.614。PubTator3 的 F1 是 0.680；配对差 −0.064 [−0.151, 0.033]。PubTator3 的关系模型训练数据包含 BioRED。 |
| 判定卡消融，BioRED Dev+Test | 判定卡降低了召回率（0.418 对 0.604）和 F1（0.551 对 0.665）。这只说明这一张卡，它同时改变了很多规则和长度。 |
| BioTriplex 锁定集，30 篇全文 | Claude Opus 5.5 的召回率从摘要的 0.366 升到全文的 0.813（25 篇有效）；F1 0.525，全部 30 篇上 BioREx 是 0.293，RENET2 是 0.361。GPT-6.1 Sol：召回率 0.858，F1 0.433（27 篇有效）。 |
| BioTriplex 新增 49 篇全文 | 事先定好的比较（失败记为空答案）显示召回率提高 +0.322 [0.133, 0.477]。35 篇有效文章上：召回率 0.829，F1 0.474，F1 比 RENET2 高 0.262 [0.148, 0.365]。这批文章是看过锁定集结果后才确定的，只算进一步验证，不算独立测试。 |
| 引用检查，每个语料 30 条 | 引用句单独足以支持的比例，BioRED 是 58.5% [45.8, 70.9]，BioTriplex 是 75.4% [58.0, 89.1]；两者都有 5.6% 被判为不支持。只有一位标注者。 |
| 任务 1：基因和疾病抽取，BioRED Test | 文档概念层 F1：基因 0.968（PubTator3 0.917），疾病 0.874（PubTator3 0.958）。PubTator3 的实体模型训练数据包含 BioRED。 |
| 任务 2：最初的 100 篇论文 | V6 完成 81 篇，报告 1,180 条 supported 关系；V1 在同样论文上是 881 条。另外 19 篇被安全机制中止，保留标注过的 V1 结果。这些论文没有金标，所以只是输出数量。 |
| V1 与 V6 | V6 在三个关系基准上召回率都更高，F1 在 BioRED 和 EU-ADR 上更高。V1 精确率相同或更高，费用更低，格式错误更少。见[对比](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/?lang=zh#v1v6)。 |
| 费用 | 按标价，含修复请求：BioTriplex 全文每篇约 $0.40，原始论文每篇约 $0.64，摘要每篇约 $0.07。 |

**局限性。**
- 安全机制让 179 篇全文中的 37 篇没有 Claude 结果（分别是 4/30、14/49 和 19/100）。我们不清楚触发原因。有效文章上的结果只适用于没被拦截的文章。失败记为空答案时，49 篇的召回率是 0.632，本地运行的 BioREx 是 0.848。
- 所有盲法检查都只有一位标注者，且了解研究假设；没有第二位标注者。
- BioTriplex 的标注不完整，测得的精确率可能低估文本上的正确性。
- F1 比较受实体范围和训练文本影响，不能理解为关系判断能力更强。
- 实体没有映射到数据库 ID（NCBI Gene / MeSH），也无法排除模型在预训练时见过这些文本。

100 篇论文的 V6 和 V1 二部图都在[这里](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/graphs/?lang=zh)。[十月的结果页](https://suzuka2788.github.io/GeneDiseaseTrace-LLM-site/results/ab/ab_summary.html)用的是判定卡和较早的评分方式。

## 仓库结构
| 路径 | 内容 |
|---|---|
| `experiments/single_shot/` | **方法本体。** 判定卡、输出格式、打包、运行、展开脚本、批量运行、冻结的 A/B 变体、与 PubTator3 的对比和页面生成；`runs/` 保存每次运行的输出。说明见 `README.md` |
| `experiments/biored/` | **基准评测。** BioRED 适配（`prepare.py`）、冻结的样本和标准答案、评分脚本（`evaluate.py`）、PubTator3 基线、结果页。说明见 `README.md` |
| `experiments/biotriplex/` | 全文基准：适配、锁定集、第二批规格、评分、配对比较、盲法检查 |
| `experiments/citation_check/` | 引用充分性：自动统计、盲标表单、评分 |
| `experiments/entity_eval/` | 任务 1：基于已有输出的基因和疾病抽取评测 |
| `experiments/safety_blocks/` | 被安全机制中止的全文的描述统计 |
| `experiments/biorex/`、`experiments/renet2_model/` | 有监督对照的本地运行（第三方代码和权重需下载，不存放在仓库） |
| `experiments/euadr/`、`experiments/renet2/` | EU-ADR 和 RENET2 基准适配 |
| `paper/` | 准备中的论文（LaTeX）和决策备忘 |
| `skill_snapshot/` | 判定规则：skill 1.3.0，含数据契约和边界案例（从另一个独立仓库导出） |
| `toolkit/` | 契约校验器和共用工具（也是早期 agent 模式的工具包） |
| `biomed/` | 各部分共用的 JATS 解析器和 PubTator3 客户端 |
| `baseline/` | 最早不用 LLM 的基线：把 PubTator3 标注对齐到全文（`run.ps1`、配置、审计和画图脚本） |
| `site_builder/` | 生成公开项目主页（`build_site.py`、中英双语主页模板、二部图页面） |
| `tests/` | 解析器、流程和基线的单元测试 |
| `archive/` | 开发过程记录：第一篇论文的人工标注、盲测、规则修订提案、一次性脚本 |

## 复现
需要 Python 3.12；跑模型需要先安装并登录 Claude CLI。仓库不包含基准文本和金标。

```bash
# BioRED：从 https://ftp.ncbi.nlm.nih.gov/pub/lu/BioRED/ 下载 BIORED.zip，放到 experiments/biored/data/ 并解压
python experiments/biored/prepare.py --all
python experiments/single_shot/batch.py --variant v6_minimal_lean_medium --run-set biored_all --papers <sample_all.json 中的编号>
python experiments/biored/evaluate.py v6_minimal_lean_medium biored_all

# BioTriplex：把 splits.json 放到 experiments/biotriplex/data/（见 experiments/biotriplex/README.md）
python experiments/biotriplex/prepare.py
python experiments/single_shot/batch.py --variant v6_minimal_lean_medium --run-set biotriplex --papers <sample_biotriplex.json 中的编号>
python experiments/biotriplex/evaluate.py v6_minimal_lean_medium biotriplex
python experiments/biotriplex/compare.py v6_minimal_lean_medium v6_minimal_lean_medium:abstract

# 第二批：规格见 experiments/biotriplex/batch2_protocol.md
python experiments/biotriplex/prepare.py --batch2

# 单元测试和公开主页
python -m unittest discover -s tests
python site_builder/build_site.py
```

100 篇运行用 `--variant v6_minimal_lean_medium --all --run-set full100`（十月用 `v1_card_lean_medium`），输入放在 `.Input/100papers/`。
与 PubTator3 的对比读取基线在 `outputs/` 中对齐好的标注；用 `.\baseline\run.ps1 -Command run` 生成（见 `baseline/README.md`）。

## 数据与致谢
BioRED：Luo et al., *Briefings in Bioinformatics* 2022（公有领域，NCBI）。BioTriplex：Collins et al., *Bioinformatics* 2026。PubTator3：NCBI。BioREx：Lai et al. 2023。RENET2：Su et al. 2021。本仓库不转载输入论文、基准文本和含全文的提示词。
