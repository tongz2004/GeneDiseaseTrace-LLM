# EU-ADR：第二个人工标注基准（2026-10-29）

目的：在 BioRED 之外，用另一个人工标注语料检验冻结的 V1（Claude Opus 5.5，判定卡未改动）。EU-ADR 的关系分为正向（PA）、推测性（SA）和负向（NA），所以还能检验 BioRED 测不到的东西：我们的 uncertain 和 explicit_negative 判断准不准。

## 数据
- **来源**：EU-ADR 语料（van Mulligen et al., *J Biomed Inform* 2012）。通过 Hugging Face 的 dataset server 以 JSON 取得（bigbio/euadr，euadr_source 配置），放在 `data/`。
- **许可证**：未注明。原始数据、生成的 JATS 文件和标准答案只留在本地，已写进 `.gitignore`，不上传。
- **范围**：Target-Disorder 子集的全部 100 篇摘要，不抽样。名单冻结在 `sample_euadr.json`。
  - 其中 64 篇至少有一对基因–疾病关系；另外 36 篇用来检验我们会不会报告不存在的关系。
- **标准答案**：只用共识标注（标记为 True，即至少两位标注者一致）。
  - 评分对象是基因–疾病配对，即 "Genes & Molecular Sequences" 与 "Diseases & Disorders" 之间的关系。SNP–疾病关系（PA 49、SA 8、NA 9）只计数、不评分，与 BioRED 的做法一致。
  - EU-ADR 是按句子标注的，我们按配对评分：同一对基因和疾病，只要任一句是 PA 就算 PA，否则看 SA，再否则为 NA。有 12 对在不同句子里标签不同。
  - 最终得到 **134 对**：PA 122、SA 9、NA 3。
- **标签含义**：PA 正向关联、SA 推测性关联、NA 负向关联（文中明确说没有关联）。依据是 EU-ADR 原论文对关系类别的描述，在 2015 年一篇关于该语料众包评估的论文里有转述。
- **适配**：`prepare.py`。
  - 文本取标题加摘要，用一个空格连接，标注的字符位置就是按这个拼接计算的。
  - 有些标注的位置偏了几个字符，按标注文字找最近的出现位置来修正；文字在摘要里找不到的提及有 34 个，予以舍弃。
  - 有的摘要在源数据中有编码错误（例如 "KrÃ¼ppel"），所以改用 difflib 对齐原文和解析后的文本。
  - 1,620 个标准答案提及全部落在正确文字上。

## 运行
```bash
python experiments/euadr/prepare.py
python experiments/single_shot/batch.py --variant v1_card_lean_medium --run-set euadr --papers <sample_euadr.json 里的 100 个编号>
python experiments/euadr/evaluate.py
```
- 100/100 通过校验，89 篇一次通过。
- 合计 $4.75（每篇 $0.048），1,841 秒。
- 没有空回答或调用失败。

## 结果
| 评分范围 | 预测 | 精确率 | 召回率 | F1 [95% CI] |
|---|---|---:|---:|---|
| PA | supported | 0.469 | 0.549 | 0.506 [0.40, 0.60] |
| PA | supported + uncertain | 0.414 | 0.680 | 0.514 [0.43, 0.59] |
| PA + SA | supported | 0.477 | 0.527 | 0.500 [0.40, 0.60] |
| PA + SA | supported + uncertain | 0.424 | 0.656 | 0.515 [0.43, 0.59] |

**我们的判定与 EU-ADR 标签对照**（每个标准配对，我们给了什么）：

| EU-ADR 标签 | supported | uncertain | explicit_negative | not_found | 实体没记录 | 合计 |
|---|---:|---:|---:|---:|---:|---:|
| PA（正向） | 67 | 16 | 4 | 10 | 25 | 122 |
| SA（推测） | 2 | 1 | 3 | 0 | 3 | 9 |
| NA（负向） | 0 | 0 | **3** | 0 | 0 | 3 |

### 为什么精确率远低于 BioRED（0.47 对 0.89）
我们判为 supported 的 128 对里，61 对与 PA 或 SA 一致。不一致的 67 对分两类：
- **37 对：两个实体 EU-ADR 都标注了，但两者之间没有共识关系。**
  - 10 对有一位标注者标了关系，但没有达到共识，所以不在标准答案里；
  - 3 对的共识关系挂在别的实体类型上；
  - 24 对完全没有标注。
  - 其中有大量明显成立的关系，例如 GJB2–听力损失、CEP290–Leber 先天性黑蒙、p53–Li-Fraumeni 综合征、NPAS2–B 细胞淋巴瘤（摘要报告 OR = 0.61）。
- **30 对：我们的基因或疾病，EU-ADR 没有作为共识实体标注。** 基因 17 对、疾病 8 对、两者都不是 5 对，例如 IL1A。

**结论**：EU-ADR 的共识标注明显不完整。这里的精确率是**下限**，不能直接与 BioRED 比较。要得到真实的精确率，需要人工复核这 67 对。我们的证据句都在 `runs/euadr/v1_card_lean_medium/` 里，复核时可以直接对照。

### 漏掉的 55 对 PA
- 16 对：我们判为 uncertain。
- 15 对：我们没把它记为基因。其中多数是染色体区段，EU-ADR 把它们归在 "Genes & Molecular Sequences" 下，例如 3q26、7q 扩增、1p 杂合性缺失，按我们的规则不算基因；也有真正的漏判，例如 IL-23。
- 10 对：判为 not_found。
- 7 对：疾病被我们记为症状。
- 4 对：判为 explicit_negative。其中 PKLR、NOS1AP 与 2 型糖尿病两对，EU-ADR 在别的句子里也标了 NA（原文是"此前的关联在丹麦人群中没有重复出来"），所以我们的判断有依据。
- 3 对：疾病没有记录。

### 负向和推测性判断
- **NA**：3 对全部判为 explicit_negative。说明我们的负向判断与专家一致，但只有 3 对，样本太小。
- **SA**：9 对里只有 1 对判为 uncertain，3 对判为 explicit_negative，2 对判为 supported。EU-ADR 的"推测性"和我们的 uncertain 并不对应：EU-ADR 把"本研究没发现关联、但值得进一步研究"之类的表述也归为推测。样本同样太小。

### 与 BioRED 的关系
- 召回率与 BioRED 相近：0.55 对 0.43；把 uncertain 也算作有关系时，0.68 对 0.58。
- 精确率差距主要来自标准答案的完整程度，而不是方法本身。要确认这一点，需要人工复核上面那 67 对。

## 与 PubTator3 对比（同样 100 篇，同一套对齐和评分规则）
`python experiments/biored/pubtator_baseline.py euadr`：取 PubTator3 的摘要版公开标注，结果写到 `runs/euadr/pubtator3/`；然后运行 `python experiments/euadr/evaluate.py pubtator3 euadr`。PubTator3 列出了关系的配对算作预测有关系。

| 评分范围 | 系统 | 精确率 | 召回率 | F1 [95% CI] |
|---|---|---:|---:|---|
| PA | 我们（supported） | **0.469** | 0.549 | 0.506 [0.40, 0.60] |
| PA | PubTator3 | 0.408 | **0.754** | 0.529 [0.45, 0.61] |
| PA + SA | 我们（supported） | **0.477** | 0.527 | 0.500 [0.40, 0.60] |
| PA + SA | PubTator3 | 0.438 | **0.771** | 0.558 [0.48, 0.64] |

- **F1**：两者相近，置信区间重叠。我们的精确率略高，PubTator3 的召回率明显更高。
- **NA（明确无关）**：3 对里，我们全部判为 explicit_negative；PubTator3 把其中 1 对判为有关系，另外 2 对没有预测。
- **SA（推测性）**：9 对 PubTator3 全部判为有关系。它只有"有关系"和"没预测"两种输出，无法表达"推测"；我们把其中 1 对判为 uncertain、3 对判为 explicit_negative、2 对判为 supported。
- **注意**：EU-ADR 的共识标注不完整，两个系统的精确率都是下限。
