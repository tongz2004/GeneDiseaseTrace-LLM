# RENET2：全文上的人工标注基准（2026-10-30）

目的：BioRED 和 EU-ADR 都只有摘要，这一轮在全文上检验冻结的 V1（Claude Opus 5.5，判定卡未改动）。

## 数据
- **来源**：RENET2（Su et al., *NAR Genomics and Bioinformatics* 2021）作者服务器上的数据包 `renet2_data_models.tar.gz`（https://www.bio8.cs.hku.hk/RENET2/ ，约 692 MB）。
  - 下载时边传边解压，只保留全文数据 `data/ft_data/` 和 `data/ft_info/`，共 22 MB。
  - 数据许可证未注明，数据、生成的文本和标准答案只留在本地，已写进 `.gitignore`。
- **文本是截断的全文**：每篇只有标题、摘要和正文的前一部分（引言加一部分结果），在文章中途截断。
  - 500 篇的中位数约 5,400 字符；抽中的 40 篇平均约 1.2 万字符，比我们那 100 篇完整论文短得多。
  - RENET2 的标注就是在这段文本上做的，所以我们也读同样的文本。
- **标签**（`labels.txt`，按基因 ID 和疾病 MeSH ID 记录配对）：
  - associated：先由他们的模型提出候选，再经专家人工确认；
  - ambiguous：专家无法确定；
  - non-associated：没有人工核对，只是三个工具都判为无关。
- **实体**：来自 PubTator Central 的自动标注，所以也带有它的识别错误。例如把 CIMT（颈动脉内膜中层厚度）和 CRS 标成了基因，把"肿瘤生长""胰岛素降低"标成了疾病。
- **抽样**：入选条件是至少有一对 associated、正文不少于 1,000 字符，500 篇里有 187 篇符合。用种子 20261028 抽 40 篇，名单冻结在 `sample_renet2_40.json`。
  - 这 40 篇的标准配对：associated 181、non-associated 83、ambiguous 1。
  - 其中 3 对 associated 的基因或疾病在截断后的文本里没有出现，任何读这段文本的人都找不到，单独列出。
- **适配**：`prepare.py`。标注位置按"标题 + 空格 + 摘要 + 空格 + 正文"计算；4,210 个标准答案提及全部落在正确文字上。

## 运行
```bash
python experiments/renet2/prepare.py
python experiments/single_shot/batch.py --variant v1_card_lean_medium --run-set renet2 --parallel 2 --papers <sample_renet2_40.json 里的 40 个编号>
python experiments/renet2/evaluate.py
```
- 40/40 通过校验，32 篇一次通过。RENET32321856 修正一次后仍有格式问题，单独重跑一次后通过，和之前几轮的处理方式相同。
- 合计 $6.36（每篇 $0.16），1,803 秒。没有空回答或调用失败。

## 结果（40 篇，以 associated 为标准答案）
| 预测 | 精确率 | 召回率 | F1 [95% CI] |
|---|---:|---:|---|
| supported | 0.318 | 0.370 | 0.342 [0.26, 0.42] |
| supported + uncertain | 0.312 | 0.464 | 0.373 [0.28, 0.46] |

**分数明显低于 BioRED，但大部分差距来自标准答案的覆盖范围：**

1. **我们多报的 152 对，几乎都是 RENET2 没有标过的配对。**
   - 105 对：两个实体 RENET2 都有，但这一对从没被他们的工具提出来，所以既没标"相关"也没标"不相关"；
   - 45 对：我们的实体不是 RENET2 的概念，例如 miRNA（let-7、miR-15/16）、"ryanodine receptor"；
   - **只有 2 对是 RENET2 明确标为 non-associated 的。** 也就是说，在 223 对 supported 里，和标准答案直接矛盾的只有约 1%。
2. **多报的大多是背景知识。** 152 对里有 136 对我们自己标为"背景"，也就是引言里复述的已知关联，例如 miR-15/16 与慢性淋巴细胞白血病、Dicer 与肺癌。RENET2 的候选是由一个在摘要上训练的模型生成的，似乎很少提出这类配对。
   - 按我们的来源标签分开算精确率：标为"本研究/两者兼有"的是 **0.53**（18/34）；标为"背景"的是 0.28（53/189）。
3. **漏掉的 114 对 associated**：
   - 54 对：我们读过这一对，判为 not_found；
   - 17 对：判为 uncertain；
   - 43 对：实体没记录。其中 27 对是疾病没记录，多为"肿瘤生长""大脑中动脉闭塞"这类宽泛或表型性的说法；9 对被我们记为症状；7 对是基因没记录，其中不少是 PubTator 的识别错误，例如把 CIMT、CRS 当成基因。
   - 把 uncertain 也算作有关系，召回率从 0.37 升到 0.46。

**结论**：
- 在这批截断的全文上，我们和 RENET2 **直接矛盾的很少**：我们判为 supported 的 223 对里，只有 2 对被 RENET2 明确判为不相关。
- 分数低的主要原因有三：RENET2 只标了他们的工具提出过的候选配对；它的实体层来自 PubTator 的自动标注；而我们会报告引言里的背景关联。
- 这里的精确率是**下限**。要得到真实数字，需要人工复核那 105 对没被标注过的配对。
- 召回率方面，54 对被我们判为 not_found，是这次最值得逐条看的地方，可以作为下一步优化规则的方向。

## 与 PubTator3 对比（同样 40 篇，同一套对齐和评分规则）
`python experiments/biored/pubtator_baseline.py renet2`：取 PubTator3 的全文版标注，这样正文里的实体也能对齐；关系预测只来自摘要，这一点已在 BioRED 上确认过。结果写到 `runs/renet2/pubtator3/`；然后运行 `python experiments/renet2/evaluate.py pubtator3 renet2`。

| 系统 | 预测的配对数 | 精确率 | 召回率 | F1 [95% CI] |
|---|---:|---:|---:|---|
| 我们（supported） | 223 | 0.318 | **0.370** | 0.342 [0.26, 0.42] |
| 我们（supported + uncertain） | 320 | 0.312 | **0.464** | 0.373 [0.28, 0.46] |
| PubTator3 | 61 | **0.623** | 0.254 | 0.361 [0.24, 0.48] |

- **F1**：两者相近。PubTator3 预测得少（只看摘要，40 篇共 61 对），所以精确率高；我们读了截断后的全文，预测多得多，找到的标准关系也更多（67 对，PubTator3 是 46 对），但多报的大多是 RENET2 没标过的背景配对，见上文。
- **和"明确不相关"的配对冲突**：在 83 对 non-associated 里，我们判为 supported 的有 4 对，PubTator3 有 2 对。
- **PubTator3 有主场优势**：RENET2 的实体本身就是 PubTator Central 标注的，所以实体对齐这一步对 PubTator3 天然有利。
- 40 篇样本小，置信区间宽。
