# R1 小鼠脑组织多组学分析

本仓库用于分析 DMA 与 CON 两组小鼠 R1 脑区的转录组、蛋白组和代谢组数据。

## 实验设计

| 组学 | DMA | CON | 分析范围 |
|---|---:|---:|---|
| 转录组 | 3 | 3 | R1 脑区 |
| 蛋白组 | 4 | 4 | R1 脑区 |
| 代谢组 | 6 | 6 | R1 脑区；另有 3 个 QC 样本 |

三种组学的样本量不同，当前资料中没有跨组学动物个体对应表。因此，同名样本不能自动视为来自同一只动物。默认开展组水平、分子映射和通路水平整合；只有补充个体对应关系后，才能进行严格的跨组学配对相关分析。

## 目录结构

```text
.
├── config/
│   └── samples.csv
├── data/
│   ├── transcriptomics/
│   ├── proteomics/
│   ├── metabolomics/
│   ├── metadata/
│   └── README.md
├── scripts/
└── results/
```

## 推荐输入矩阵

- 转录组差异分析：`data/transcriptomics/r1_frontal_cortex_gene_counts.csv`
- 转录组可视化和联合分析：`data/transcriptomics/r1_vst_expression.csv`
- 蛋白组分析：`data/proteomics/protein_expression_log2.xlsx` 的“数据矩阵”工作表
- 代谢组分析：`data/metabolomics/metabolite_expression.xlsx` 的“数据矩阵”工作表

`all_regions_vst_expression.csv` 保留了 R1–R4 全部转录组样本，用于脑区特异性比较；R1 三组学联合分析应使用已经提取好的 `r1_vst_expression.csv`。

## 分析注意事项

- 转录组原始 count 用于 DESeq2 等计数模型，不能把 VST 值作为 DESeq2 输入。
- 蛋白组现有差异列表按名义 `P < 0.05` 和 FC 阈值筛选，但多重校正后没有 `q < 0.05` 的蛋白。
- 代谢组现有差异列表也没有 `q < 0.05` 的代谢物，现有结果应作为探索性发现。
- 原始 FASTQ、质谱 RAW 文件和 tar 包体积较大，未复制到本仓库。
- 当前按原始文件命名将 R1 解释为额叶（frontal cortex）。

数据来源、工作表说明和校验值见 `data/README.md`。
