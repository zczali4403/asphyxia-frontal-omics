# 机械性窒息对小鼠额叶影响的多组学分析

本仓库用于分析机械性窒息组（DMA）与对照组（CON）小鼠额叶的转录组、蛋白组和代谢组数据。额叶在原始数据中编号为 R1。

## 实验设计

| 组学 | DMA | CON | 分析范围 |
|---|---:|---:|---|
| 转录组 | 3 | 3 | 额叶（R1） |
| 蛋白组 | 4 | 4 | 额叶（R1） |
| 代谢组 | 6 | 6 | 额叶（R1）；另有 3 个 QC 样本 |

同一实验组内编号相同的样本来自同一只动物。转录组样本名中的 `CON_R1_1`、`DMA_R1_1` 分别对应统一动物编号 `CON-1`、`DMA-1`，其余编号以此类推。

- 三组学完整配对：每组 3 只动物，即 CON-1～3 和 DMA-1～3。
- 蛋白组与代谢组配对：每组 4 只动物，即 CON-1～4 和 DMA-1～4。
- 代谢组另有每组 2 只动物，即 CON-5～6 和 DMA-5～6。

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

- 转录组差异分析：`data/transcriptomics/frontal_cortex_gene_counts.csv`
- 转录组可视化和联合分析：`data/transcriptomics/frontal_cortex_vst_expression.csv`
- 蛋白组分析：`data/proteomics/protein_expression_log2.xlsx` 的“数据矩阵”工作表
- 代谢组分析：`data/metabolomics/metabolite_expression.xlsx` 的“数据矩阵”工作表

`all_regions_vst_expression.csv` 保留了 R1–R4 全部转录组样本，用于脑区特异性比较；额叶三组学联合分析应使用已经提取好的 `frontal_cortex_vst_expression.csv`。

## 分析注意事项

- 转录组原始 count 用于 DESeq2 等计数模型，不能把 VST 值作为 DESeq2 输入。
- 蛋白组现有差异列表按名义 `P < 0.05` 和 FC 阈值筛选，但多重校正后没有 `q < 0.05` 的蛋白。
- 代谢组现有差异列表也没有 `q < 0.05` 的代谢物，现有结果应作为探索性发现。
- 原始 FASTQ、质谱 RAW 文件和 tar 包体积较大，未复制到本仓库。
- 元数据同时保留原始脑区编号 `R1` 和标准名称 `frontal_cortex`。
- 个体水平三组学分析应以 6 只完整配对动物为准；蛋白组—代谢组配对分析可以使用 8 只动物。

数据来源、工作表说明和校验值见 `data/README.md`。
