# 工作站流程与成熟项目参考

本次改进依据用户上传报告的正文章节和软件方法，并参考以下实际可访问的项目及其源码/使用文档。报告 HTML 引用的图片、外部 Result 数据未上传，因此没有把报告里的生物学结果当作本地重新分析结果。用户报告文件、合同资料、样本清单和客户数据不随源码发布。

| 项目 | 本工作站采用或借鉴的部分 | 说明 |
|---|---|---|
| [你的 16S 工作站](https://github.com/kangweihou522-cmyk/16s-) | Windows/WSL 浏览器操作、分阶段分析、分组与图表编辑 | 参考本地服务与操作流程；转录组使用独立定量和统计模型 |
| [nf-core/rnaseq](https://github.com/nf-core/rnaseq) | 样本表、文库方向、预处理、定量/QC、追溯 | 参考其工作流设计；本工作站不是完整 nf-core 分支，也未集成其 BAM 全流程 |
| [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance) | 全样本共同建模、对比组合、协变量、报告 | 不把多个组拆成独立统计，以保留适当的离散度估计 |
| [MultiQC](https://github.com/MultiQC/MultiQC) | 汇总 fastp / Salmon | 实际运行 MultiQC 1.27.1；导入计数不伪造 reads QC |
| [DESeq2](https://github.com/thelovelab/DESeq2) | 负二项模型、Wald 检验、normal LFC shrinkage | 实际使用 Bioconductor 包；未校正与 BH 方案同时保留 |
| [tximport](https://github.com/mikelove/tximport) | Salmon 到基因、长度 offset | TPM 不输入 DESeq2 |
| [fgsea](https://github.com/ctlab/fgsea) | multilevel GSEA、leading edge | 实际使用 1.24.0 / Bioconductor 3.16 稳定源码 |
| [GSVA](https://github.com/rcastelo/GSVA) | ssGSEA 单样本通路活性 | 实际使用 1.46；同模型设计下用 limma 比较通路分数 |
| [WGCNA](https://github.com/cran/WGCNA) | signed 网络、soft threshold、模块/eigengene/kME | 实际使用 1.72.1；小样本/重复对象会跳过 |
| [decoupleR](https://github.com/saezlab/decoupleR) | 基于小鼠先验网络的 ULM 活性 | 实际使用 2.4.0；必须提供有效的小鼠 source/target/weight 网络 |
| [DoRothEA](https://github.com/saezlab/dorothea) | 小鼠调控网络可用来源 | 数据资源需要另备；不自动把人类网络当成小鼠网络 |
| [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler) / [enrichplot](https://github.com/YuLab-SMU/enrichplot) | 富集背景、方向、gene ratio、富集曲线展示 | 借鉴展示和分析组织；ORA 当前是明确记录背景的超几何检验/BH，并非直接调用 clusterProfiler |
| [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) | 聚类、分组注释 | 借鉴布局；实际导出为 SciPy/Matplotlib，交互为 Plotly |
| [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano) | 点边框、标签、指定基因标注、避让 | 借鉴显示思路；导出使用 adjustText，交互标注可拖动 |
| [ideal](https://github.com/federicomarini/ideal) / [pcaExplorer](https://github.com/federicomarini/pcaExplorer) | 先诊断样本、再进入基因/通路探索 | 借鉴工作站信息组织与交互，不直接复制其 Shiny 界面 |
| [msigdbr](https://github.com/igordot/msigdbr) | Hallmark 等基因集的获取与物种管理 | 使用已获取且有权使用的小鼠 GMT；跨物种转换必须说明 |

## 按报告章节组织的主线

1. 样本与实验设计：CSV/TSV/文本xls、Windows路径转换、比较组合、协变量、参考版本、预检查。
2. 测序与定量质控：fastp 原始/过滤 reads、Q20/Q30、GC、接头及每碱基曲线；Salmon 映射率/方向；MultiQC。
3. 表达与样本诊断：文库大小、检测基因数、表达分布（VST / log2 normalized / log2 TPM）、PCA、Pearson 相关、样本距离和离散度。
4. 差异分析：DESeq2 共同模型、可选 LFC normal 收缩、双筛选方案计数、上/下调清单、交互火山图、MA、聚类热图、单基因表达。
5. GO/KEGG/GMT：按 all/up/down 分开 ORA；明确 eligible universe 与 BH 范围；fgsea multilevel、NES/log2err/leading edge、富集曲线。
6. 通路与关键基因：ssGSEA 单样本活性、limma 协变量调整比较、leading-edge 与 DEG 联合筛选、ID/Symbol/Entrez 基因查询、探索性表型关联。
7. 共表达与调控：WGCNA 模块/trait/kME/边表和 soft-threshold 诊断；decoupleR ULM 先验调控活性及同设计的 limma 比较。
8. 报告与导出：自包含离线 HTML；SVG/PDF/PNG、绘图数据和设置；RDS、日志、版本、输入记录和结果 ZIP。

## 两套差异筛选

每个比较同时生成：

- **nominal**：`pvalue < 0.05` 且 `abs(log2FoldChange) > 0`。
- **strict**：`padj < 0.05` 且 `abs(log2FoldChange) > 1`。

这里的比较符号确实为严格小于/大于。两个方案共用同一个 DESeq2 模型与同一 LFC 列（是否收缩由参数指定），分别保留 complete differential、DEG all/up/down 与 ORA。工作站先展示两种方案的数量，再选择 active profile。切换不会删除另一方案，也不重拟合模型；会更新当前清单、ORA 和相关图表，并保存 `active-profile.json`。显示阈值编辑只影响图，不替换 active profile 的统计。

nominal 未控制多重检验的 FDR，适用于明确标注的探索；本工作站允许选择，不通过自动放宽阈值制造结果。strict 没有基因时明确显示为空。手动指定非 DEG 基因的热图属于自选基因探索。GSEA、ssGSEA、WGCNA、调控活性使用全部保留表达基因，不先按任何 DEG 阈值筛选。

## 可视化交互

火山/MA：点带边框，悬停显示基因，点击显示 ID、Symbol、LFC、P/校正P、平均表达；可输入 ID/Symbol 定位高亮并加标签；可拖动标签，保存后同步到 SVG/PDF/PNG 和报告。导出标签自动避让，字体、点大小/透明度、阈值、配色、轴范围和自动/指定标签都可调整。

热图：变量基因/差异基因、自选基因、聚类排序、Symbol/ID、样本分组色条和点击基因探索。PCA：主成分选择、分组配色、悬停样本。富集：按方向筛选，点大小表示命中数量，颜色表示校正P；点击通路查看 GSEA。GSEA：运行富集曲线、命中位置与导出图的 Wald rank track。高级 HTML 图形使用本地 vendored Plotly，页面无需 CDN 联网。

模块–表型表当前为未调整协变量的 Pearson 关联，并对表格检验作 BH；单基因–数字表型为探索性 Spearman，未做全基因多重检验。都不能作为因果证据。WGCNA 至少15个独立样本、推荐20+；存在重复 subject 时当前模块跳过，不以重复样本数抬高关联自由度。soft-threshold 需要同时考察 slope、signed R² 与连接度；默认 power=6 仅为起始参数。

## 需要 BAM 的报告扩展

用户报告中的 HISAT2 基因组比对、外显子/内含子分布、IGV、StringTie 新转录本、rMATS 剪接和 GATK/SnpEff 变异分析，不能从 Salmon quant.sf 推算。它们未被本工作站主线执行。可以另用 nf-core/rnaseq 生成基因组 BAM，再运行对应成熟模块；小鼠基因组/品系/注释版本需匹配。没有把这些模块显示成“已完成”。蛋白互作也需要 STRING 等真实互动资源，不把表达相关性当作蛋白互作。

## 2026-10-27 方法与图形修订

GO背景按BP/CC/MF分开；Salmon可选匹配基因组full-decoy；WGCNA提供自动/手动power与模块功能富集；新增调整后基因–表型关联、候选证据表、leading-edge热图和通路重叠网络。具体假设、校正范围、输入限制和旧任务兼容见 [方法审查](ANALYSIS_REVIEW.md) 与 [教程第27章](USER_GUIDE.md)。
