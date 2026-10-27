# 分析与可视化实现审查（2026-10-27）

此记录对应本次源码修订。逐段审查主线代码并使用本地合成数据验证；不代表用户真实样本已通过质量检查，也不把演示注释当作小鼠生物学验证。完整原理、操作、术语与参考在 [使用教程](USER_GUIDE.md)。

| 环节 | 对照依据与审查结论 | 本次处理 / 边界 |
|---|---|---|
| 样本与设计 | DESeq2 官方 vignette 的 full-rank、paired design、biological replicates | 保留满秩/残差自由度校验；增加单水平分类协变量和 SE/PE 文库一致性检查；全链保留样本/基因前导零，阻止重复FASTQ/quant文件伪作多个样本。重复动物需要显式 subject 与合适模型。 |
| FASTQ 预处理 | fastp 官方项目 | 保留 fastp；补足原始 Q20/Q30/GC 记录、质量/碱基曲线 R1/R2 与过滤阶段。counts/quant 导入不伪造 reads 质控。 |
| Salmon | Salmon / nf-core/rnaseq 的 decoy-aware selective alignment | 新增可选匹配基因组 FASTA，流式生成 gentrome + decoys.txt，显式传给 Salmon。单端显式提供可调片段长度均值/标准差（默认200/80，不能误作read长度）；未提供genome时继续转录本索引并提示限制；小鼠全基因组索引需另评估内存。 |
| tximport | 官方 summarizeToGene 接口和 vignette | 同时规范化 tx2gene 转录本版本号和管道符后缀；量化导入校验有效长度/丰度/一致转录本列表；保留标准 DESeqDataSetFromTximport 长度处理。不能用 TPM/FPKM 替代原始计数。 |
| 表达与 PCA | DESeq2 VST / plotPCA | 保留 VST、最高变异500基因、center=TRUE/scale=FALSE；旧任务 PCA 后备路径统一最高变异规则。相关图按实际 1-r 距离聚类，距离图按已有距离聚类。 |
| DESeq2 | 官方 Wald test、independent filtering、LFC shrinkage | 保留两种严格边界筛选及同一模型；保存缺失 P、独立过滤、排序 ties 数量。P=0.05 / |LFC|=1 不满足严格 < / >。筛 LFC 不是对 |LFC|>1 的正式阈值假设检验。normal 收缩不改原始 Wald stat。 |
| GO / GMT ORA | 超几何分布、BH、clusterProfiler 的 universe 思路 | **修正 GO 背景按 BP/CC/MF 分开**；每方向/ontology 全部可测试集合 BH；GMT 同方向统一 BH。先排除非有限 Wald，再按平均表达选择重复ID代表行，避免无效行遮蔽有效行。 |
| GSEA | fgseaMultilevel 官方源代码与输出定义 | 全部保留、有限 Wald 排序，独立于 DEG 方案；保持 NES/ES/BH/log2err/leadingEdge。三层曲线展示 running ES、命中条带和 Wald stat。相同 stat 顺序确定，但 ties 仍影响解释。 |
| ssGSEA / ULM | GSVA、decoupleR、limma | 保持 VST 与原设计组间检验；活性热图增加逐行Z或原分数选择。不同通路/调控子的绝对分数不具有统一强度单位。 |
| WGCNA | WGCNA 方法及 pickSoftThreshold / blockwiseModules | 新任务 UI 默认自动选首个 signed R²≥0.85 且 slope<0 的 power；无候选明确跳过并保留诊断，也可手动。此阈值是启发式，不是网络质量认证。grey 排除出生物模块数/hub；保留15独立样本门槛与重复subject跳过。 |
| 表型关联 | 带协变量的线性残差/偏相关与 t 分布；BH | 新增 VST 基因和连续表型对完整设计矩阵（含分组）的残差 Pearson；df=n-rank(X)-1，≥8独立样本且df≥3；全基因×可检验表型统一BH。混杂/不足/重复subject明确跳过。保留原始单基因 Spearman 为探索性描述并说明渐近P的局限；重复subject不计算普通通路相关P，基因卡片只给描述系数。 |
| 信息整合 | enrichplot 的关系网络/核心富集可视化思路 | 新增 leading-edge VST 热图、Jaccard 通路重叠网络、DEG跨比较重叠、候选证据表、WGCNA hub/子网络及模块GO/GMT富集（网络输入基因为背景，跨模块×术语按ontology统一BH）。网络不伪称PPI，交集不伪称独立验证，无任意加权“生物标志物概率”。 |
| 预览与导出 | Plotly、Matplotlib、EnhancedVolcano 的交互与标注思路 | 同一设置与选中数据驱动两种引擎；方形火山图，比例居中，非写入预览，配色预设，独立行列聚类、色阶、表达尺度、通路样式、数值与图例。保留SVG/PDF/PNG/数据/设置/报告。布局由两种引擎分别排版，不保证像素相同。 |
| 兼容和复现 | 保存原始结果、方法与版本 | 旧任务不自动重算统计；新统计输出需新任务。图表错误单独记录 plot-warnings.json，统计表保留；ZIP排除原始reads、量化、索引和新建的大型reference目录。 |

## 实际查阅的官方 / 成熟源码

查阅日期为 2026-10-27；以下链接会随上游开发更新，工作站使用的软件版本以任务 `R-session.txt` 与容器为准。这些项目是方法或界面参考，不意味着工作站直接执行 nf-core/clusterProfiler/EnhancedVolcano。

- [DESeq2 官方 vignette](https://github.com/thelovelab/DESeq2/blob/devel/vignettes/DESeq2.Rmd)：实验设计、归一化、Wald、缺失 P、独立过滤、收缩和阈值检验。
- [tximport 官方源码](https://github.com/thelovelab/tximport/blob/master/R/tximport.R)：版本号处理、转录本到基因及长度入口。
- [nf-core/rnaseq 官方 usage](https://github.com/nf-core/rnaseq/blob/master/docs/usage.md)：提供 genome 时 full-decoy 索引的推荐与不提供时的边界。
- [Salmon](https://github.com/COMBINE-lab/salmon)：Patro et al., Nature Methods 2017, DOI 10.1038/nmeth.4197；Srivastava et al., Genome Biology 2020, DOI 10.1186/s13059-020-02151-8。
- [fgseaMultilevel](https://github.com/alserglab/fgsea/blob/master/R/fgseaMultilevel.R)：多层自适应估计、最小集合、NES、误差和BH。
- [WGCNA 官方 R 实现镜像](https://github.com/cran/WGCNA/blob/master/R/Functions.R)：signed network、软阈值与模块；Langfelder & Horvath, BMC Bioinformatics 2008, DOI 10.1186/1471-2105-9-559。
- [GSVA](https://github.com/rcastelo/GSVA)：Hänzelmann et al., BMC Bioinformatics 2013, DOI 10.1186/1471-2105-14-7。本工作站兼容镜像内 GSVA 1.x 接口，未声称已迁移 GSVA 2.x 新接口。
- [decoupleR ULM](https://github.com/saezlab/decoupleR/blob/master/R/statistic-ulm.R)：靶基因网络的单变量线性模型活性。
- [enrichplot](https://github.com/YuLab-SMU/enrichplot)、[EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano/blob/master/R/EnhancedVolcano.R)、[ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)：借鉴交互问题、网络关系和注释可读性，实际绘图仍为本项目 Plotly/Matplotlib。

## 真实数据仍需人工判断

测序深度与 RNA 完整性是否满足实验问题；样本是否独立；批次、组织组成和处理是否混杂；注释/参考版本是否匹配；小样本结果的不确定性；过度解释 nominal P；模型残差/离群点、GSEA ties 和网络稳定性。自动流程不能替代这些判断。当前仍未内置 BAM 剪接/变异分析、细胞类型去卷积、多组织随机效应、多时间点交互模型或独立队列验证。
