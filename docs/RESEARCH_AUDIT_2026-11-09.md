# 转录组工作站：上游检索、方法与交互审查

审查日期：2026-11-09。被审查版本：`0379ec893f16c61f862ccd83ea1af8592dd2f3d0`。对象：小鼠 bulk RNA-seq 工作站的当前源码、统计设计、界面和交付方式。

本次交付是审查结论与改进规格。文中的“建议”不代表功能已经上线；当前功能以“现状”及对应代码为准。同步操作教程在 [USER_GUIDE.md](USER_GUIDE.md) 第29章，网页和 Word 从该原稿生成。

后续状态：本报告保留`0379ec8`基线审查。增强版已经落实的功能与仍有边界的方向见 [教程第30章](USER_GUIDE.md) 和 [验证记录](../VALIDATION.md)，不能把本报告的历史缺口当作当前状态。

## 1. 判断与证据范围

**常规独立样本的两组、多组两两比较，现有主线有合理方法基础；对于复杂设计、跨比较解释、混杂控制与初学者操作，仍需明显完善。** 不宜仅因已经拥有很多图，就称其覆盖了所有研究设计。

- 方法链是 fastp → Salmon → tximport → DESeq2 → GO/GMT ORA、fgsea → ssGSEA/limma、ULM、可选 WGCNA。Salmon/tximport/DESeq2 的组合有公开方法论文与官方工作流支持。
- 现有三类主输入、两个固定 DEG 方案、单步模块、显式结果接力、清空新任务输入、FASTQ 成功步骤复用，以及图表点击/定位/导出，值得保留。
- 主要缺口是实验设计容易填错、整体多组检验缺失、一些表型关联没有调整混杂、参考与注释匹配缺少可见诊断、多组结果缺乏统一比较视图，以及手工文件准备仍偏多。

检索以 GitHub 官方项目的 README、vignette、usage、源码、CITATION 为主，辅以其中可核对的论文元数据。也浏览了 GitHub `rna-seq` topic 和 `bulk RNA-seq analysis` 仓库检索页作为发现入口。没有取得同一口径的实时文献被引数据库，**不提供“引用次数排名”或虚构 GitHub 星数**；经典方法论文、社区维护工作流和界面参考分开使用。部分参考为综述或教学工作流，fgsea 文献为预印本，不能全部统称高分 SCI 原始研究。

本次实际核对的上游文件 URL、抓取内容 SHA256、源码证据位置和检查结果保存在 [检索与验证记录](research/RNASEQ_AUDIT_2026-11-09.json)。SHA256 是本次获取内容的指纹，不等于上游发布版本；上游开发分支后续可能变化。

本次执行了56项现有 Python 回归测试、2/3/4组的合成输入与预检查测试、桌面/移动宽度的真实 Chromium 操作。没有用本次输入检查冒充新的真实 FASTQ 全链统计验证；较早的实际 R/Salmon 合成运行记录见 [VALIDATION.md](../VALIDATION.md)。没有读取或发布用户真实测序数据。

## 2. 值得重点参考的成熟体系

| 项目 / 官方材料 | 成熟依据或论文 | 对本工作站最有价值的内容 |
|---|---|---|
| [nf-core/rnaseq](https://github.com/nf-core/rnaseq)，[usage](https://github.com/nf-core/rnaseq/blob/master/docs/usage.md)，[output](https://github.com/nf-core/rnaseq/blob/master/docs/output.md) | 社区策展工作流；nf-core 框架论文，Nature Biotechnology 2020，[DOI](https://doi.org/10.1038/s41587-020-0439-x) | 样本表、跨 lane 合并、链特异性判断、参考管理、QC 汇总、执行与版本记录。该框架论文不等同于对本工作站的验证；nf-core/rnaseq 本身也不负责完整差异统计。 |
| [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance)，[usage](https://github.com/nf-core/differentialabundance/blob/master/docs/usage.md) | nf-core 下游统计工作流；本次文档包含复杂公式、contrast 与 DREAM | 将样本信息、模型和比较分开定义；多比较报告、交互公式、重复测量入口。不能因为它现在支持 DREAM 就声称本工作站已经支持。 |
| [DESeq2](https://github.com/thelovelab/DESeq2)，[官方 vignette](https://github.com/thelovelab/DESeq2/blob/devel/vignettes/DESeq2.Rmd) | Love et al., Genome Biology 2014，[DOI](https://doi.org/10.1186/s13059-014-0550-8) | 配对、混杂、Wald/LRT、交互、Cook 距离、独立过滤、LFC 收缩与正式效应阈值检验。 |
| [tximport](https://github.com/thelovelab/tximport)；[rnaseqGene](https://github.com/mikelove/rnaseqGene) | Soneson et al., F1000Research 2015，[DOI](https://doi.org/10.12688/f1000research.7563.2)；Bioconductor 教学工作流 | Salmon 长度处理、基因级推断、检查与可复现教学。保留正确 tximport 入口，不把 TPM 取整后做 DESeq2。 |
| [Salmon](https://github.com/COMBINE-lab/salmon) | Patro et al., Nature Methods 2017，[DOI](https://doi.org/10.1038/nmeth.4197) | 偏倚校正、转录本参考、映射与文库诊断；结合 nf-core 使用 full-decoy 的建议。 |
| [iDEP 当前代码](https://github.com/gexijin/idepGolem) | Ge et al., BMC Bioinformatics 2018，[DOI](https://doi.org/10.1186/s12859-018-2486-6) | 矩阵导入、模块化引导、表达探索到通路的连续操作。旧 `gexijin/iDEP`、`iDEP-SDSU/idep` 已指向新代码，应避免照抄旧入口。其许可条件也需单独遵守。 |
| [DEBrowser](https://github.com/UMMS-Biocore/debrowser) | Kucukural et al., BMC Genomics 2019，[DOI](https://doi.org/10.1186/s12864-018-5362-x) | 过滤前后对照、交互结果浏览、样本信息与可下载结果。成熟项目的每条旧操作不一定适合本工作站：其 README 中对非整数取整的做法，不作为导入 TPM 的依据。 |
| [iSEE](https://github.com/iSEE/iSEE) | Bioconductor 交互探索工具；F1000Research 2018，[DOI](https://doi.org/10.12688/f1000research.14966.1) | 火山图、热图、样本和表格之间传递选择；套索、元数据着色与会话状态。它也面向单细胞，但其 SummarizedExperiment 交互思路适用于 bulk。 |
| [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)、[InteractiveComplexHeatmap](https://github.com/jokergoo/InteractiveComplexHeatmap) | Gu et al., Bioinformatics 2016，[DOI](https://doi.org/10.1093/bioinformatics/btw313) | 多层样本注释、分组切片、相邻统计面板、热图选区到基因明细。实际工作站仍用 Plotly/Matplotlib，参考其设计不要求整体替换引擎。 |
| [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)、[enrichplot](https://github.com/YuLab-SMU/enrichplot) | Wu et al., The Innovation 2021，[DOI](https://doi.org/10.1016/j.xinn.2021.100141)；Xu et al., Nature Protocols 2024，[DOI](https://doi.org/10.1038/s41596-024-01020-z) | compareCluster 式多组通路对照、GO 语义去冗余、基因–通路联动与可解释图例。当前 ORA 为自有超几何实现，不能写成已经运行 clusterProfiler。 |
| [fgsea](https://github.com/alserglab/fgsea)、[GSVA](https://github.com/rcastelo/GSVA)、[decoupleR](https://github.com/saezlab/decoupleR) | fgsea 方法预印本；GSVA 2013，[DOI](https://doi.org/10.1186/1471-2105-14-7)；decoupleR 2022，[DOI](https://doi.org/10.1093/bioadv/vbac016) | 完整排序、leading-edge、富集不确定性、单样本活性、调控网络覆盖与分数含义。GSVA 包中使用 ssGSEA 不等于同时执行所有 GSVA 方法。 |
| [ExploreModelMatrix](https://github.com/csoneson/ExploreModelMatrix) | Bioconductor 模型矩阵可视化工具 | 把“参数是否可估计、调整了谁、比较了什么”直接展示给用户。 |
| [variancePartition / dream](https://github.com/GabrielHoffman/variancePartition) | Hoffman & Roussos, Bioinformatics 2021卷期，[DOI](https://doi.org/10.1093/bioinformatics/btaa687) | 同一动物重复测量与随机效应；不能用多次观测数量代替独立动物数。 |

这些来源是不同层面的参考。工作站可以保留轻量 Salmon 路线，同时吸收 nf-core 的 QC 和记录方式；不需要把完整 STAR/BAM 分支作为每个任务的默认前提。用户的32GB物理内存也不能保证任意参考和全量 BAM 分析均能运行。

## 3. 是否支持两组及两组以上

### 3.1 当前确实支持的范围

完整分析或独立差异分析：一份样本表可以有2组以上；所有入选样本共同拟合 `~ 协变量 + group`，然后提取指定的处理组/对照组 Wald 对比。默认不拆成多个独立两组模型。这符合 DESeq2 常见多水平分组用法。

| 研究问题 | 当前状态 | 具体边界 |
|---|---|---|
| 两个独立组，例如 Control 与 Treatment | 支持 | 每组至少2个重复，建议至少3；最低门槛不是统计效能保证。 |
| 三组：Control、Model、Treatment | 支持 | 手工添加 Model/Control、Treatment/Control、Treatment/Model。不会因有3组就自动生成3个比较。 |
| 四组或更多的预先定义两两比较 | 支持 | 每任务1–30个比较。4组全两两为6个；8组为28个；9组全两两为36个，超过当前单任务上限。 |
| batch、性别等分类因素，年龄等连续因素的调整 | 支持加性协变量 | 连续因素须声明 numeric；满秩和残差自由度检查保留。分组与批次完全混杂时拒绝。 |
| 简单、可估计的配对设计 | 有条件支持 | 需显式加入 categorical subject；仅有 subject 列不会自动进入 DE 模型。 |
| 三组“任何一组不同”的总体检验 | 尚未实现 | 两两 Wald 不等于全局检验；建议增加 full `~ batch + group` 对 reduced `~ batch` 的 LRT。 |
| 基因型×处理、差中之差 | 尚未实现正式交互对比 | 可把组合编码为 group 得到组合间两两比较，但不能直接回答处理效应在不同基因型间是否不同。 |
| 时间趋势、非线性剂量趋势、复杂纵向/随机效应 | 尚未实现 | 加入一个普通 numeric 协变量不等于完成专门趋势或纵向检验；需单独的模型、假设和验证。 |
| 独立通路/TF活性任务输入多个组 | 支持同类两两比较 | 转换后表达与样本表，ssGSEA/ULM 后用相同加性设计 limma。 |
| 独立富集一次导入多个比较的外部差异表 | 当前不支持批量导入 | 一次一个比较的完整差异表；完整主流程可以输出多个比较的富集。 |
| 表达诊断、WGCNA 单步任务 | 接受多组样本 | 不自动生成基因差异比较；WGCNA还要求至少15个独立样本，推荐20以上。 |

源码依据：[validation.py](../validation.py#L96)、[analyze.R](../scripts/analyze.R#L11)、[analyze.R 对比](../scripts/analyze.R#L43)、[单步模块](../scripts/module_inputs.py#L62)。

### 3.2 本次实际验证

合成整数计数，交叉分布的两个 batch，每组3个样本：

| 输入 | 提交比较数 | 校验和 `/api/preflight` |
|---|---|---|
| 2组、6样本 | 1 | 通过，HTTP200 |
| 3组、9样本 | 3 | 通过，HTTP200；Chromium 实际填入并提交预检查通过 |
| 4组、12样本 | 6 | 通过，HTTP200 |
| subject列存在、未选subject协变量 | 1 | 当前会接受，实际协变量为空；这是误用风险 |
| 同一简单配对数据，明确选择subject | 1 | 通过，设计包含subject |
| batch完全等同group | 1 | 拒绝，提示混杂/不可估计 |

### 3.3 多组FDR必须说清

每个 DESeq2 对比中的 padj 是该对比内的基因多重检验校正，**不是所有对比×所有基因整体的FDR保证**。当前 ORA 又按方向、GO ontology/GMT 分别校正；GSEA、活性和模块关联也有各自检验家族。这些局部校正通常合理，但必须在结果页和导出中写出范围。

建议保留标准对比内 padj；若用户确实需要跨全部比较的发现控制，可提供定义明确的全局/分层选项，另存列和方法。不能对现成 padj 再简单做一次 BH 冒充全局校正，也不能悄悄改写旧结果。

## 4. 优先解决的统计与分析细节

### P0：降低会改变结论的误用风险

**A. 配对/重复测量引导。** 目前只默认识别多水平 batch，不默认加入 subject。重复subject会让WGCNA、普通关联跳过，但不会自动使DESeq2变为配对。应加入“独立动物 / 简单配对 / 重复测量”设计选择、对象数与观测数预览，并在有重复subject而模型未处理时阻止静默开跑、要求明确选择处理方式。简单配对可生成 `~ subject + group`；其他情况先检查是否可估计，不能无条件自动加subject。依据：DESeq2 配对设计与 dream 重复测量文档。

**B. 表型关联的混杂口径统一。** 当前 [基因–表型](../discovery.py#L30) 对完整设计残差化；[通路–表型](../scripts/mining.R#L70) 和 [模块–表型](../scripts/mining.R#L111) 仍是原始 Pearson，并做各自BH。BH不消除batch、group、年龄的混杂。建议并列提供原始相关与预先定义协变量调整关联，展示协变量、独立样本n、自由度、效应、区间和校正范围。是否调整group取决于问题：总体共变与组内关联不是同一估计目标，应由用户选择并标明。模块构建是否做批次调整是另一个决定，不能因DE模型有batch就声称网络已去批次。

**C. 注释、ORA背景和不可检验行可见化。** 当前有正确的ontology覆盖限定、重复ID代表和映射失败处理，但缺少“输入→匹配tx2gene→保留基因→GO/Symbol/Entrez/GMT命中”的覆盖率漏斗。应输出未映射名单、版本/物种、实际背景与独立过滤/Cook原因；低覆盖不能静默解释为没有生物学通路。

ORA现以有限stat定义背景，缺stat时用有限pvalue。DESeq2官方 [results源码](https://github.com/thelovelab/DESeq2/blob/devel/R/results.R) 对Cook异常行可将pvalue设NA而保留stat，因此两种口径不必相同。建议明确比较“表达/有限stat背景”和“当前检验可用P背景”的差别，结合独立过滤规则制定并记录适用策略、做敏感性核对，避免未经说明地把不同方案的分母混在一起。不能简单用DEG自身作背景。

**D. 对比与统计语义。** 每张图/表显示“处理/对照”、DEG方案、原始或收缩LFC、统计模型、校正家族。保留用户要求的两套规则：`pvalue<0.05 & |LFC|>0` 与 `padj<0.05 & |LFC|>1`。前者保持探索性标签，不因候选数少而把它自动改成确认性结论。绘图颜色阈值与真实下游名单阈值分开标明，防止画面变色被误认为重算分析。

### P1：完善常规多组研究与解释

**E. 多组总体检验与比较生成器。** 增加“全部处理对同一对照”“所有两两”“自定义”选择，自动生成方向一致且不重复的比较名。增加有协变量的LRT总检验；LRT回答是否存在组差异，不自带一个通用方向或统一log2FC，应返回单独表、再联动两两效应。不能将无方向LRT统计量直接作为有方向GSEA的排序。

**F. 诊断应能解释异常来源。** 现有PCA、文库大小、相关、离散度、P分布和cooks_max应保留；增加按batch/subject/性别切换PCA颜色/形状、样本级Cook异常基因数、实际model matrix可视化、不可检验行原因。低表达默认规则本身不是错误，但可提供基于实验组大小/设计的过滤建议和过滤前后敏感性摘要；不能为增加DEG而自动调参或自动删样本。

**G. 效应估计与检验分离。** 当前可选normal收缩，保留raw LFC，但替换展示LFC/lfcSE，原Wald stat/P不变。建议明列raw与shrunk LFC/SE、适用区间与估计来源；将apeglm/ashr作为经过对比适配测试后的选项，不直接替换任意contrast。`padj<0.05`再筛`|LFC|>1`不等于正式检验`|真实LFC|>1`；正式阈值检验应另设模式和输出，不更改用户固定的两套历史方案。

**H. 质控不止Q30和映射率。** 当前有fastp/MultiQC、Salmon映射率、expected_format和compatible_ratio；建议强化实际文库推断/兼容性、片段长度、跨样本异常汇总和适用的污染/rRNA筛查。映射率60%的提示是启发式，不能当统一质量合格证。没有BAM时不能从Salmon结果推算基因体覆盖、内外显子比例或剪接质量；确有需求时提供可选BAM/nf-core结果导入。PCR重复率也不能单独作为普通RNA-seq按坐标去重的理由；UMI文库需专门分支。

**I. 稳定的统计检查点与模块隔离。** FASTQ成功步骤已经能复用；目前完整R主脚本在富集/高级挖掘之后才保存dds.rds与最终summary。建议在DE拟合后先落盘可追溯检查点，将富集、活性、网络分成独立状态，某个可选模块失败仍能浏览已完成DE。依赖输入、模型、包版本的缓存失效规则必须明确；不能把任何旧dds都拿来重用。

## 5. 可视化：已有能力与应该升级的方向

当前已有近方形火山图（默认7×6.5英寸）、悬停/点击基因、输入ID/Symbol定位、拖动标注、图形预览、单图参数保存、配色/尺寸/色阶/聚类调节，以及SVG/PDF/PNG/数据/设置导出。这些应继续作为基础，而不是重新包装成“尚待实现”的功能。

| 图形或操作 | 当前不足 | 建议形成的交互 |
|---|---|---|
| 火山图 / MA | 能查单基因；多比较和多基因持续选择不足 | 套索形成候选列表→热图/表达图联动；多比较分面、共享坐标范围；标注过密提示、锁定标签和撤销；NA/极端P显示数量及截断说明。 |
| 表达热图 | 基因、排序、行Z和色阶可调；多层注释与选区解释不足 | 顶部同时标group/batch/sex/subject；组内聚类、分组切片、可见树状图；选行/列后更新样本/基因明细。Z-score保留明确含义，不能解释为绝对表达量。 |
| PCA / 样本相关 | 主要按group着色 | 颜色、形状、标签映射到可选元数据；PC选择；点击样本联动QC和样本表。展示调整前/用于展示的协变量调整视图时明确方法，避免删掉研究效应。 |
| 基因表达 | 有箱线/小提琴、各样本点、配对连线 | 小样本优先展示全部点；多基因×多组分面与统一顺序；可选模型调整效应及区间。原始表达均值不是协变量调整后LFC，不能从图上组均值比反推模型结果。 |
| 多比较总体结果 | 有双方案DEG数和Jaccard重叠 | 同一方案、同方向的UpSet；基因×比较LFC热图、效应森林图；通路×比较NES/活性效应矩阵；交叉比较时同时可见各自FDR口径。 |
| GO/GSEA | 有点/条/棒棒糖、曲线、leading-edge热图、重叠网络 | 可选择GO语义去冗余、通路分组与比较矩阵；点击通路后基因列表/样本表达/候选清单一致更新；去冗余仅压缩展示，不隐瞒原始检验全集。 |
| WGCNA / 表型 | 有软阈值、模块、hub、有限边子网与相关热图 | 模块树与表型热图联动、eigengene逐样本分布、原始/调整关联切换；将模块稳定性和preservation证据独立展示。hub与高kME不能称为已验证关键调控基因。 |
| 统一编辑和成图 | 主要逐图调节与下载 | 项目级字体、色盲友好色板、组别顺序、导出尺寸模板；一键批量套用可兼容图；保存/恢复可命名绘图会话；图注自动带模型、n、阈值与单位。 |

画布宽高比问题已在当前版得到基础修正；仍需用长中文通路名、密集标签、20+样本和多比较做排版压力验证。交互预览与静态导出使用不同引擎，应验证标签、选集、坐标和图例一致，不能承诺像素完全相同。可视化不需要追求3D或装饰性网络。

## 6. 输入与输出是否简单明了

**判断：熟悉数据文件的用户可以使用；对新手仍偏繁琐。** 当前浏览器支持本机路径和后端文件浏览，不等于已经具备拖入目录自动建样本、Excel列映射或参考资源向导。

实际Chromium检查中，3组计数任务设置页在1500×1000视口下总高度约3531px。390px宽度未发生整页横向溢出，但比较输入框仅约92px宽，长比较名需要横向查看；未开启的WGCNA设置、counts入口不需要的文库类型等仍占据认知空间。这些是具体的可用性问题，不能用“页面能打开”替代评价。

### 6.1 建议的输入流程

1. **选择目标和已有材料。** 完整FASTQ、已有Salmon、原始counts、独立富集、独立表达/活性/网络。每个入口只显示必需和适用参数，高级项默认折叠。
2. **读取数据并建立样本表。** 选择本地FASTQ目录后按R1/R2及lane规则生成待确认草稿；预览配对、重复ID、缺文件与lane合并。CSV/TSV继续支持，增加真实xlsx读取和列映射、可编辑样本表及模板导出。真正的二进制xls/xlsx不能冒充文本表；当前只兼容文本“.xls”。
3. **选择研究设计。** 独立/配对/复杂设计；分组和协变量用列选择器，连续/分类类型可确认；显示独立动物数、各组n、批次交叉表、设计矩阵和所问假设。
4. **选择比较和参考。** 自动生成比较、正负方向预览；参考资源按assembly/release管理，将transcripts、GTF/tx2gene、可选genome作为有版本的组合。初次准备完成后可复用，不要求每次重新找文件。
5. **启动前概览。** 必需项、可选项、将跳过的模块、参考匹配/注释覆盖、磁盘空间及内存预算。可以提出资源建议，但不能未测量就保证用户32GB能完成特定全基因组索引。

数据留在本机，FASTQ使用路径引用即可，不应让新手为了分析先在浏览器上传几十GB文件。参考和基因集可离线导入；下载受限不能阻止本地已有资源的使用。

### 6.2 建议的结果与导出流程

- 保留技术文件的同时，在界面默认显示“先看什么”：质量/样本 → 对比概览 → DEG双方案 → 通路 → 候选/表型 → 报告。
- 输出卡片按用途分为“完整统计表、候选名单、表达矩阵、富集/活性、图形、方法与运行记录”。显示中文解释、文件名、单位和能否作为下游输入；原始文件列表放高级区。
- 比较名、方向、active profile和分析阶段写进文件名/索引与方法页；根目录首个比较兼容副本应有清楚标识，避免被当作所有比较汇总。
- 表格导出增加多sheet xlsx和跨比较汇总；基因ID强制文本、保留完整数值；Excel显著性格式只是展示，不改变数据。
- 增加“分析结果包”和“可复现项目备份”区别。当前ZIP明确排除qc/quant/index/reference，**不是**包含全部输入和缓存的完整备份。未来完整备份应展示所需空间并允许路径清单方式引用大型原始数据。
- 报告方法应按实际入口生成。当前完整counts/salmon任务使用通用主线方法段，仍可能展示“Salmon selective alignment”措辞；建议按实际执行步骤细分，不把导入描述成重新运行。
- 路径与内容指纹、软件版本、镜像digest、数据库时间/物种/ID来源都应可下载。当前基础镜像digest和vendor校验值得保留，但apt软件和间接Python依赖并未全部锁定，重新构建不保证字节一致；更新缓存也应纳入工具/参考版本。

## 7. 高级信息挖掘的合理扩展顺序

| 方向 | 推荐参考 | 适用前提与输出 |
|---|---|---|
| 跨比较通路效应与核心基因 | clusterProfiler/enrichplot、fgsea、iSEE | 优先完成。提供NES/效应矩阵、共享leading-edge、基因集重叠解释；同源数据库和同一模型口径。 |
| 通路足迹与TF活性增强 | [PROGENy](https://github.com/saezlab/progeny)、decoupleR | 将基于响应基因的足迹和ssGSEA的成员打分区分；核对小鼠先验来源/同源转换、覆盖、权重与版本。PROGENy项目说明原始人类模型和小鼠扩展，[Nature Communications论文](https://doi.org/10.1038/s41467-017-02391-6)。 |
| WGCNA稳定性/独立队列保留 | WGCNA | 增加样本离群/批次敏感性、可配置稳健相关、重抽样稳定性与模块preservation；每一步需要足够独立样本。小样本时不靠放宽门槛伪造稳定网络。 |
| 组织细胞组成解释 | [MuSiC](https://github.com/xuranw/MuSiC) | [Nature Communications 2019](https://doi.org/10.1038/s41467-018-08023-x)；需要匹配组织、物种和质量合格的单细胞参考。组成估计帮助解释bulk变化，但分数不自动等于真实细胞比例，不能默认套人类免疫参考到小鼠。 |
| 交互、时间、重复测量 | DESeq2、dream、ExploreModelMatrix | 先补实验设计UI，再增加经验证的公式/contrast与模型诊断；输出效应和区间，不能仅多画几组小提琴图。 |
| 富集稳健性对照 | [limma camera](https://github.com/bioc/limma/blob/RELEASE_3_22/man/camera.Rd) | 相邻方法可帮助评估基因间相关带来的差异，但需要合适的表达/voom模型、设计和统计假设，不能把DESeq2的padj当作camera原始输入。不是要求所有方法必须显著才允许结论。 |
| 剪接/转录本使用等 | 专门的BAM或转录本统计工作流 | 另立分析分支与输入规范。现有gene-level Salmon/DESeq2输出不足以直接提供剪接、变异或新转录本结论。 |

不优先增加小样本机器学习“诊断模型”、无独立验证的风险评分或自动“核心致病基因”排序。若后续确有预测任务，应另做样本量、嵌套交叉验证、数据泄漏防范和独立验证设计。

## 8. 建议实施与验收顺序

| 顺序 | 交付 | 可检查的验收标准 |
|---|---|---|
| 第一批 | 设计向导、subject误用防护、原始/调整表型口径、覆盖率与方法标识 | 配对漏填不能静默运行；完全混杂明确拒绝；已知混杂的合成数据可区分原始与调整关联；比较方向/方案/单位始终可见。 |
| 第二批 | 多组比较生成、LRT、模型/NA诊断、DE检查点 | 3组/4组自动生成3/6个无重复对比；LRT与官方DESeq2逐项对照；可选模块失败不丢失已完成DE；模型或软件改变正确使缓存失效。 |
| 第三批 | 跨图选择、多层热图、UpSet、跨比较基因/通路效应图、批量样式与导出 | 同一选择在表格、热图和基因卡片一致；不同方案不会混合；长标签/密集样本不溢出；SVG/PDF/PNG/数据遵循同一配置。 |
| 第四批 | 参考/资源管理、xlsx与样本向导、PROGENy、网络稳定性、可选细胞组成 | 离线资源可用；来源/物种/版本/ID覆盖明确；没有适用输入时跳过而不造结果；高级分析有独立验证数据与清楚边界。 |

所有更新继续保持两套用户指定的DEG筛选，不覆盖历史统计；网页教程与Word同步。实施前后用相同合成基准与官方函数核对数值，再用适用的公开真实数据做外部复核；公开数据下载失败不应阻断程序开发或本地已有输入的分析。
