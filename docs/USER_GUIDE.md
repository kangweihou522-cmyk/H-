# 小鼠 bulk RNA-seq 工作站完整使用说明、方法原理与参考

说明书版本：2026.10.26

适用程序：小鼠 bulk RNA-seq 工作站 2.0，Windows + WSL2 Ubuntu 24.04 或原生 Ubuntu。本文对应当前仓库实际实现；每个任务的真实设计、软件版本与完成状态，以该任务保存的参数、R-session、日志和报告为准。说明书编写日期：2026-10-26。

本说明既是操作手册，也是分析方法和结果解读索引。示例样本名、数值、基因和路径仅用于说明格式，不代表实验发现。工作站不会把演示注释或合成数据当作真实小鼠生物学证据。

## 1. 工作站能完成什么

本工作站面向小鼠组织、细胞或其他混合样本的 bulk RNA-seq。常规主线是：样本与实验设计 → 原始读段预处理 → Salmon 转录本定量 → tximport 基因汇总 → DESeq2 差异分析 → GO/KEGG/其他 GMT 富集 → 单样本通路和调控活性 → 共表达与表型探索 → 图表编辑、报告及结果归档。

分析不要求先下载官方测试数据。可以直接使用自己的 FASTQ、已有 Salmon 定量或原始基因计数；合成演示仅用于验证程序是否能执行。首次安装依赖和获取参考文件通常需要联网，准备好镜像、参考、真实注释及基因集后，核心统计和网页操作可离线完成。

| 输入入口 | 从哪一步开始 | 能提供的结果 |
|---|---|---|
| FASTQ | fastp 预处理、Salmon 定量 | 测序/定量质控、表达、差异、富集和满足条件的挖掘模块 |
| Salmon quant.sf | tximport 基因汇总 | 保留下来的定量元数据、表达与后续统计；不重做原始 reads 质控 |
| 原始基因计数 | DESeq2 建模 | 表达与后续统计；没有原始 reads 质控和转录本长度信息 |

界面按“样本与实验设计、测序与定量质控、表达与样本诊断、差异基因、GO/KEGG/GSEA、通路与关键基因、共表达与表型、报告与导出”组织，最后提供本使用说明。按顺序查看结果比直接寻找显著基因更可靠：输入或样本质量有问题时，后面的精美图表不能纠正前面的偏差。

## 2. 本地部署、启动与保存位置

### 2.1 Windows、Ubuntu、Docker 分别负责什么

Windows 提供文件管理和浏览器；WSL2 Ubuntu 承载 Linux 分析环境；Docker 将 Salmon、R 及相关分析包放在相对固定的计算环境中。网页是控制界面，真正的分析在本机后台执行。仅打开 Ubuntu 窗口不会自动启动网站；关闭浏览器也不等于停止计算。

对于已有 Windows + WSL2 的电脑，推荐继续使用该组合，不需要为了本工作站改成纯 Ubuntu。Java 和 Nextflow 可以保留，但本工作站目前由本地任务管理器调度，并不要求通过 Nextflow 执行。

### 2.2 首次部署与每天打开

首次下载 GitHub 源码 ZIP 后全部解压，在解压目录双击 `部署到D盘.bat`。默认目标为 `D:\转录组分析工作站`，默认网页端口为 `8767`。完成部署后，每天双击目标目录里的 `启动工作站.bat` 或桌面快捷方式，即可调用 Ubuntu 后台并打开浏览器。

也可以打开 Ubuntu，执行部署时生成的命令：

```bash
bash ~/启动转录组工作站.sh
```

如果同名脚本已被用户修改，安装器可能保留旧脚本并生成 `~/启动D盘转录组工作站.sh`，以部署输出的实际路径为准。采用 Windows D 盘入口时，浏览器地址通常是 `http://localhost:8767`；直接在 Ubuntu 源码目录运行 `bash launch.sh` 时默认端口为 `8765`，不要混淆两个运行位置。

分析期间保持工作站终端运行、电脑不休眠；使用 Docker Desktop 时保持其引擎运行，使用 Ubuntu Docker Engine 时保证其服务可用。正常关闭先确认任务结束，再在工作站终端按 Ctrl+C。

### 2.3 哪些东西在 D 盘

| 位置 | 保存内容 |
|---|---|
| `D:\转录组分析工作站\程序\版本摘要\` | 各版本源码；更新通常生成新目录 |
| `原始数据\` | 用户自行放入的 FASTQ、样本表、计数等；不会自动搬移已有文件 |
| `参考文件\` | 用户准备的 FASTA、GTF、tx2gene、注释和 GMT |
| `分析结果\任务ID\` | 每次分析的数据、图表、报告、参数及日志 |
| `导出报告\` | 可手动选择的导出目标；浏览器默认下载目录不被强制改写 |
| `配置\workstation.json` | 生效的程序目录、发行版、端口和 Python 环境路径 |
| `日志\server.log` | 后台服务日志 |
| `日志\wsl-launcher.log` | Windows 到 WSL 的调用与错误记录 |

Python 虚拟环境保存在 Ubuntu 文件系统，Docker 镜像由原来使用的 Docker 管理。源码和结果位于 D 盘，不代表 Ubuntu 虚拟磁盘和 Docker 存储也被自动搬到 D 盘。大型 FASTQ 经 `/mnt/d` 访问可能比 Ubuntu 内部文件系统慢；有空间和性能需求时，可把原始数据放到 Ubuntu 目录并在样本表中填写实际路径。

### 2.4 网络与代理的已知情况

Windows 能下载文件，不代表 Ubuntu 或 Docker 守护进程使用了相同代理。在 WSL NAT 模式下，Windows 的 `127.0.0.1` 代理地址不能直接当作 Ubuntu 的本地代理。`wsl: 检测到 localhost 代理配置，但未镜像到 WSL` 是代理提示，需结合随后命令的退出码和具体网络报错判断。

如果使用 Ubuntu 原生 Docker Engine，构建第一步连接 `registry-1.docker.io/v2/` 超时，而 Python 包已安装成功，可按部署文档检查引擎代理。现有辅助脚本为：

```bash
sudo python3 scripts/configure_wsl_docker_proxy.py --port 10808
```

仅在实际使用对应 HTTP 代理端口时运行。脚本先检查通过 Windows 网关访问代理和 Docker Hub，再检查活动容器，备份并合并 Docker 配置，验证、重启和核对生效情况，最后拉取锁定摘要的基础镜像。它适用于本机 Ubuntu Docker Engine，不适用于 Docker Desktop 的代理设置；有运行容器时会要求先正常结束相关任务。WSL 网关重启后可能变化，需要时重新运行。详细步骤见仓库 `docs/WINDOWS_D_DRIVE.md`。

## 3. 开始前先明确实验问题

先写清楚组织、品系、性别、年龄、处理、时间、建库方式和比较方向。小鼠不同品系、参考组装及注释版本可能影响可定量基因和映射情况，不能只记录“mouse”。记录组装如 GRCm38/GRCm39、注释来源与 release，并保证参考资源相互匹配。

程序至少要求两个组、总样本数至少 4、每组至少 2 个生物学重复，推荐每组至少 3 个独立重复。两个重复只满足最低计算门槛，不能替代足够的实验设计。技术重复、同一样本不同测序 lane，以及同一动物重复取样，都不应直接当作新的独立生物学重复。

bulk 数据反映混合细胞的平均信号：免疫细胞比例变化、组织取样差异或 RNA 质量变化都可能改变表达。差异不能自动解释为某一种细胞内的调控改变。DESeq2 的常用归一化依赖表达组成的可用信息；存在全局 RNA 产量大幅变化时，应结合实验设计和外部标准，不从常规归一化结果推断绝对转录总量。

对于同一样本不同 lane，先按 R1/R2 配对正确合并原始文件，再建立一个样本条目。配对或重复测量实验应在样本表明确 `subject`；普通配对差异分析可把 subject 作为协变量，但更复杂的纵向或随机效应模型不是当前内置模型。

## 4. 准备样本表和三类输入

### 4.1 通用样本表规则

样本表为带表头的 UTF-8 CSV 或 TSV，必须包含 `sample`、`group`。样本名唯一，以字母或数字开头，最多 80 个字符，只使用字母、数字、点、下划线和短横线。分组名不能为空，界面填写的组名必须与表内完全一致。

输入框接受 Ubuntu 路径，也可填写 `D:\data\...` 或本地 file URI，后两者会转换为 `/mnt/d/...`。路径必须指向当前 Ubuntu 能读取的本机文件；`file:///...` 不是上传文件的机制。真实二进制 Excel 工作簿需先另存 CSV/TSV；公司命名为 `.xls` 的纯文本 TSV 可直接读取。

FASTQ 样本表示例：

```csv
sample,group,fastq1,fastq2,batch,subject,disease_score
C1,Control,D:/data/C1_R1.fastq.gz,D:/data/C1_R2.fastq.gz,B1,M01,1.2
C2,Control,D:/data/C2_R1.fastq.gz,D:/data/C2_R2.fastq.gz,B2,M02,1.0
C3,Control,D:/data/C3_R1.fastq.gz,D:/data/C3_R2.fastq.gz,B1,M03,1.1
T1,Treatment,D:/data/T1_R1.fastq.gz,D:/data/T1_R2.fastq.gz,B1,M04,3.1
T2,Treatment,D:/data/T2_R1.fastq.gz,D:/data/T2_R2.fastq.gz,B2,M05,2.8
T3,Treatment,D:/data/T3_R1.fastq.gz,D:/data/T3_R2.fastq.gz,B2,M06,3.2
```

这里的 `subject` 是独立动物标识，不代表应把每个唯一 subject 都加入这个非配对示例的模型；把每个样本独有的标识当作协变量会消耗全部自由度。样本表可以保存描述性列，只有被选中的协变量才进入 DESeq2 模型。

### 4.2 FASTQ 入口

选择“FASTQ → fastp → Salmon”，样本表需有 `fastq1`；双端增加且完整填写 `fastq2`。同一任务不能混合单端和双端。准备匹配版本的小鼠转录本 FASTA 及 tx2gene。FASTA 应为转录本序列，例如同版本 Ensembl cDNA 加 ncRNA；不能把完整基因组 FASTA 填到转录本输入框。

从匹配 GTF 生成 tx2gene 的命令为：

```bash
.venv/bin/python scripts/prepare_reference.py /data/mouse.gtf.gz /data/tx2gene.tsv
```

对于 D 盘部署，命令中的 Python 应替换为配置记录的实际虚拟环境解释器，或使用已有可运行该脚本的 Python 3。tx2gene 为带表头的两列 TSV：

```tsv
transcript	gene
ENSMUST00000000001	ENSMUSG00000000001
ENSMUST00000000002	ENSMUSG00000000001
```

转录本去掉版本后必须唯一；转录本对应多个不同基因时不能直接用于汇总。tximport 忽略转录本版本后缀，基因注释映射也去除基因版本后缀，但这不解决不同参考 release 的实质差异。

### 4.3 已有 Salmon 入口

样本表包含 `sample,group,quant`，`quant` 指向每个样本的 `quant.sf`，另提供相匹配的 tx2gene。标准字段为 `Name、Length、EffectiveLength、TPM、NumReads`；Name 不重复，NumReads 为非负有限值。

最好保留 quant.sf 旁边的 `aux_info/meta_info.json`、`lib_format_counts.json` 和 `cmd_info.json`，工作站可复制相应元数据。只有 quant.sf 时仍可汇总表达，但缺失的映射率或方向诊断不会被补造，也不会重新产生 fastp 原始质量报告。

### 4.4 原始计数入口

样本表只需基础设计列，另选原始基因计数矩阵。第一列为唯一基因 ID，其他样本列至少包含所有入选的 `sample`，矩阵内容必须是非负整数，不能含 NA、无穷值或总计数为零的样本。示例：

```csv
gene,C1,C2,C3,T1,T2,T3
ENSMUSG00000000001,110,95,120,300,280,310
ENSMUSG00000000002,40,55,48,42,50,45
```

不要把 TPM、FPKM、RPKM、log 表达或已经归一化的数值四舍五入后输入 DESeq2。整数检查只检查格式，不会识别被人为取整的错误数据；数据来源仍需用户确认。原始计数本身缺少转录本长度，当前工作站不会从它反推 TPM。

## 5. 配置实验设计并通过启动检查

进入“样本与实验设计”，先填写路径并点击“读取样本”，检查分组、数量和是否误选文件。默认保留全部样本；如需排除，在此明确记录。排除后会重新检查每组重复数和设计是否可估计，不会因为某样本在 PCA 中偏离而自动删除。

每个比较填写唯一比较名、对照组和处理组，最多 30 个比较。所有样本共同拟合一个 DESeq2 模型，分别提取所需比较，方向始终是“处理组 / 对照组”。例如 `Treatment_vs_Control` 中 log2FC 为正表示处理组高于对照组，反向比较会翻转方向。

协变量输入逗号分隔的样本列名，例如 `batch,subject`；连续变量写为 `age:numeric`。当样本表存在多水平 batch 且未指定其他设置时，程序默认加入 batch；显式输入 `none` 表示不调整协变量。应以预检查返回的实际协变量和保存的设计为准。

| 设计例子 | 含义与适用情况 |
|---|---|
| `~ group` | 无其他需调整因素的独立两组或多组比较 |
| `~ batch + group` | 批次与组别存在可区分的信息，可估计组别条件效应 |
| `~ subject + group` | 配对对象在不同组均有观测，且设计有残差自由度 |
| `~ age + group` | age 被明确声明为 numeric，估计年龄调整后的组别效应 |

若所有对照都在批次 B1、所有处理都在 B2，则批次与处理完全混杂，任何软件都不能仅靠统计把两者分开。设计矩阵不满秩、选中协变量含空值、数值协变量无变化、没有残差自由度时，程序拒绝启动。当前模型是加性设计，不提供任意交互项、随机效应或专门时间序列模型。

在“资源与启动检查”点击“检查样本与参数”，核对实际样本、比较、协变量、可用输出空间和提示，再点击“启动完整分析”。默认 4 线程、容器内存上限 24GB，可在界面配置 1–8 线程、4–28GB；32GB 电脑还需为 Windows、WSL 和浏览器保留内存。任务串行执行，提交后可以切换项目查看状态。

## 6. fastp：接头和质量预处理

FASTQ 默认使用 fastp，质量碱基阈值为 Phred 20，保留 read 最短长度为 30；可在界面调整。双端使用接头检测，输出清洁 R1/R2、JSON 和 HTML。对于已经规范预处理的 reads，可以选择跳过 fastp，报告会标记该步骤未运行。

Phred 质量与错误概率的关系是 Q = -10 log10(Perror)，因此 Q20 对应约 1% 单碱基错误概率，Q30 对应约 0.1%。Q20/Q30 比例用于质量概览，不等于整体比对准确率，也不能保证某个基因定量可靠。不要把公司针对特定读长设置的最短 150bp 参数直接用于 100bp 数据。

建议依次查看原始与保留 reads、保留比例、接头量、每碱基质量和 A/T/G/C 曲线、GC 分布，再看样本间一致性。双端 fastp 汇总常按 read 端计数，而 Salmon 可按 fragment 计数；比较不同工具的数量前先确认单位。

明显低质量、接头多、保留比例很低或 GC 异常时，先核对建库和测序信息、原始文件配对、物种及污染可能性。当前工作站不内置完整污染分类或 rRNA 去除诊断，不应仅凭一条曲线确定污染来源。

## 7. Salmon 与 tximport：从转录本到基因

Salmon 将 reads/fragment 与转录本参考进行匹配并估计转录本丰度，工作站启用 selective alignment、序列偏好校正和 GC 偏好校正。每个任务建立转录本索引，随后逐样本量化。当前索引未集成完整基因组 decoy；与使用 decoy-aware 索引或基因组比对的其他流程比较时应记录这一差别。[1–3]

文库类型默认 `A` 自动判断；双端支持 IU/ISR/ISF，单端支持 U/SR/SF。链特异性应根据建库试剂说明和 Salmon 实际方向记录共同确认，不能只凭是否双端判断。查看映射率、兼容方向比例、处理 fragments、估计数量与样本间差异。

程序发现 Salmon 映射率低于 60% 会提出检查提示。这是定位问题的起点，不是所有实验通用的合格线。先检查转录本参考物种/组装/release、R1/R2 匹配、文库方向、数据预处理和潜在污染，再决定是否重新量化。

Salmon 的 fragment mapping rate 不是 HISAT2/STAR 的“基因组 unique mapping rate”，也不能由此推导外显子/内含子比例、基因体覆盖或 IGV 轨迹。Salmon 主线不会生成用于这些分析的基因组 BAM。

tximport 根据 tx2gene 把转录本估计数量汇总到基因，并把长度相关信息交给 DESeq2 的标准入口。工作站保存 TPM 用于表达展示，但统计仍通过 tximport 的估计计数和相应归一化信息建立模型；不能把 TPM 矩阵当作原始计数替换这个过程。

## 8. 表达诊断：先检查样本结构

低表达过滤默认要求至少 2 个样本 count ≥ 10；可根据实验和比较设计调整。过滤后少于 10 个基因时停止。这个过滤定义在全体保留样本上，不是只保留差异显著基因。

DESeq2 根据计数分布估计归一化因素和离散度。工作站保存 normalized count，以及采用 `blind=FALSE` 的方差稳定变换 VST。归一化计数用于同一基因的样本表达比较，VST 用于 PCA、距离、相关性和后续挖掘；两者都不是绝对分子数。不同基因间仅比较 normalized count 大小会受到基因长度等因素影响。[4]

PCA 默认取 VST 中变异最高的最多 500 个基因，对样本进行中心化、不逐基因标准差缩放的主成分分析。PC1/PC2 百分比是对应主成分解释的变异比例；分组接近或分离是总体表达结构，不是差异检验的 P 值。

| 图表或指标 | 主要看什么 | 不宜单独得出的结论 |
|---|---|---|
| 文库大小/检测基因数 | 深度、检测范围及离群样本 | 文库更大就一定生物学表达更高 |
| 表达分布 | VST 或 log2(normalized+1)/log2(TPM+1) 的样本分布 | 箱线图相似就消除了批次效应 |
| PCA | 组别、批次、对象和潜在离群结构 | 离群样本必须删除；组间不分离就没有差异 |
| VST Pearson 相关 | 样本整体模式的一致程度 | 相关高就不存在真实组间效应 |
| VST 样本距离 | 多基因表达距离及聚类结构 | 聚在一起就证明样本标签正确 |
| 离散度/Cook 最大值 | 表达均值、变异与异常影响 | 单项数值可以替代完整统计诊断 |

VST、normalized 和 TPM 是不同尺度，图注应写明所用尺度。仅 Salmon/FASTQ 入口保存 TPM。VST 使用了设计信息，不是“已经去除 batch 的表达矩阵”；批次调整体现在模型中，PCA 仍可能看到批次结构。

## 9. DESeq2：模型、效应与结果列

DESeq2 用负二项广义线性模型描述基因计数及其样本间变异，估计归一化信息、离散度和设计中的组别/协变量效应，再以 Wald 检验评估指定比较。批次或年龄调整后的 log2FC 是模型条件效应，不一定等于结果表简单分组均值之比。[4]

| 列 | 含义 |
|---|---|
| `gene` / `symbol` | 原始基因 ID / 真实注释得到的符号；缺失映射保留为空 |
| `baseMean` | 归一化计数的平均表达水平 |
| `log2FoldChange` | 当前使用的处理/对照 log2 倍数；是否收缩以任务设置为准 |
| `log2FoldChange_raw` | 保留的未收缩 LFC |
| `lfcSE` | 当前 LFC 的标准误，选收缩时随之更新 |
| `stat` | 原始 Wald 统计量；用于全基因 GSEA 排名 |
| `pvalue` / `padj` | 单基因原始 P 值 / DESeq2 的 BH 校正 P 值 |
| `mean_control` / `mean_treatment` | 两组 normalized count 的简单平均，属于描述统计 |
| `significant` / `direction` | 当前双方案筛选定义下是否入选及 up/down/not_significant |

默认不做额外 LFC 收缩，也可选择 DESeq2 normal 先验收缩。收缩可缓和低计数和不稳定大效应，便于展示和排序，但不是增加显著性证据；本工作站保留原始 LFC，原始 Wald stat 仍用于 GSEA。选择收缩可能改变以 LFC 阈值定义的 DEG 数量，方法部分必须记录。

`padj=NA` 不应当成 1 或 0。低表达独立过滤、异常计数处理或不可估计情况都可能产生 NA，需结合基因均值、诊断和 DESeq2 日志判断。当前筛选不把缺失 P 值认定为显著。

## 10. 两种差异筛选方案与后续选择

每个比较同时计算、保存以下两套结果，比较符号均为严格小于/大于：

| 方案 | 显著性要求 | 效应要求 | 使用定位 |
|---|---|---|---|
| 方案一 nominal | `pvalue < 0.05` | `abs(log2FoldChange) > 0` | 探索性结果，与部分公司报告筛选口径接近 |
| 方案二 strict | `padj < 0.05` | `abs(log2FoldChange) > 1` | BH 多重检验校正并加入效应大小门槛 |

log2FC 为 1 对应约 2 倍，-1 对应约 0.5 倍；由于条件是 `>1`，恰好 ±1 不入选 strict。`>0` 几乎不设置实际效应大小门槛，但恰好 0 仍排除。

“差异基因”页面先显示两种方案上调、下调和总数，再选择后续使用哪一种。nominal 不控制跨基因检验的 FDR，不能因为 DEG 更多就说证据更强；即使富集阶段再做 BH，也不会补救前一步未控制的基因选择错误。[5]

| 切换 active profile 会更新 | 保持全基因分析，不随之重新计算 |
|---|---|
| 当前 DEG 清单和数量、完整差异表的 significant/direction | 同一 DESeq2 拟合、原始统计量 |
| all/up/down 的 GO 与 GMT ORA | Wald stat 排名的 GO/GMT GSEA |
| 火山图、MA、差异热图和相应富集展示 | ssGSEA 活性、调控活性、WGCNA 及原有表型关联 |

两个方案共用模型及所选 LFC 列，另一方案不会因切换而删除。strict 没有 DEG 时，应如实报告；可以选择 nominal 做明确标注的探索，或查看不依赖 DEG 阈值的 GSEA，而不是把非显著基因改称 DEG。

图表编辑器里的 P 列、P 阈值和 LFC 显示门槛只控制当前图的显示，不会改写 active profile 的 DEG 和 ORA。需要改变下游方案，应使用两套方案卡片上的选择操作。

## 11. 准备真实小鼠注释与基因集

### 11.1 注释与 ID 映射

默认数据库路径使用小鼠 `org.Mm.eg.db`，输入基因按 ENSEMBL 映射到 Symbol、Entrez、GOALL 和 ONTOLOGYALL；GOALL 包含祖先术语。程序会去掉基因版本后缀后匹配。Symbol 或 Entrez 为主的计数 ID 不应直接假定能走这个 ENSEMBL 数据库入口，需要先可靠转换，或准备以实际输入 ID 为 gene 的离线注释并核对基因集 ID。[6]

离线注释 TSV 必须包含下面五列表头，一个 gene–GO 关系一行，同一基因可占多行：

```tsv
gene	symbol	entrez	go	ontology
ENSMUSG00000000001	Example1	12345	GO:0008150	BP
ENSMUSG00000000001	Example1	12345	GO:0005575	CC
```

此例仅演示列格式，不是可用于真实分析的基因关联。ontology 只能为 BP、CC、MF 或空，GO 编号为 `GO:` 后接七位数字。自定义注释应包含需要的祖先 GO 关系；当前导入不会自动把仅有的直接 GO 关联扩展成 GOALL。

有网络且满足依赖条件时，可从源码目录执行 `bash scripts/install_mouse_annotation.sh` 安装真实小鼠注释，也可从已有可信注释包导出。没有真实 GO 注释时会明确提示 GO 未运行，而不是产生模拟结果。ID 映射率低时先修正 ID/物种/版本，再解释富集。

### 11.2 GMT：KEGG、Reactome、Hallmark 与自定义集合

GMT 每行是一个基因集，使用 TAB 分隔：`通路ID<TAB>描述<TAB>基因1<TAB>基因2...`。通路 ID 不重复，文件至少有一个有效集合。界面选择 ENSEMBL、ENTREZID 或 SYMBOL，必须与 GMT 成员一致；不同集合混用 ID 类型会造成不完整匹配。

小鼠 KEGG 物种代码为 mmu。现有获取命令为：

```bash
.venv/bin/python scripts/fetch_kegg.py /data/mouse_kegg.gmt
```

输出 GMT 和对应 `.source.json`，记录 mmu、ENTREZID、来源、获取时间和 SHA256；界面选 Entrez ID。网络不可用时可以导入已有合法 GMT，KEGG 下载不是制作或打开工作站的前提。遵守 KEGG、MSigDB、Reactome 和其他资源的使用条件，保存资源版本/日期。[7,8]

人类通路或调控网络不能仅把 species 标签改成 mouse 后使用。跨物种转换要记录同源关系版本、映射数量、一对多处理和丢失基因，再确认小鼠 ID 一致。无任何测试基因匹配的 GMT 会报错，避免输出没有实际基因交集的结果。

## 12. GO/KEGG ORA：差异列表是否过度代表某功能

ORA 是对选出的 DEG 做过度代表分析。对每个集合，记 N 为合格背景基因数、n 为背景内被选基因数、K 为该集合在背景内的基因数、k 为选中且落入该集合的数量，通过超几何分布计算观察到至少 k 个命中的尾概率。

本工作站使用本次低表达过滤后、具有有限 Wald stat、ID 可用且能映射到对应集合的基因作为候选背景，再与各数据库集合的成员取交集。它不是“整个小鼠基因组”背景。一个 ID 对应多个表达行时优先保留平均表达较高的一行，不按最显著 P 值挑行。不同输入、映射和参考版本会改变背景及统计结果。

GO 分 BP（生物过程）、CC（细胞组分）、MF（分子功能）；GMT 对应输入的 KEGG/Reactome/Hallmark 等集合。每种数据库都分别生成 all、up、down 三个方向。up/down 表示成员基因的差异方向，不等于整个通路已被功能性激活或抑制。

默认基因集大小范围为 15–500，以实际可测试背景中的集合大小为准，可在界面调整。每个方向内按 GO ontology 对所有可测试集合进行 BH，GMT 在其方向内统一校正；零命中的可测试集合仍参与 ORA 的检验集合。

| 字段 | 解读 |
|---|---|
| overlap / genes | 命中数量及实际成员 ID |
| size / background / selected | 集合在背景内大小 / 背景总数 / 选中总数 |
| gene_ratio | overlap / selected；不是对全部表达基因的比例 |
| rich_factor | overlap / size；与 gene_ratio 分母不同 |
| pvalue / padj | 超几何检验原始 P / 该 ORA 校正范围内 BH P |
| direction / ontology | all/up/down，以及 BP/CC/MF/GMT |

阅读时先确认背景、ID 和方向，再看 padj、命中数量和具体基因，不只看最上方几个通路名。GO 术语具有层级关系，显著术语常共享成员；多个相近术语并不代表多项独立生物学证据。ORA 当前为本项目明确实现的超几何检验，不是直接调用 clusterProfiler。[6,9]

## 13. GSEA：不先截断 DEG 的排序富集

GSEA 把可用基因按原始 Wald stat 从高到低排序，检验某基因集是否更集中在排序一端。它可以发现多个成员各有适中变化、但 DEG 截断后不明显的协调信号。因此本工作站的 GSEA 不先按 nominal 或 strict 筛 DEG，切换 DEG 方案不会重跑其排序分析。[10,11]

实际计算使用 fgsea multilevel，默认集合大小 15–500、随机种子 42、初步排列 `nPermSimple=1000`、精度下限 `eps=1e-10`、单进程。界面的“fgsea 初步排列次数”不是传统固定排列次数 GSEA 的完整定义；增大它不保证每个集合都能计算出更小 P 值。

| 字段 | 解读 |
|---|---|
| ES | 基于排序命中位置的运行富集分数 |
| NES | 对集合大小等因素标准化后的富集分数 |
| NES 正/负 | 相对富集于处理高于对照的一端 / 对照高于处理的一端 |
| padj | fgsea 此次调用结果的 BH 校正值 |
| log2err | P 值估计不确定性指标之一，需结合日志理解 |
| leadingEdge | 对富集峰贡献突出的核心成员，不等于已验证关键驱动基因 |

GO GSEA 的 BH 范围是该次 GO fgsea 调用的可测试集合，并不是 ORA 那样按 BP/CC/MF 分别校正；GMT GSEA 是另一调用。不要把不同检验集合的 padj 当成可直接互换的数值。样本量、排序信号与基因集结构会影响结果。

富集曲线显示运行 ES、命中位置及 Wald rank track；点击结果表通路可查看。统计达到精度下限、正负分布不平衡或无法可靠估计时，日志会提示，NA 保持 NA，不填造 P 值。单条曲线外观很好看并不替代多重检验、重复样本和成员表达检查。

## 14. GSVA 包的 ssGSEA：每个样本的通路活性

勾选“GSVA ssGSEA 通路活性 + limma 比较”并提供匹配 GMT 后，工作站对 VST 表达计算 ssGSEA 单样本分数，采用 GSVA 包的 `method='ssgsea'` 和归一化选项。这里使用的是 ssGSEA 方法，不是 GSVA 包中另一个同名 GSVA 算法。[12]

分数表示样本表达对基因集的相对支持程度，不是酶活、代谢通量或直接的绝对通路活性测量。集合大小、成员、预处理和物种映射都会影响分数，不宜跨不一致数据库版本直接比较绝对数值。

组间比较用 limma 对分数建立线性模型，使用与 DESeq2 相同的协变量和比较方向，再进行经验贝叶斯检验。结果中的 `logFC` 是通路分数差，不是“通路基因表达的 log2 倍数”；显著性通常看 `adj.P.Val`。[13]

还可输出通路–表型 Pearson 相关及表内 BH 校正，这个相关步骤没有像上述 limma 比较一样调整协变量。把差异活性和表型相关视作两个不同问题，不把未调整的相关结果称为批次调整后的效应。

## 15. 调控因子活性：decoupleR ULM

提供真实、已转换为小鼠 ID 的 `source,target,weight` TSV 后，可按先验靶基因表达估计调控因子活动。source 通常为转录因子或其他调控者，target 为基因，weight 表示有符号的调控关系权重。每个 source–target 对唯一，权重为有限非零数；在当前表达数据中至少检测到 5 个靶基因的 source 才进入计算。[14,15]

工作站使用 decoupleR ULM，在 VST 表达和先验网络上估计逐样本分数，再用同设计的 limma 比较各组。网络可来自有权使用的 DoRothEA mouse 等资源，需记录来源、版本、证据等级和 target ID 类型。

调控活性由靶基因模式推断，不等于 source 本身 mRNA 的升降，也不直接测量蛋白量、磷酸化或 DNA 结合。网络缺失或每个调控者匹配靶基因不足时，会标记跳过，不生成假的调控信号。候选调控者仍需结合实验和其他证据验证。

## 16. WGCNA 共表达模块与表型关联

WGCNA 研究样本间共同变化的基因模块。当前要求至少 15 个独立样本，推荐 20 个以上；如果样本表有重复 subject，模块会跳过，不能把同一动物的多个观测当作独立样本增加自由度。每组 3 个、总计 6 个样本仍可做差异和通路探索，但不满足当前 WGCNA 条件。[16]

默认选取 VST 变异最高的 3000 个基因，可调 100–6000；使用 signed 网络和 signed TOM，模块最小大小 30、每块最多 3000 基因、模块合并阈值 0.25、随机种子 42。先输出软阈值诊断，再使用所填 power 建网络，默认 power=6 仅是初始参数，程序不会自动以最佳拟合替换它。

软阈值诊断同时看 signed scale-free 拟合、slope 和平均连接度；不要只追求一个很高的 R²，导致网络连接过稀。小样本、明显批次结构或强离群可能制造不稳定模块。希望比较其他 power 或筛选基因数，应新建并完整记录一个分析任务。

| 结果 | 含义与使用 |
|---|---|
| 模块颜色/编号 | 共表达簇标识；颜色不是生物功能或稳定跨项目标签 |
| eigengene | 模块表达主要变化轴的概括，不是所有基因的算术均值 |
| kME | 基因与对应模块 eigengene 的相关型成员度 |
| module–trait | 模块 eigengene 与表型的当前未调整 Pearson 相关及表内 BH |
| hub 边表 | 候选高成员度子集的共表达邻接展示，不是完整网络或蛋白互作 |

当前网络边表取绝对 kME 较高的最多 120 个基因，保留 signed 邻接权重大于 0.2 的边，并按权重最多输出 500 条。这个可视化子集不是“所有模块的完整调控网络”。灰色模块通常代表未明确分配的基因，不应像稳定模块那样直接作强解释。

分类表型按水平转成指示变量，数值表型保留其数值，group 自动纳入表型矩阵。模块–表型 Pearson 当前未调整批次、性别等协变量，分类的互补指示变量也不构成彼此独立的生物学验证。相关、共表达和网络中心性都不证明因果；复杂重复测量和混杂需进一步建模。

## 17. 关键基因探索：把统计结果落回每个样本

可在“通路与关键基因”输入 Ensembl ID、Symbol 或有真实映射的 Entrez ID 查询。匹配不区分大小写；基因必须仍在当前过滤保留的表达结果中。若一个查询匹配多行，应检查候选记录和注释，确认当前显示的实际 gene ID，而不是只依赖符号名称。

查看各重复的 normalized count、VST，以及有保存时的 TPM，同时核对 LFC、pvalue、padj 和组内变异。若总体差异由一个样本驱动，或某组几乎全为零，应比对原始计数、Cook 诊断和实验信息；不要只展示最符合预期的组均值。

候选基因建议结合：同方向、足够效应的差异证据；可信通路的 leading-edge；跨重复一致的表达；适用时模块成员度和调控先验；与实验表型或已有文献的一致性。DEG、leading-edge 和 hub 的交集只是候选优先级工具，不会自动转化为已验证的关键靶点。

基因与样本表数值列的关联使用探索性 Spearman。它未调整协变量，也未对“全基因 × 全表型”的搜索作整体校正，P 值不能直接用于高维因果结论。正式证据需要独立样本、预先定义的验证问题，以及合适的 qPCR、蛋白、功能干预等实验。

## 18. 图表交互、编辑与导出

### 18.1 火山图定位和标注

火山图横轴为当前使用的 log2FC，纵轴为 `-log10(pvalue)` 或 `-log10(padj)`，以轴标签为准。悬停查看基因与数值；点击感兴趣的点可进入基因信息；在“输入基因 ID / Symbol”中填写目标并点击“定位并标注”，可高亮该点。

可以拖动标签位置、调整自动标签数量或指定标签，点击“重绘并保存”保存设置和导出图。极小 P 值绘图时有数值下限处理，以避免无穷坐标；图上的上限并不代表原始统计 P 值被改写。缺失 P 值的基因可能不显示在当前火山图中。

### 18.2 不同图的用途

PCA 可选择存在且不同的主成分，并设置分组配色；热图可按当前 DEG 选择，也可手动填写基因，使用 VST 的逐行 Z-score 展示相对模式。行 Z-score 不能用来比较两个基因的绝对表达高低。严格方案无 DEG 时，不会自动换成全基因差异热图；手动基因热图应标为自选探索。

热图聚类可切换，当前层次聚类使用 average linkage；样本相关和距离图需结合其实际矩阵含义解读。GO/KEGG 点图的大小可代表命中数量，颜色对应校正 P 值，横轴显示 gene ratio；GSEA 条形图显示 NES。通路活性和 TF 活性热图显示的是各自分数，不要把它们的色标当作 log2FC。

### 18.3 编辑参数和可交付格式

通用编辑包含标题、坐标名称、宽高、字号、配色、显示数量等；火山图等还提供点大小/透明度、标签、范围和显示阈值。默认图宽 9、高 6 英寸、字号 11、300 dpi；界面允许的范围由程序校验。选择颜色时兼顾组别含义和色觉可辨性，整个报告保持同一组同一颜色。

SVG/PDF 适合矢量排版，PNG 适合常规展示；同时保留每图数据和 `.settings.json`，便于追溯图从何来。SVG/PDF 在外部软件编辑后要保存独立副本，不应把外部美化结果当作工作站自动生成的原始导出。重绘会更新本任务报告中的相应图。

图注至少说明物种/组织、组别、样本量、比较方向、表达尺度、检验和校正、筛选阈值及重要预处理。图形美观不能替代统计设计，也不应通过裁剪坐标或隐藏不一致重复误导读者。

## 19. 输出目录、报告和复现记录

D 盘部署结果位于 `分析结果/任务ID/`；直接源码运行时通常为 `runs/任务ID/`，也可由环境变量选择输出位置。网页“报告与导出”显示实际目录，不要仅凭旧版源码所在位置猜测结果保存处。

| 文件或目录 | 作用 |
|---|---|
| `config.json`、`samples.csv` | 实际参数、保留样本、协变量、比较和排除记录 |
| `input-provenance.json` | 部分输入的路径、大小、修改时间及源码版本记录 |
| `status.json`、`steps.json`、`log.txt`、`execution.jsonl` | 状态、分步信息、工具输出及执行命令 |
| `qc/样本/` | fastp 清洁 reads、JSON/HTML 和完成标记 |
| `index/`、`quant/样本/` | Salmon 索引、quant.sf 及对应元数据 |
| `raw-counts.csv`、`normalized.csv`、`vst.csv`、`tpm.csv` | 各尺度表达；TPM 仅相应入口存在 |
| `pca-scores.csv`、`pca-variance.csv`、`sample-qc.csv` | 样本表达与主成分诊断 |
| `gene-annotation.csv`、`gene-sets.json` | 本次实际映射及可用集合 |
| `contrasts/比较名/profiles/nominal或strict/` | 每个比较两套预计算差异及 ORA 结果 |
| `contrasts/比较名/active-profile.json` | 当前下游方案记录 |
| `gsea.csv`、`pathway_gsea.csv` 及排名表 | 当前比较全基因 GSEA 结果与排名 |
| `pathway-activity.csv`、`tf-activity.csv`、`wgcna-*.csv` | 满足输入和样本条件时的挖掘结果 |
| `mining-status.json`、`analysis-summary.json` | 模块执行/跳过理由、设计及统计汇总 |
| `dds.rds`、`R-session.txt` | DESeq2 对象和真实 R 会话版本 |
| `*.svg`、`*.pdf`、`*.png`、`*.settings.json` | 图、绘图设置及对应数据导出 |
| `report.html`、`analysis-results.zip` | 自包含图像的离线报告与结果包 |

`input-provenance.json` 不是所有大型输入的完整内容 SHA256 清单；续跑的输入标记主要依赖路径、文件大小和修改时间。正式发表或长期归档时，应另行保存完整原始数据校验值、参考/GMT来源及版本、镜像标识与样本采集记录。

结果 ZIP 当前排除 `qc/`、`quant/`、`index/` 目录，也不包含外部原始 FASTQ。因此 ZIP 并不是整个实验的全部数据备份；fastp 独立报告、清洁 reads、完整 Salmon 目录和索引如需留存，要单独备份。报告内嵌导出图，可以离线阅读，但不会自动包含尚未执行模块的结果。

## 20. 重试、取消、更新与数据保留

任务失败时先看“报告与导出”的执行日志，找到第一个实际失败的工具及其错误，而不是只看最后一行的通用退出码。确认输入路径、参考和网络后，可对失败或取消任务点击“重试”。已完成且输入文件路径/大小/修改时间标记一致的 fastp、Salmon 索引和定量可复用；统计和图表会按任务配置继续执行。

重试用于恢复原任务，不是随意更改模型后沿用旧结果。需要改变样本设计、过滤、参考、收缩、基因集或 WGCNA 参数时，建议创建新分析，以保留可比较的参数与结果。切换已有两套 DEG 方案不需要重新拟合模型。

点击“取消任务”会停止对应计算进程和本工作站对应的分析容器，保存已有文件。服务重启后未完成任务会标记失败，便于显式重试，不会宣称中断任务已经完成。删除或移动正在使用的输入、结果或程序目录会破坏任务和可追溯性。

更新前完成或正常取消任务，再使用 D 盘 `更新工作站.bat`。更新采用独立版本目录准备，成功后切换配置，保留旧版本和结果；活动任务或用户修改过的关键入口可能阻止更新，需按提示处理。源码模式可在确认本地修改和数据已保留后，使用仓库支持的快进更新和安装步骤。下载 Word 使用说明不会自动更新已经部署的程序。

## 21. 常见问题与处理顺序

| 现象 | 优先检查和处理 |
|---|---|
| 浏览器打不开 | 后台是否启动、终端是否退出、实际端口是8767还是8765、WSL本地转发 |
| 页面仍是旧版 | 旧进程是否仍占端口、更新是否生效、浏览器强制刷新、程序目录与版本 |
| 容器镜像找不到 | Docker服务是否可用、镜像是否安装完成，不仅看Python是否安装成功 |
| Docker Hub超时 | 引擎本身的代理和网络，与Windows浏览器代理分别确认 |
| 样本表读取失败 | UTF-8、CSV/TSV、列名、实际文件路径、是否为真正二进制Excel |
| 设计混杂/自由度不足 | 组别与批次分布、协变量数量、是否把每个样本唯一ID当协变量 |
| Salmon映射率低 | 物种、参考版本、转录本FASTA、双端配对、方向与污染可能 |
| strict没有DEG | 先查设计、重复变异、效应和检验功效；如做nominal探索须标注未控FDR |
| GO未运行/ID查不到 | 真实注释是否准备、基因ID是否为相应版本、映射率和低表达过滤 |
| GMT零匹配 | 物种及ENSEMBL/ENTREZID/SYMBOL选择是否一致 |
| WGCNA跳过 | 独立样本是否达到15、是否有重复subject，而不是简单增加技术重复 |
| 图为空或标签缺失 | 实际结果是否为空、P值是否NA、显示过滤、标签ID、重绘保存是否完成 |

排查时保留原始报错、模块状态和参数。不要为绕过错误而关闭 TLS/校验、伪造注释、强制把 TPM 当计数、自动删除离群样本或把缺失统计值改成显著值。

## 22. 建议的结果阅读和验证路线

先确认样本及比较方向，阅读 fastp/Salmon/MultiQC 和样本诊断；再核对模型、低表达过滤、实际注释覆盖和两套 DEG 数量。选择一个明确标注的方案作为当前 DEG/ORA 主线，保留另一套用于敏感性对照。

随后联合查看 GO/GMT ORA、全基因 GSEA 和成员表达。若 ssGSEA 或 ULM 提供额外信息，检查其分数差和同设计 limma 结果，而不是把不同统计量放在一张表中直接按数值大小排名。WGCNA 仅在独立样本和质量前提满足时使用。

最后选择少量有一致证据的候选基因/通路，回到每个生物学重复检查；记录阴性和不一致结果，准备独立实验验证。报告应区分“统计差异”“富集关联”“活性推断”“共表达关联”和“已验证机制”，不要把探索级结果写成因果结论。

## 23. 当前未实现的分析边界

本主线不从 Salmon 结果推算基因组 BAM，不内置 STAR/HISAT2 完整基因组比对、基因体覆盖、外显子/内含子分布、IGV 轨迹、新转录本组装、rMATS 可变剪接或 GATK/SnpEff 变异注释。需要这些结果时，应使用独立成熟流程并保留 BAM、对应参考及专门质量控制。

单细胞分析、空间转录组、差异转录本使用/剪接、细胞成分去卷积、正式中介/因果推断、任意随机效应和复杂时间交互模型也不是当前内置能力。可将本工作站的结果作为进一步分析输入，但需要各方法独立的数据要求、统计设计和验证。

STRING 蛋白互作网络与 WGCNA 共表达网络含义不同；本工作站的共表达边不能标成 STRING 实验互作。尚未运行或所需资源缺失的模块，应在报告中如实写“未运行/跳过”，而非填入示例结果。

## 24. 方法论文与成熟项目参考

以下论文支持分析原理；实际软件版本以任务的 `R-session.txt`、定量元数据和执行日志为准。以下资源不等于工作站已经实现其中所有功能。

[1] Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 2018. DOI: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)。项目：[OpenGene/fastp](https://github.com/OpenGene/fastp)。本工作站实际运行 FASTQ 预处理。

[2] Patro R, et al. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 2017. DOI: [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)。项目：[COMBINE-lab/salmon](https://github.com/COMBINE-lab/salmon)。本工作站实际运行转录本定量。

[3] Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research, 2015. DOI: [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)。项目：[mikelove/tximport](https://github.com/mikelove/tximport)。本工作站实际使用基因汇总与标准 DESeq2 入口。

[4] Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 2014. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)。项目：[thelovelab/DESeq2](https://github.com/thelovelab/DESeq2)。本工作站实际使用负二项模型、Wald 检验与可选 normal LFC 收缩。

[5] Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society B, 1995. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)。支持本工作站所用 BH 多重检验校正原理。

[6] The Gene Ontology Consortium. The Gene Ontology resource: enriching a GOld mine. Nucleic Acids Research, 2021. DOI: [10.1093/nar/gkaa1113](https://doi.org/10.1093/nar/gkaa1113)。官网：[geneontology.org](https://geneontology.org/)。注释需与物种和版本匹配。

[7] Kanehisa M, et al. KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Research, 2023. DOI: [10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963)。官网：[KEGG](https://www.kegg.jp/)。工作站导入或获取小鼠 GMT，不包含全部 KEGG 可视化和商业资源。

[8] Liberzon A, et al. The Molecular Signatures Database hallmark gene set collection. Cell Systems, 2015. DOI: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)。资源：[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/)。小鼠或可靠同源转换和授权条件需另外确认。

[9] Wu T, et al. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation, 2021. DOI: [10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141)。项目：[YuLab-SMU/clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)、[enrichplot](https://github.com/YuLab-SMU/enrichplot)。借鉴富集组织和展示，当前 ORA 使用本项目实现的超几何检验。

[10] Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 2005. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)。支持 GSEA 排序富集原理。

[11] Korotkevich G, et al. Fast gene set enrichment analysis. bioRxiv, 2021 revision. DOI: [10.1101/060012](https://doi.org/10.1101/060012)。项目：[ctlab/fgsea](https://github.com/ctlab/fgsea)。这是预印本，工作站实际使用 fgsea multilevel；应同时记录所用包版本。

[12] Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics, 2013. DOI: [10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7)。项目：[rcastelo/GSVA](https://github.com/rcastelo/GSVA)。本工作站实际使用该包的 ssGSEA 功能。

[13] Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research, 2015. DOI: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007)。本工作站用 limma 比较通路/调控活性分数，不以 limma 替代主线 DESeq2 基因计数模型。

[14] Badia-i-Mompel P, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances, 2022. DOI: [10.1093/bioadv/vbac016](https://doi.org/10.1093/bioadv/vbac016)。项目：[saezlab/decoupleR](https://github.com/saezlab/decoupleR)。本工作站实际使用 ULM，不代表已运行其所有推断方法。

[15] Garcia-Alonso L, et al. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research, 2019. DOI: [10.1101/gr.240663.118](https://doi.org/10.1101/gr.240663.118)。项目：[saezlab/dorothea](https://github.com/saezlab/dorothea)。论文主要为人类资源，不能直接替代本工作站需要的小鼠网络；使用的 mouse 导出或同源转换需另记录。

[16] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 2008. DOI: [10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)。项目：[cran/WGCNA](https://github.com/cran/WGCNA)。本工作站实际运行满足样本前提的 signed 网络。

[17] Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 2016. DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)。项目：[MultiQC/MultiQC](https://github.com/MultiQC/MultiQC)。本工作站 FASTQ 流程实际运行汇总。

[18] [nf-core/rnaseq](https://github.com/nf-core/rnaseq) 与 [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance)：参考样本表、定量/质控、设计、比较和可追溯报告组织；本工作站不是这些流程的完整分支，也未声称执行了其完整 BAM 流程。

[19] [EnhancedVolcano](https://github.com/kevinblighe/EnhancedVolcano)、[ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)、[ideal](https://github.com/federicomarini/ideal)、[pcaExplorer](https://github.com/federicomarini/pcaExplorer)：参考图形布局、标签和诊断优先的交互思路；本工作站实际交互使用本地 Plotly，静态导出使用 Matplotlib、SciPy 与 adjustText。

[20] [msigdbr](https://github.com/igordot/msigdbr)、[Reactome](https://reactome.org/)：可作为通路资源来源；需自行准备有权使用且 ID/物种一致的 GMT，当前不自动保证所下载资源与样本版本匹配。

[21] [kangweihou522-cmyk/16s-](https://github.com/kangweihou522-cmyk/16s-)：参考 Windows/WSL 本地网页工作站、分步操作与图表编辑方式；转录组主线使用独立定量和计数统计模型。

第三方源代码、稳定源码包及许可范围另见 `docs/THIRD_PARTY.md`。复用软件方法时应引用实际所用方法和数据资源，而不是仅引用本说明书或一个 GitHub 地址。

## 25. 简明术语表

| 术语 | 含义 |
|---|---|
| bulk RNA-seq | 多细胞混合样本的转录组测序；平均信号受细胞组成影响 |
| biological replicate | 独立生物学对象形成的重复，不等于重复测序 |
| transcript / gene | 转录本 / 基因；一个基因可有多个转录本 |
| count / TPM | 计数或估计数量 / 长度与文库规模相关的丰度尺度，分析用途不同 |
| normalized / VST | DESeq2归一化计数 / 适用于诊断和探索的方差稳定变换 |
| LFC | log2 fold change，必须同时记录比较方向 |
| pvalue / padj / FDR | 原始P值 / 校正P值 / 预期错误发现比例的控制概念 |
| DEG / ORA | 差异表达基因 / 选中列表的过度代表分析 |
| GSEA / NES / leading-edge | 全排序富集 / 标准化富集分数 / 富集峰核心成员 |
| ssGSEA / ULM | 单样本基因集评分 / 基于先验权重的单变量线性模型活性推断 |
| eigengene / kME | 模块主要变化轴 / 基因与模块轴的成员度相关 |
| confounding | 两个因素在现有设计中不可区分，不能靠简单批次校正修复 |
| provenance | 输入、参数、版本和执行来源记录，不自动等于完整数据备份 |

## 26. 说明书、网页与 Word 的同步规则

`docs/USER_GUIDE.md` 是唯一规范内容来源。后续修改分析方法、输入约束、阈值、图表含义、输出或部署操作时，应在同一次变更中更新对应章节和版本说明，再重新生成网页说明及 Word 文档，并运行同步检查。

网页和 Word 是同一规范源的发布形式，不应分别手工编辑成互相冲突的版本。生成物应随代码一起进入 GitHub；用户本机则需更新工作站版本并重启后台，才会看到新的使用说明。已下载的旧 Word、已经生成的历史实验报告和旧任务参数不会因 GitHub 上的说明书更新自动变化，应保留其版本以便追溯。
