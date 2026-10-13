# 小鼠 bulk RNA-seq 本地工作站

中文本地网页，界面布局参考 [kangweihou522-cmyk/16s-](https://github.com/kangweihou522-cmyk/16s-) 的深绿色分步工作站。分析实现使用 Salmon、Bioconductor tximport / DESeq2 / GO.db / org.Mm.eg.db，提供 FASTQ、现成 Salmon 定量、原始计数三种入口。所有样本数据与结果留在本机。

## 在你现有 Windows + WSL2 Ubuntu 24.04 上部署

Ubuntu 足以完成部署和分析。Windows 仅提供 WSL2 和浏览器；也可在原生 Ubuntu 上运行。已有 Java 21 和 Nextflow 25.10.4 不需要改动。本版由 Python 工作站顺序调度 Salmon/R，部署不依赖 Nextflow。

1. Docker Desktop → Settings → Resources → WSL Integration，启用 Ubuntu-24.04。若在原生 Ubuntu 使用 Docker Engine，确认当前用户可以运行 `docker info`。
2. 推荐把项目和 FASTQ 放在 WSL 的 `/home/用户名/rnaseq/`，不放 `/mnt/c`。将本仓库源码或提供的源码压缩包放入该目录。
3. 在 Ubuntu 终端执行：

```bash
sudo apt update
sudo apt install -y python3 python3-venv git
cd ~/rnaseq/H-                       # 换成你的实际源码目录
bash install.sh                     # Python 依赖 + Salmon/R Docker 镜像
bash scripts/install_mouse_annotation.sh  # 安装真实小鼠 GO、Symbol、Entrez 注释
bash launch.sh
```

在本机浏览器打开 `http://localhost:8765`。也可双击 `启动工作站.bat`；若发行版名称不同，设置 `WSL_DISTRO`，可用 `wsl -l -q` 查看名称。终端 Ctrl+C 停止网页服务。关闭服务时避免关闭尚在分析的任务；重启后运行中的任务会标记为中断，可重试。

初次安装需要访问 Docker Hub、Debian 软件源、PyPI、CRAN 和 Bioconductor；完成镜像和参考数据准备后，分析可离线运行。镜像中的 Debian 包签名、HTTPS 与 R 包 MD5 检查保持开启。默认 Docker 内运行 Salmon 1.10.1、R 4.2、DESeq2 1.38.3、tximport 1.26；精确 R 依赖版本保存在每个任务的 `R-session.txt`。主机已安装 Salmon/R 时优先使用主机命令，主机 R 需要相同依赖。

32GB 内存建议 4 线程，Docker 容器限制 24GB，任务顺序运行；界面允许 1–8 线程。WSL2 可在 Windows 用户目录 `%UserProfile%\.wslconfig` 设置：

```ini
[wsl2]
memory=28GB
processors=4
swap=8GB
```

修改后先结束分析，再在 PowerShell 运行 `wsl --shutdown` 并重新启动 WSL/Docker。此设置是建议，不是安装前提。实际需要的磁盘空间取决于 FASTQ 体积；每个任务的索引、定量、RDS 和图表会保留，预留原始数据以外的足够空间。Docker Engine 服务和 WSL 文件系统必须可读写。

## 准备数据

### 从 FASTQ 开始

准备匹配版本的 Ensembl 小鼠 **cDNA + ncRNA 转录本** FASTA（如 GRCm39）及 GTF。可用你已下载的参考文件，不要求先下载任何官方测试 reads。FASTA 转录本标识符必须与 GTF 的 `transcript_id` 对应；不要混合不同 release、品系或参考组装。将 cDNA 和 ncRNA FASTA 合并后输入工作站；不能把基因组 FASTA 当作转录本 FASTA。

从 GTF 建立转录本与基因映射：

```bash
.venv/bin/python scripts/prepare_reference.py /data/Mus_musculus.gtf.gz /data/tx2gene.tsv
```

`tx2gene.tsv` 是两列、含表头的 TSV。转录本版本号会在 tximport 时忽略，基因版本号会在注释时去除。若你的 FASTA 标识符格式不是 Ensembl，先提供对应映射，不能凭字符串猜测基因。

样本表 CSV，每组至少 2 个独立生物学重复，推荐 ≥3。生物学重复不是同一个样本的不同测序 lane；同一样本多 lane 请事先合并匹配的 R1/R2 gzip 文件。单端去掉 `fastq2` 列；不混用单双端。

```csv
sample,group,fastq1,fastq2,batch
C1,Control,/data/C1_R1.fastq.gz,/data/C1_R2.fastq.gz,B1
C2,Control,/data/C2_R1.fastq.gz,/data/C2_R2.fastq.gz,B2
T1,Treatment,/data/T1_R1.fastq.gz,/data/T1_R2.fastq.gz,B1
T2,Treatment,/data/T2_R1.fastq.gz,/data/T2_R2.fastq.gz,B2
```

路径使用 Ubuntu 绝对路径，例如 Windows `D:\data` 对应 `/mnt/d/data`。样本名使用字母、数字、点、下划线或短横线。`batch` 可省略；存在多个批次时设计为 `~ batch + group`，否则为 `~ group`。分组与批次完全混杂会被拒绝。

Salmon 文库类型 `A` 为自动判断；双端常见 `IU`（非链特异）、`ISR`（反向链特异）、`ISF`（正向链特异）；单端 `U`、`SR`、`SF`。应以建库试剂说明和 Salmon `lib_format_counts.json` 为依据。流程使用 `--validateMappings --seqBias --gcBias`。本版不自动进行接头剪切或去除 rRNA；需要时先用 fastp 等工具处理并保留 QC 报告。Salmon 输出保存完整定量、映射率与文库方向统计。使用转录本索引；尚未提供基因组 decoy-aware 索引入口。

### 导入已有 Salmon 或计数矩阵

Salmon 模式使用样本表 `sample,group,quant`，quant 为各样本 `quant.sf` 路径，并提供 tx2gene。定量文件旁若存在 `aux_info/meta_info.json`、`lib_format_counts.json`、`cmd_info.json` 会一起保留。

计数矩阵 CSV 首列为唯一基因 ID，其余列与 sample 名匹配，必须为非负整数**原始基因计数**，不是 TPM、FPKM、归一化计数。工作站拒绝小数与负数；无法从整数外观证明它就是原始计数，输入来源需你确认。GO 默认识别 Ensembl 小鼠基因 ID。

## GO、KEGG、Reactome 与 GSEA

- GO：默认使用 `org.Mm.eg.db` 的 GOALL 注释，包含祖先术语，按 BP/CC/MF 分别进行 BH 校正。背景为本次保留且具有 GO 注释的表达基因。
- KEGG：可导入你有权使用的小鼠 GMT，或运行下面脚本从 KEGG REST 建立。下载不属于工作站启动前提；请遵守 KEGG 数据使用条款。
- Reactome / MSigDB / 自定义通路：导入 GMT，明确 Ensembl、Symbol 或 Entrez 类型。每行 `通路ID<TAB>描述<TAB>基因1<TAB>基因2…`。人类基因集需要事先可靠的同源转换，不能直接当作小鼠通路使用。
- ORA：使用超几何检验，对所有满足大小条件的集合进行 BH 校正（包括无命中集合），输入为满足 padj 和 LFC 阈值的全部差异基因，未分开上/下调。
- GSEA：按 DESeq2 的 Wald stat 排序，权重为 |stat|，基因集合大小 15–500；固定随机种子 42，基因集合标签排列，按同符号零分布计算 NES 与经验 P，再做 BH 校正。默认 1000 次，允许 100–10000；P 值分辨率受排列次数限制。这个实现采用 preranked GSEA 的运行和统计定义，**不是 fgsea 的 multilevel 算法，也不是原始 GSEA 的 phenotype permutation FDR**。小样本解释需谨慎。

```bash
.venv/bin/python scripts/fetch_kegg.py /data/mouse_kegg.gmt
```

界面 GMT ID 类型选 `ENTREZID`。脚本保存来源时间及 SHA256。KEGG REST 网络不可用时，导入已有 GMT 即可。没有 GMT 时仅运行 GO。

离线注释表可替代 org.Mm.eg.db：带表头 `gene,symbol,entrez,go,ontology` 的 **TSV**，同一基因多个 GO 术语需要多行，ontology 为 BP/CC/MF。应包含传播后的 GO 祖先注释。可从已安装的包导出：

```bash
docker run --rm -v "$PWD:/project" mouse-rnaseq:1 \
  Rscript /project/scripts/export_mouse_annotation.R /project/annotation.tsv
```

如果既无 org.Mm.eg.db 又没有离线注释，工作站会明确显示“GO 未运行”，保留 DESeq2 和可映射的 GMT 结果，不会把空 GO 结果称作分析成功。Entrez/Symbol GMT 还需要对应的映射；零映射将报告错误。合成验证注释仅用于测试，不能用于实际数据。

## 编辑、导出与追溯

PCA、火山图、差异基因热图、GO 富集、GO GSEA、GMT/KEGG 富集和 GSEA 可逐图修改标题、轴标题、颜色、大小、字号、PNG DPI、展示数量及火山图阈值。每次重绘输出 SVG（文字可编辑）、PDF（矢量）、PNG，并保留 `.settings.json` 与 `.plot-data.csv`。显示阈值修改不重算 DESeq2/富集；修改统计阈值需创建新分析任务。

PCA 使用全体保留基因的 VST 值；热图选择 padj 最小的基因，对每行 VST 做 Z-score。标题中可输入中文，但 Matplotlib 所在主机需安装中文字体（Ubuntu 可安装 `fonts-noto-cjk`；字体可进一步在系统 Matplotlib 配置中设置）。

结果位于 `runs/任务ID/`，含输入样本表/配置、执行命令 JSONL、日志、`differential.csv`、`normalized.csv`、`vst.csv`、`dds.rds`、GO/GSEA/通路 CSV、Salmon QC、图表与 R 版本记录。输入 FASTQ/参考 FASTA 不会复制。只在本机 loopback 监听，不要直接暴露到公网。

## 验证

无需下载官方测试数据。合成数据由本地随机种子生成，含 6 个双端样本、240 条转录本、预设上/下调基因，以及明确标注为人工的 GO/GMT 集合。运行真实 Salmon / tximport / DESeq2 与富集运算，不使用伪造定量或假差异结果。

```bash
.venv/bin/pip install pytest
.venv/bin/python -m pytest -q
.venv/bin/python scripts/make_smoke_data.py "$PWD/runs/smoke_inputs"
bash launch.sh
# 另一个终端提交 runs/smoke_inputs/request.json 到本地 /api/jobs
# 等待 completed 后：
.venv/bin/python scripts/verify_smoke.py "$PWD/runs/返回的任务ID"
```

验证结果与已知限制见 [VALIDATION.md](VALIDATION.md)。真实小鼠 org.Mm.eg.db 下载与 KEGG REST 实时下载单独记录；它们失败不影响使用离线注释/GMT验证程序。当前版本没有原始 reads FastQC/fastp 页面、多因素交互模型、成对设计、自动同源转换或完整 nf-core/rnaseq 等全流程功能。
