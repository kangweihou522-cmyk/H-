# v2 验证记录

本次在云端 Linux、Python 3.12、Docker 上实际执行。目标为 Windows + WSL2 Ubuntu 24.04（32GB），原生 Ubuntu 同样支持；没有启动 Windows/WSL2 VM，因此本地配置与 bat 脚本检查不能视为 Windows 实机测试。用户报告只用于需求与方法参考，未执行其 JS，未发布其文件或真实数据。

## 安装与回归

- `bash install.sh` 实际完成并再次执行成功，构建 `mouse-rnaseq:2-core` / `mouse-rnaseq:2`。镜像中的 DESeq2、tximport、GO.db、jsonlite、fgsea 1.24.0、GSVA 1.46、limma、WGCNA 1.72.1、decoupleR 2.4.0 实际加载通过。稳定源码和本地 Plotly 的 SHA256 校验通过，签名/TLS 保持开启。
- 服务实际启动、停止、重新启动；健康检查、任务读取、文件下载、错误响应通过。
- **24 项 pytest**：原始计数约束、样本名/路径保护、组重复数、批次混杂、多组比较、数值协变量、样本排除后复查、Windows/file URI 转换、TSV 格式 xls、GMT、Origin/Host 保护、PCA/导出、GSEA 加权运行曲线、极小 P/NA 表格、缺失 Entrez 条件下基因查询、全 NA padj 时表达图、R 特有 WGCNA 颜色导出。
- Python、JavaScript、Shell 语法与 Git 空白检查通过。

## 真实工具端到端（本地合成输入）

| 入口/场景 | 实际执行 | 检查结果 |
|---|---|---|
| 6 个双端 FASTQ，240 转录本 | fastp、Salmon 建索引/定量、tximport、DESeq2、GO/GMT ORA、fgsea multilevel、GSVA ssGSEA、MultiQC、图表/报告 | 6 个样本定量完成，预设 40 上调/40 下调方向正确，padj<0.05；生成 fastp 原始/过滤 reads、Q20/Q30、碱基曲线和真实 MultiQC HTML |
| 同一合成计数入口 | DESeq2、GO/GMT/fgsea、表达诊断、ssGSEA | 完成；不伪造 FASTQ 质控/TPM |
| 已有 Salmon 导入 | tximport、DESeq2、双方案、ORA/GSEA、ssGSEA、导出 | nominal 88（43 上/45 下），strict 80（40 上/40 下）；小样本 WGCNA 明确跳过 |
| 18 样本、600 基因、3 组、批次协变量、2 个比较、数值表型、小鼠 ID 的人工 TF 网络 | DESeq2 normal LFC 收缩、双方案、fgsea、ssGSEA+limma、WGCNA、decoupleR ULM+limma | 3 个通路活性集合、3 个共表达模块、3 个调控因子/600 个目标；所有完成的活性矩阵均为有限数值 |

18 样本结果中：

| 比较 | nominal 总数（上/下） | strict 总数（上/下） |
|---|---:|---:|
| Treatment / Control | 134（107/27） | 0（0/0） |
| Rescue / Treatment | 275（200/75） | 1（1/0） |

通过独立 Python 条件计算核对所有 profile 的 DEG 清单、方向和计数。strict 为零时热图显示为空，不使用未显著基因冒充 DEG。两种方案均保留完整结果，切换不会改变全基因 GSEA 的 running score。

用 SciPy 从导出的 eigengene 和 trait 矩阵重新计算 WGCNA 的每个 Pearson r/P，与 R 输出逐项一致。该合成数据的 disease_score 关联并未达到 BH 显著，报告如实保留。另在真实 R 容器中验证：6 个独立样本跳过 WGCNA；36 次测量来自18个重复 subject 也跳过，避免用重复测量抬高独立样本数。

## 浏览器和图表

system Chromium + Playwright 实际验证桌面 1440/1500px、手机 390px，无横向页面溢出或 JavaScript 异常：

- 8 步导航、输入模式切换、合成演示样本读取、项目/比较选择、分页/搜索表格、基因查询。
- 同时显示两套 DEG 数量；实际点击方案按钮，更新 active profile、清单和火山图 P 轴；全基因 GSEA 曲线保持相同。
- 用真实鼠标点击已渲染火山图点，显示对应基因信息；按 Symbol 定位/高亮，实际鼠标拖动标签，重绘后位置持久保存。另用 Plotly click/relayout 事件检查处理逻辑；改标题、导出 SVG 校验基因标签/标题/设置。
- GSEA 表格通路点击进入曲线；真实浏览器下载 ZIP，校验归档和内嵌 SVG 的离线 report.html。
- WGCNA soft-threshold SVG/PDF/PNG/数据、模块/表型和调控图表实际生成。

## 已修复的实际问题

Docker 代理别名解析使用主机解析地址传入 build --add-host，保留现有 Docker 认证。Pandas 3 的只读/整数 PCA 数组改为浮点副本。R 多协变量列表读取和 WGCNA 版本 API 修正；WGCNA/BiocGenerics 的 cor 方法冲突修正；R 特有模块颜色转换到可用绘图颜色。缺失 Entrez 不再令数字 ID 查询失效。严格方案没有 DEG 时不回退到全基因差异热图。

## 复现

先安装开发检查依赖（无需下载官方测试 reads）：

```bash
.venv/bin/pip install pytest==9.1.1 playwright
.venv/bin/python -m pytest -q
docker run --rm -v "$PWD:/project:ro" mouse-rnaseq:2 Rscript /project/tests/check_mining.R
.venv/bin/python scripts/make_smoke_data.py runs/smoke_inputs
.venv/bin/python scripts/make_mining_data.py runs/mining_inputs
```

在本机已经启动的工作站，可填写生成的路径，或直接提交各目录的 request.json。需要保存返回的任务 ID：

```bash
curl -H 'Content-Type: application/json' --data-binary @runs/mining_inputs/request.json http://127.0.0.1:8765/api/jobs
```

完成后将 TASK_ID 替换为实际 ID：

```bash
.venv/bin/python scripts/verify_v2.py runs/TASK_ID
.venv/bin/python tests/browser_check.py --chromium /usr/bin/chromium
.venv/bin/python tests/browser_results.py --chromium /usr/bin/chromium --job SIX_SAMPLE_TASK_ID
```

浏览器检查需要系统 Chromium 或已安装的 Playwright 浏览器；不作为运行工作站的前提。browser_results 只用于6样本合成 fixture，会编辑其图表和 active profile。18样本 fixture 输入数值与本次检查一致；人工调控网络/GO/GMT 仅证明运算与已知效应恢复，不代表真实小鼠生物学关系。

## 验证边界

- 云端网络策略拒绝 Bioconductor 域名，**真实 org.Mm.eg.db 安装/映射路径未在本云环境验证**。真实 GO 可在本地运行注释安装脚本，或导入真实离线 TSV；缺失时明确报告 GO 未运行。人工 GO 注释的真实统计运算已验证。
- KEGG REST 在线下载未验证；GMT 运算已验证。真实 KEGG/Reactome/Hallmark 必须另备匹配小鼠 ID、来源清晰且有权使用的集合。
- 低次数合成 fgsea 检查中有统计分布不平衡/eps 下限提示；不可靠的 P/NES 保留 NA，未伪造为显著。真实分析默认初步排列1000，应查看日志和 log2err。
- 未运行真实实验/官方 benchmark，不从合成 fixture 宣称真实生物学结论；需确认样本设计、参考/注释版本和建库方向。
- 当前 Salmon 索引未包含基因组 decoy。基因组区域分布、IGV、新转录本、剪接、变异和 STRING 蛋白互作未由本主线执行，详见 docs/WORKFLOW.md。
- module–trait 为未调整协变量的 Pearson+BH，单基因–表型为探索性 Spearman。WGCNA 不适用于本实现跳过的小样本/重复测量情形。
- 本地服务需持续运行；云环境配置草稿不等于已发布快照，Docker/网页进程不会跨任务自动保留。

`runs/`、依赖环境、镜像、客户报告与数据均不进入源码提交。

## D 盘部署入口验证（2026-10-14）

- 24项Python检查包含新增7项：准备/激活事务、重复部署、保留数据/旧版本、保护修改过的启动文件、并发版本切换、中文/空格/引号与shell特殊字符、Ubuntu入口跟随当前版本、其他服务/PID重用保护、未结束分析拒绝更新。
- 使用官方 PowerShell 7.4.6 运行时（官方发布SHA256已核对）解析Windows脚本；13项WSL选择、原生参数引号、UNC路径、现有服务、端口归属与新窗口/浏览器调用检查通过。这些wsl.exe/浏览器调用被模拟，没有宣称Windows实机通过。另实际运行.NET ZIP解压，验证合法源码与拒绝路径穿越，下载调用模拟。
- 在Linux的中文/空格目录模拟D盘层级，实际安装独立于程序目录的Python环境、复用Docker构建、激活版本、健康/目录检查、启动、只停止本工作站空闲进程、保留结果并切换新版本通过。
- 从该部署路径运行6样本合成FASTQ的真实fastp、Salmon、tximport、DESeq2、GO/GMT ORA、fgsea、ssGSEA与图表报告；最终nominal88、strict80，6份quant.sf齐全。实际MultiQC由外部venv/bin/multiqc执行，报告HTML和结果ZIP完整性通过；WGCNA小样本明确跳过。
- 完整验证发现并修复Rscript在--file中将空格编码成~+~的兼容问题：应用向容器传递准确的RNASEQ_SCRIPT_DIR，直接R调用保留解码后备；复用已完成的定量重试统计/图表后完成。pytest指定正式tests目录，避免部署副本被重复收集。
- Chromium验证新部署结果、双方案卡片、工作站/结果目录显示、MultiQC入口、手机布局与无JavaScript异常通过。
- **未运行Windows桌面或WSL2 VM**：COM桌面快捷方式、真实wsl.exe参数桥接、Windows localhost转发、Docker Desktop UI和首次sudo交互，需要在你的电脑执行部署入口后确认。没有直接写入用户电脑D盘。核心程序/结果的Linux路径验证不能替代Windows实机部署。

## 首次部署启动器修复（2026-10-14）

- 针对用户在「部署到D盘.bat」遇到代理警告和退出码126，确认并复现旧启动器将WSL标准错误混入路径的缺陷。此复现确认输出污染，没有在Windows上重现完整的`/usr/bin/w`失败链。
- WSL调用使用ProcessStartInfo和明确的CRT参数引用，标准输出/标准错误并行独立读取，兼容UTF-16LE发行版列表。路径必须为单行绝对路径；先检查源码脚本，再以`/bin/bash --`调用。失败包含步骤、参数及D盘诊断日志。
- 15项PowerShell逻辑检查通过；另以真实Linux可执行文件替代wsl.exe，验证中文/空格/空参数/引号/尾部反斜杠/字面shell字符、代理警告与路径隔离、UTF-16LE列表、20万字符双管道、126退出码、诊断记录及非捕获调用继承stdin。执行方法：安装PowerShell后`bash tests/check_windows_native.sh`；这是Linux进程传输测试，仍不构成Windows/WSL2实机验证。

## 多个WSL入口导致Start异常（2026-10-23）

- 实际复现：PATH中存在两个wsl.exe时，`Get-Command -CommandType Application`返回两项，旧版将`.Source`数组转成一个含两条路径的FileName；Process.Start报找不到文件。选择单个可执行文件后相同场景通过。Windows优先系统入口，32位宿主优先Sysnative；后台启动使用相同解析函数，启动异常也记录文件路径。
- PowerShell检查新增双入口回归场景。仍在Linux运行PowerShell及可执行WSL替身，没有Windows/WSL2虚拟机；未宣称Windows系统路径分支已在实机验证。
- 参考16s-的Ubuntu直接启动方式，增加可选`部署到D盘.sh`，复用已有准备/安装/激活事务；Ubuntu快捷脚本直接运行Bash，不再转回PowerShell。验证配置激活、安装失败保留旧版本和数据，以及存在故障Windows命令时Ubuntu启动仍工作。日常Windows双击入口继续保留。

## WSL选项与输出编码修复（2026-10-24）

- 用户实机直接调用确认：Ubuntu-24.04为WSL2，`/bin/echo UBUNTU_OK`成功，Docker服务版本29.1.3。当前127错误属于启动器参数传递；这些检查没有证明后续下载或全部依赖安装完成。
- 对照Microsoft/WSL提交`4d85ce8dcbcf8091a66da89fd777d0b0cc96d685`：`src/windows/common/helpers.cpp`中的ParseArgument保留引号；WslClient.cpp直接匹配管理选项，只有识别裸`--exec`后才以CommandLineToArgvW解析剩余内容。旧版`"--list" "--quiet"`会落入默认Shell，产生用户报告的命令不存在。现在简单参数不加引号，需要引号的exec路径及参数继续使用CRT规则。
- WslClient.cpp支持子进程环境`WSL_UTF8=1`，用于统一WSL本身与Linux输出为UTF-8。保留旧版UTF-16列表后备，父进程环境不变。已从用户乱码还原出`/bin/bash: line 1: --list: command not found`。
- 17项PowerShell逻辑检查、33项Linux真实进程/管道检查通过，新增原始Windows命令行格式、UTF-8作用域及WSL警告与Linux错误同时解码检查。先前Linux argv替身会先去掉引号，无法发现WSL的原始命令行解析问题；新增原始字符串检查补齐这一缺口。没有Windows/WSL2 VM，实机完整部署仍以用户执行为准。

## WSL原生Docker代理配置辅助（2026-10-24）

- 用户实机确认Python依赖已安装，Ubuntu Docker Engine 29.1.3未配置代理；Windows代理为127.0.0.1:10808，构建在Docker Hub连接超时。该网络配置与启动器修复分开处理。
- `scripts/configure_wsl_docker_proxy.py`在用户本机探测WSL默认网关、验证HTTPS代理连通，保留Docker其他设置，检查活动容器，验证/备份配置，重启及核对生效值，失败恢复，最后拉取锁定摘要的基础镜像。
- 6项聚焦检查通过：保留已有设置和备份、活动容器拒绝修改、无效配置拒绝修改、重启失败恢复、代理未生效恢复、新建配置失败移除。Docker调用为替身，未修改云端或用户的Docker配置，未宣称已验证用户代理连通或本机镜像下载。
