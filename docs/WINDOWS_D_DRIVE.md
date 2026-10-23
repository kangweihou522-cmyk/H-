# D 盘转录组分析工作站

目标电脑：Windows + WSL2 Ubuntu 24.04、32GB、已安装 Docker。云端不能直接写入你电脑 D 盘；请在你的电脑执行以下安装。无需先下载官方测试数据。

## 首次部署

1. 从 GitHub 下载最新源码 ZIP，**全部解压**到任意文件夹；不要只从 ZIP 内运行单个 bat。
2. 双击解压目录里的 **部署到D盘.bat**。程序默认创建 **D:\转录组分析工作站**。
3. 脚本自动选择 Ubuntu-24.04（否则 Ubuntu）、检查 Docker，并在需要时打开 Docker Desktop。如果提示 WSL Integration 未开启，在 Docker Desktop → Settings → Resources → WSL Integration 中启用相应 Ubuntu，再重试。Ubuntu Docker Engine 也支持已有服务启动。
4. 初次缺少 Python/venv/中文字体时会调用 Ubuntu sudo 安装；密码仅输入你电脑的 Ubuntu 终端，不发给他人。随后安装核心依赖与 Docker 分析镜像，保留签名、SHA256 和 TLS 检查。
5. 成功后创建 **转录组分析工作站** 桌面快捷方式，启动服务并打开浏览器。初次准备依赖需要联网；日常启动不重复下载安装。

如果失败，窗口会保留报错，可双击检查环境.bat。没有 Docker 时脚本显示官方安装地址，不会悄悄创建第二套 Docker 或修改现有 WSL 发行版。

### 部署提示代理警告或退出码126

如果旧版本在首次部署时显示 `wsl: 检测到 localhost 代理配置，但未镜像到 WSL`，随后出现 `/usr/bin/w: cannot execute binary file`：下载最新完整源码ZIP，解压到一个新文件夹，再双击其中的 **部署到D盘.bat**。无需删除已经创建的D盘工作站目录；原始数据、结果和旧版本会保留。首次部署尚未完成时，不必先运行依赖已安装配置的检查入口。

旧启动器将标准错误合并到标准输出，WSL启动警告可能混入 `wslpath` 的返回路径或JSON。新版分开读取两条输出、明确引用每个原生参数、检查路径及脚本存在性，并显示具体失败阶段。代理警告本身不等于Ubuntu启动失败；即使路径读取正常，NAT模式下Ubuntu仍可能无法使用Windows的localhost代理，需要根据实际下载报错处理网络。

如果失败命令显示为 `"--list" "--quiet"`，退出码127，或乱码还原后为 `--list: command not found`，也是旧启动器的参数问题。WSL管理选项需要保持裸参数，不能把每一项都加引号；新版仅对需要引用的路径/参数加引号。启动器子进程设置`WSL_UTF8=1`，使WSL警告与Linux错误使用相同编码。可先在Windows PowerShell直接运行 `wsl.exe --list --verbose`、`wsl.exe -d Ubuntu-24.04 --exec /bin/echo UBUNTU_OK`、`wsl.exe -d Ubuntu-24.04 --exec docker info --format '{{.ServerVersion}}'`，确认本机这三层是否工作，再对启动器做定点修复。

新版调用诊断保存在 **D:\转录组分析工作站\日志\wsl-launcher.log**，包含失败步骤、参数、退出码和标准错误。若仍失败，可提供失败窗口最后几行及此日志末尾相应记录（分享前隐藏私人路径/信息）。我们没有在Windows实机上复现该退出码，因此修复已确认的输出污染缺陷后，仍需以你电脑的重试结果确认。

## 每天打开

双击桌面 **转录组分析工作站**，或 `D:\转录组分析工作站\启动工作站.bat`。启动器自动检查 Docker、识别本工作站是否已经运行、检查服务归属并打开浏览器，不会把其他程序的页面当成工作站。

默认网页 `http://localhost:8767`，与 16S 常用的8765分开。请保持弹出的 Ubuntu 工作站终端运行；关闭时按 Ctrl+C。下一次启动后可从当前分析项目选择已有结果。

也可以在 Ubuntu 输入：

```bash
bash ~/启动转录组工作站.sh
```

Windows双击会调用Ubuntu后台，通常不用先手动打开Ubuntu窗口。选择Ubuntu方式时，打开Ubuntu后还需执行上面的启动命令，单纯打开终端不会启动网站。浏览器页面与后台是两部分：分析期间保持工作站终端和Docker运行，电脑不要休眠。

### 可选：直接在Ubuntu首次部署

Windows双击是正常部署入口。若需要完全从Ubuntu操作，在解压后的源码目录执行 `bash 部署到D盘.sh`。它使用同一套版本目录、依赖安装及配置激活流程，默认目标仍为 `/mnt/d/转录组分析工作站`，端口8767。无需通过PowerShell启动进程，部署成功后直接运行后台；以后执行 `bash ~/启动转录组工作站.sh`。源码所在目录可以通过Ubuntu的 `cd` 命令进入，也可以在Windows资源管理器中查看对应D盘位置。

如果以前已有不同内容的Ubuntu快捷脚本，新脚本会保存在 `~/启动D盘转录组工作站.sh`，使用部署输出列出的路径。此备用方式需要Ubuntu能访问Docker；使用Docker Desktop时应先打开它并启用WSL Integration。

若原来已有你修改过的同名脚本，会保留它，改写到 `~/启动D盘转录组工作站.sh`；如果两个都存在不同内容，则不覆盖，可使用 Windows 入口。

## 固定位置

| 目录/文件 | 用途 |
|---|---|
| `D:\转录组分析工作站\程序\版本摘要\` | 工作站源码；更新时创建新版本，保留旧版本及修改 |
| `原始数据\` | 放你准备导入的FASTQ、计数、样本表；现有原始文件不会自动移动 |
| `参考文件\` | 放匹配版本的转录本FASTA、GTF/tx2gene、注释TSV、GMT及调控网络 |
| `分析结果\任务ID\` | 自动保存全部任务结果、图表、参数、RDS、报告和日志 |
| `导出报告\` | 下载图表/报告时可手动选择的保存位置；浏览器默认下载目录不会被改动 |
| `配置\workstation.json` | 实际程序路径、Ubuntu发行版、端口与运行环境位置 |
| `日志\server.log` | 工作站服务日志；部署/更新源码归档在日志/updates |
| `日志\wsl-launcher.log` | Windows启动器的WSL调用与错误诊断 |
| `启动工作站.bat` / `检查环境.bat` / `更新工作站.bat` | 日常启动、诊断与保留数据的更新 |

网页“报告与导出”也会显示工作站目录和实际结果目录。Windows D盘在Ubuntu中为 `/mnt/d`。

**程序与分析数据位于 D 盘；Python 环境放在 Ubuntu 文件系统内以避免跨文件系统依赖问题，Docker 镜像由已有 Docker 管理。安装脚本不会迁移 Ubuntu 的虚拟磁盘或 Docker 的存储位置。** 这些环境占用哪个盘取决于你已有 WSL/Docker 的设置；若目前在C盘，整个计算环境不会因本脚本而自动搬到D盘。大型FASTQ通过/mnt/d访问可能慢于Ubuntu文件系统，是否另放原始数据由你决定。

## 更新和数据保留

结束分析后双击 `更新工作站.bat`。它通过 HTTPS 获取 GitHub main 源码，记录归档SHA256，检查ZIP越界路径，在独立版本目录准备依赖，成功后切换并重启属于本工作站的空闲服务。

存在运行/排队任务时拒绝更新，先在网页完成或取消任务。不覆盖原始数据和分析结果、不删除旧程序版本，也不覆盖你修改过的根目录启动文件。部署失败时当前配置不会切到未安装成功的版本；可查看错误后重试。旧任务的新模块仍需创建新分析。

如果已有旧工作站结果，先完成/取消旧任务，在旧终端 Ctrl+C 停止服务。确认目标任务ID目录没有同名结果后，再将旧 `runs/` 中的任务目录**复制**到新的 `分析结果/`，不要删除旧副本。脚本不会扫描全盘、自动迁移或重写旧数据。旧任务config中的输入绝对路径仍需保持可用。

## 自定义

默认不改 `.wslconfig`、WSL发行版、Docker配置、代理/凭据和Ubuntu系统账户。其他位置/端口/发行版可在解压源码的 PowerShell 中显式设置，例如：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\windows\Workstation.ps1 -Action Install -InstallRoot 'D:\转录组分析工作站' -Distro Ubuntu-24.04 -Port 8768
```

ExecutionPolicy仅对该次PowerShell进程生效，不修改全局策略。不要移动正在运行的工作站目录。

真实 GO 需要 org.Mm.eg.db 或真实离线 TSV；KEGG/Reactome/Hallmark 等需要匹配小鼠 ID 的 GMT；调控活性需要真实小鼠网络。它们不作为创建目录、安装和打开工作站的前提。真实注释的联网安装方法见主README；D盘部署的Python在配置的venv_wsl路径，可从检查环境输出/配置文件查看。

## 验证范围

见 VALIDATION.md。Python部署事务、保留用户数据、中文/空格/引号路径、进程归属、外部Python环境、真实分析链在Linux执行验证；PowerShell由官方校验过的运行时解析并用模拟的WSL/浏览器调用检查。没有 Windows 桌面/WSL2 VM，COM桌面快捷方式、真实wsl.exe桥接和Windows本地浏览器转发必须以你电脑的实际运行结果为准。
