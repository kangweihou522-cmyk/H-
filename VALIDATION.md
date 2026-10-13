# 验证记录

本次验证在云端 Linux、Python 3.12、Docker 上实际执行。目标部署为 Windows + WSL2 Ubuntu 24.04（32GB），原生 Ubuntu 同样支持；没有实际启动 Windows 或 WSL2 VM，因此 bat、WSL 配置以文档与脚本检查为证据，不能视为 Windows 实机验证。

## 已通过

- `bash install.sh`：实际安装并再次执行成功，复用 Docker 构建缓存；容器内 DESeq2、tximport、GO.db、jsonlite 加载通过。所有签名/校验/TLS 检查保持开启。
- `bash launch.sh`：实际启动，停止本次服务后重启，健康请求与已有任务读取通过。
- 8 项 pytest 测试：原始计数约束、样本名/路径保护、批次混杂、每组重复数、有效批次设计、GMT 格式、跨站/Host 保护、PCA 数据回归及 SVG/PDF/PNG/数据导出。
- 实际浏览器（system Chromium / Playwright）：桌面1440px和手机390px，无横向溢出；切换输入模式、填入演示数据；结果7张图、修改标题/颜色、SVG内容验证及PDF下载通过，无 JavaScript 异常。
- 本地合成 FASTQ 端到端：真实 Salmon 建索引和6个双端样本定量（240转录本，全部映射率100%），tximport、DESeq2、GO ORA、GO GSEA、自定义GMT ORA/GSEA、7类图表生成通过。预设40个上调与40个下调均显著（padj<0.05），方向正确；GMT GSEA上调 NES约2.497，下调 NES约−2.846。
- 原始计数入口：实际完成240基因差异分析与所有GO/GMT计算；上调方向正确。
- 已有Salmon入口：实际导入6个quant.sf，240基因差异结果与FASTQ入口一致。
- 重试：PCA兼容问题修复后，复用完成的Salmon定量重新执行R分析和绘图，任务达到completed。
- Python、JavaScript与Shell语法检查通过。

合成fixture仅证明运算、任务管理、输入/输出与已知效应恢复正常。GO成员关系和GMT集合是人工构造，未用于推断真实小鼠通路；验证没有把官方测试reads下载成功作为前提。

## 已诊断、修复的问题

云端 Docker 构建无法解析代理别名：复用主机解析的代理地址到build --add-host，正常保留签名校验；脚本不会保存代理凭据。Debian不提供org.Mm.eg.db包：小鼠真实注释改为Bioconductor安装或用户提供离线TSV。pandas 3中的只读数组及整数数组PCA中心化错误：显式创建浮点可写副本并加入回归检查。R requireNamespace参数名改为quietly。

## 未验证与限制

- 云端网络策略拒绝Bioconductor域名，未实际安装真实org.Mm.eg.db，也未验证该包的实际基因映射路径；离线注释输入与GO运算已验证。提供了 `scripts/install_mouse_annotation.sh` 和 `scripts/export_mouse_annotation.R`，需要在允许联网的本地执行或导入既有真实注释。无GO注释时界面明确报告GO未运行。
- KEGG REST实时下载未运行；下载器使用HTTPS与明确的mmu来源，记录日期/哈希。GMT计算路径已用本地人工集合验证；不能把它称作真实KEGG生物学验证。
- 未运行真实实验样本或官方benchmark，不承诺生物学结论正确；样本设计、参考版本、注释版本、建库方向与技术QC需要使用者确认。
- 不自动剪接头、不提供decoy genome索引和多因素交互/配对模型。GSEA为基因集合标签排列与BH校正，不等同于fgsea multilevel/phenotype permutation FDR。
- 本地服务需持续运行；Docker容器和任务队列进程不会随云环境发布快照保留。启动说明已另存于环境配置。

验证输出默认保存在当前机器的 `runs/`，为忽略文件；源码压缩包不包含原始/合成reads、依赖环境、Docker镜像或用户数据。
