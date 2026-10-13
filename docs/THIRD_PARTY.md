# 第三方组件与来源

- `vendor/fgsea-1.24.0.tar.gz`：ctlab/fgsea，稳定 RELEASE_3_16，commit `ac74ccd935c15623b8584caa791835aec514144b`，MIT，许可证随源归档保留。
- `vendor/decoupleR-2.4.0.tar.gz`：Bioconductor 官方 GitHub mirror bioc/decoupleR，RELEASE_3_16，commit `284793542e6c4d5c9ecc72a35a0390ff1a554e77`，GPL-3 + LICENSE，许可证与完整包源码随归档保留。
- `static/vendor/plotly.min.js`：plotly.js-dist-min 2.35.3，MIT；文件头包含版权与许可说明。
- 下载使用 HTTPS；固定 SHA256 保存在 `vendor/SHA256SUMS`，安装前校验。包源完整保留，以便审查与离线构建。
- Docker 中 DESeq2、GSVA、WGCNA、limma、tximport 等来自 Debian 签名软件源；各自许可证、源码与系统依赖遵守上游条款。R 会在结果中写入 sessionInfo。
- 客户报告只用作需求与方法参考，没有随仓库发布，也没有复制其图片、JS/CSS 或结果数据。
