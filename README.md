# NHANES neutrophil count and depressive symptoms — final code archive

研究：**Neutrophil count and depressive symptoms in U.S. adults: a cross-sectional analysis of attenuation after socioeconomic and clinical adjustment**。NHANES 2005–2018，7 个周期，主要样本 33,174 人。主设计为全成人复杂抽样设计后的子总体分析，合并 MEC 权重 WTMEC2YR/7，保留分层与 PSU，nest=TRUE。

本包以第七阶段安全修复代码为唯一代码基准，第五阶段锁定结果为数值基准，并纳入第八阶段限定运行验证记录。本次仅整理代码、文档与聚合参考结果，**未重新运行科研模型，未公开发布；已准备限定原创代码范围的 MIT LICENSE，等待最后发布批准**。

## 目录与使用入口

- `scripts/nhanes/`：队列、全成人设计、Model1–3、联合计数/Wald、症状、调整衰减、插补、敏感性、睡眠及图表。
- `scripts/mr/`：冻结 MR/MVMR、诊断、外部本地共定位及明确 UNVALIDATED 的上游在线流程。
- `scripts/lib/`：共同路径/设备安全层与全成人子总体设计函数。
- `config/storage_policy.json`：强制保护源码的父目录。此归档位于研究总目录内，统计工作区必须在研究总目录之外。
- `environment/`：观测软件版本、依赖与种子；不是完全锁定的跨平台安装环境。
- `supplementary/validated/`：当前 S5/S8 与 Figure3 的锁定聚合结果；`supplementary/reference/` 为历史数值核对基准，不能据此替代修订结果。
- `reference_results/README.md`：当前参考结果与历史比较文件的使用规则，避免重复复制同一套数据。
- `figures/`：锁定 Figure3 PDF（PNG未纳入此次公开候选，原因见下文），属于结果图而非投稿手稿。
- `provenance/workflow.csv`：逐步输入—代码—输出—论文用途；`docs/SCRIPT_CATALOG.md`：每个脚本的作用及正式执行顺序。
- `docs/VALIDATION_STATUS.md`、`docs/RELEASE_NOTES.md`：验证边界及本次整理检查结论。

## 数据获取与输入

从 CDC NHANES 公共数据页面按原周期下载 XPT 和相应说明，入口：https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/ 。本包不含参与者数据，也不自动下载。

保存为 `2005-2006/Demographics/demo_d.xpt`、`2005-2006/Examination/bmx_d.xpt`、`2005-2006/Laboratory/cbc_d.xpt`、`2005-2006/Questionnaire/dpq_d.xpt` 等周期/组件目录。七周期后缀为 D–J，文件匹配不区分大小写。核心 DEMO/BMX/CBC/SMQ/DPQ/MCQ/DIQ/BPQ/ALQ/RXQ_RX；额外 COT/COTNAL、CRP/HSCRP、SLQ 和无周期后缀 RXQ_DRUG 字典。睡眠模块另需全部七周期的 MCQ/KIQ_U/PAQ/SLQ，共 28 个真实文件。周期特例及完整对象契约见 `data/README.md`。

MR/共定位输入与 NHANES 分开：OpenGWAS 暴露 ieu-b-34、ieu-b-32、ukb-d-30080_irnt；PGC MDD2025 指定不含 23andMe/UKBB 的结果及本地 VCF/索引。外部输入须自行取得授权、确认版本和 SHA256，绝不放入本代码包。数据权利见 `DATA_SOURCES_AND_RIGHTS.md`，其历史核对日期不等于本次重新确认许可。

## 软件与安装

已观测 R 4.5.2、macOS arm64，Python 3.12.14（普通启动器要求 Python ≥3.11），survey 4.5、mice 3.19.0、jsonlite 2.0.0、TwoSampleMR 0.6.29、MVMR 0.4.2、coloc 5.2.3。完整观测清单见 `environment/package_versions.csv`；安装说明见 `environment/INSTALL.md`。先建立独立软件环境并人工核对包版本，再执行 `Rscript --vanilla environment/check_versions.R environment/package_versions.csv`。检查脚本只显示版本。

Figure1/2 原绘图后端使用 macOS Quartz；并行 Bootstrap 使用 mclapply 的 4 个 fork worker；Windows/Linux 未验证。H10 需要 shasum -a256。ieugwasr 1.1.0.9000 的开发提交未保留，不能声称环境完全锁定。

## 初始化：源码目录不能作为输出目录

以下均为占位绝对路径，请替换为自己的真实路径；`/separate/new_session` 必须**尚不存在**，且与研究总目录、源码、历史候选包、全部输入目录互不重叠。归档策略会自动保护源码父目录。外部研究者只需声明自己实际存在的输入/保护目录，无须建立阶段5–8目录；本地作者应额外声明原研究及历史目录。

```sh
python3 scripts/init_workspace.py --workspace /separate/new_session \
  --protected-root /data/public_nhanes \
  --raw-nhanes /data/public_nhanes
```

布局：`new_session/inputs/gwas_cache/` 是只读导入副本；`new_session/output/work/`、`output/logs/`、`output/figures/`、`output/tables_revised/` 才允许写入。私有 `run_config.json` 不应提交。未初始化拒绝执行；符号链接、祖先/子目录重叠、写回输入、覆盖已有输出均拒绝。绘图设备初始化失败即停止，不回退 `Rplots.pdf`。不要手写 marker、伪造配置、重复使用已有工作区来绕过拒绝。

## 运行顺序

先做无模型检查：

```sh
python3 scripts/run_job.py scripts/nhanes/05_joint_splines.R --workspace /separate/new_session --check-only
Rscript --vanilla scripts/nhanes/05_joint_splines.R --workspace /separate/new_session --check-only
python3 scripts/run_pipeline.py --workspace /separate/new_session --frozen-imputations --plan
```

分析顺序为：公开模块队列 → 全成人设计 → 扩展协变量 → 完整分析对象 → 主模型/插补或冻结对象分支 → 症状及敏感性 → 描述/表格/Figure1–2 → 衰减/AUC/Shapley/异质性 → 种族/睡眠 → 联合样条 → 完整三计数性别交互 → **S5/S8 最终导出** → Figure3。准确脚本名称与依赖顺序见 `docs/SCRIPT_CATALOG.md` 和 `provenance/workflow.csv`，不要按文件数字简单排序；例如 27 在 26 前。

只有在安排新的科学运行时才执行以下完整流水线：

```sh
# 冻结分支：初始化时额外提供 --import-cache imputations_standardized.rds=/authorized/imputations_standardized.rds
python3 scripts/run_pipeline.py --workspace /separate/new_session --frozen-imputations
# 或另建独立新工作区，重新执行完整 20×20 插补；两条命令不能在同一输出目录接续运行。
python3 scripts/run_pipeline.py --workspace /separate/another_new_session
```

冻结分支需要合法 mice 对象与原样本行对齐，不能由普通用户凭空重建。单独运行下游脚本还需 analytic_v3.rds、全成人设计、缩放参数及对应上游输出；缺失即停止。第八阶段为验证限定模块导入过这些已验证前置对象，不代表完整上游已重新执行。单独运行说明与 MR 示例见 `docs/RUN_GUIDE.md`。

## 输出与论文对应

Table1/2 及主联合模型/Wald → `work/export/`、`work/imputation/`；症状/阈值/临床敏感性 → 各同名 work 子目录；协变量衰减 → `work/covariates/`、`work/attenuation/`；AUC、Shapley、规格 → `work/discrimination/`、`work/shapley/`、`work/heterogeneity/`；最终性别标量/联合检验 → `work/sex/`、`work/sex_total/`；Figure3/S8 → `work/spline/`；S5/S8 显示 CSV → `output/tables_revised/`；MVMR/S6 → `work/mr_export/`；共定位/S7 → `work/coloc/`。逐项映射见 workflow。

四条 Figure3 曲线均联合调整淋巴细胞和血小板。当前性别探索性模型区分 N×sex 单项 P=0.455 与三项交互联合 P=0.170；不采用旧受限交互结果生成当前 S5。README 中的显示值不替代最高精度参考 CSV。

## 已验证范围、随机性与限制

- 第五阶段：分组件重建 NHANES、主模型及多类敏感性/Bootstrap/规格结果；主插补模型使用冻结对象。它不是单条命令从原始数据到新插补再到全部论文结果的端到端证明。
- 第七阶段：安全层、静态模型一致性和 98 项安全测试。
- 第八阶段：MVMR b/SE/P 9 项一致；睡眠 28 输入匹配、40 次模型拟合、7 表 11,305 字段一致；S5/S8 与精度共 219 显示项一致；9 类启动器场景通过。S5/S8 执行了 pipeline 的真实最终导出语句，整条 pipeline 未执行。
- 新插补设定 20 份 ×20 迭代、种子 20261008；AUC 500 次/20261013；Shapley 200 次及衰减 1,000 次/20261014。详见 `environment/seeds.csv`。本次未重跑任何科学模型。
- 历史随机 MR 种子未知，median/mode 的重拟合 SE/P 可有差异；不以新随机值替换锁定稿件值。在线提取、LD clumping、历史 LD 快照继续 **UNVALIDATED**。本地冻结 MR/共定位与在线重建必须分开。
- TwoSampleMR 0.6.29 加载时会访问 GitHub 公告/版本文件；第八阶段未调用 OpenGWAS 数据接口，但不能称严格无网络复现。

## 许可与公开

MIT LICENSE 已按 `Copyright (c) 2026 Wenjin Wang` 准备。适用范围仅为你有权授权的原创源码及配置：`scripts/**`、`tests/**`、`environment/check_versions.R`、`config/storage_policy.json`。依赖包、图件、聚合/参考 CSV 和外部数据不由此重新授权；其他权利与条款保持原样。来源核查和作者对 MR 原创性的确认见 `docs/LICENSE_REVIEW.md`。

本包无个体级数据、私人 RDS/.RData、GWAS 输入、文献库、投稿手稿或凭据。公开范围依照作者接受的147文件/PDF-only方案，图件及聚合结果不默认采用MIT。尚未发布 GitHub/Zenodo。

## 样本量与发布边界

|分析|样本量|
|---|---:|
|全成人设计基底|42,143|
|主要分析/主MI分支/完整三计数性别交互|33,174|
|Figure3 Model1/2共同完整样本|29,415|
|Figure3 Model3完整样本|29,360|
|睡眠Model3/4共同样本|29,024|

这些样本量属于不同分析，不能互相替代；其余敏感性分析以各输出模型注册表为准。完整原始数据→新插补→全部结果的单次执行没有得到独立验证。

公开候选仅保留Figure3 PDF与绘图代码/聚合坐标；PNG因内嵌第三方ICC profile暂不捆绑。PDF科研内容未修改，也未重新拟合或绘图。代码许可与图件、结果、外部数据的权利范围分开；MIT文本已准备，实际发布仍需最后批准。
