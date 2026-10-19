# 打开模型与复核结果

## 环境

- Windows 与 Ansys Zemax OpticStudio 24R1；运行 API 需要相应许可。
- 原项目记录使用 ZOSPy 2.1.5；依赖见根目录 `requirements.txt`。本次发布仅做离线完整性与语法检查，没有重跑 Zemax。
- 材料使用 CDGM/NHG 目录。请在本机安装匹配的目录版本；仓库不附带厂商 AGF 文件或软件安装包。
- 导出的脚本可读取 `OPTICSTUDIO_DIR`，默认 `D:\Zemax`。`match_patent_glasses.py` 可读取 `ZEMAX_GLASSCAT_DIR`，默认当前用户 Documents/Zemax/Glasscat。

## 查看与下载

完整下载仓库，或解压 `packages/codex-4-zemax-simulation-bundle.zip`，保留目录层级。先阅读 `design_35mm_f28/README.md` 与设计状态记录。

在 OpticStudio 中打开 `candidate_v31_thermal_probe_review.zos` 查看基线；打开 `candidate_v35_native_fourgroup.zos` 查看四组试验。v35 的配置顺序为 500、200、800、近似无穷远、300、400、2000 mm；请以模型配置编辑器及相应 validation_configs JSON 为准。

v35 保存的评价函数保留最后 2 m 的 5 变量焦位目标。完整联合优化目标需查看父模型及相应脚本，不能仅从该 MFE 推断全部配置均已验收。

## 离线校验

```powershell
python tools/verify_bundle.py
```

该命令只校验归档文件哈希、Python 语法、JSON 可读性、README/指南本地链接及 ZIP 完整性，不导入仿真模块，不连接 Zemax。

## 独立复核示例

在另一份工作副本中运行，避免覆盖随包提供的历史证据。安装依赖后，可由脚本启动独立 ZOS-API 会话：

```powershell
$env:OPTICSTUDIO_DIR = 'D:\Zemax'
python design_35mm_f28/native_audit_multiconfig.py candidate_v35_native_fourgroup --curves --verify-schema --output-suffix _local_recheck --extra-curve-configs 3 4 5 6 7
```

脚本读取 `_schema.json`、`_validation_configs.json` 和 `.zos`，另存带 `_local_recheck` 后缀的审计及曲线。先检查许可、玻璃解析和配置顺序，再与逐配置证据来源比较。

部分历史脚本会立即执行或以固定候选名输出，且优化脚本通常拒绝覆盖已有候选。请先阅读脚本及参数；复核已有模型无需重跑全部历史优化。

基线温度证据见 `candidate_v31_thermal_probe_thermal_compensated_native.json` 与对应 WD500/WD200 `.zos`。其主动调焦和垫圈假设不能直接转用到新模型。
