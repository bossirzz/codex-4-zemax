# 当前设计入口

最新：[设计状态与验收记录](设计状态与验收记录.md)（2026-10-16；新增v33/v34联合优化及v35四组调焦试验）。

**审查模型：candidate_v31_thermal_probe_review.zos。**常温500/200/800已原生验收；500/200冷热补偿见最新记录。完整需求尚未全部通过，300/2m离轴、公差等仍待完成。

**最新试验模型：[candidate_v35_native_fourgroup.zos](candidate_v35_native_fourgroup.zos)。**通过ZOSPy实际调用Zemax内置DLS，先完成68变量及78变量的两批联合优化，再固定v34处方，分别对300/400/2000 mm执行5变量、10轮焦位优化，测试L6独立浮动的四组方案。

300 mm最差离轴180 lp/mm从本轮共同对照v32焦位版的0.12112提高到0.13751；400 mm从0.18614提高到0.19830。2 m新增浮动试验未改善最差离轴，已恢复v34焦位，最终约0.03026。300/400/2m按0.20作扩展工程比较，尚未全部达到。500/200 mm物理配置与v34相同，主MTF通过；新四组机构、温度、公差和连续调焦仍未验证，未替换v31基线。

见[v35试验与复核记录](candidate_v35_native_fourgroup_review.md)、[高分辨率物距MTF图](candidate_v35_native_fourgroup_range_review.png)、[v34联合优化记录](candidate_v34_native_thickness_review.md)。逐配置原生证据来源见v35记录；500/200/800及无穷远沿用物理配置相同的v34证据。

前序v32同样通过ZOSPy实际调用Zemax内置DLS；模型：[candidate_v32_native_contrast_c.zos](candidate_v32_native_contrast_c.zos)。500/200 mm主MTF、800 mm工程检查和独立几何复核通过；300/400/2000 mm高频像质仍弱，新参数尚未重新验收温度性能。

此前焦位候选：[candidate_v32_native_focus_refine.zos](candidate_v32_native_focus_refine.zos)，见[焦位修正版记录](candidate_v32_native_focus_refine_review.md)。三个弱物距各完成10轮、4变量的原生DLS；400 mm中心225从0.132提高到0.324，2 m中心从0.237提高到0.313，但离轴仍未达到扩展工程目标。500/200/800 mm及无穷远物理配置与v32c相同，保存文件已独立核对；v31仍保留作为已验证基线。

见[原生优化复核记录](candidate_v32_native_contrast_c_review.md)、[高分辨率MTF比较图](candidate_v32_native_contrast_c_native_review.png)和[原生运行日志](candidate_v32_native_contrast_c.log)。v31仍保留为此前已验证的基线。

v29及其800配置保留作历史对照；v30多物距候选只是预测，未替换已验证模型。
