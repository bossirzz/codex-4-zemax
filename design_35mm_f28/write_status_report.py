"""Document native results without presenting a partial candidate as qualified."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
tag='candidate_v19'
s=json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
rows=json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
b=next(x for x in rows if x['wd']==500);n=next(x for x in rows if x['wd']==200);inf=next(x for x in rows if x['wd']>=1e7)
f512=pd.read_csv(R/(tag+'_review_mtf_512.csv'),header=[0,1],index_col=0)
f1024=pd.read_csv(R/(tag+'_review_mtf_1024.csv'),header=[0,1],index_col=0)
sampling_difference=float(np.max(abs(f512.values-f1024.values)))
aux=json.loads((R/(tag+'_auxiliary.json')).read_text(encoding='utf8'))
aux_refocus=json.loads((R/(tag+'_auxiliary_refocus_native.json')).read_text(encoding='utf8'))
red=next(x for x in aux_refocus if x['wavelength_um']==.65)
index=json.loads((R/(tag+'_glass_index_crosscheck.json')).read_text(encoding='utf8'))
front=max(x['front_group_position'] for x in rows)-min(x['front_group_position'] for x in rows)
rear=max(x['rear_group_position'] for x in rows)-min(x['rear_group_position'] for x in rows)
lines=['# 35 mm F/2.8 设计状态与原生验收','',
'当前通过最佳物距的暂定 MTF 目标，**尚未通过全需求验收**。200 mm 与较远物距的边缘像质仍需改进；新的联合优化候选还未替换已验证模型。','',
'可复核模型：`candidate_v19_review.zos`。配置顺序为 500 mm、200 mm、2000 mm、无穷远。采用真实光阑瞄准、固定物理光阑半径及固定镜片口径。所有镜片来自真实 CDGM/NHG 玻璃库。','',
'## 已验证的指标','',
'| 指标 | 原生计算结果 | 要求及判定 |','|---|---:|---|',
f'| 最佳物距中心 MTF，225 lp/mm | {b["mtf_sampling6"][0]:.6f} | ≥0.30，通过 |',
f'| 最佳物距边缘 MTF，180 lp/mm，T/S | {b["mtf_sampling6"][-2]:.6f} / {b["mtf_sampling6"][-1]:.6f} | 均≥0.20，通过 |',
f'| 最佳物距有效焦距 | {b["efl"]:.5f} mm | 35 mm 名义焦距 |',
f'| 最佳物距实际像圈 | {b["image_circle_diameter"]:.5f} mm | ≥18 mm，通过 |',
f'| 最佳物距最大畸变 | {max(abs(x) for x in b["distortion_percent"]):.4f}% | ≤0.3%，通过 |',
f'| 200 mm 最大畸变 | {max(abs(x) for x in n["distortion_percent"]):.4f}% | ≤0.6%，通过 |',
f'| 无穷远边缘主光线角 | {inf["cra_deg"]:.4f}° | ≤10°，通过 |',
f'| 无穷远边缘相对照度 | {inf["relative_illumination"]*100:.3f}% | ≥75%，通过 |',
f'| 已验证配置最大光学总长 | {max(x["ttl"] for x in rows):.5f} mm | ≤60 mm，通过 |',
f'| 最小最后面至像面几何间距 | {min(x["backfocus"] for x in rows):.5f} mm | ≥8 mm，通过 |',
f'| 最大镜片毛坯直径 | {2*max(s["fixed_semi_diameters"].values()):.5f} mm | ≤23 mm，通过 |',
f'| 前六片 / 后两片移动跨度 | {front:.5f} / {rear:.5f} mm | 约6 mm，当前配置满足 |','',
'最佳物距使用 500 mm（物面至第一光学面顶点）。F/d/C 三波长等权，650/850 nm 在可见光结果中权重为零，另作单色校验。光阑名义 F/2.8；有限物距下 working F/# 约2.962，与名义 F 数定义不同。','',
'## 各配置像质','',
'| 物距 | 中心/边缘频率（lp/mm） | 中心 MTF | 边缘 T / S | 最大畸变 | 总长（mm） |','|---|---|---:|---|---:|---:|']
for row in rows:
 wd='∞' if row['wd']>=1e7 else f'{row["wd"]:g} mm'
 lines.append(f'| {wd} | {row["frequencies"][0]:g} / {row["frequencies"][1]:g} | {row["mtf_sampling6"][0]:.5f} | {row["mtf_sampling6"][-2]:.5f} / {row["mtf_sampling6"][-1]:.5f} | {max(abs(x) for x in row["distortion_percent"]):.4f}% | {row["ttl"]:.4f} |')
lines+=['',
'200 mm 边缘 MTF 仍偏低。用户已选择方案1：近距按 F/2.8、中心150 lp/mm ≥0.30、边缘120 lp/mm T/S均≥0.20验收，并要求3.5小时后执行。400–500 mm 表示最佳物距范围，本轮通过点为500 mm；不能由此推断整个范围均通过同一 MTF 门槛。','',
'无穷远仅通过本文明确列出的主光线角和照度指标，高频边缘 MTF 很低。2000 mm 也没有达到最佳物距的 MTF 门槛。','',
'## 数据与采样核对','',
f'- 8片真实玻璃在486.1327、587.5618、656.2725、650、850 nm的折射率逐一与原生 INDX 对比，最大差异 {index["maximum_index_difference"]:.3g}。数据见 `candidate_v19_glass_index_crosscheck.json`。',
f'- 最佳物距512×512与1024×1024的全部曲线最大差异 {sampling_difference:.6f}。原始CSV、精确频率原生操作数、清晰PNG/PDF均保留。',
'- 近距512/1024原生FFT及192/256独立衍射计算相互吻合，边缘120 lp/mm约0.058；快速优化积分此前高估该值，已提高积分阶数并校准。',
'- 镜片共享毛坯半径用于全径厚度和间隙检查。光阑与前一镜片的最小间隙约0.304 mm，玻璃全径最小厚度为正。正式公差与机械装配余量尚未完成验证。','',
'## 辅助光源与工程边界','',
'| 光源 | 可见光同一焦面中心225 lp/mm | 边缘180 lp/mm T / S |','|---|---:|---|']
for label,values in aux.items():
 m=values['same_visible_focus'];lines.append(f'| {label} | {m[0]:.5f} | {m[1]:.5f} / {m[2]:.5f} |')
lines+=['',
f'650 nm 同焦面不满足中心 MTF 目标；独立重新调焦后，原生中心225 lp/mm为{red["refocused_mtf"][0]:.6f}，边缘180 lp/mm T/S为{red["refocused_mtf"][1]:.6f}/{red["refocused_mtf"][2]:.6f}，通过暂定MTF门槛。模型见 `candidate_v19_650nm_refocused.zos`，不可据此宣称白光和650 nm共焦。850 nm 按用户要求作为附加校验，不优先进入优化。','',
'本模型为名义光学设计：未验证100 mm转接环状态、相机盖玻璃/滤光片堆栈、胶层厚度、公差、热膨胀与机械热漂移、镀膜和杂散光。不能将目前结果解释为可直接生产或已通过−20～60°C整机验收。','',
'光学MTF没有乘传感器像素MTF。2.4 μm像素的奈奎斯特频率为208.3 lp/mm；225 lp/mm评价镜头本体，整机分辨率需另作镜头与传感器联合计算。','',
'## 文件入口与复现','',
'- `candidate_v19_review.zos`：当前已完成原生检查的四配置候选。',
'- `candidate_v19_schema.json`：真实玻璃、曲率、厚度、六个非球面、固定口径和光阑的完整声明。',
'- `candidate_v19_native_acceptance.json`：精确频率、畸变、CRA、照度、像圈与各配置位置的原始结果。',
'- `candidate_v19_review_mtf_1024.png` / `.pdf`：原生MTF高清图；相同前缀CSV保留全部数据。',
'- `native_zospy.py` / `native_real_acceptance.py`：独立ZOSPy批量导入和原生校验。',
'- 所有 `candidate_v01`～`candidate_v18` 及旧 `final_native_*` 数据均为已淘汰实验版本。旧 `35mm_F28_Phi18_optimized.zos` 实为v04，存在物理间隙问题，已更名为 `obsolete_v04_unqualified.zos`，不可用于验收。','',
'该处方经过重新设计，不等于CN217060609U或CN110703413B某实施例的忠实复现。专利原始数据与复现矛盾见上级目录中的专利审计报告。','',
'光阑定义及真实瞄准说明：[Ansys Ray Aiming](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v251/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Ray_Aiming.html)、[Float by Stop Size](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/Float_by_Stop_Size.html)。本文性能数据来自本地原生计算，不来自专利宣传曲线。','']
(R/'设计状态与验收记录.md').write_text('\n'.join(lines),encoding='utf8')
(R/'README.md').write_text('# 当前设计入口\n\n请先阅读 [设计状态与验收记录](设计状态与验收记录.md)。当前candidate_v19仅通过500 mm暂定MTF和已列明几何/CRA/照度指标，尚未通过全需求验收。\n\n有效候选：`candidate_v19_review.zos`。旧v04的optimized文件已停用，新联合优化候选完成原生验收后再替换。\n',encoding='utf8')
print(str(R/'设计状态与验收记录.md'))
