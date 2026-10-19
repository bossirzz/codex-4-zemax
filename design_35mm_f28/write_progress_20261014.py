"""Refresh the project entry using verified native data, keeping history."""
import json
import shutil
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent
original = R / '设计状态与验收记录.md'
archive = R / '设计状态与验收记录_v19历史.md'
if original.exists() and not archive.exists():
    shutil.copy2(original, archive)
nominal = json.loads((R / 'candidate_v29_native_acceptance.json').read_text())
extra = json.loads((R / 'candidate_v29_800_native.json').read_text())[0]
thermal = json.loads((R / 'candidate_v29_thermal_compensated_native.json').read_text())['cases']
s = json.loads((R / 'candidate_v29_schema.json').read_text())
lines = [
    '# 35 mm F/2.8 当前设计状态（2026-10-14）', '',
    '当前可复核模型为 **candidate_v29_plus800_review.zos**，包含500、200、2000、无穷远、800 mm共5个配置。镜片形状和玻璃保持v29，新增800 mm三组对焦姿态已经独立原生复核。**完整需求尚未全部通过。**', '',
    '所有下表数值来自Zemax原生操作数；可见光F/d/C等权、20°C、裸像面、名义F/2.8、真实光阑瞄准、实际CDGM/NHG玻璃库和固定实体镜片口径。无相机盖玻璃数据，按用户要求采用裸像面。', '',
    '## 常温成像', '',
    '| 物距mm | 中心/离轴频率lp/mm | 中心MTF | 边缘T/S | 全离轴最小MTF | 最大畸变 | 总长mm |',
    '|---:|---|---:|---|---:|---:|---:|',
]
for row in sorted(nominal + [extra], key=lambda x: x['wd']):
    mtf = row.get('mtf_sampling6', row.get('native_mtf'))
    label = '∞' if row['wd'] >= 1e7 else f'{row["wd"]:g}'
    fc, fe = row['frequencies']
    lines.append(f'| {label} | {fc:g}/{fe:g} | {mtf[0]:.5f} | {mtf[-2]:.5f}/{mtf[-1]:.5f} | {min(mtf[2:]):.5f} | {max(abs(x) for x in row["distortion_percent"]):.4f}% | {row["ttl"]:.4f} |')
lines += ['',
    '500 mm按中心225 lp/mm≥0.30、全部检查离轴180 lp/mm≥0.20通过；200 mm按中心150≥0.30、全部检查离轴120≥0.20通过。800 mm也达到与最佳物距相同的MTF门槛，但其他物距并未另外约定用户验收阈值。2m/无穷远的高频离轴像质仍较低，不能声称整个物距范围达到20MP目标。', '',
    '## 温度补偿原生验收', '',
    '以下采用**23 ppm/K垫圈假设和主动三组温度补偿**；这属于假定装配模型，未验证实际制造装配。每个case按名义F/2.8检查。', '',
    '| 物距mm | 温度°C | 中心MTF | 边缘T/S | 全离轴最小MTF | MTF判定 | 最大畸变 |',
    '|---:|---:|---:|---|---:|---|---:|',
]
for row in thermal:
    mtf = row['native_mtf']
    if mtf is None:
        continue
    lines.append(f'| {row["wd"]:g} | {row["temperature_c"]:g} | {mtf[0]:.5f} | {mtf[-2]:.5f}/{mtf[-1]:.5f} | {min(mtf[2:]):.5f} | {"通过" if row["sampled_mtf_passed"] else "未通过"} | {max(abs(x) for x in row["distortion_percent"]):.4f}% |')
lines += ['',
    '500 mm冷态广起点与0.5μm扫焦后的最弱预测仍只有约0.20157，而原生与快速模型存在约0.004偏差，因此不能通过该预测宣称冷态合格。需要镜面联合优化。', '',
    '## 已核验几何与辅助指标', '',
    f'- 名义最佳物距焦距35.045 mm；当前模型最大实体镜片直径{2 * max(s["fixed_semi_diameters"].values()):.4f} mm。',
    '- 常温5配置最大总长58.0512 mm；当前原生温度配置最大总长58.2591 mm。',
    '- 当前常温/温度/800 mm姿态合并后三组最大行程约5.8521/5.5288/5.8601 mm。尚未验证连续对焦轨迹各中间位置。',
    '- 常温实际像圈均≥18 mm、相对照度均≥75%；无穷远CRA约9.9608°。冷/热无穷远CRA约9.9789°/9.9423°。',
    '- 最佳物距原生SystemData的近轴后焦距12.00442 mm；它与最后光学面至传感器的几何距离不同。', '',
    '## 仍需完成', '',
    '- 300 mm和2m离轴高频像质，以及300～800 mm连续对焦性能。',
    '- 500 mm／−20°C补偿后MTF。',
    '- 色边的最终核验、650/850 nm在最新模型上的辅助对焦校验、100 mm接圈配置、公差和杂散光。',
    '- 多物距与温态镜面候选目前属于数值预测，只有通过独立原生验收后才替换当前模型。', '',
    '## 证据文件', '',
    '- `candidate_v29_native_acceptance.json`：常温4配置原生操作数。',
    '- `candidate_v29_800_native.json`：800 mm独立原生核验。',
    '- `candidate_v29_plus800_review_provenance.json`：新增配置来源。',
    '- `candidate_v29_thermal_compensated_native.json`：原生补偿温态结果。',
    '- `optimization_progress_20261014_native.png/.pdf`：原生1024曲线及温态MTF图。',
    '- `candidate_v29_best_near_comparison.png/.pdf`：500/200原生1024曲线。',
    '- `far_focus_diagnostic_v29_结论.md`：冻结镜片的远距数值诊断，预测与原生明确区分。',
    '- `candidate_v29_far_frozen_native.json`：2m三组调焦候选未采纳的原生证据。',
    '- `设计状态与验收记录_v19历史.md`：旧v19记录，仅供历史追溯。', '',
]
original.write_text('\n'.join(lines), encoding='utf8')
(R / 'README.md').write_text(
    '# 当前设计入口\n\n'
    '最新记录：[设计状态与验收记录](设计状态与验收记录.md)（2026-10-14）。\n\n'
    '可复核模型：**candidate_v29_plus800_review.zos**。500/200 mm常温MTF已通过，新增800 mm配置已原生复核；2m高频离轴和500 mm冷态尚未解决，完整需求仍未通过。\n\n'
    '当前多物距/温态镜面优化候选须通过独立原生验收后才能替换。旧v19及obsolete_v04为历史文件。\n', encoding='utf8')
print(original)
