"""Build the current project report after native warm and thermal checks finish."""
import json
import shutil
from pathlib import Path

import numpy as np

R = Path(__file__).resolve().parent
TAG = 'candidate_v31_thermal_probe'
warm = json.loads((R / (TAG + '_native_acceptance.json')).read_text())
thermal = json.loads((R / (TAG + '_thermal_compensated_native.json')).read_text())['cases']
schema = json.loads((R / (TAG + '_schema.json')).read_text())
required = {(500., -20.), (500., 60.), (200., -20.), (200., 60.)}
if not required.issubset({(q['wd'], q['temperature_c']) for q in thermal}):
    raise RuntimeError('Thermal validation incomplete; report not replaced.')
if not (R / (TAG + '_review.zos')).exists():
    raise RuntimeError('Review model is missing.')
checks = {}
for row in warm:
    wd = row['wd']
    values = row['mtf_sampling6']
    if wd in [500, 200, 800]:
        checks[f'warm_{wd:g}_mtf'] = min(values[:2]) >= .3 and min(values[2:]) >= .2
    checks[f'warm_{wd:g}_geometry'] = row['ttl'] <= 60 and row['image_circle_diameter'] >= 18 and row['relative_illumination'] >= .75
    if wd in [500, 200]:
        checks[f'warm_{wd:g}_distortion'] = max(abs(x) for x in row['distortion_percent']) <= (.3 if wd == 500 else .6)
    if wd >= 1e7:
        checks['warm_infinity_cra'] = row['cra_deg'] <= 10
for row in thermal:
    key = f'thermal_{row["wd"]:g}_{row["temperature_c"]:g}'
    if row['native_mtf'] is not None:
        checks[key + '_mtf'] = row['sampled_mtf_passed']
        checks[key + '_distortion'] = max(abs(x) for x in row['distortion_percent']) <= (.3 if row['wd'] == 500 else .6)
    else:
        checks[key + '_cra'] = row['cra_deg'] <= 10
    checks[key + '_geometry'] = row['ttl_mm'] <= 60 and row['image_circle_mm'] >= 18 and row['relative_illumination'] >= .75
positions = []
for row in warm:
    t = np.array(schema['thicknesses'])
    t[5], t[10], t[-1] = row['iris_gap'], row['gap'], row['backfocus']
    positions.append([sum(t), sum(t[6:]), sum(t[11:])])
positions += [q['group_positions_from_sensor_mm'] for q in thermal]
strokes = np.ptp(positions, axis=0)
checks['tested_group_strokes'] = bool(np.max(strokes) <= 6.01)
summary = dict(candidate=TAG, checks=checks, checked_subset_passed=all(checks.values()),
               complete_requirement_acceptance=False,
               temperature_compensation='Active three-group compensation; assumed thin-shell spacers at 23 ppm/K',
               group_strokes_mm=strokes.tolist(),
               missing_work=['300/2000 mm and continuous working-distance performance', '850 nm refocus and colour-fringe verification',
                             '100 mm adapter', 'manufacturing tolerances and stray light', 'thermal performance over all distances'])
if not summary['checked_subset_passed']:
    raise RuntimeError('At least one required subset check failed; report not replaced.')
(R / (TAG + '_acceptance_summary.json')).write_text(json.dumps(summary, indent=2), encoding='utf8')
current = R / '设计状态与验收记录.md'
archive = R / '设计状态与验收记录_v29_20261014.md'
if current.exists() and not archive.exists():
    shutil.copy2(current, archive)
lines = ['# 35 mm F/2.8 当前设计状态（2026-10-14，v31）', '',
    '**当前审查模型：candidate_v31_thermal_probe_review.zos。**常温500/200/800 mm的MTF已原生复核；500/200 mm冷、热补偿结果见下表。**尚未完成完整需求验收。**', '',
    '条件：裸像面、F/d/C等权、名义F/2.8、实际CDGM/NHG玻璃、真实光阑瞄准及固定镜片口径。四个常温配置顺序：500、200、800 mm、无穷远。所有数值来自Zemax原生计算；预测文件中的INCOMPLETE_CHECKPOINT不代表原生判定，工程安全余量与用户实际门槛分别记录。', '',
    '## 常温原生验收', '',
    '| 物距mm | 中心/离轴频率lp/mm | 中心MTF | 边缘T/S | 全离轴最低MTF | 最大畸变 | 总长mm |',
    '|---:|---|---:|---|---:|---:|---:|']
for row in warm:
    if row['wd'] >= 1e7:
        continue
    v = row['mtf_sampling6']
    fc, fe = row['frequencies']
    lines.append(f'| {row["wd"]:g} | {fc:g}/{fe:g} | {v[0]:.5f} | {v[-2]:.5f}/{v[-1]:.5f} | {min(v[2:]):.5f} | {max(abs(x) for x in row["distortion_percent"]):.4f}% | {row["ttl"]:.4f} |')
lines += ['', '500/200沿用用户确认的中心≥0.30、离轴≥0.20门槛；800为与最佳物距同频率的工程对照，未另行约定用户阈值。所有四个检查视场的T/S均保留。500和200原生512/1024曲线最大差分别约0.000437和0.000500。', '',
    '## 原生温度补偿', '',
    '**采用23 ppm/K垫圈假设和主动三组补偿。**按新镜片参数重新生成Make Thermal模型，未直接使用旧热态几何；未验证实际制造装配，也未验证全部物距的温度性能。', '',
    '| 物距mm | 温度°C | 中心MTF | 边缘T/S | 全离轴最低MTF | MTF判定 | 最大畸变 |',
    '|---:|---:|---:|---|---:|---|---:|']
for row in thermal:
    v = row['native_mtf']
    if v is None:
        continue
    lines.append(f'| {row["wd"]:g} | {row["temperature_c"]:g} | {v[0]:.5f} | {v[-2]:.5f}/{v[-1]:.5f} | {min(v[2:]):.5f} | {"通过" if row["sampled_mtf_passed"] else "未通过"} | {max(abs(x) for x in row["distortion_percent"]):.4f}% |')
inf = next(q for q in warm if q['wd'] >= 1e7)
lines += ['', f'与v29相比，500 mm冷态边缘从0.19661/0.19800提高到本表结果，冷态像质问题在该检查配置已解决。', '',
    '## 几何与其它检查', '',
    f'- 最佳物距焦距约{warm[0]["efl"]:.5f} mm；初始原生SystemData近轴后焦距约12.00 mm（区别于最后面至像面的几何距离）。',
    f'- 本轮原生配置最大光学总长{max([q["ttl"] for q in warm] + [q["ttl_mm"] for q in thermal]):.4f} mm；最大实体镜片直径{2*max(schema["fixed_semi_diameters"].values()):.4f} mm。',
    f'- 已检查配置组三组最大行程{strokes[0]:.4f}/{strokes[1]:.4f}/{strokes[2]:.4f} mm；连续轨迹尚未验收。',
    f'- 无穷远常温CRA {inf["cra_deg"]:.5f}°、相对照度{100*inf["relative_illumination"]:.3f}%；冷热无穷远角度/照度已独立检查。',
    '- 八种玻璃的5波长原生折射率与v29完全一致；本轮调整6个非球面系数、L5/L6厚度和焦位。', '',
    '## 未完成事项', '',
    '- 300 mm、2 m高频离轴及连续工作物距表现。多物距小范围联合优化尚未找到满意远距解，不将旧2m数据作为新v31原生验收。',
    '- 色边最终检查、650/850 nm辅助对焦、100 mm接圈、制造公差和杂散光；全部物距的温度性能。',
    '- 本模型属于原生验证的阶段候选，不能作为已通过全需求的生产设计。', '',
    '## 审查文件', '',
    f'- `{TAG}_review.zos`：四物距常温审查模型。',
    f'- `{TAG}_native_acceptance.json`：常温操作数、粗扫及0.25μm细扫焦。',
    f'- `{TAG}_thermal_compensated_native.json`：新的原生热模型验收结果。',
    f'- `{TAG}_thermal_compensated_native_WD500.zos` / `_WD200.zos`：20/−20/60°C温度配置。',
    f'- `{TAG}_best_near_comparison.png/.pdf`：原生1024高分辨率曲线。',
    f'- `optimization_progress_20261014_{TAG}_native.png/.pdf`：500/800曲线与温度补偿汇总图。',
    f'- `{TAG}_acceptance_summary.json`：已检查子集及未完成项。',
    '- `candidate_v29_plus800_review.zos`及v29历史记录：已保留，便于回退和对照。', '']
current.write_text('\n'.join(lines), encoding='utf8')
(R / 'README.md').write_text('# 当前设计入口\n\n'
    '最新：[设计状态与验收记录](设计状态与验收记录.md)（2026-10-14，v31）。\n\n'
    '**审查模型：candidate_v31_thermal_probe_review.zos。**常温500/200/800已原生验收；500/200冷热补偿见最新记录。完整需求尚未全部通过，300/2m离轴、公差等仍待完成。\n\n'
    'v29及其800配置保留作历史对照；v30多物距候选只是预测，未替换已验证模型。\n', encoding='utf8')
print(json.dumps(summary, ensure_ascii=False))
