"""Combine fresh focus audits with explicitly identified unchanged parent evidence."""
import copy
import json
from pathlib import Path

R = Path(__file__).resolve().parent
TAG = 'candidate_v35_native_fourgroup'
PARENT = 'candidate_v34_native_thickness'


def read(tag, suffix):
    return json.loads((R / (tag + suffix + '.json')).read_text())


schema = read(TAG, '_schema')
parent_schema = read(PARENT, '_schema')
poses = read(TAG, '_validation_configs')
parent_poses = read(PARENT, '_validation_configs')
checked = read(TAG, '_native_schema_check')
geometry = read(TAG, '_dense_geometry_audit')
fresh = read(TAG, '_native_post_audit')
parent = read(PARENT, '_native_post_audit')
runs = read(TAG, '_runs')
keys = ['curvatures', 'thicknesses', 'aspheres', 'glasses', 'media',
        'fixed_semi_diameters', 'fixed_stop_radius', 'fields']
common_unchanged = all(schema[k] == parent_schema[k] for k in keys)
unchanged = []
for index in range(4):
    a, b = poses[index], parent_poses[index]
    assert all(a[k] == value for k, value in b.items())
    assert a['extra_thicknesses'] == {'9': parent_schema['thicknesses'][8]}
    unchanged.append(index+1)
assert common_unchanged and checked['native_prescription_verified']
assert sorted(q['configuration'] for q in fresh) == [5, 6, 7]
assert len(runs) == 3 and all(q['succeeded'] and q['error'] == 'None' for q in runs)

cases = []
for cfg in range(1, 8):
    inherited = cfg in unchanged
    source = parent if inherited else fresh
    item = copy.deepcopy(next(q for q in source if q['configuration'] == cfg))
    item['evidence_source'] = (PARENT if inherited else TAG) + '_native_post_audit.json'
    item['native_values_inherited_from_unchanged_parent_configuration'] = inherited
    item['focus_group_starts'] = schema['focus_group_starts']
    shape = next(q for q in geometry['cases'] if q['wd'] == item['wd'])
    item['group_positions_from_sensor_mm'] = shape['group_positions_from_sensor_mm']
    cases.append(item)

primary = all(next(q for q in cases if q['wd'] == wd)['primary_mtf_passed'] for wd in [500, 200])
structure = all(geometry[key] for key in ['spacing_passed', 'length_passed',
    'lens_diameter_passed', 'geometric_back_clearance_passed', 'paraxial_bfl_passed', 'stroke_passed'])
infinity = next(q for q in cases if q['true_infinity'])
subset = primary and structure and infinity['cra_deg'] <= 10
subset = subset and min(q['relative_illumination'] for q in cases) >= .75
subset = subset and min(q['image_circle_diameter'] for q in cases) >= 18
subset = subset and all(max(abs(x) for x in next(q for q in cases if q['wd'] == wd)['distortion_percent']) <= bound
                        for wd, bound in [(500, .3), (200, .6)])
evidence = dict(candidate=TAG, parent=PARENT, actual_native_optimizer=True, runs=runs,
    common_optical_prescription_unchanged=common_unchanged,
    unchanged_parent_configurations=unchanged, native_saved_model_schema_check=checked,
    cases=cases, geometry=geometry, room_temperature_primary_mtf_passed=primary,
    room_temperature_checked_subset_passed=bool(subset),
    thermal_result_carried_over=False, complete_requirements_passed=False,
    mechanism_validated=False, previous_validated_baseline='candidate_v31_thermal_probe_review.zos')
(R / (TAG + '_review_evidence.json')).write_text(json.dumps(evidence, indent=2), encoding='utf8')

lines = ['# Zemax 四组浮动调焦试验（2026-10-16）', '',
    '**ZOSPy 实际调用 OpticStudio 内置 DLS，固定 v34 镜片处方，增加 L6 独立浮动。完整需求尚未通过。**', '',
    '分组：L1–L3及光阑 / L4–L5 / L6 / L7–L8。运动机构和连续轨迹尚未验证；v31仍保留为此前已验证基线。', '',
    '每个试验物距采用5个变量：三个原有空气间隔、光阑半径及新增加的L5–L6空气间隔。'
    '玻璃、曲率、非球面系数和材料厚度固定；500/200/800 mm及无穷远的物理配置固定。', '',
    '## 实际原生运行', '',
    '| 物距 mm | 变量/目标 | 请求轮次 | 耗时秒 | 到达时限 | API成功 | 采用焦位 |',
    '|---:|---|---:|---:|---|---|---|']
for q in runs:
    lines.append(f'| {q["wd"]} | {q["variables"]}/{q["targets"]} | {q["requested_cycles"]} | '
                 f'{q["elapsed_seconds"]:.1f} | {q["time_limit_reached"]} | {q["succeeded"]} | {q["selected"]} |')
lines += ['', '未到时限且API正常结束的批次完成所设轮次；采用判据为原生高采样最差余量改善且几何检查通过。'
          '评价函数下降不等于全部MTF通过。', '',
    '## 保存文件独立复核', '',
    '300/400/2000 mm从新保存文件重新打开，执行原生1024采样及完整MTF曲线。'
    '500/200/800 mm与无穷远继承v34原生证据；公共处方和这四个配置的全部焦位数值相同，'
    '新保存文件处方及七个配置焦位已通过原生核对。每个配置的证据来源单列在JSON中。', '',
    '| 物距 mm | 中心/离轴 lp/mm | 中心MTF | 边缘T/S | 最差离轴 | v34最差离轴 | 名义F/# |',
    '|---:|---|---:|---|---:|---:|---:|']
for q in sorted(cases, key=lambda q:q['wd']):
    old = next(v for v in parent if v['wd'] == q['wd'])
    mtf = q['mtf_sampling6']
    label = '∞' if q['true_infinity'] else f'{q["wd"]:g}'
    fc, fe = q['frequencies']
    lines.append(f'| {label} | {fc:g}/{fe:g} | {mtf[0]:.5f} | {mtf[-2]:.5f}/{mtf[-1]:.5f} | '
                 f'{min(mtf[2:]):.5f} | {min(old["mtf_sampling6"][2:]):.5f} | {q["nominal_fno"]:.5f} |')
lines += ['', f'500/200 mm主MTF门槛及中间场检查：**{"通过" if primary else "未通过"}**。'
    '300/400/800/2000 mm按中心225、离轴180 lp/mm作扩展工程比较，0.20不是新增用户验收条件。', '',
    f'独立501半径位置几何检查：{"通过" if structure else "未全部通过"}；最大总长'
    f'{max(q["ttl_mm"] for q in geometry["cases"]):.5f} mm；最大镜片直径{geometry["maximum_lens_diameter_mm"]:.5f} mm。', '',
    '四组最大行程 mm：' + ' / '.join(f'{x:.5f}' for x in geometry['group_strokes_mm']) + '。', '',
    '新增四组模型仍需温度、连续调焦、机构、制造公差、色边和杂散光验证，不能沿用v31的热补偿结论。', '',
    '## 文件', '',
    f'- [{TAG}.zos]({TAG}.zos)：7配置试验模型，MFE保留最后2 m的5变量焦位目标。',
    f'- [{TAG}_runs.json]({TAG}_runs.json)：原生DLS过程与取舍。',
    f'- [{TAG}_review_evidence.json]({TAG}_review_evidence.json)：复核结果及逐配置证据来源。',
    f'- [{TAG}_range_review.png]({TAG}_range_review.png)：高分辨率原生MTF曲线。',
    f'- [{PARENT}_review.md]({PARENT}_review.md)：父模型78变量联合优化及主验收证据。', '']
(R / (TAG + '_review.md')).write_text('\n'.join(lines), encoding='utf8')
print(json.dumps(dict(primary_mtf_passed=primary, checked_subset_passed=bool(subset),
                     common_unchanged=common_unchanged, strokes=geometry['group_strokes_mm'])))
