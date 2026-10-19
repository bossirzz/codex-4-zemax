"""Summarize native optimization evidence without promoting unchecked candidates."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent


def run(tag, reference_evidence=None):
    load = lambda suffix: json.loads((R / (tag + suffix + '.json')).read_text())
    info = load('_run_summary')
    cases = load('_native_post_audit')
    geometry = load('_dense_geometry_audit')
    baseline = json.loads((R / 'candidate_v31_thermal_probe_native_acceptance.json').read_text())
    primary_pass = all(next(q for q in cases if q['wd'] == wd)['primary_mtf_passed'] for wd in [500, 200])
    structure_pass = all(geometry[key] for key in ['spacing_passed', 'length_passed',
        'lens_diameter_passed', 'geometric_back_clearance_passed', 'paraxial_bfl_passed', 'stroke_passed'])
    infinity = next(q for q in cases if q['true_infinity'])
    distortion_pass = all(max(abs(v) for v in next(q for q in cases if q['wd'] == wd)['distortion_percent']) <= limit
                          for wd, limit in [(500, .3), (200, .6)])
    illumination_pass = min(q['relative_illumination'] for q in cases) >= .75
    image_circle_pass = min(q['image_circle_diameter'] for q in cases) >= 18.
    convergence = {}
    for wd, suffix in [(500, '_review'), (200, '_near')]:
        low = pd.read_csv(R / (tag + suffix + '_mtf_512.csv'), header=[0, 1], index_col=0)
        high = pd.read_csv(R / (tag + suffix + '_mtf_1024.csv'), header=[0, 1], index_col=0)
        x = np.asarray(high.index, float)
        comparison = np.column_stack([np.interp(x, np.asarray(low.index, float), low.iloc[:, i])
                                      for i in range(low.shape[1])])
        convergence[str(wd)] = float(np.max(np.abs(np.asarray(high)-comparison)))
    summary = dict(tag=tag, actual_native_optimizer=True, native_run=info,
        room_temperature_primary_mtf_passed=primary_pass, geometry_passed=structure_pass,
        primary_distortion_passed=distortion_pass, checked_image_circle_passed=image_circle_pass,
        infinity_cra_passed=infinity['cra_deg'] <= 10.,
        checked_relative_illumination_passed=illumination_pass,
        room_temperature_checked_subset_passed=all([primary_pass, structure_pass, distortion_pass,
            illumination_pass, image_circle_pass, infinity['cra_deg'] <= 10.]),
        sampling_512_1024_max_difference=convergence,
        previous_validated_baseline='candidate_v31_thermal_probe_review.zos',
        thermal_result_carried_over=False, complete_requirements_passed=False,
        promotion='experimental candidate; v31 retained as previously validated baseline')
    (R / (tag + '_acceptance_summary.json')).write_text(json.dumps(summary, indent=2), encoding='utf8')
    reduction = 100*(1-info['current_merit']/info['initial_merit'])
    lines = [f'# Zemax 原生优化复核（2026-10-16，{tag}）', '',
        '**本次通过 ZOSPy 调用 OpticStudio 内置 DLS，搜索步骤由 Zemax 执行。**', '',
        f"- 运行状态：{info['status']}；API Succeeded={info['succeeded']}；错误信息：{info['error']}。",
        f"- {info['variables']} 个变量，{info['targets']} 个目标，{info['cores']} 核；耗时 {info['elapsed_seconds']:.1f} 秒。",
        f"- 本次直接输入模型：`{info.get('source', '见运行记录')}`。",
        f"- 同一评价函数从 {info['initial_merit']:.10f} 降到 {info['current_merit']:.10f}，下降 {reduction:.2f}%。",
        '- 各轮权重不同，轮与轮之间不能用评价函数绝对值比较像质。',
        '- 条件：20°C、裸像面、F/d/C 等权、真实玻璃及光阑瞄准；复核时不另行调焦。', '',
        '## 独立原生 MTF', '',
        '| 物距 mm | 中心/离轴 lp/mm | 中心 | 边缘 T/S | 全离轴最低 | 名义 F/# |',
        '|---:|---|---:|---|---:|---:|']
    for q in sorted(cases, key=lambda q:q['wd']):
        mtf = q['mtf_sampling6']
        wd = '∞' if q['true_infinity'] else f"{q['wd']:g}"
        fc, fe = q['frequencies']
        lines.append(f'| {wd} | {fc:g}/{fe:g} | {mtf[0]:.5f} | {mtf[-2]:.5f}/{mtf[-1]:.5f} | {min(mtf[2:]):.5f} | {q["nominal_fno"]:.5f} |')
    lines += ['', '用户确认的主验收值：500 mm 中心225≥0.30、边缘180≥0.20；200 mm 中心150≥0.30、边缘120≥0.20。'
        '中间视场另按≥0.20检查；300/400/800/2000 mm所加优化目标是工程目标，不能当作新增的用户验收条件。', '',
        f'常温主验收及中间场检查：**{"通过" if primary_pass else "未通过"}**。', '',
        f'512/1024曲线最大差：500 mm {convergence["500"]:.6f}；200 mm {convergence["200"]:.6f}。', '',
        '## 与已验证基线比较', '',
        '| 物距 mm | v31 中心/全离轴最低 | 本轮中心/全离轴最低 |',
        '|---:|---|---|']
    for wd in [500, 200, 800]:
        a = next(q for q in baseline if q['wd'] == wd)['mtf_sampling6']
        b = next(q for q in cases if q['wd'] == wd)['mtf_sampling6']
        lines.append(f'| {wd} | {a[0]:.5f}/{min(a[2:]):.5f} | {b[0]:.5f}/{min(b[2:]):.5f} |')
    if reference_evidence:
        reference = json.loads((R / reference_evidence).read_text())
        lines += ['', '## 与共同对照模型比较', '',
            f'对照：`{reference["candidate"]}`。它不一定是本次直接输入；直接输入见上方运行记录。表内各物距沿用相同的空间频率定义；以本轮实际名义光圈为准。', '',
            '| 物距 mm | 对照中心/最差离轴 | 本轮中心/最差离轴 |', '|---:|---|---|']
        comparison = []
        for wd in [200, 300, 400, 500, 800, 2000]:
            a = next(q for q in reference['cases'] if q['wd'] == wd)['mtf_sampling6']
            b = next(q for q in cases if q['wd'] == wd)['mtf_sampling6']
            lines.append(f'| {wd} | {a[0]:.5f}/{min(a[2:]):.5f} | {b[0]:.5f}/{min(b[2:]):.5f} |')
            comparison.append(dict(wd=wd, source_centre=a[0], source_worst_offaxis=min(a[2:]),
                final_centre=b[0], final_worst_offaxis=min(b[2:])))
        summary['reference_candidate_comparison'] = comparison
        summary['reference_candidate'] = reference['candidate']
        (R / (tag + '_acceptance_summary.json')).write_text(json.dumps(summary, indent=2), encoding='utf8')
    lines += ['', '## 独立几何检查', '',
        f'- 采用每个间隙501个半径位置检查；结构检查{"通过" if structure_pass else "未全部通过"}。',
        f'- 最大镜片直径 {geometry["maximum_lens_diameter_mm"]:.5f} mm；最大光学总长 {max(q["ttl_mm"] for q in geometry["cases"]):.5f} mm。',
        '- 各组行程 mm：' + ' / '.join(f'{x:.5f}' for x in geometry['group_strokes_mm']) + '。',
        f'- 无穷远原生CRA {infinity["cra_deg"]:.6f}°（上限10°）；全部检查配置相对照度最低 {100*min(q["relative_illumination"] for q in cases):.3f}%。',
        '- 主验收点原生最大畸变：' + '；'.join(
            f'{wd} mm {max(abs(v) for v in next(q for q in cases if q["wd"] == wd)["distortion_percent"]):.5f}%'
            for wd in [500, 200]) + '。',
        f'- 全部检查配置原生最小像圈直径 {min(q["image_circle_diameter"] for q in cases):.5f} mm。',
        f'- 全口径间隙/材料厚度检查：{geometry["spacing_passed"]}；详情见 `{tag}_dense_geometry_audit.json`。', '',
        '## 判定与文件', '',
        '**保留 v31 为此前已验证基线；本轮作为原生优化候选保存。完整需求尚未全部通过。**'
        '原生 v31 的热补偿验收不能沿用到本轮；连续物距、温度、公差、色边和杂散光仍需检查。', '',
        f'- `{tag}.zos`：优化器实际保存的7配置模型，包含评价函数和变量。',
        f'- `{tag}_run_summary.json` / `.log`：原生优化运行证据。',
        f'- `{tag}_native_post_audit.json`：优化结果的高采样、精确频率独立复核。',
        f'- `{tag}_native_system.txt`：原生系统报告。',
        f'- `{tag}_native_review.png/.pdf`：原生1024曲线及与v31的比较。',
        '', '首批 candidate_v32_native_contrast_b 虽然原生评价函数下降14.94%，但主验收MTF下降、L5边缘材料厚度约0.430 mm低于0.45 mm工程下限，未采用。'
        '后续批次强化了MTF和径向厚度约束。无穷远使用已验证可被OpticStudio识别的API表示，避免IEEE Infinity使离轴光线无效。', '']
    if info['status'] == 'cancelled_at_time_limit':
        lines.insert(6, '- 本批按设定时限结束，保存当前迭代结果；没有声称完成全部预设轮次或已经收敛。')
    target = R / (tag + '_review.md')
    target.write_text('\n'.join(lines), encoding='utf8')
    print(target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    parser.add_argument('--reference-evidence')
    args = parser.parse_args()
    run(args.tag, args.reference_evidence)
