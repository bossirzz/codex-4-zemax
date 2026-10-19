"""Combine native metrics with declared aperture and full-radius geometry checks."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import numpy as np
import pandas as pd
import optical_engine as e
import fast_otf_optimize as f
from optimize_best_conjugate import load_glasses

R = Path(__file__).resolve().parent

def run(tag):
    s = json.loads((R/(tag+'_schema.json')).read_text(encoding='utf8'))
    rows = json.loads((R/(tag+'_native_acceptance.json')).read_text(encoding='utf8'))
    load_glasses([(g['catalog'], g['name']) for g in s['glasses']])
    surfaces = f.SURFS + [2, 9]
    c = np.r_[s['curvatures'], np.zeros(24)]
    for a in s['aspheres']:
        j = 15 + 3*surfaces.index(a['surface'])
        c[j:j+3] = a['coefficients']
    t = np.array(s['thicknesses'], float)
    caps = np.array([s['fixed_semi_diameters'].get(str(i+1), 10.) for i in range(15)])
    tests = []
    physical = []
    def check(name, value, passed, requirement):
        tests.append(dict(name=name, value=value, passed=bool(passed), requirement=requirement))
    for row in rows:
        wd = row['wd']
        label = '∞' if wd >= 1e7 else f'{wd:g} mm'
        mtf = row['mtf_sampling6']
        if wd in (500, 200):
            check(label+'中心MTF', min(mtf[:2]), min(mtf[:2]) >= .30,
                  '225 lp/mm ≥0.30' if wd == 500 else '150 lp/mm ≥0.30')
            check(label+'边缘MTF T/S', mtf[-2:], min(mtf[-2:]) >= .20,
                  '180 lp/mm 均≥0.20' if wd == 500 else '120 lp/mm 均≥0.20')
            distortion = max(abs(x) for x in row['distortion_percent'])
            limit = .3 if wd == 500 else .6
            check(label+'畸变', distortion, distortion <= limit, f'≤{limit}%')
        check(label+'测试光圈', row['nominal_fno'], abs(row['nominal_fno']-2.8) <= .02,
              'ISFN F/2.8，数值核对容差±0.02')
        check(label+'真实瞄准', row['ray_aiming'], row['ray_aiming'] == 'Real', 'Real')
        check(label+'像圈', row['image_circle_diameter'], row['image_circle_diameter'] >= 18., '≥18 mm')
        check(label+'总长', row['ttl'], row['ttl'] <= 60., '≤60 mm')
        check(label+'后部几何间距', row['backfocus'], row['backfocus'] >= 8., '≥8 mm')
        if wd >= 1e7:
            check('无穷远CRA', row['cra_deg'], row['cra_deg'] <= 10., '≤10°')
            check('无穷远相对照度', row['relative_illumination'], row['relative_illumination'] >= .75, '≥75%')
        tt = t.copy()
        tt[5], tt[10], tt[-1] = row['iris_gap'], row['gap'], row['backfocus']
        clear = f.blank_clearance(c, tt, caps) + np.where(e.MEDIA[:14] >= 0, .05, 0.)
        _, ray_radii, _, invalid = f.quantities(c, tt, wd=wd,
            fields=np.array([0., 3., 6.3, 9.08]), fixed_stop_radius=row['stop_radius'], low_frequency=False)
        physical.append(dict(wd=wd, minimum_clearance_margin=float(clear.min()),
                             minimum_ray_radius_margin=float((caps-ray_radii).min()), invalid_rays=int(invalid)))
        check(label+'全径厚度及间隙', float(clear.min()), clear.min() >= -1e-5,
              '玻璃≥0.45 mm，空气≥0.10 mm，光阑相邻间隙≥0.30 mm；33个径向检查点')
        check(label+'采样光线通过固定口径', float((caps-ray_radii).min()),
              (caps-ray_radii).min() >= -1e-4 and invalid == 0, '采样光线无口径截断、无无效光线')
    best = next(x for x in rows if x['wd'] == 500)
    check('500 mm名义焦距', best['efl'], abs(best['efl']-35.) <= .05,
          '35 mm；本轮数值核对容差±0.05 mm')
    diameter = 2*max(s['fixed_semi_diameters'].values())
    check('最大镜片毛坯直径', diameter, diameter <= 23., '≤23 mm')
    positions = np.array([x['group_positions'] for x in rows])
    strokes = np.ptp(positions, axis=0)
    check('三组移动跨度', strokes.tolist(), np.max(strokes) <= 6.01, '约6 mm；本轮按≤6.01 mm核对')
    a = pd.read_csv(R/(tag+'_review_mtf_512.csv'), header=[0, 1], index_col=0)
    b = pd.read_csv(R/(tag+'_review_mtf_1024.csv'), header=[0, 1], index_col=0)
    sampling_difference = float(np.max(abs(a.values-b.values)))
    check('512/1024 MTF采样收敛', sampling_difference, sampling_difference <= .005, '全部曲线最大差异≤0.005')
    passed = all(x['passed'] for x in tests)
    unverified = ['400–500 mm整个最佳物距范围的逐点MTF', '100 mm转接状态',
                  '传感器盖玻璃及滤光片堆栈', '−20～60°C热性能', '制造与装配公差',
                  '镀膜、杂散光及色边', '650 nm独立调焦核验；850 nm仅为附加指标']
    result = dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                  checked_targets_passed=passed, full_requirements_passed=False,
                  native_model=str(R/(tag+'_review.zos')), tests=tests, physical_checks=physical,
                  unverified_requirements=unverified,
                  aperture_note='各配置明确设置名义F/2.8；光阑半径若变化，以MCE SDIA控制并记录，工作F数另外列出。')
    (R/(tag+'_scheduled_result.json')).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf8')
    lines = ['# 定时续跑结果', '',
             '本轮已检查项目全部通过；完整需求尚未验收。' if passed else '本轮存在未通过项目，候选不能认定为符合需求。', '',
             f'原生模型：`{tag}_review.zos`；原始数据：`{tag}_native_acceptance.json`。', '',
             result['aperture_note'], '', '| 项目 | 结果 | 要求 | 判定 |', '|---|---|---|---|']
    for x in tests:
        value = json.dumps(x['value'], ensure_ascii=False) if isinstance(x['value'], list) else str(x['value'])
        lines.append(f'| {x["name"]} | {value} | {x["requirement"]} | '+('通过' if x['passed'] else '**未通过**')+' |')
    lines += ['', '以下项目仍待继续验证：', ''] + ['- '+x for x in unverified]
    lines += ['', '2000 mm及无穷远的高频MTF见原始数据；本轮未替这两种物距另设用户未规定的MTF门槛。',
              '本模型经过重新设计，不是专利实施例的忠实复现。', '']
    (R/(tag+'_scheduled_result.md')).write_text('\n'.join(lines), encoding='utf8')
    print(json.dumps(dict(checked_targets_passed=passed, failed=[x['name'] for x in tests if not x['passed']]), ensure_ascii=False), flush=True)

if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('tag')
    q = a.parse_args()
    run(q.tag)
