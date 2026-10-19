"""Publish only native curves and native temperature acceptance results."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--tag', default='candidate_v29')
args = parser.parse_args()
tag = args.tag
font = FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
fig = plt.figure(figsize=(14, 10))
grid = fig.add_gridspec(2, 2, height_ratios=[1.25, .85])
colors = ['#202124', '#4285f4', '#15966e', '#d23d35']
nominal = json.loads((ROOT / (tag + '_native_acceptance.json')).read_text())
row800 = next((q for q in nominal if q['wd'] == 800), None)
if row800 is None:
    row800 = json.loads((ROOT / (tag + '_800_native.json')).read_text())[0]
values800 = row800.get('mtf_sampling6', row800.get('native_mtf'))
curve800 = ('candidate_v29_plus800_review_wd800_mtf_1024.csv' if tag == 'candidate_v29'
            else tag + '_wd800_mtf_1024.csv')

for j, (wd, filename, values) in enumerate([
    (500, tag + '_review_mtf_1024.csv', next(x for x in nominal if x['wd'] == 500)['mtf_sampling6']),
    (800, curve800, values800),
]):
    ax = fig.add_subplot(grid[0, j])
    df = pd.read_csv(ROOT / filename, header=[0, 1], index_col=0)
    for field, h in enumerate([0., 3., 6.3, 9.08]):
        for direction, style in enumerate(['-', '--']):
            ax.plot(np.asarray(df.index, float), df.iloc[:, 2 * field + direction], style,
                    color=colors[field], lw=1.6, label=f'{h:g} mm ' + ['T', 'S'][direction])
    ax.scatter([225, 180, 180], [values[0], values[-2], values[-1]],
               c=[colors[0], colors[3], '#e7891d'], s=45, zorder=5)
    ax.text(.035, .045, f'中心225：{values[0]:.4f}\n边缘180：T {values[-2]:.4f} / S {values[-1]:.4f}',
            transform=ax.transAxes, fontproperties=font, fontsize=11,
            bbox=dict(facecolor='white', edgecolor='#ddd', alpha=.9))
    ax.set_title(f'物距 {wd} mm · 20°C · 原生1024×1024', fontproperties=font, fontsize=13)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1)
    ax.set_xlabel('空间频率（lp/mm）', fontproperties=font)
    ax.set_ylabel('光学 MTF', fontproperties=font)
    ax.grid(alpha=.2)
    ax.legend(prop=font, loc='upper right', fontsize=9)

ax = fig.add_subplot(grid[1, :])
thermal = json.loads((ROOT / (tag + '_thermal_compensated_native.json')).read_text())['cases']
bars = []
for wd in [500, 200]:
    for temp in [-20, 20, 60]:
        if temp == 20:
            mtf = next(x for x in nominal if x['wd'] == wd)['mtf_sampling6']
        else:
            mtf = next(x for x in thermal if x['wd'] == wd and x['temperature_c'] == temp)['native_mtf']
        bars.append((wd, temp, min(mtf[2:])))
x = np.arange(len(bars))
heights = [b[2] for b in bars]
ax.bar(x, heights, color=['#d23d35' if h < .2 else '#15966e' for h in heights], width=.6)
ax.axhline(.2, color='#555', linestyle='--', label='离轴门槛 0.20')
for i, value in enumerate(heights):
    ax.text(i, value + .001, f'{value:.4f}', ha='center', fontsize=11)
ax.set_xticks(x, [f'{wd} mm\n{temp:+g}°C' for wd, temp, _ in bars], fontproperties=font)
ax.set_ylim(.17, .235)
ax.set_ylabel('所有检查离轴视场的最低 MTF', fontproperties=font)
ax.set_title('原生温度补偿验收：500 mm @180 lp/mm；200 mm @120 lp/mm',
             fontproperties=font, fontsize=13)
ax.grid(axis='y', alpha=.2)
ax.legend(prop=font, loc='upper right')
fig.suptitle('本轮确认结果 · 名义F/2.8 · F/d/C等权 · 裸像面', fontproperties=font, fontsize=16)
thermal_status = ('500/200 mm检查冷热配置均通过。' if all(q['sampled_mtf_passed'] for q in thermal if q['native_mtf'] is not None)
                  else '500 mm冷态未通过。')
fig.text(.5, .025, '温度结果采用23 ppm/K垫圈假设及主动三组补偿；' + thermal_status + '完整需求尚未全部验收。',
         ha='center', fontproperties=font, fontsize=10)
fig.tight_layout(rect=[0, .06, 1, .95])
target = ROOT / ('optimization_progress_20261014_native.png' if tag == 'candidate_v29'
                 else f'optimization_progress_20261014_{tag}_native.png')
fig.savefig(target, dpi=220)
fig.savefig(target.with_suffix('.pdf'))
print(target)
