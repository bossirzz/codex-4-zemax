"""Plot native FFT curves and exact acceptance operands against the v31 baseline."""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd

R = Path(__file__).resolve().parent


def run(tag):
    current = json.loads((R / (tag + '_native_post_audit.json')).read_text())
    baseline = json.loads((R / 'candidate_v31_thermal_probe_native_acceptance.json').read_text())
    font = FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
    colors = ['#202124', '#3079c1', '#13866c', '#c43d42']
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for row, (wd, suffix) in enumerate([(500, '_review'), (200, '_near')]):
        new = next(q for q in current if q['wd'] == wd)
        old = next(q for q in baseline if q['wd'] == wd)
        df = pd.read_csv(R / (tag + suffix + '_mtf_1024.csv'), header=[0, 1], index_col=0)
        frequency = np.asarray(df.index, float)
        for field, height in enumerate([0., 3., 6.3, 9.08]):
            for k, style in enumerate(['-', '--']):
                axes[row, 0].plot(frequency, df.iloc[:, 2*field+k], style, c=colors[field],
                    lw=1.6, label=f'{height:g} mm ' + ['T', 'S'][k])
        axes[row, 0].set(xlim=(0, 300), ylim=(0, 1))
        axes[row, 0].set_title(f'物距 {wd} mm · 原生 1024×1024 全视场曲线', fontproperties=font)
        axes[row, 0].set_xlabel('空间频率（lp/mm）', fontproperties=font)
        axes[row, 0].legend(ncol=2, prop=font, fontsize=9)
        x = [0., 3., 6.3, 9.08]
        for k, label in enumerate(['T', 'S']):
            axes[row, 1].plot(x, old['mtf_sampling6'][k::2], '--', c='#999999', lw=1.3,
                marker=['o', 's'][k], label='v31 ' + label)
            axes[row, 1].plot(x, new['mtf_sampling6'][k::2], '-', c=['#287caa', '#c43d42'][k],
                lw=2, marker=['o', 's'][k], label='本轮原生优化 ' + label)
        axes[row, 1].scatter(x, [.3, .2, .2, .2], marker='x', c='#222', s=80, label='验收及中间场检查下限')
        verdict = '通过' if new['primary_mtf_passed'] else '未通过'
        axes[row, 1].set_title(f'精确验收频率 · {verdict}', fontproperties=font)
        axes[row, 1].set(xlim=(-.4, 9.5), ylim=(.1, .5), xticks=x)
        axes[row, 1].set_xlabel('近轴像高（mm）', fontproperties=font)
        fc, fe = new['frequencies']
        axes[row, 1].text(.03, .97,
            f'中心 {fc:g} lp/mm：{new["mtf_sampling6"][0]:.4f}\n'
            f'离轴 {fe:g} lp/mm 最低：{min(new["mtf_sampling6"][2:]):.4f}\n'
            f'名义光圈 F/{new["nominal_fno"]:.4f}',
            va='top', transform=axes[row, 1].transAxes, fontproperties=font, fontsize=10)
        axes[row, 1].legend(prop=font, loc='upper right', fontsize=9)
        for ax in axes[row]:
            ax.grid(alpha=.2)
            ax.set_ylabel('光学 MTF', fontproperties=font)
    fig.suptitle('Zemax 内置 DLS 优化后的独立复核 · F/d/C 等权 · 裸像面', fontproperties=font, fontsize=16)
    fig.text(.5, .02, '仅展示常温主验收点；冷热、连续物距、公差等不能沿用旧版本的验证结论。',
             ha='center', fontproperties=font, fontsize=10)
    fig.tight_layout(rect=[0., .04, 1., .95])
    target = R / (tag + '_native_review.png')
    fig.savefig(target, dpi=220)
    fig.savefig(R / (tag + '_native_review.pdf'))
    print(target)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    run(parser.parse_args().tag)
