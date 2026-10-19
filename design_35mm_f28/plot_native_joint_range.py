"""Native FFT curves and exact field operands for the engineering range checks."""
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
    reference = json.loads((R / 'candidate_v32_native_focus_refine_review_evidence.json').read_text())['cases']
    font = FontProperties(fname=r'C:\Windows\Fonts\msyh.ttc')
    colors = ['#202124', '#3079c1', '#13866c', '#c43d42']
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    for row, wd in enumerate([300, 400, 2000]):
        new = next(q for q in current if q['wd'] == wd)
        old = next(q for q in reference if q['wd'] == wd)
        df = pd.read_csv(R / (tag + f'_WD{wd}_mtf_1024.csv'), header=[0, 1], index_col=0)
        for field, height in enumerate([0., 3., 6.3, 9.08]):
            for direction, style in enumerate(['-', '--']):
                axes[row, 0].plot(np.asarray(df.index, float), df.iloc[:, 2*field+direction],
                    style, c=colors[field], lw=1.5, label=f'{height:g} mm '+['T', 'S'][direction])
        axes[row, 0].set(xlim=(0, 300), ylim=(0, 1))
        axes[row, 0].set_title(f'物距 {wd} mm · 原生1024×1024曲线 · F/{new["nominal_fno"]:.4f}', fontproperties=font)
        axes[row, 0].set_xlabel('空间频率（lp/mm）', fontproperties=font)
        axes[row, 0].legend(ncol=2, prop=font, loc='upper right', fontsize=9)
        x = [0., 3., 6.3, 9.08]
        for k, name in enumerate(['T', 'S']):
            axes[row, 1].plot(x, old['mtf_sampling6'][k::2], '--', c='#999999', lw=1.4,
                marker=['o', 's'][k], label='优化前 '+name)
            axes[row, 1].plot(x, new['mtf_sampling6'][k::2], '-', c=['#287caa', '#c43d42'][k],
                lw=2, marker=['o', 's'][k], label='本轮 '+name)
        axes[row, 1].scatter(x, [.3, .2, .2, .2], c='#222', marker='x', s=70, label='扩展工程对照')
        axes[row, 1].set(xlim=(-.4, 9.5), ylim=(0, .46), xticks=x)
        axes[row, 1].set_title(f'中心225 / 离轴180 lp/mm · 最差离轴 {min(new["mtf_sampling6"][2:]):.4f}', fontproperties=font)
        axes[row, 1].set_xlabel('近轴像高（mm）', fontproperties=font)
        axes[row, 1].legend(prop=font, loc='upper right', fontsize=9)
        for ax in axes[row]:
            ax.set_ylabel('光学 MTF', fontproperties=font)
            ax.grid(alpha=.2)
    title = 'Zemax原生四组调焦试验' if 'fourgroup' in tag else 'Zemax原生联合优化'
    fig.suptitle(title + ' · 主工作段与远距 · F/d/C等权 · 裸像面', fontproperties=font, fontsize=16)
    fig.text(.5, .015, '本图为扩展工程检查；500/200 mm主验收另行核对。优化前数据来自已保存的v32焦位修正版。',
        ha='center', fontproperties=font, fontsize=10)
    fig.tight_layout(rect=[0, .035, 1, .96])
    path = R / (tag + '_range_review.png')
    fig.savefig(path, dpi=220)
    fig.savefig(R / (tag + '_range_review.pdf'))
    print(path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag')
    run(parser.parse_args().tag)
