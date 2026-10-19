import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

OUT = Path(__file__).with_name('zospy_image_results')
path = OUT / 'report.json'
report = json.loads(path.read_text(encoding='utf-8'))
checks = []
for name in ('cross_section.png', 'cross_section_960x640.png', 'viewer3d.png'):
    image_path = OUT / name
    if not image_path.exists():
        continue
    with Image.open(image_path) as image:
        image.load()
        pixels = np.asarray(image.convert('RGB'))
        nonwhite = int(np.any(pixels < 245, axis=-1).sum())
        checks.append({'file': name, 'width': image.width, 'height': image.height,
                       'nonwhite_pixels': nonwhite, 'decoded': True,
                       'dimensions_match_960x640': image.size == (960, 640)})
report['image_checks'] = checks
report['image_size_convention'] = {'CrossSection': '(height, width)', 'Viewer3D': '(width, height)'}

curves = {}
for p in OUT.glob('mtf_*.json'):
    size = int(p.stem.split('_')[1])
    data = json.loads(p.read_text(encoding='utf-8'))['data']
    curves[size] = (np.array(data['index'], dtype=float), np.array(data['data'], dtype=float))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
colors = plt.get_cmap('viridis')(np.linspace(0.05, 0.9, len(curves)))
for (size, (frequency, values)), color in zip(sorted(curves.items()), colors):
    for ax in axes:
        ax.plot(frequency, values[:, 2], label=f'{size} x {size}', color=color, linewidth=1.6)
for ax in axes:
    ax.set_xlabel('Spatial frequency (cycles/mm)')
    ax.set_ylabel('MTF (on-axis tangential)')
    ax.set_ylim(-0.02, 1.02)
    ax.grid(alpha=0.25)
axes[0].set_xlim(0, 100)
axes[0].set_title('Full frequency range')
axes[1].set_xlim(0, 5)
axes[1].set_title('Low-frequency detail')
axes[1].legend(fontsize=8)
fig.suptitle('Zemax FFT MTF: sampling comparison', fontsize=14)
fig.savefig(OUT / 'mtf_sampling_comparison.png', dpi=160)
plt.close(fig)

mtf_tests = [t for t in report['tests'] if t['name'].startswith('fft_mtf_') and t['status'] == 'passed']
if mtf_tests:
    latest = max(mtf_tests, key=lambda t: t['details']['sampling'])['details']
    report['numerical_accuracy_check'] = {
        'highest_sampling': latest['sampling'], 'messages': latest['messages'],
        'compared_with_sampling': latest.get('compared_with_sampling'),
        'max_absolute_difference': latest.get('max_absolute_difference'),
        'mean_absolute_difference': latest.get('mean_absolute_difference'),
        'convergence_established': False if latest['messages'] else None,
        'note': 'Successful API execution and valid image export do not establish optical numerical accuracy.'}
path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'image_checks': checks, 'accuracy': report.get('numerical_accuracy_check'),
                  'lens_parameters_unchanged': report.get('lens_parameters_unchanged'),
                  'plot': str(OUT / 'mtf_sampling_comparison.png')}, ensure_ascii=True, indent=2))
