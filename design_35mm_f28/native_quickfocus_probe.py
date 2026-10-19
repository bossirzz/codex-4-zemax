"""Check native Quick Focus at weak conjugates on a disposable saved copy."""
import os
import json
from pathlib import Path
import sys

import zospy as zp

R = Path(__file__).resolve().parent
SOURCE = 'candidate_v32_native_contrast_c'
TAG = 'candidate_v32_native_quickfocus_probe'
sys.stdout.reconfigure(encoding='utf8')
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = tool = None
rows = []
reference = json.loads((R / (SOURCE + '_native_post_audit.json')).read_text())
poses = json.loads((R / (SOURCE + '_validation_configs.json')).read_text())
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / (SOURCE + '.zos'), saveifneeded=False)
    typ = zp.constants.Editors.MFE.MeritOperandType
    controller = next(o.MCE.GetOperandAt(i) for i in range(1, o.MCE.NumberOfOperands+1)
                      if str(o.MCE.GetOperandAt(i).Type) == 'THIC' and int(o.MCE.GetOperandAt(i).Param1) == 15)

    def mtf(sampling=4):
        return [float(o.MFE.GetOperandValue(getattr(typ, kind), sampling, 0,
                     field, 225. if field == 1 else 180., 1, 0, 0, 0))
                for field in range(1, 5) for kind in ['MTFT', 'MTFS']]

    def score(values):
        return min(min(values[:2])-.3, min(values[2:])-.2)

    for config, wd in [(5, 300), (6, 400), (7, 2000)]:
        o.MCE.SetCurrentConfiguration(config)
        cell = controller.GetOperandCell(config)
        original = float(cell.DoubleValue)
        before = mtf()
        tool = o.Tools.OpenQuickFocus()
        tool.Criterion = zp.constants.Tools.General.QuickFocusCriterion.SpotSizeRadial
        tool.UseCentroid = True
        tool.RunAndWaitForCompletion()
        succeeded, error = bool(tool.Succeeded), str(tool.ErrorMessage)
        tool.Close()
        tool = None
        focused = float(o.LDE.GetSurfaceAt(15).Thickness)
        cell.DoubleValue = focused
        o.MCE.SetCurrentConfiguration(1)
        o.MCE.SetCurrentConfiguration(config)
        after = mtf()
        selected = succeeded and score(after) > score(before)
        after_high = mtf(6) if selected else None
        before_high = next(q['mtf_sampling6'] for q in reference if q['wd'] == wd)
        if selected:
            selected = score(after_high) > score(before_high)
        cell.DoubleValue = focused if selected else original
        poses[config-1]['backfocus'] = focused if selected else original
        rows.append(dict(wd=wd, configuration=config, tool='OpticStudio QuickFocus SpotSizeRadial Centroid',
            succeeded=succeeded, error=error, source_backfocus_mm=original, probed_backfocus_mm=focused,
            native_mtf_sampling4_before=before, native_mtf_sampling4_after=after,
            native_mtf_sampling6_before=before_high, native_mtf_sampling6_after=after_high,
            selected_in_disposable_probe=selected,
            high_sampling_acceptance_required=True, complete_requirements_passed=False))
        (R / (TAG + '.json')).write_text(json.dumps(rows, indent=2), encoding='utf8')
        print('QUICKFOCUS_PROBE', json.dumps(rows[-1]), flush=True)
    o.MCE.SetCurrentConfiguration(1)
    o.save_as(R / (TAG + '.zos'))
    (R / (TAG + '_schema.json')).write_text((R / (SOURCE + '_schema.json')).read_text(), encoding='utf8')
    (R / (TAG + '_validation_configs.json')).write_text(json.dumps(poses, indent=2), encoding='utf8')
finally:
    if tool is not None:
        tool.Close()
    if o is not None:
        o.close(saveifneeded=False)
    if z.Application is not None:
        z.disconnect()
