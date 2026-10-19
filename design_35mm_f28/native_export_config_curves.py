"""Export native curves from a verified configuration of a separate system copy."""
import os
import argparse
import json
import sys
from pathlib import Path

import zospy as zp

R = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')


def run(tag, config, expected_wd):
    z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    o = None
    try:
        p = z.connect('standalone', message_logging=False)
        o = p.copy_system()
        o.load(R / (tag + '_review.zos'), saveifneeded=False)
        o.MCE.SetCurrentConfiguration(config)
        wd = float(o.LDE.GetSurfaceAt(0).Thickness)
        if abs(wd - expected_wd) > 1e-9:
            raise ValueError(f'Configuration has WD {wd}, expected {expected_wd}')
        a = zp.analyses.mtf.FFTMTF(sampling='1024x1024', maximum_frequency=300,
                                 field='All', wavelength='All')
        result = a.run(o, oncomplete='Close')
        prefix = tag + f'_wd{wd:g}_mtf_1024'
        result.data.to_csv(R / (prefix + '.csv'))
        (R / (prefix + '.json')).write_text(result.to_json(), encoding='utf8')
        metadata = dict(candidate=tag, configuration=config, wd=wd, sampling=1024,
                        messages=[str(x) for x in result.messages or []])
        (R / (prefix + '_metadata.json')).write_text(json.dumps(metadata, indent=2), encoding='utf8')
        print('CONFIG_CURVE_EXPORTED', json.dumps(metadata), flush=True)
    finally:
        if o is not None:
            o.close(saveifneeded=False)
        if z.Application is not None:
            z.disconnect()


if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('tag')
    a.add_argument('--configuration', type=int, required=True)
    a.add_argument('--wd', type=float, required=True)
    q = a.parse_args()
    run(q.tag, q.configuration, q.wd)
