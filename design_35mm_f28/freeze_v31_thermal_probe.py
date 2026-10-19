"""Freeze a consistent numerical checkpoint for independent native review."""
import argparse
import json
import time
from pathlib import Path

R=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output-tag',default='candidate_v31_thermal_probe')
args=parser.parse_args()
source='candidate_v31_thermal_joint'
suffixes=['_schema.json','_validation_configs.json','_focus_targets.json',
          '_cold_reference_targets.json','_prediction.json','_delta_parameters.npy']
paths=[R/(source+s) for s in suffixes]
for attempt in range(20):
    before=[p.stat().st_mtime_ns for p in paths]
    raw=[p.read_bytes() for p in paths]
    after=[p.stat().st_mtime_ns for p in paths]
    if before==after and after[-1]>=max(after[:-1]):
        for s,data in zip(suffixes,raw):
            (R/(args.output_tag+s)).write_bytes(data)
        prediction=json.loads(raw[4])
        prediction['frozen_from']=source
        prediction['native_validation_required']=True
        (R/(args.output_tag+'_prediction.json')).write_text(
            json.dumps(prediction,ensure_ascii=False,indent=2),encoding='utf8')
        print(json.dumps(dict(tag=args.output_tag,status=prediction['status'],
            min_constraint=prediction['minimum_constraint'],cases=[dict(
                wd=q['wd'],temperature_c=q['temperature_c'],mtf=q.get('predicted_mtf'))
                for q in prediction['cases']])))
        break
    time.sleep(.03)
else:
    raise RuntimeError('Checkpoint changed during snapshot; no files were copied')
