import json
from pathlib import Path
import numpy as np
import optical_engine as e
from optimize_best_conjugate import load_glasses
R=Path(__file__).resolve().parent
s=json.loads((R/'candidate_v19_schema.json').read_text(encoding='utf8'))
load_glasses([(g['catalog'],g['name']) for g in s['glasses']])
rows=json.loads((R/'candidate_v19_native_indices.json').read_text(encoding='utf8'))
ns=e.indices(np.array([.4861327,.5875618,.6562725,.65,.85]));differences=[]
for row in rows:
 expected=ns[:,row['surface']-1];row['independent_indices']=expected.tolist()
 row['max_abs_difference']=float(np.max(abs(expected-np.array(row['indices']))));differences.append(row['max_abs_difference'])
result=dict(maximum_index_difference=max(differences),rows=rows)
(R/'candidate_v19_glass_index_crosscheck.json').write_text(json.dumps(result,indent=2),encoding='utf8')
print('MAXIMUM_INDEX_DIFFERENCE',max(differences))
