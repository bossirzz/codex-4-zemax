from pathlib import Path
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from match_patent_glasses import parse_agf,enrich,CATALOG_DIR
glasses=[]
for g in parse_agf(CATALOG_DIR/'CDGM.AGF')[0]:
 if g['status_code']==1 and 1.78<=g['nd_nominal']<=1.86 and g['vd_nominal']>=40:
  gg=enrich(g)
  if gg['covers_486_to_850nm']:glasses.append(dict(name=gg['name'],nd=gg['nd'],vd=gg['vd']))
print(json.dumps(sorted(glasses,key=lambda g:abs(g['nd']-1.834))[:15],indent=2))
