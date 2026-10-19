import os
import json,sys
from pathlib import Path
import zospy as zp
R=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));o=None
try:
 p=z.connect('standalone',message_logging=False);o=p.copy_system();o.load(R/'candidate_v19.zos',saveifneeded=False)
 typ=zp.constants.Editors.MFE.MeritOperandType;rows=[]
 for i in [1,3,4,7,8,10,12,14]:
  rows.append(dict(surface=i,material=str(o.LDE.GetSurfaceAt(i).Material),indices=[float(o.MFE.GetOperandValue(typ.INDX,i,w,0,0,0,0,0,0)) for w in range(1,6)]))
 (R/'candidate_v19_native_indices.json').write_text(json.dumps(rows,indent=2),encoding='utf8');print('INDICES',json.dumps(rows),flush=True)
 cfg=next(x for x in json.loads((R/'candidate_v19_real_focus_targets.json').read_text(encoding='utf8')) if x['wd']==200)
 o.LDE.GetSurfaceAt(0).Thickness=200.;o.LDE.GetSurfaceAt(11).Thickness=cfg['gap'];o.LDE.GetSurfaceAt(15).Thickness=cfg['backfocus']
 for sample in [512,1024]:
  a=zp.analyses.mtf.FFTMTF(sampling=f'{sample}x{sample}',maximum_frequency=200,field='All',wavelength='All')
  rr=a.run(o,oncomplete='Close');rr.data.to_csv(R/f'candidate_v19_near_native_{sample}.csv');print('NEAR_FFT',sample,rr.data.iloc[(abs(rr.data.index-120)).argmin()].tolist(),flush=True)
finally:
 if o is not None:o.close(saveifneeded=False)
 if z.Application is not None:z.disconnect()
