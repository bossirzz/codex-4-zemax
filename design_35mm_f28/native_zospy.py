"""Build a complete declared design in an isolated OpticStudio system.

Does not edit the user's existing lens; keeps an independent session backup.
"""
import os
import argparse
from datetime import datetime
import json
import math
from pathlib import Path
import sys
import numpy as np
import zospy as zp

ROOT=Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf-8')

def snap(oss):
    return dict(file=str(oss.SystemFile),needs_save=bool(oss.NeedsSave),
        rows=[[str(getattr(oss.LDE.GetSurfaceAt(i),k)) for k in ['Radius','Thickness','Material','Comment']]
              for i in range(oss.LDE.NumberOfSurfaces)])

def run(schema_name,sampling=512,connection='extension'):
    schema=json.loads((ROOT/schema_name).read_text(encoding='utf8'))
    tag=Path(schema_name).stem.replace('_schema','')
    target=ROOT/(tag+'.zos')
    z=zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'));isolated=None
    try:
        primary=z.connect(connection,**({'instance_number':1} if connection=='extension' else {}),message_logging=False)
        before=snap(primary)
        isolated=primary.copy_system()
        if isolated.SystemID==primary.SystemID:raise RuntimeError('Independent system required')
        backup=ROOT/f'user_session_backup_{datetime.now():%Y%m%d_%H%M%S}.zos'
        isolated.save_as(backup)
        isolated.new(saveifneeded=False)
        isolated.SystemName='35 mm F2.8 | optimized actual glass | WD500 | candidate'
        cats=isolated.SystemData.MaterialCatalogs
        for cat in sorted({g['catalog'] for g in schema['glasses']}):
            cats.AddCatalog(cat)
        while isolated.LDE.NumberOfSurfaces<len(schema['curvatures'])+2:
            isolated.LDE.InsertNewSurfaceAt(isolated.LDE.NumberOfSurfaces-1)
        isolated.LDE.GetSurfaceAt(0).Thickness=schema['wd']
        for i,(cv,t,m) in enumerate(zip(schema['curvatures'],schema['thicknesses'],schema['media']),1):
            row=isolated.LDE.GetSurfaceAt(i)
            row.Radius=1/cv if cv else math.inf
            row.Thickness=t
            row.Material='' if m<0 else schema['glasses'][m]['name']
            row.Comment='STOP' if i==schema['stop_index']+1 else f'S{i:02d}'
            if i==schema['stop_index']+1:row.IsStop=True
        for asp in schema.get('aspheres',[]):
            row=isolated.LDE.GetSurfaceAt(asp['surface'])
            row.ChangeType(row.GetSurfaceTypeSettings(zp.constants.Editors.LDE.SurfaceType.EvenAspheric))
            for number,value in zip([2,3,4],asp['coefficients']):
                row.GetSurfaceCell(getattr(zp.constants.Editors.LDE.SurfaceColumn,f'Par{number}')).DoubleValue=value
        isolated.LDE.GetSurfaceAt(isolated.LDE.NumberOfSurfaces-1).Comment='IMAGE Phi18 mm'
        isolated.LDE.GetSurfaceAt(isolated.LDE.NumberOfSurfaces-1).Thickness=0.
        ap=isolated.SystemData.Aperture
        ap.ApertureType=zp.constants.SystemData.ZemaxApertureType.ImageSpaceFNum
        ap.ApertureValue=2.8
        waves=isolated.SystemData.Wavelengths
        while waves.NumberOfWavelengths>1:waves.RemoveWavelength(waves.NumberOfWavelengths)
        for i,(w,weight) in enumerate([(.4861327,1),(.5875618,1),(.6562725,1),(.65,0),(.85,0)],1):
            wave=waves.GetWavelength(1) if i==1 else waves.AddWavelength(w,weight)
            wave.Wavelength=w;wave.Weight=weight
        waves.GetWavelength(2).MakePrimary()
        fields=isolated.SystemData.Fields
        fields.SetFieldType(zp.constants.SystemData.FieldType.ParaxialImageHeight)
        while fields.NumberOfFields>1:fields.DeleteFieldAt(fields.NumberOfFields)
        for i,y in enumerate(schema['fields'],1):
            if i==1:fields.GetField(1).Y=y
            else:fields.AddField(0,y,1)
        if 'ray_aiming' in schema:
            import System
            aiming=isolated.SystemData.RayAiming
            aiming.RayAiming=System.Enum.Parse(aiming.RayAiming.GetType(),schema['ray_aiming'])
        for surface,cap in schema.get('fixed_semi_diameters',{}).items():
            row=isolated.LDE.GetSurfaceAt(int(surface))
            row.SemiDiameterCell.SetSolveData(row.SemiDiameterCell.CreateSolveType(zp.constants.Editors.SolveType.Fixed))
            row.SemiDiameter=float(cap)
        if 'fixed_stop_radius' in schema:
            ap.ApertureType=zp.constants.SystemData.ZemaxApertureType.FloatByStopSize
            row=isolated.LDE.GetSurfaceAt(schema['stop_index']+1)
            row.SemiDiameterCell.SetSolveData(row.SemiDiameterCell.CreateSolveType(zp.constants.Editors.SolveType.Fixed))
            row.SemiDiameter=float(schema['fixed_stop_radius'])
        isolated.save_as(target)
        isolated.load(target,saveifneeded=False)
        rows=[]
        for i in range(isolated.LDE.NumberOfSurfaces):
            row=isolated.LDE.GetSurfaceAt(i)
            rows.append(dict(surface=i,radius=str(row.Radius),thickness=float(row.Thickness),material=str(row.Material),semi_diameter=float(row.SemiDiameter),stop=bool(row.IsStop)))
        print(json.dumps(dict(stage='built',target=str(target),backup=str(backup),rows=rows),ensure_ascii=False),flush=True)
        analysis=zp.analyses.mtf.FFTMTF(sampling=f'{sampling}x{sampling}',maximum_frequency=300,field='All',wavelength='All')
        mtf=analysis.run(isolated,oncomplete='Close')
        ROOT.joinpath(tag+'_native_mtf.json').write_text(mtf.to_json(),encoding='utf8')
        if mtf.data is not None:
            mtf.data.to_csv(ROOT/(tag+'_native_mtf.csv'))
            print(mtf.data.loc[mtf.data.index.to_series().sub(225).abs().nsmallest(1).index].to_string(),flush=True)
            print(mtf.data.loc[mtf.data.index.to_series().sub(180).abs().nsmallest(1).index].to_string(),flush=True)
        print('MESSAGES',str(mtf.messages),flush=True)
        system=zp.analyses.new_analysis(isolated,zp.constants.Analysis.AnalysisIDM.SystemData)
        try:
            system.ApplyAndWaitForCompletion()
            report=system.get_text_output(str(ROOT/(tag+'_native_system.txt')),z.get_txtfile_encoding())
            print('\n'.join(s for s in report.splitlines() if any(k in s for k in ['有效焦距','后焦距','总长','F/#','光阑半径','Effective Focal','Total Track','Back Focal'])),flush=True)
        finally:system.Close()
        if snap(primary)!=before:raise RuntimeError('Original system changed unexpectedly')
        ROOT.joinpath(tag+'_native_build.json').write_text(json.dumps(dict(status='built_and_analyzed',sampling=sampling,target=str(target),backup=str(backup),rows=rows,original_unchanged=True,messages=[str(m) for m in mtf.messages or []]),ensure_ascii=False,indent=2),encoding='utf8')
    finally:
        if isolated is not None:isolated.close(saveifneeded=False)
        if z.Application is not None:z.disconnect()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('schema');ap.add_argument('--sampling',type=int,default=512)
    ap.add_argument('--connection',choices=['extension','standalone'],default='extension')
    a=ap.parse_args();run(a.schema,a.sampling,a.connection)
