"""Import CN110703413B example 3 as an explicitly provisional ZOSPy model."""
import os
import argparse
import csv
import json
import math
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'zemax_import' / 'CN110703413B_example3'
OUT.mkdir(parents=True, exist_ok=True)
TARGET = OUT / 'CN110703413B_example3_modelglass.zos'

# Original Table 7: face label, radius, following thickness, Nd, Vd.
SOURCE = [
    ('L1R1', 23.35, 3.133, 1.77, 49.1),
    ('L1R2', 59.04, .14, None, None),
    ('L2R1', 10.81, 4.03, 1.59, 68.5),
    ('L2R2', 44.65, .23, None, None),
    ('L3R1', 35.96, .79, 1.57, 41.5),
    ('L3R2', 7.83, 8.76, None, None),
    ('L4R1', -8.35, 1.45, 1.72, 28.3),
    ('L4R2', -21.31, 0, None, None),
    ('L5R1', -21.31, 1.97, 1.80, 46.7),
    ('L5R2', -13.35, .14, None, None),
    ('L6R1', -1159, 2.13, 1.77, 49.1),
    ('L6R2', -15.21, 3.13, None, None),
    ('L7R1', 114.24, 1.11, 1.76, 26),
    ('L7R2', 24.8, 3.84, None, None),
    ('L8R1', 32.62, 4.065, 1.77, 49.1),
    ('L8R2', 242.95, 14.65, None, None),
]
CONFIGS = [
    {'name': 'WD500_source', 'object_mm': 500, 'gap_mm': 3.13},
    {'name': 'WD150_source', 'object_mm': 150, 'gap_mm': 7.71},
    {'name': 'WD300_source', 'object_mm': 300, 'gap_mm': 4.88},
    {'name': 'WD1000_source', 'object_mm': 1000, 'gap_mm': 1.88},
    {'name': 'Infinity_source', 'object_mm': 'Infinity', 'gap_mm': .79},
]
WAVES = [(0.4861327, 1), (0.5875618, 1), (0.6562725, 1), (0.650, 0), (0.850, 0)]
FIELDS = [0, 4.5, 6.3, 9.0]


def emit(stage, **data):
    print(json.dumps({'stage': stage, **data}, ensure_ascii=True, default=str), flush=True)


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str), encoding='utf-8')


def compile_prescription():
    """Resolve model glasses and compile shared cemented interfaces once."""
    rows = []
    for label, radius, thickness, nd, vd in SOURCE:
        if label == 'L4R2':
            continue
        if label == 'L5R1':
            label = 'L4R2_L5R1_cemented'
        if label == 'L3R2':
            thickness = 4.38
        rows.append({'comment': label, 'radius': radius, 'thickness': thickness,
                     'glass': None if nd is None else {'kind': 'model', 'nd': nd, 'vd': vd, 'dpgf': 0}})
        if label == 'L3R2':
            rows.append({'comment': 'STOP_assumed_midgap', 'radius': 'Infinity',
                         'thickness': 4.38, 'glass': None, 'stop': True})
    assert len(rows) == 16
    assert math.isclose(sum(row['thickness'] for row in rows), 49.568)
    assert sum(row['glass'] is not None for row in rows) == 8
    assert rows[11]['comment'] == 'L6R2'
    assert rows[8]['comment'] == 'L4R2_L5R1_cemented'
    return rows


def prepare():
    rows = compile_prescription()
    schema = {
        'patent': 'CN110703413B', 'example': 3,
        'source': 'https://patents.google.com/patent/CN110703413B/zh',
        'source_pdf': str(ROOT / 'patent_research_35mm' / 'CN110703413B.pdf'),
        'units': 'mm', 'sequential': True, 'source_table': 'Table 7 (PDF page 8)',
        'surfaces': rows, 'configurations': CONFIGS,
        'aperture': {'type': 'ImageSpaceFNum', 'value': 2.8},
        'fields': {'type': 'ParaxialImageHeight', 'y_mm': FIELDS, 'origin': 'user requirement Phi18mm'},
        'wavelengths_um_weights': WAVES, 'primary_wavelength_index': 2,
        'assumptions': [
            'Exact catalog glass names and partial dispersion are not disclosed. Model glasses reproduce published rounded Nd/Vd; dPgF=0 is provisional.',
            'Stop is between L3 and L4; exact axial position is not disclosed. 8.76mm gap is provisionally split 4.38mm/4.38mm.',
            'Object distance is interpreted from the object plane to the first surface; patent WD reference is not specified.',
            'Fields are evaluation settings for user Phi18mm requirement; coverage is not validated.',
            '650nm and 850nm have zero polychromatic weights but remain selectable for monochromatic analyses; NIR glass dispersion is not validated.',
            'Semi-diameters use Zemax automatic sizing because the patent does not supply clear apertures; no 23mm diameter compliance is asserted.',
            'Source spacings and fixed 14.65mm last-surface-to-image distance are retained; no autofocus or optimization is applied.',
            'Rounded Table 7 track is 49.568mm at WD500, whereas Table 8 nominal TTL is 50mm. The prescription is not stretched to force 50mm.',
        ],
    }
    write_json(OUT / 'prescription.json', schema)
    glasses = [{'lens': label[:2], 'resolution': 'ModelGlass', 'nd': nd, 'vd': vd,
                'dpgf': 0, 'catalog_name': None}
               for label, radius, thickness, nd, vd in SOURCE if nd is not None]
    write_json(OUT / 'glass_resolution.json', glasses)
    with (OUT / 'source_table7.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['face', 'radius_mm', 'following_thickness_mm', 'Nd', 'Vd'])
        writer.writerows(SOURCE)
    emit('validated_schema', refracting_elements=8, sequential_rows=len(rows),
         configurations=len(CONFIGS), track_mm=49.568, glass_resolution='published Nd/Vd model glasses')
    return schema


def snapshot(oss):
    return {'file': str(oss.SystemFile), 'needs_save': bool(oss.NeedsSave),
            'surfaces': [[str(getattr(oss.LDE.GetSurfaceAt(i), key))
                          for key in ['Radius', 'Thickness', 'Material', 'Conic', 'Comment']]
                         for i in range(oss.LDE.NumberOfSurfaces)],
            'configuration': int(oss.MCE.CurrentConfiguration)}


def validate_loaded(oss, schema, zp):
    """Compare every imported geometric and glass value after a disk round trip."""
    oss.MCE.SetCurrentConfiguration(1)
    assert int(oss.LDE.NumberOfSurfaces) == 18
    assert int(oss.MCE.NumberOfConfigurations) == len(CONFIGS)
    records = []
    for i, expected in enumerate(schema['surfaces'], 1):
        row = oss.LDE.GetSurfaceAt(i)
        wanted_r = math.inf if expected['radius'] == 'Infinity' else expected['radius']
        assert math.isclose(float(row.Radius), wanted_r, rel_tol=0, abs_tol=1e-9), (i, 'radius')
        assert math.isclose(float(row.Thickness), expected['thickness'], abs_tol=1e-9), (i, 'thickness')
        assert bool(row.IsStop) == bool(expected.get('stop', False)), (i, 'stop')
        if expected['glass'] is not None:
            solve = row.MaterialCell.GetSolveData()._S_MaterialModel
            assert math.isclose(float(solve.IndexNd), expected['glass']['nd'], abs_tol=1e-9)
            assert math.isclose(float(solve.AbbeVd), expected['glass']['vd'], abs_tol=1e-9)
        else:
            assert str(row.Material) == '', (i, 'expected air')
        records.append({'surface': i, 'comment': str(row.Comment), 'radius': str(row.Radius),
                        'thickness_mm': float(row.Thickness), 'material': str(row.Material),
                        'semi_diameter_mm': float(row.SemiDiameter), 'stop': bool(row.IsStop)})
    for i, expected in enumerate(WAVES, 1):
        wave = oss.SystemData.Wavelengths.GetWavelength(i)
        assert math.isclose(float(wave.Wavelength), expected[0], abs_tol=1e-10)
        assert math.isclose(float(wave.Weight), expected[1], abs_tol=1e-10)
    assert int(oss.SystemData.Wavelengths.NumberOfWavelengths) == len(WAVES)
    assert str(oss.SystemData.Aperture.ApertureType) == 'ImageSpaceFNum'
    assert math.isclose(float(oss.SystemData.Aperture.ApertureValue), 2.8)
    assert str(oss.SystemData.Fields.GetFieldType()) == 'ParaxialImageHeight'
    assert [float(oss.SystemData.Fields.GetField(i).Y) for i in range(1, 5)] == FIELDS
    configurations = []
    for i, config in enumerate(CONFIGS, 1):
        oss.MCE.SetCurrentConfiguration(i)
        wanted_obj = math.inf if config['object_mm'] == 'Infinity' else config['object_mm']
        assert math.isclose(float(oss.LDE.GetSurfaceAt(0).Thickness), wanted_obj, abs_tol=1e-8)
        assert math.isclose(float(oss.LDE.GetSurfaceAt(12).Thickness), config['gap_mm'], abs_tol=1e-8)
        assert math.isclose(float(oss.LDE.GetSurfaceAt(16).Thickness), 14.65, abs_tol=1e-8)
        configurations.append({'number': i, **config,
                               'track_mm': sum(float(oss.LDE.GetSurfaceAt(k).Thickness) for k in range(1, 17))})
    oss.MCE.SetCurrentConfiguration(1)
    return {'status': 'passed', 'surface_count': 18, 'glass_count': 8,
            'rows': records, 'configurations': configurations,
            'primary_configuration': 'WD500_source', 'disk_round_trip': True}


def build(schema):
    import zospy as zp
    import zospy_zh as zh
    if TARGET.exists():
        raise FileExistsError(f'Refusing to overwrite existing model: {TARGET}')
    zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
    isolated = None
    report = {'status': 'running', 'target': str(TARGET), 'assumptions': schema['assumptions']}
    try:
        primary = zos.connect('extension', instance_number=0, message_logging=False)
        before = snapshot(primary)
        write_json(OUT / 'previous_session_snapshot.json', before)
        emit('connected', file=before['file'], needs_save=before['needs_save'])
        isolated = primary.copy_system()
        if int(isolated.SystemID) == int(primary.SystemID):
            raise RuntimeError('CopySystem did not create an independent optical system')
        backup = OUT / f'previous_session_backup_{datetime.now():%Y%m%d_%H%M%S}.zos'
        isolated.save_as(backup)
        assert backup.exists() and backup.stat().st_size > 0
        report['previous_session_backup'] = str(backup)
        emit('backed_up', backup=str(backup))
        isolated.new(saveifneeded=False)
        if str(isolated.Mode) != 'Sequential':
            isolated.MakeSequential()
        isolated.SystemName = 'CN110703413B Ex3 | provisional ModelGlass | WD500'
        # New sequential systems contain OBJ, STOP and IMA. Extend to compiled row count.
        while int(isolated.LDE.NumberOfSurfaces) < len(schema['surfaces']) + 2:
            isolated.LDE.InsertNewSurfaceAt(int(isolated.LDE.NumberOfSurfaces) - 1)
        isolated.LDE.GetSurfaceAt(0).Thickness = 500
        for i, data in enumerate(schema['surfaces'], 1):
            row = isolated.LDE.GetSurfaceAt(i)
            row.Comment = data['comment']
            row.Radius = math.inf if data['radius'] == 'Infinity' else data['radius']
            row.Thickness = data['thickness']
            glass = data['glass']
            if glass:
                zp.solvers.material_model(row.MaterialCell, refractive_index=glass['nd'],
                                          abbe_number=glass['vd'], partial_dispersion=glass['dpgf'])
            else:
                row.Material = ''
            if data.get('stop'):
                row.IsStop = True
        isolated.LDE.GetSurfaceAt(17).Comment = 'IMAGE_requirement_Phi18mm_unvalidated'
        isolated.SystemData.Aperture.ApertureType = zp.constants.SystemData.ZemaxApertureType.ImageSpaceFNum
        isolated.SystemData.Aperture.ApertureValue = 2.8
        waves = isolated.SystemData.Wavelengths
        while int(waves.NumberOfWavelengths) > 1:
            waves.RemoveWavelength(int(waves.NumberOfWavelengths))
        waves.GetWavelength(1).Wavelength = WAVES[0][0]
        waves.GetWavelength(1).Weight = WAVES[0][1]
        for wavelength, weight in WAVES[1:]:
            waves.AddWavelength(wavelength, weight)
        waves.GetWavelength(2).MakePrimary()
        fields = isolated.SystemData.Fields
        fields.SetFieldType(zp.constants.SystemData.FieldType.ParaxialImageHeight)
        while int(fields.NumberOfFields) > 1:
            fields.DeleteFieldAt(int(fields.NumberOfFields))
        fields.GetField(1).X = 0
        fields.GetField(1).Y = FIELDS[0]
        for y in FIELDS[1:]:
            fields.AddField(0, y, 1)
        mce = isolated.MCE
        for _ in CONFIGS[1:]:
            mce.AddConfiguration(False)
        # Source configurations vary object distance and only the L6-L7 air gap.
        for operand_number, surf, values in [
            (1, 0, [math.inf if c['object_mm'] == 'Infinity' else c['object_mm'] for c in CONFIGS]),
            (2, 12, [c['gap_mm'] for c in CONFIGS]),
        ]:
            if int(mce.NumberOfOperands) < operand_number:
                mce.AddOperand()
            operand = mce.GetOperandAt(operand_number)
            operand.ChangeType(zp.constants.Editors.MCE.MultiConfigOperandType.THIC)
            operand.Param1 = surf
            for config_index, value in enumerate(values, 1):
                operand.GetOperandCell(config_index).DoubleValue = value
        mce.SetCurrentConfiguration(1)
        emit('built', surfaces=int(isolated.LDE.NumberOfSurfaces), configurations=int(mce.NumberOfConfigurations))
        report.update(validate_loaded(isolated, schema, zp))
        # Require a basic real and paraxial ray trace, without changing focus or prescription.
        ray = zh.SingleRayTrace(hx=0, hy=0, px=0, py=.5, wavelength=2)
        try:
            ray_result = ray.run(isolated, oncomplete='Close')
            (OUT / 'ray_trace.json').write_text(ray_result.to_json(), encoding='utf-8')
            report['ray_trace'] = {'status': 'completed', 'messages': [str(m) for m in ray_result.messages or []]}
        finally:
            if ray.analysis is not None:
                ray.analysis.Close()
        system = zp.analyses.new_analysis(isolated, zp.constants.Analysis.AnalysisIDM.SystemData)
        try:
            system.ApplyAndWaitForCompletion()
            text = system.get_text_output(str(OUT / 'system_data_native.txt'), zos.get_txtfile_encoding())
            if not text.strip():
                raise RuntimeError('Native SystemData returned no report')
            report['system_data_metrics'] = [line.strip() for line in text.splitlines()
                                             if any(label in line for label in [
                                                 '有效焦距', '后焦距', '总长', 'F/#', '光阑半径',
                                                 'Effective Focal Length', 'Back Focal Length', 'Total Track'])]
        finally:
            system.Close()
        isolated.save_as(TARGET)
        assert TARGET.exists() and TARGET.stat().st_size > 0
        isolated.load(TARGET, saveifneeded=False)
        report.update(validate_loaded(isolated, schema, zp))
        isolated.save()
        if snapshot(primary) != before:
            raise RuntimeError('The existing primary lens changed during isolated import; saved model is ready but will not replace the current session.')
        # User requested import into Zemax. The original unsaved system is already backed up.
        primary.load(TARGET, saveifneeded=False)
        report.update(validate_loaded(primary, schema, zp))
        primary.MCE.SetCurrentConfiguration(1)
        primary.save()
        report['loaded_in_open_zemax'] = str(primary.SystemFile)
        report['status'] = 'passed'
        write_json(OUT / 'import_report.json', report)
        emit('completed', target=str(TARGET), backup=str(backup), readback='passed', loaded_in_open_zemax=True)
    except Exception as exc:
        report['status'] = 'failed'
        report['error'] = {'type': type(exc).__name__, 'message': str(exc)}
        write_json(OUT / 'import_report.json', report)
        emit('failed', **report['error'])
        raise
    finally:
        if isolated is not None:
            isolated.close(saveifneeded=False)
        if zos.Application is not None:
            zos.disconnect()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', action='store_true')
    args = parser.parse_args()
    schema = prepare()
    if args.build:
        build(schema)
