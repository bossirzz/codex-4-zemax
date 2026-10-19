"""Diagnose native operand ranges on an isolated prepared system."""
import os
import json
import sys
from pathlib import Path
import zospy as zp

R = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding='utf8')
z = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
o = None
try:
    p = z.connect('standalone', message_logging=False)
    o = p.copy_system()
    o.load(R / 'candidate_v32_native_contrast_prepared.zos', saveifneeded=False)
    mf = o.MFE
    out = dict(mfe_members=[n for n in dir(mf) if 'Error' in n or 'Valid' in n], rows=[])
    for number in [1, 4, 7, 8, 9, 10, 11, 12, 860, 861, 862]:
        row = mf.GetOperandAt(number)
        cells = []
        for i in range(2, 10):
            cell = row.GetCellAt(i)
            dtype = str(cell.DataType)
            value = cell.IntegerValue if dtype == 'Integer' else cell.DoubleValue if dtype == 'Double' else cell.StringValue
            cells.append(dict(column=i, type=dtype, value=str(value)))
        out['rows'].append(dict(index=number, type=str(row.Type), cells=cells,
                               error_members=[n for n in dir(row) if 'Error' in n or 'Valid' in n]))
    # Restrict the automatically generated physical boundaries to actual lens surfaces.
    for i in range(1, mf.NumberOfOperands + 1):
        row = mf.GetOperandAt(i)
        if str(row.Type) in ['MNCA', 'MXCA', 'MNEA', 'MNCG', 'MXCG', 'MNEG']:
            row.GetCellAt(2).IntegerValue = 1
            row.GetCellAt(3).IntegerValue = 15
    out['corrected_merit'] = float(mf.CalculateMeritFunction())
    out['error_count'] = str(getattr(mf, 'NumberOfErrors', 'unavailable'))
    out['corrected_bounds'] = [dict(index=i, type=str(mf.GetOperandAt(i).Type),
        value=str(mf.GetOperandAt(i).Value)) for i in [7, 8, 9, 10, 11, 12]]
    out['invalid_rows'] = [dict(index=i,type=str(mf.GetOperandAt(i).Type))
        for i in range(1,mf.NumberOfOperands+1) if not mf.GetOperandAt(i).IsValidRow]
    # Isolate the generated objective from explicit native constraints.
    while mf.NumberOfOperands >= 860:
        mf.RemoveOperandAt(mf.NumberOfOperands)
    out['generated_only_merit'] = float(mf.CalculateMeritFunction())
    audit = json.loads((R / 'candidate_v32_native_contrast_merit_audit.json').read_text())
    out['explicit_prefix_tests'] = []
    for spec in audit['explicit_operands']:
        row = mf.AddOperand()
        row.ChangeType(getattr(zp.constants.Editors.MFE.MeritOperandType, spec['type']))
        for j, value in enumerate(spec['parameters'], 2):
            cell = row.GetCellAt(j)
            if str(cell.DataType) == 'Integer':
                cell.IntegerValue = int(value)
            else:
                cell.DoubleValue = float(value)
        row.Target, row.Weight = spec['target'], spec['weight']
        merit = float(mf.CalculateMeritFunction())
        result = dict(index=mf.NumberOfOperands, type=spec['type'], merit=merit, value=str(row.Value))
        out['explicit_prefix_tests'].append(result)
        if merit >= 1e8 or merit != merit:
            out['first_explicit_failure'] = result
            break
    while mf.NumberOfOperands >= 860:
        mf.RemoveOperandAt(mf.NumberOfOperands)
    # Keep only the first configuration's first Contrast block.
    while mf.NumberOfOperands > 143:
        mf.RemoveOperandAt(mf.NumberOfOperands)
    out['first_block_merit'] = float(mf.CalculateMeritFunction())
    out['first_me_pairs'] = [dict(index=i,type=str(mf.GetOperandAt(i).Type),
        value=str(mf.GetOperandAt(i).Value),weight=str(mf.GetOperandAt(i).Weight))
        for i in range(1,mf.NumberOfOperands+1) if str(mf.GetOperandAt(i).Type) in ['MECS','MECT']][:10]
    (R / 'native_merit_problem_inspection.json').write_text(json.dumps(out, indent=2), encoding='utf8')
    print(json.dumps(out), flush=True)
finally:
    if o is not None:
        o.close(saveifneeded=False)
    if z.Application is not None:
        z.disconnect()
