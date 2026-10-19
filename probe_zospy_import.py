import os
import json
import zospy as zp

zos = zp.ZOS(opticstudio_directory=os.environ.get('OPTICSTUDIO_DIR', r'D:\Zemax'))
try:
    oss = zos.connect('extension', instance_number=0, message_logging=False)
    def members(obj):
        return sorted({str(m.Name) for m in obj.GetType().GetMembers()})
    result = {
        'zospy': zp.__version__, 'file': str(oss.SystemFile),
        'needs_save': bool(oss.NeedsSave), 'surfaces': int(oss.LDE.NumberOfSurfaces),
        'mode': str(zos.Application.Mode), 'license': str(zos.Application.LicenseStatus),
        'api_valid': bool(zos.Application.IsValidLicenseForAPI),
        'system_members': members(oss._System),
        'mce_members': members(oss.MCE),
        'mce_row_members': members(oss.MCE.GetOperandAt(1)),
        'fields_members': members(oss.SystemData.Fields),
        'waves_members': members(oss.SystemData.Wavelengths),
        'wave_row_members': members(oss.SystemData.Wavelengths.GetWavelength(1)),
        'aperture_members': members(oss.SystemData.Aperture),
        'field_types': str(zp.constants.SystemData.FieldType),
        'aperture_types': str(zp.constants.SystemData.ZemaxApertureType),
        'mce_types': str(zp.constants.Editors.MCE.MultiConfigOperandType),
    }
    print(json.dumps(result, indent=2, ensure_ascii=True), flush=True)
finally:
    if zos.Application is not None:
        zos.disconnect()
