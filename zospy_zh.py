"""Opt-in Simplified Chinese report adapters for ZOSPy 2.1.5 / OpticStudio 24R1.

Use these three wrappers instead of their zospy.analyses counterparts. The
installed package and OpticStudio language are unchanged. Unknown formats still
raise parser errors; no missing data or numeric values are fabricated.
"""
import re
from types import SimpleNamespace

from zospy.analyses.parsers.base import load_grammar, parse
from zospy.analyses.raysandspots.single_ray_trace import SingleRayTrace as _SingleRayTrace
from zospy.analyses.raysandspots.ray_fan import RayFan as _RayFan
from zospy.analyses.reports.system_data import SystemData as _SystemData


COMMON = {'文件': 'File', '题目': 'Title', '日期': 'Date', '单位': 'Units'}
RAY = {
    **COMMON, '波长': 'Wavelength', '坐标': 'Coordinates',
    '归一化X视场坐标(Hx)': 'Normalized X Field Coord (Hx)',
    '归一化Y视场坐标(Hy)': 'Normalized Y Field Coord (Hy)',
    '归一化X光瞳坐标(Px)': 'Normalized X Pupil Coord (Px)',
    '归一化Y光瞳坐标(Py)': 'Normalized Y Pupil Coord (Py)',
    '实际光线追迹数据': 'Real Ray Trace Data',
    '近轴光线追迹数据': 'Paraxial Ray Trace Data',
}
SYSTEM = {
    **COMMON, '常用透镜数据': 'GENERAL LENS DATA', '面': 'Surfaces', '光阑': 'Stop',
    '系统孔径': 'System Aperture', '半口径快速计算': 'Fast Semi-Diameters',
    '视场无偏振': 'Field Unpolarized', '将膜层相位转化为等效几何光线': 'Convert thin film phase to ray equivalent',
    'J/E 转化方法': 'J/E Conversion Method', '玻璃库': 'Glass Catalogs', '光线瞄准': 'Ray Aiming',
    '切趾法': 'Apodization', 'OPD参考': 'Reference OPD', '近轴光线': 'Paraxial Rays Setting',
    'F/#计算': 'Method to Compute F/#', '惠更斯积分计算': 'Method to Compute Huygens Integral',
    '显示坐标断点': 'Print Coordinate Breaks', '多线程': 'Multi-Threading', 'OPD以2π取模': 'OPD Modulo 2 Pi',
    '温度 (℃)': 'Temperature (C)', '压强 (ATM)': 'Pressure (ATM)',
    '折射率数据与环境匹配': 'Adjust Index Data To Environment', '有效焦距': 'Effective Focal Length',
    '后焦距': 'Back Focal Length', '总长': 'Total Track', '像方空间 F/#': 'Image Space F/#',
    '近轴处理 F/#': 'Paraxial Working F/#', '工作F/#': 'Working F/#',
    '像方空间 NA': 'Image Space NA', '物方空间 NA': 'Object Space NA', '光阑半径': 'Stop Radius',
    '近轴成像高度': 'Paraxial Image Height', '近轴放大率': 'Paraxial Magnification',
    '入瞳直径': 'Entrance Pupil Diameter', '入瞳位置': 'Entrance Pupil Position',
    '出瞳直径': 'Exit Pupil Diameter', '出瞳位置': 'Exit Pupil Position', '视场类型': 'Field Type',
    '最大径向视场': 'Maximum Radial Field', '主波长': 'Primary Wavelength [µm]',
    '角放大率': 'Angular Magnification', '透镜单位': 'Lens Units', '光源单位': 'Source Units',
    '分析单位': 'Analysis Units', '无焦模式单位': 'Afocal Mode Units', 'MTF 单位': 'MTF Units',
    '计算数据保存于Session文件': 'Include Calculated Data in Session File',
    '视场': 'Fields', '波长': 'Wavelengths', '预测坐标 ABCD 矩阵': 'Predicted coordinate ABCD matrix',
}
BOOL_KEYS = {
    'Fast Semi-Diameters', 'Field Unpolarized', 'Convert thin film phase to ray equivalent',
    'Ray Aiming', 'Print Coordinate Breaks', 'Multi-Threading', 'OPD Modulo 2 Pi',
    'Adjust Index Data To Environment', 'Include Calculated Data in Session File',
}


def _number_tokens(value):
    # The upstream grammar misses scientific notation with an integer mantissa.
    return re.sub(r'(?<![\w.])([+-]?\d+)([eE][+-]?\d+)(?!\w)', r'\1.0\2', value)


def normalize_text(text, grammar):
    mapping = {'single_ray_trace': RAY, 'system_data': SYSTEM}[grammar]
    if not any(key in text for key in ('光线追迹数据', '常用透镜数据')):
        return text
    lines = []
    for line in text.splitlines():
        stripped = line.lstrip('\ufeff').strip()
        if not stripped:
            lines.append('')
            continue
        match = re.match(r'^\s*([^:：]+?)\s*[:：](.*)$', line.lstrip('\ufeff'))
        if match and (key := match[1].strip()) in mapping:
            name, value = mapping[key], match[2]
            if grammar == 'single_ray_trace' and key == '坐标':
                value = {'局部': 'Local', '全局': 'Global'}.get(value.strip(), value)
            if grammar == 'system_data':
                if name in BOOL_KEYS:
                    value = {'关': 'Off', '开': 'On'}.get(value.strip(), value)
                elif key == '切趾法':
                    value = value.replace('因子', 'factor')
                    # Only the categorical apodization name changes; preserve its number.
                    value = value.replace('均匀', 'Uniform').replace('高斯', 'Gaussian')
                elif key == '有效焦距':
                    value = value.replace('(在系统温度和压强的空气中)', '(in air at system temperature and pressure)')
                    value = value.replace('(在像方空间)', '(in image space)')
                elif key == '主波长':
                    value = re.sub(r'\s*[µμ]m\s*$', '', value)
                if key not in {'文件', '题目', '日期'}:
                    value = _number_tokens(value)
            line = f'{name}: {value}'
        elif grammar == 'single_ray_trace' and stripped.startswith('表面\t'):
            header = stripped
            for zh, en in {'表面': 'Surf', '坐标': 'coordinate', '余弦': 'cosine', '法线': 'normal',
                           '角': 'Angle in', '路径长度': 'Path length', '注释': 'Comment'}.items():
                header = header.replace(zh, en)
            line = header
        elif grammar == 'single_ray_trace' and stripped.startswith('物面'):
            line = re.sub(r'^\s*物面\s+无限\s+无限\s+无限', 'OBJ infinity infinity infinity', line)
        elif grammar == 'system_data' and stripped == '渐晕因子':
            line = 'Vignetting Factors'
        elif grammar == 'system_data' and stripped.startswith('#'):
            line = line.replace('X-角度', 'X-Value').replace('Y-角度', 'Y-Value')
            line = line.replace('X-高度', 'X-Value').replace('Y-高度', 'Y-Value')
            line = line.replace('权重', 'Weight').replace('值', 'Value')
        lines.append(line)
    return '\n'.join(lines) + '\n'


def normalize_fan_description(description):
    match = re.fullmatch(r'(子午|弧矢)光扇[，,]\s*视场号\s*(\d+)\s*=\s*(.+)', description)
    if not match:
        return description
    direction = {'子午': 'Tangential', '弧矢': 'Sagittal'}[match[1]]
    # Coordinate values and unit strings are kept verbatim.
    return f'{direction} fan, field number {match[2]} = {match[3]}'


class _ChineseText:
    def parse_output(self, grammar, transformer, result_type=None):
        parsed = parse(normalize_text(self.get_text_output(), grammar), load_grammar(grammar), transformer)
        return parsed if result_type is None else result_type(**parsed)


class SingleRayTrace(_ChineseText, _SingleRayTrace, analysis_type='RayTrace', needs_text_output_file=True):
    """SingleRayTrace with Chinese text labels normalized before validation."""


class SystemData(_ChineseText, _SystemData, analysis_type='SystemData', needs_text_output_file=True):
    """SystemData with Chinese text labels normalized before validation."""


class RayFan(_RayFan, analysis_type='RayFan'):
    """RayFan using native numeric arrays and normalized series descriptions."""

    def get_data_series(self):
        results = self.analysis.Results
        series = []
        for i in range(results.NumberOfDataSeries):
            ds = results.GetDataSeries(i)
            series.append(SimpleNamespace(Description=normalize_fan_description(str(ds.Description)),
                                          XData=ds.XData, XLabel=ds.XLabel,
                                          SeriesLabels=ds.SeriesLabels, YData=ds.YData))
        proxy = SimpleNamespace(analysis=SimpleNamespace(Results=SimpleNamespace(
            NumberOfDataSeries=len(series), GetDataSeries=lambda i: series[i])))
        return _RayFan.get_data_series(proxy)
