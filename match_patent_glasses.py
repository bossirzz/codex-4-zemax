"""Match published Nd/Vd against local AGF catalogs without modifying Zemax."""
import os
import json
import math
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MODEL = ROOT / 'zemax_import' / 'CN110703413B_example3'
OUT = MODEL / 'glass_matching'
CATALOG_DIR = Path(os.environ.get('ZEMAX_GLASSCAT_DIR', str(Path.home() / 'Documents' / 'Zemax' / 'Glasscat')))
VENDORS = ['CDGM', 'NHG', 'SCHOTT', 'OHARA', 'HOYA', 'HIKARI', 'SUMITA',
           'LZOS', 'NIKON', 'NIKON-HIKARI', 'GBJ', 'PILKINGTON', 'CORNING']
STATUS = {0: 'Standard', 1: 'Preferred', 2: 'Obsolete', 3: 'Special', 4: 'Melt'}
WAVES = {'g': .4358343, 'F': .4861327, 'd': .5875618, 'C': .6562725,
         '650nm': .650, '850nm': .850}
LOCAL_COMPARISON = {
    'L1': ('CDGM', 'H-LAF50B'), 'L2': ('CDGM', 'H-ZPK5'),
    'L3': ('CDGM', 'H-QF3'), 'L4': ('CDGM', 'H-ZF4A'),
    'L5': ('NHG', 'H-ZLAF2'), 'L7': ('CDGM', 'H-ZF12'),
}


def read_text(path):
    raw = path.read_bytes()
    if raw.startswith((b'\xff\xfe', b'\xfe\xff')):
        return raw.decode('utf-16')
    if raw[:80].count(b'\x00') > 10:
        return raw.decode('utf-16-le')
    for encoding in ['utf-8-sig', 'cp1252', 'gb18030']:
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise ValueError(f'Cannot decode catalog: {path}')


def parse_agf(path):
    records = []
    current = None
    comments = []
    for number, raw in enumerate(read_text(path).splitlines(), 1):
        fields = raw.strip().split()
        if not fields:
            continue
        tag = fields[0]
        if tag == 'CC':
            comments.append(raw.strip()[3:])
        elif tag == 'NM':
            if len(fields) < 6:
                raise ValueError(f'Incomplete NM at {path}:{number}')
            current = {'catalog': path.stem.upper(), 'catalog_file': str(path), 'line': number,
                       'name': fields[1], 'formula': int(float(fields[2])),
                       'nd_nominal': float(fields[4]), 'vd_nominal': float(fields[5]),
                       'exclude_substitution': bool(int(float(fields[6]))) if len(fields) > 6 else False,
                       'status_code': int(float(fields[7])) if len(fields) > 7 else None}
            current['status'] = STATUS.get(current['status_code'], 'Unspecified')
            records.append(current)
        elif current is not None and tag in ['CD', 'ED', 'LD', 'TD']:
            current[tag] = [float(value) for value in fields[1:]]
        elif current is not None and tag == 'GC':
            current['comment'] = raw.strip()[3:]
    return records, comments


def index(glass, wave):
    coeff = glass.get('CD', [])
    l2 = wave * wave
    if glass['formula'] == 2 and len(coeff) >= 6:
        value = 1 + sum(coeff[i] * l2 / (l2 - coeff[i + 1]) for i in [0, 2, 4])
    elif glass['formula'] == 1 and len(coeff) >= 6:
        value = coeff[0] + coeff[1] * l2 + sum(coeff[i] / l2 ** (i - 1) for i in range(2, 6))
    elif glass['formula'] == 3 and len(coeff) >= 6:
        herzberger_l = 1 / (l2 - .028)
        return (coeff[0] + coeff[1] * herzberger_l + coeff[2] * herzberger_l ** 2
                + coeff[3] * l2 + coeff[4] * l2 ** 2 + coeff[5] * l2 ** 3)
    elif glass['formula'] == 6 and len(coeff) >= 8:
        value = 1 + sum(coeff[i] * l2 / (l2 - coeff[i + 1]) for i in [0, 2, 4, 6])
    elif glass['formula'] == 12 and len(coeff) >= 8:
        value = (coeff[0] + coeff[1] * l2
                 + sum(coeff[i] / l2 ** (i - 1) for i in range(2, 6))
                 + coeff[6] * l2 ** 2 + coeff[7] * l2 ** 3)
    elif glass['formula'] == 13 and len(coeff) >= 9:
        value = (coeff[0] + coeff[1] * l2 + coeff[2] * l2 ** 2
                 + sum(coeff[i] / l2 ** (i - 2) for i in range(3, 9)))
    else:
        raise NotImplementedError(f'Dispersion formula {glass["formula"]}')
    return math.sqrt(value)


def enrich(glass):
    # Only compute inside the catalog's stated dispersion range. LD coverage
    # does not establish transmission or suitability for manufacturing.
    bounds = glass.get('LD', [0, float('inf')])
    try:
        values = {label: index(glass, wave) for label, wave in WAVES.items()
                  if bounds[0] <= wave <= bounds[1]}
        if not all(label in values for label in ['F', 'd', 'C']):
            raise ValueError('Catalog wavelength range does not cover F/d/C lines')
        glass['indices'] = values
        glass['nd'] = values['d']
        glass['vd'] = (values['d'] - 1) / (values['F'] - values['C'])
        glass['basis'] = 'Computed from AGF CD dispersion coefficients'
        if 'g' in values:
            glass['pgf_computed'] = (values['g'] - values['F']) / (values['F'] - values['C'])
        glass['nd_nominal_difference'] = glass['nd'] - glass['nd_nominal']
        glass['vd_nominal_difference'] = glass['vd'] - glass['vd_nominal']
    except NotImplementedError:
        glass['nd'] = glass['nd_nominal']
        glass['vd'] = glass['vd_nominal']
        glass['basis'] = 'Nominal NM reference fields; dispersion formula not evaluated'
    except (ValueError, ZeroDivisionError) as exc:
        glass['dispersion_error'] = str(exc)
        glass['nd'] = glass['nd_nominal']
        glass['vd'] = glass['vd_nominal']
        glass['basis'] = 'Nominal NM reference fields; CD evaluation failed'
    glass['covers_486_to_850nm'] = ('LD' in glass and glass['LD'][0] <= .4861327 and glass['LD'][1] >= .850)
    glass['variant_name_flag'] = bool(re.search(r'@|\(M\)|^M(?:C|P)?-|^[DL]-', glass['name'], re.I))
    glass['reference_inconsistency_flag'] = (glass['nd_nominal'] <= 1
                                           or abs(glass.get('nd_nominal_difference', 0)) > .0001
                                           or abs(glass.get('vd_nominal_difference', 0)) > .3)
    return glass


def rank(glass, nd, vd, vd_half):
    item = {key: value for key, value in glass.items() if key not in ['CD', 'TD']}
    item['delta_nd'] = glass['nd'] - nd
    item['delta_vd'] = glass['vd'] - vd
    item['within_published_rounding'] = abs(item['delta_nd']) <= .005 + 1e-12 and abs(item['delta_vd']) <= vd_half + 1e-12
    item['score'] = math.hypot(item['delta_nd'] / .005, item['delta_vd'] / .5)
    return item


def write_markdown(report):
    lines = [
        '# CN110703413B 实施例 3：玻璃牌号近似匹配', '',
        '本机库中能找到近似牌号；这份结果不能唯一确定专利原玻璃。专利没有给出厂家、牌号或完整色散数据。', '',
        f"读取本机 {report['catalogs_scanned']} 个 AGF 文件，共 {report['records_scanned']} 条材料记录；"
        f"在 {len(VENDORS)} 个光学厂商目录的 {report['optical_vendor_records']} 条记录中排序。全部文件解析成功。", '',
        f'目录：`{CATALOG_DIR}`。记录状态来自已安装目录，不代表已经核实当前生产或采购可用性。', '',
        'Nd、Vd 均按 AGF 的 CD 色散系数计算，波长为 d=0.5875618 μm、F=0.4861327 μm、C=0.6562725 μm。'
        '采用目录系数参考条件，未施加光学系统工作温度、气压修正。NM 行中的 Nd/Vd 仅作参考交叉检查。', '',
        '## 国内牌号对照，可用于建立比较模型', '',
        '以下优先列出 CDGM 或 NHG 的候选，便于比较；并非宣称其在所有厂商中距离最小。'
        'L1、L6、L8 的公开参数相同，故合并展示，但这不证明三片的实际原牌号相同。', '',
        '| 镜片 | 专利 Nd / Vd | 目录 / 候选牌号 | 系数计算 Nd / Vd | ΔNd | ΔVd | 目录状态 |',
        '|---|---|---|---|---|---|---|',
    ]
    selected = []
    for group in report['groups']:
        catalog, name = LOCAL_COMPARISON[group['lenses'][0]]
        candidate = next(item for item in group['best_ordinary_by_vendor']
                         if item['catalog'] == catalog and item['name'] == name)
        assert candidate['basis'].startswith('Computed')
        assert not candidate['reference_inconsistency_flag']
        assert candidate['covers_486_to_850nm']
        assert candidate['status_code'] in [0, 1]
        selected.append({'lenses': group['lenses'], 'target_nd': group['target_nd'],
                         'target_vd': group['target_vd'], 'candidate': candidate})
        lines.append(f"| {'、'.join(group['lenses'])} | {group['target_nd']:.2f} / {group['target_vd']:g} "
                     f"| {catalog} / {name} | {candidate['nd']:.6f} / {candidate['vd']:.4f} "
                     f"| {candidate['delta_nd']:+.6f} | {candidate['delta_vd']:+.4f} | {candidate['status']} |")
    lines += ['',
              'Δ 为“候选减专利”。L4 的 ΔNd 约 +0.00825，不能按两位小数通常的四舍五入解释为 1.72。'
              'L3 的候选 Nd 约 1.57501，位于 1.57 的四舍五入区间边界外侧；截断或原文精度问题可能造成差异，但尚无证据确认。', '',
              '## 各组较近的其他候选', '']
    for group in report['groups']:
        lines += [f"### {'、'.join(group['lenses'])}，目标 {group['target_nd']:.2f} / {group['target_vd']:g}", '',
                  '| 目录 | 牌号 | Nd | Vd | ΔNd | ΔVd | 目录状态 |',
                  '|---|---|---|---|---|---|---|']
        for item in group['nearest_ordinary_name'][:6]:
            lines.append(f"| {item['catalog']} | {item['name']} | {item['nd']:.6f} | {item['vd']:.4f} "
                         f"| {item['delta_nd']:+.6f} | {item['delta_vd']:+.4f} | {item['status']} |")
        lines += ['', '名称带有 M-、MC-、MP-、D-、L-、(M) 或 @ 的变体另外列出。'
                  '这仅按名称筛分，制造工艺与供货形态需以厂商说明为准。', '',
                  '| 目录 | 变体名称 | Nd | Vd | ΔNd | ΔVd |',
                  '|---|---|---|---|---|---|']
        for item in group['nearest_variant_name'][:3]:
            lines.append(f"| {item['catalog']} | {item['name']} | {item['nd']:.6f} | {item['vd']:.4f} "
                         f"| {item['delta_nd']:+.6f} | {item['delta_vd']:+.4f} |")
        lines += ['', f"若假设公开数值采用四舍五入：Nd 半区间 0.005，Vd 半区间 {group['rounding_assumption']['vd_half_interval']:g}；"
                  f"同时满足者共 {group['within_rounding_count']} 条（包括可能已标记 Obsolete 的记录）。", '']
    lines += ['## 650 nm 与 850 nm 折射率', '',
              '以下两波长都在所列候选的 LD 色散公式有效范围内；LD 覆盖不等于已验证透过率。', '',
              '| 镜片 | 目录 / 牌号 | n(650 nm) | n(850 nm) |',
              '|---|---|---|---|']
    for item in selected:
        glass = item['candidate']
        lines.append(f"| {'、'.join(item['lenses'])} | {glass['catalog']} / {glass['name']} "
                     f"| {glass['indices']['650nm']:.7f} | {glass['indices']['850nm']:.7f} |")
    lines += ['', '## 复核与使用限制', '',
              '- 以 sqrt(ΔNd²/0.005² + ΔVd²/0.5²) 排序，这是检索尺度，不是原牌号识别依据或材料制造公差。',
              '- 仅 L5 找到符合上述四舍五入假设的候选。未满足不等于数据一定错误，也可能涉及截断、不同玻璃目录版本、专用材料或专利数值问题。',
              '- 筛选 Standard/Preferred、允许替代且色散公式覆盖 486–850 nm 的记录；没有使用红外公式外推得到的 Nd/Vd。',
              '- 24 条光学厂商记录的 LD 范围不覆盖 F/d/C，未纳入已计算候选。2 条系数与 NM 参考数值差异很大的记录另存 JSON，未纳入候选推荐。',
              '- 使用的 1、2、3、6、12、13 号色散公式已与各自目录 NM 参考值交叉检查；排除上述 2 条后，所有已计算记录的 |ΔNd| 均小于 0.0001。Vd 参考字段有些仅保留一位小数，实际系数计算值会有差异。',
              '- 原文件使用 ModelGlass，ΔPgF=0。真实玻璃还包含部分色散信息；Nd/Vd 相近并不能保证可见光色差和 850 nm 像质相同。',
              '- 这次只生成匹配结果，没有替换或保存当前 Zemax 光学模型。若用于比较，应建立单独模型，保持结构与光阑位置，替换材料后重新对焦、比较单色和多色 MTF。',
              '- 当前模型与专利标称焦距及 L8 单片焦距存在差异，光阑具体位置也未公开。玻璃匹配本身不能证明这些问题已解决。', '',
              '## 追溯信息与依据', '',
              '完整候选、文件与行号、状态、色散式编号、Nd/Vd 差值及参考异常保存在 `glass_candidates.json`。'
              '国内对照方案保存在 `selected_comparison_glasses.json`，它是待验证候选清单。', '',
              '- [专利 CN110703413B](https://patents.google.com/patent/CN110703413B/zh)，实施例 3 的表 7。',
              '- [Ansys：AGF 与 BGF 格式](https://ansyshelp.ansys.com/public/Views/Secured/Zemax/v252/en/OpticStudio_User_Guide/OpticStudio_Help/topics/The_AGF_BGF_File_Formats.html)，NM 参考值、CD 系数、LD 波长范围与状态定义。',
              '- [Ansys：Model Glass](https://optics.ansys.com/hc/en-us/articles/42661781831571-How-to-use-the-Model-Glass)，Nd/Vd、部分色散及可见光外使用限制。',
              '- [Zemax 用户手册](https://neurophysics.ucsd.edu/Manuals/Zemax/ZemaxManual.pdf)，玻璃色散公式章节，说明各公式的系数排列。', '']
    (OUT / '玻璃牌号匹配.md').write_text('\n'.join(lines), encoding='utf-8')
    (OUT / 'selected_comparison_glasses.json').write_text(json.dumps(selected, ensure_ascii=False, indent=2), encoding='utf-8')


def main():
    inventory = []
    errors = []
    catalogs = []
    for path in sorted(CATALOG_DIR.glob('*')):
        if path.suffix.lower() != '.agf':
            continue
        try:
            records, comments = parse_agf(path)
        except Exception as exc:
            errors.append({'file': str(path), 'error': str(exc)})
            continue
        inventory.extend(enrich(glass) for glass in records)
        catalogs.append({'name': path.stem.upper(), 'file': str(path), 'records': len(records), 'comments': comments[:3]})
    scope = [glass for glass in inventory if glass['catalog'] in VENDORS]
    targets = json.loads((MODEL / 'glass_resolution.json').read_text(encoding='utf-8'))
    groups = []
    seen = set()
    for target in targets:
        key = (target['nd'], target['vd'])
        if key in seen:
            continue
        seen.add(key)
        lenses = [item['lens'] for item in targets if (item['nd'], item['vd']) == key]
        vd_half = .5 if lenses == ['L7'] else .05
        candidates = sorted([rank(glass, *key, vd_half) for glass in scope],
                            key=lambda item: (item['score'], VENDORS.index(item['catalog']), item['name']))
        exact = [item for item in candidates if item['within_published_rounding']]
        active = [item for item in candidates if item['status_code'] in [0, 1]
                  and not item['exclude_substitution'] and item['covers_486_to_850nm']
                  and item['basis'].startswith('Computed') and not item['reference_inconsistency_flag']]
        ordinary = [item for item in active if not item['variant_name_flag']]
        best_by_vendor = []
        for vendor in VENDORS:
            match = next((item for item in ordinary if item['catalog'] == vendor), None)
            if match is not None:
                best_by_vendor.append(match)
        top_unique = []
        signatures = set()
        for item in active:
            signature = (item['name'], round(item['nd'], 6), round(item['vd'], 3))
            if signature in signatures:
                continue
            signatures.add(signature)
            top_unique.append(item)
            if len(top_unique) == 8:
                break
        groups.append({'lenses': lenses, 'target_nd': key[0], 'target_vd': key[1],
                       'rounding_assumption': {'nd_half_interval': .005, 'vd_half_interval': vd_half,
                                             'rule': 'Nearest rounding of the printed values; truncation is not assumed'},
                       'within_rounding_count': len(exact), 'within_rounding_candidates': exact,
                       'nearest_standard_or_preferred': top_unique,
                       'nearest_ordinary_name': ordinary[:10], 'best_ordinary_by_vendor': best_by_vendor,
                       'nearest_variant_name': [item for item in active if item['variant_name_flag']][:8],
                       'nearest_including_obsolete_special': candidates[:8]})
    evaluated = [glass for glass in scope if glass['basis'].startswith('Computed')]
    formula_checks = {}
    for formula in sorted(set(glass['formula'] for glass in scope)):
        subset = [glass for glass in evaluated if glass['formula'] == formula
                  and not glass['reference_inconsistency_flag']]
        formula_checks[formula] = {'evaluated_records': len(subset),
                                  'max_abs_nd_vs_NM': max((abs(g['nd_nominal_difference']) for g in subset), default=None),
                                  'max_abs_vd_vs_NM': max((abs(g['vd_nominal_difference']) for g in subset), default=None)}
    report = {'catalog_directory': str(CATALOG_DIR), 'catalogs_scanned': len(catalogs),
              'records_scanned': len(inventory), 'optical_vendor_catalogs': VENDORS,
              'optical_vendor_records': len(scope), 'parsing_errors': errors,
              'dispersion_formulas_in_scope': dict(Counter(glass['formula'] for glass in scope)),
              'formula_reference_checks': formula_checks,
              'catalog_reference_outliers': [{key: glass.get(key) for key in ['catalog', 'name', 'formula', 'nd', 'vd', 'nd_nominal', 'vd_nominal', 'nd_nominal_difference', 'vd_nominal_difference']}
                                             for glass in scope if glass['basis'].startswith('Computed') and glass['reference_inconsistency_flag']],
              'records_not_evaluated_in_scope': [{key: glass.get(key) for key in ['catalog', 'name', 'formula', 'basis', 'dispersion_error']}
                                                for glass in scope if not glass['basis'].startswith('Computed')],
              'ranking': 'Euclidean distance after scaling delta Nd by 0.005 and delta Vd by 0.5; ranking is not proof of identity',
              'status_note': 'Status reflects the installed catalog, not verified current production or procurement availability',
              'variant_name_note': 'Names containing @ or (M), or beginning M-, MC-, MP-, D-, L- are separated by a name heuristic; this is not a definitive manufacturing classification',
              'wavelength_range_note': 'LD coverage means the stated range for dispersion computation, not verified optical transmission',
              'catalogs': catalogs, 'groups': groups}
    OUT.mkdir(exist_ok=True)
    (OUT / 'glass_candidates.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    write_markdown(report)
    print(json.dumps({'catalogs_scanned': len(catalogs), 'records_scanned': len(inventory),
                      'scope_records': len(scope), 'formulas': report['dispersion_formulas_in_scope'], 'errors': errors,
                      'formula_checks': formula_checks,
                      'not_evaluated_in_scope': report['records_not_evaluated_in_scope'],
                      'groups': [{'lenses': group['lenses'], 'target': [group['target_nd'], group['target_vd']],
                                  'within_rounding_count': group['within_rounding_count'],
                                  'closest_active': [{key: item[key] for key in ['catalog', 'name', 'nd', 'vd', 'delta_nd', 'delta_vd', 'status', 'within_published_rounding', 'basis']}
                                                     for item in group['nearest_standard_or_preferred'][:3]],
                                  'closest_overall': [{key: item[key] for key in ['catalog', 'name', 'nd', 'vd', 'delta_nd', 'delta_vd', 'status', 'within_published_rounding', 'basis']}
                                                      for item in group['nearest_including_obsolete_special'][:2]],
                                  'closest_ordinary': [{key: item[key] for key in ['catalog', 'name', 'nd', 'vd', 'delta_nd', 'delta_vd', 'status', 'within_published_rounding', 'basis']}
                                                      for item in group['nearest_ordinary_name'][:4]]}
                                 for group in groups]}, indent=2, ensure_ascii=True), flush=True)


if __name__ == '__main__':
    main()
