#!/usr/bin/env python3
"""
check_fastq.py — folder-aware filename & FASTQ integrity checker

Usage:
  python3 check_fastq.py "/path/to/folder" [more folders or files...]

Exit code:
  0 if everything passes
  1 if any problem (empty read is FATAL)
"""
import argparse
import gzip
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Union

_R_TOKEN_RE = re.compile(r'(?<![A-Za-z0-9])R([12])(?![A-Za-z0-9])', re.IGNORECASE)
_FASTQ_EXT_RE = re.compile(r'\.(fastq|fq)(\.gz)?$', re.IGNORECASE)


def check_filename(p: Path) -> List[str]:
    issues: List[str] = []
    name = p.name

    if '(' in name or ')' in name:
        issues.append('filename contains parentheses')

    tokens = [m.group(1) for m in _R_TOKEN_RE.finditer(name)]
    c1 = sum(1 for t in tokens if t == '1')
    c2 = sum(1 for t in tokens if t == '2')

    if c1 > 1:
        issues.append(f"more than one 'R1' token (found {c1})")
    if c2 > 1:
        issues.append(f"more than one 'R2' token (found {c2})")
    if c1 >= 1 and c2 >= 1:
        issues.append('contains both R1 and R2 tokens (ambiguous)')

    return issues


def _open_text(p: Path):
    if p.suffix == '.gz':
        return gzip.open(p, 'rt', encoding='utf-8', errors='replace')
    return open(p, 'rt', encoding='utf-8', errors='replace')


def check_fastq_content(p: Path, max_examples: int = 10) -> Dict[str, Union[int, List[str]]]:
    issues: List[str] = []
    n_records = 0
    n_len_mismatch = 0
    n_format_err = 0

    with _open_text(p) as fh:
        while True:
            header = fh.readline()
            if not header:
                break  # EOF

            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            n_records += 1

            if not (seq and plus and qual):
                issues.append(f'record {n_records}: file ends mid-record')
                n_format_err += 1
                break

            header_clean = header.strip()

            if not header.startswith('@'):
                n_format_err += 1
                if len(issues) < max_examples:
                    issues.append(f"record {n_records}: header invalid -> {header[:30]}")

            if not plus.startswith('+'):
                n_format_err += 1
                if len(issues) < max_examples:
                    issues.append(f"record {n_records}: plus line invalid -> {plus[:30]}")

            s = seq.rstrip('\r\n')
            q = qual.rstrip('\r\n')

            # ✅ FATAL EMPTY READ CHECK WITH HEADER NAME
            if len(s) == 0 or len(q) == 0:
                raise RuntimeError(
                    f'EMPTY READ DETECTED in file {p}\n'
                    f'Header: {header_clean}'
                )

            if len(s) != len(q):
                n_len_mismatch += 1
                if len(issues) < max_examples:
                    issues.append(
                        f'record {n_records}: sequence length {len(s)} != quality length {len(q)}; header: {header_clean}'
                    )

    return {
        'records': n_records,
        'length_mismatches': n_len_mismatch,
        'format_errors': n_format_err,
        'examples': issues,
    }


def _strip_ext(name: str) -> str:
    return _FASTQ_EXT_RE.sub('', name)


def pairing_key(p: Path) -> Tuple[str, str]:
    base = _strip_ext(p.name)
    m = _R_TOKEN_RE.search(base)
    role = 'UNK'
    key = base
    if m:
        role = f'R{m.group(1)}'
        key = _R_TOKEN_RE.sub('R?', base, count=1)
    return key, role


def iter_fastq_files(root: Path, recursive: bool = False):
    if root.is_file():
        if _FASTQ_EXT_RE.search(root.name):
            yield root
        return

    if not root.is_dir():
        return

    files = root.rglob('*') if recursive else root.iterdir()
    for p in files:
        if p.is_file() and _FASTQ_EXT_RE.search(p.name):
            yield p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('roots', nargs='+')
    ap.add_argument('--recursive', action='store_true')
    ap.add_argument('--max-examples', type=int, default=10)
    args = ap.parse_args()

    roots = [Path(r) for r in args.roots]
    all_files: List[Path] = []

    for r in roots:
        if not r.exists():
            print(f'ERROR: path does not exist -> {r}')
            sys.exit(1)

        all_files.extend(list(iter_fastq_files(r, args.recursive)))

    if not all_files:
        print('No FASTQ files found.')
        sys.exit(1)

    groups: Dict[str, dict] = {}
    for p in sorted(all_files):
        key, role = pairing_key(p)
        groups.setdefault(key, {'R1': None, 'R2': None, 'UNK': []})
        if role in ('R1', 'R2'):
            if groups[key][role] is not None:
                print(f'WARNING: duplicate {role} in group {key}')
            groups[key][role] = p
        else:
            groups[key]['UNK'].append(p)

    any_problem = False

    for key, slot in groups.items():
        R1, R2, UNK = slot['R1'], slot['R2'], slot['UNK']
        print(f'== Pair group: {key} ==')

        per_file_results = {}

        for p in [x for x in (R1, R2) if x]:
            print(f'--> {p}')

            for msg in check_filename(p):
                print(f'  [name] {msg}')
                any_problem = True

            try:
                res = check_fastq_content(p, args.max_examples)
                per_file_results[p] = res
            except RuntimeError as e:
                print(f'\n[FATAL ERROR]\n{e}\n')
                sys.exit(1)

            print(f'  records: {res["records"]}')
            print(f'  length mismatches: {res["length_mismatches"]}')
            print(f'  format errors: {res["format_errors"]}')

            if res['examples']:
                any_problem = True
                for ex in res['examples']:
                    print(f'    - {ex}')
            else:
                print('  [fastq] OK')

        if R1 and R2 and R1 in per_file_results and R2 in per_file_results:
            c1 = per_file_results[R1]['records']
            c2 = per_file_results[R2]['records']
            if c1 != c2:
                print(f'  [pair] ERROR: R1={c1}, R2={c2}')
                any_problem = True

        print('')

    sys.exit(1 if any_problem else 0)


if __name__ == "__main__":
    main()
