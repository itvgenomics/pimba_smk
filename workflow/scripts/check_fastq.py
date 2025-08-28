#!/usr/bin/env python3
"""
check_fastq.py — folder-aware filename & FASTQ integrity checker

Usage:
  python3 check_fastq.py "/path/to/folder" [more folders or files...]

What it does:
  • Scans folders for FASTQ/FASTQ.GZ files (non-recursive by default; use --recursive).
  • Auto-groups paired files by standalone 'R1'/'R2' tokens in the filename.
  • Checks filenames:
      - Flags parentheses in names.
      - Flags more than one 'R1' or more than one 'R2' token.
      - Flags names that contain both 'R1' and 'R2' tokens.
  • Checks FASTQ integrity per file:
      - Each record must be 4 lines (@, seq, +, qual).
      - Sequence length must equal quality length.
  • Checks pairing integrity:
      - If both R1 and R2 exist for a pair, they must have the SAME number of records.

Exit code: 0 if everything passes; 1 if any problem (including unmatched pairs or count mismatches).
"""
import argparse
import gzip
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Union

# Detect R1/R2 tokens that stand alone (e.g., '_R1_' or '_R2.'), avoiding 'R10' or 'XR1A'.
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
    # Open as text, transparently handling .gz
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
                break  # EOF cleanly between records
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            n_records += 1

            if not (seq and plus and qual):
                issues.append(f'record {n_records}: file ends mid-record (not multiple of 4 lines)')
                n_format_err += 1
                break

            if not header.startswith('@'):
                if len(issues) < max_examples:
                    issues.append(f"record {n_records}: header does not start with '@' (got: {header[:30].rstrip()})")
                n_format_err += 1
            if not plus.startswith('+'):
                if len(issues) < max_examples:
                    issues.append(f"record {n_records}: third line does not start with '+' (got: {plus[:30].rstrip()})")
                n_format_err += 1

            s = seq.rstrip('\r\n')
            q = qual.rstrip('\r\n')
            if len(s) != len(q):
                n_len_mismatch += 1
                if len(issues) < max_examples:
                    hdr = header.strip()
                    issues.append(f'record {n_records}: sequence length {len(s)} != quality length {len(q)}; header: {hdr}')

    return {
        'records': n_records,
        'length_mismatches': n_len_mismatch,
        'format_errors': n_format_err,
        'examples': issues,
    }

def _strip_ext(name: str) -> str:
    # Remove .fastq(.gz) or .fq(.gz) from the end to get a base for pairing
    return _FASTQ_EXT_RE.sub('', name)

def pairing_key(p: Path) -> Tuple[str, str]:
    """
    Returns a tuple (key, role) where role is 'R1', 'R2', or 'UNK'.
    The key normalizes the first standalone R[12] token to 'R?' to match mates.
    """
    base = _strip_ext(p.name)
    m = _R_TOKEN_RE.search(base)
    role = 'UNK'
    key = base
    if m:
        role = 'R{}'.format(m.group(1))
        # Replace ONLY the first occurrence for a stable key
        key = _R_TOKEN_RE.sub('R?', base, count=1)
    return (key, role)

def iter_fastq_files(root: Path, recursive: bool = False):
    if root.is_file():
        if _FASTQ_EXT_RE.search(root.name):
            yield root
        return
    if not root.is_dir():
        return
    if recursive:
        for p in root.rglob('*'):
            if p.is_file() and _FASTQ_EXT_RE.search(p.name):
                yield p
    else:
        for p in root.iterdir():
            if p.is_file() and _FASTQ_EXT_RE.search(p.name):
                yield p

def main():
    ap = argparse.ArgumentParser(
        description='Scan folders for paired FASTQ files and check filenames, base/quality length, and pair record-count equality.'
    )
    ap.add_argument('roots', nargs='+', help='Folders and/or files to scan.')
    ap.add_argument('--recursive', action='store_true', help='Recurse into subdirectories.')
    ap.add_argument('--max-examples', type=int, default=10,
                    help='Max number of example issues to show per file (default: 10).')
    args = ap.parse_args()

    roots = [Path(r) for r in args.roots]
    all_files: List[Path] = []
    for r in roots:
        if not r.exists():
            print('ERROR: path does not exist -> {}'.format(r))
            sys.exit(1)
        all_files.extend(list(iter_fastq_files(r, recursive=args.recursive)))

    if not all_files:
        print('No FASTQ files found.')
        sys.exit(1)

    # Group by pairing key
    groups: Dict[str, dict] = {}
    for p in sorted(all_files):
        key, role = pairing_key(p)
        if key not in groups:
            groups[key] = {'R1': None, 'R2': None, 'UNK': []}
        if role in ('R1', 'R2'):
            if groups[key][role] is not None:
                print("WARNING: multiple files mapped to {} for key '{}': {} and {}".format(role, key, groups[key][role], p))
            groups[key][role] = p
        else:
            groups[key]['UNK'].append(p)

    any_problem = False

    # Process each group
    for key, slot in groups.items():
        R1 = slot['R1']
        R2 = slot['R2']
        UNK = slot['UNK']

        if R1 is None and R2 is None and not UNK:
            continue

        print('== Pair group: {} =='.format(key))
        if R1 and R2:
            print('  Found pair:')
            print('    R1: {}'.format(R1))
            print('    R2: {}'.format(R2))
        elif R1 and not R2:
            print('  UNMATCHED: found R1 but no R2')
            print('    R1: {}'.format(R1))
            any_problem = True
        elif R2 and not R1:
            print('  UNMATCHED: found R2 but no R1')
            print('    R2: {}'.format(R2))
            any_problem = True

        # Analyze R1/R2 first (store results to compare counts)
        per_file_results: Dict[Path, Dict[str, Union[int, List[str]]]] = {}
        for p in [x for x in (R1, R2) if x is not None]:
            print('--> {}'.format(p))
            # 1) Filename checks
            name_issues = check_filename(p)
            if name_issues:
                any_problem = True
                for msg in name_issues:
                    print('  [name] {}'.format(msg))
            else:
                print('  [name] OK')

            # 2) FASTQ content checks
            try:
                res = check_fastq_content(p, max_examples=args.max_examples)
                per_file_results[p] = res
            except Exception as e:
                print('  [fastq] ERROR reading file: {}'.format(e))
                any_problem = True
                continue

            print("  [fastq] records: {}".format(res['records']))
            print("  [fastq] length mismatches: {}".format(res['length_mismatches']))
            print("  [fastq] format errors: {}".format(res['format_errors']))

            if res['examples']:
                any_problem = True
                for ex in res['examples']:
                    print('    - {}'.format(ex))
            else:
                print('  [fastq] OK')

        # Pair-level record-count equality check
        if R1 and R2 and (R1 in per_file_results) and (R2 in per_file_results):
            r1c = per_file_results[R1]['records']  # type: ignore[index]
            r2c = per_file_results[R2]['records']  # type: ignore[index]
            print('  [pair] record counts -> R1: {}, R2: {}'.format(r1c, r2c))
            if r1c != r2c:
                any_problem = True
                print('  [pair] ERROR: R1 and R2 record counts differ for key {}.'.format(key))

        # Files with unknown role
        for p in UNK:
            print('--> {}'.format(p))
            name_issues = check_filename(p)
            if name_issues:
                any_problem = True
                for msg in name_issues:
                    print('  [name] {}'.format(msg))
            else:
                print('  [name] OK')

            try:
                res = check_fastq_content(p, max_examples=args.max_examples)
            except Exception as e:
                print('  [fastq] ERROR reading file: {}'.format(e))
                any_problem = True
                continue

            print("  [fastq] records: {}".format(res['records']))
            print("  [fastq] length mismatches: {}".format(res['length_mismatches']))
            print("  [fastq] format errors: {}".format(res['format_errors']))
            if res['examples']:
                any_problem = True
                for ex in res['examples']:
                    print('    - {}'.format(ex))
            else:
                print('  [fastq] OK')
            print("  NOTE: file doesn't contain a standalone R1/R2 token -> {}".format(p.name))

        print('')

    sys.exit(1 if any_problem else 0)

if __name__ == "__main__":
    main()
