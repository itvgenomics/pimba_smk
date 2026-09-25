#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import re
from pathlib import Path

import pandas as pd


# =========================================================
# Command-line arguments
# =========================================================

def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description=(
            'Create a single HTML report from PIMBA Plot outputs, '
            'optionally including PIMBA-Curate and Random Forest results.'
        )
    )

    parser.add_argument(
        '--results-dir',
        type=Path,
        default=Path("results/03-plot"),
        help='PIMBA plotting results directory. Default: results/03-plot'
    )

    parser.add_argument(
        '--output',
        type=Path,
        default=None,
        help='Output HTML file. Default: <results-dir>/PIMBA_report.html'
    )

    parser.add_argument(
        '--title',
        default='PIMBA Analysis Report',
        help='Report title.'
    )

    parser.add_argument(
        '--ordination-marker-size',
        type=float,
        default=16,
        help=(
            'Marker size for interactive NMDS Bray-Curtis and '
            'PCoA unweighted UniFrac plots. Default: 16.'
        )
    )

    parser.add_argument(
        '--include-curate',
        action='store_true',
        help='Include PIMBA-Curate results in the HTML report.'
    )

    parser.add_argument(
        '--curate-xlsx',
        type=Path,
        default=None,
        help=(
            'Path to a specific PIMBA-Curate Excel table. '
            'Only used with --include-curate. If omitted, the script '
            'searches <results-root>/02-curate/*/*.xlsx.'
        )
    )

    parser.add_argument(
        '--include-rf',
        action='store_true',
        help='Include Random Forest results in the HTML report.'
    )

    parser.add_argument(
        '--rf-dir',
        type=Path,
        default=None,
        help=(
            'Path to the Random Forest results directory. '
            'Only used with --include-rf. If omitted, the script '
            'tries <results-dir>/random_forest.'
        )
    )

    return parser.parse_args()


# =========================================================
# Resolve PIMBA-Curate Excel table
# =========================================================

def resolve_curate_xlsx(
    results_dir: Path,
    curate_xlsx: Path | None
) -> Path:

    # Explicit path supplied
    if curate_xlsx is not None:

        curate_xlsx = curate_xlsx.resolve()

        if not curate_xlsx.is_file():

            raise FileNotFoundError(
                'PIMBA-Curate Excel table not found:\n'
                f'{curate_xlsx}\n\n'
                'A valid --curate-xlsx parameter is required.'
            )

        print(
            f'PIMBA-Curate Excel table: '
            f'{curate_xlsx}'
        )

        return curate_xlsx

    # Automatic search
    curate_dir = (
        results_dir.parent
        / "02-curate"
    )

    if not curate_dir.is_dir():

        raise FileNotFoundError(
            'PIMBA-Curate directory was not found:\n'
            f'{curate_dir}\n\n'
            'Use --curate-xlsx to specify the Excel table explicitly.'
        )

    xlsx_files = sorted(
        path.resolve()
        for path in curate_dir.glob("*/*.xlsx")
        if path.is_file()
    )

    if not xlsx_files:

        raise FileNotFoundError(
            'No PIMBA-Curate Excel table was found under:\n'
            f'{curate_dir}/*/*.xlsx\n\n'
            'Use --curate-xlsx to specify the Excel table explicitly.'
        )

    if len(xlsx_files) == 1:

        print(
            f'PIMBA-Curate Excel table found: '
            f'{xlsx_files[0]}'
        )

        return xlsx_files[0]

    # Prefer the standard Curate filename
    preferred_name = (
        'tax_assignments__to_validate_singleSeq.xlsx'
    )

    preferred_files = [
        path
        for path in xlsx_files
        if path.name == preferred_name
    ]

    if len(preferred_files) == 1:

        print(
            f'PIMBA-Curate Excel table found: '
            f'{preferred_files[0]}'
        )

        return preferred_files[0]

    file_list = '\n'.join(
        f'  - {path}'
        for path in xlsx_files
    )

    raise RuntimeError(
        'Multiple PIMBA-Curate Excel tables were found and '
        'the script cannot determine which one should be used.\n\n'
        f'Files found:\n{file_list}\n\n'
        'The --curate-xlsx parameter is required. Example:\n\n'
        '--curate-xlsx "results/02-curate/folder_name/table.xlsx"'
    )


# =========================================================
# Resolve Random Forest directory
# =========================================================

def resolve_rf_dir(
    results_dir: Path,
    rf_dir: Path | None
) -> Path:

    # Explicit path supplied
    if rf_dir is not None:

        rf_dir = rf_dir.resolve()

        if not rf_dir.is_dir():

            raise FileNotFoundError(
                'Random Forest directory not found:\n'
                f'{rf_dir}\n\n'
                'A valid --rf-dir parameter is required.'
            )

        print(
            f'Random Forest directory: '
            f'{rf_dir}'
        )

        return rf_dir

    # Default location
    default_rf_dir = (
        results_dir
        / "random_forest"
    )

    if not default_rf_dir.is_dir():

        raise FileNotFoundError(
            'Default Random Forest directory was not found:\n'
            f'{default_rf_dir}\n\n'
            'Use --rf-dir to specify the Random Forest '
            'directory explicitly.'
        )

    default_rf_dir = (
        default_rf_dir.resolve()
    )

    print(
        f'Random Forest directory found: '
        f'{default_rf_dir}'
    )

    return default_rf_dir


# =========================================================
# Human-readable file titles
# =========================================================

def natural_title(
    path: Path
) -> str:

    name = path.stem

    replacements = {
        'rf_': '',
        '_': ' ',
        'braycurtis': 'Bray-Curtis',
        'unifrac': 'UniFrac',
        'top20': 'Top 20',
    }

    for old, new in replacements.items():

        name = name.replace(
            old,
            new
        )

    preserve = {
        'PCA',
        'PCOA',
        'NMDS',
        'MDS',
        'ASV',
        'ASVS',
        'RF',
    }

    words = []

    for word in name.split():

        if word.upper() in preserve:

            if word.upper() == 'PCOA':

                words.append(
                    'PCoA'
                )

            elif word.upper() == 'ASVS':

                words.append(
                    'ASVs'
                )

            else:

                words.append(
                    word.upper()
                )

        elif word.lower() == 'bray-curtis':

            words.append(
                'Bray-Curtis'
            )

        elif word.lower() == 'unifrac':

            words.append(
                'UniFrac'
            )

        else:

            words.append(
                word.capitalize()
            )

    return ' '.join(
        words
    )


# =========================================================
# SVG rendering
# =========================================================

def svg_to_html(
    path: Path
) -> str:

    svg = path.read_text(
        encoding='utf-8',
        errors='replace'
    )

    svg = svg.replace(
        '<?xml version="1.0" encoding="utf-8" standalone="no"?>',
        ''
    )

    svg = svg.replace(
        '<?xml version="1.0" encoding="utf-8"?>',
        ''
    )

    return (
        '<div class="figure-card">'
        f'<h3>{html.escape(natural_title(path))}</h3>'
        f'<div class="svg-wrap">{svg}</div>'
        f'<div class="file-name">{html.escape(path.name)}</div>'
        '</div>'
    )


# =========================================================
# Interactive HTML rendering
# =========================================================

def interactive_html_to_html(
    path: Path,
    ordination_marker_size: float = 16
) -> str:

    raw = path.read_text(
        encoding='utf-8',
        errors='replace'
    )

    lower = raw.lower()

    body_start = lower.find(
        '<body'
    )

    body_end = lower.rfind(
        '</body>'
    )

    if body_start != -1:

        body_start = raw.find(
            '>',
            body_start
        )

    if (
        body_start != -1
        and body_end != -1
        and body_end > body_start
    ):

        embedded = raw[
            body_start + 1:
            body_end
        ]

    else:

        embedded = raw

    ordination_files = {
        'NMDS_bray.html',
        'PCoA_unweighted_unifrac.html',
    }

    marker_script = ''

    if path.name in ordination_files:

        marker_size = float(
            ordination_marker_size
        )

        marker_script = f"""
<script>
(function() {{

    const script = document.currentScript;

    const card = script
        ? script.closest('.interactive-card')
        : null;

    function resizeMarkers() {{

        if (
            !card
            || typeof Plotly === 'undefined'
        ) {{

            return false;

        }}

        const plots = card.querySelectorAll(
            '.plotly-graph-div'
        );

        if (!plots.length) {{

            return false;

        }}

        plots.forEach(plot => {{

            try {{

                Plotly.restyle(
                    plot,
                    {{
                        'marker.size': {marker_size}
                    }}
                );

            }} catch (error) {{

                console.warn(
                    'Could not resize Plotly markers:',
                    error
                );

            }}

        }});

        return true;

    }}

    if (!resizeMarkers()) {{

        let attempts = 0;

        const timer = setInterval(
            () => {{

                attempts += 1;

                if (
                    resizeMarkers()
                    || attempts >= 40
                ) {{

                    clearInterval(
                        timer
                    );

                }}

            }},
            250
        );

    }}

}})();
</script>
"""

    return (
        '<div class="figure-card interactive-card">'
        f'<h3>{html.escape(natural_title(path))}</h3>'
        '<div class="interactive-wrap">'
        f'{embedded}'
        f'{marker_script}'
        '</div>'
        f'<div class="file-name">{html.escape(path.name)}</div>'
        '</div>'
    )


# =========================================================
# Krona rendering
# =========================================================

def krona_html_to_html(
    path: Path
) -> str:

    raw = path.read_text(
        encoding='utf-8',
        errors='replace'
    )

    srcdoc = html.escape(
        raw,
        quote=True
    )

    return (
        '<div class="figure-card krona-card">'
        '<h3>Krona — Class, Order and Family</h3>'
        '<div class="krona-wrap">'
        f'<iframe '
        f'class="krona-frame" '
        f'srcdoc="{srcdoc}" '
        f'loading="lazy" '
        f'allowfullscreen>'
        '</iframe>'
        '</div>'
        f'<div class="file-name">{html.escape(path.name)}</div>'
        '</div>'
    )


# =========================================================
# TSV rendering
# =========================================================

def tsv_to_html(
    path: Path
) -> str:

    with path.open(
        'r',
        encoding='utf-8',
        errors='replace',
        newline=''
    ) as handle:

        rows = list(
            csv.reader(
                handle,
                delimiter='\t'
            )
        )

    if not rows:

        body = (
            '<p class="empty">'
            'Empty table.'
            '</p>'
        )

    else:

        header = rows[0]

        data = rows[1:]

        head_html = ''.join(
            f'<th>{html.escape(cell)}</th>'
            for cell in header
        )

        row_html = []

        for row in data:

            cells = row + [''] * max(
                0,
                len(header) - len(row)
            )

            row_html.append(
                '<tr>'
                + ''.join(
                    f'<td>{html.escape(cell)}</td>'
                    for cell
                    in cells[:len(header)]
                )
                + '</tr>'
            )

        body = (
            '<div class="table-wrap">'
            '<table>'
            f'<thead><tr>{head_html}</tr></thead>'
            f'<tbody>{"".join(row_html)}</tbody>'
            '</table>'
            '</div>'
        )

    return (
        '<div class="table-card">'
        f'<h3>{html.escape(natural_title(path))}</h3>'
        f'{body}'
        f'<div class="file-name">{html.escape(path.name)}</div>'
        '</div>'
    )


# =========================================================
# XLSX rendering for PIMBA-Curate
# =========================================================

def xlsx_to_html(
    path: Path
) -> str:

    workbook = pd.read_excel(
        path,
        sheet_name=None,
        dtype=object
    )

    if not workbook:

        body = (
            '<p class="empty">'
            'No worksheets were found in the Excel file.'
            '</p>'
        )

    else:

        sheet_blocks = []

        for sheet_index, (
            sheet_name,
            df
        ) in enumerate(
            workbook.items()
        ):

            # Replace missing values with empty strings.
            # Using .where() avoids the pandas FutureWarning
            # generated by .fillna('') on object dtype columns.
            df = df.where(
                pd.notna(df),
                ''
            )

            if (
                df.empty
                and len(df.columns) == 0
            ):

                sheet_blocks.append(
                    '<div class="excel-sheet">'
                    f'<h4>'
                    f'Worksheet: '
                    f'{html.escape(str(sheet_name))}'
                    f'</h4>'
                    '<p class="empty">'
                    'Empty worksheet.'
                    '</p>'
                    '</div>'
                )

                continue

            columns = [
                str(col)
                for col in df.columns
            ]

            column_lookup = {
                col.lower(): idx
                for idx, col
                in enumerate(columns)
            }

            pid_index = (
                column_lookup.get(
                    'pid'
                )
            )

            species_index = (
                column_lookup.get(
                    'species'
                )
            )

            if species_index is not None:

                sample_indices = list(
                    range(
                        species_index + 1,
                        len(columns)
                    )
                )

            else:

                sample_indices = [
                    idx
                    for idx, col
                    in enumerate(columns)
                    if re.match(
                        r'^(SRR|ERR|DRR|ITV)',
                        col,
                        flags=re.IGNORECASE
                    )
                ]

            table_id = (
                f'curate-table-{sheet_index}'
            )

            pid_filter_id = (
                f'pid-filter-{sheet_index}'
            )

            pid_operator_id = (
                f'pid-operator-{sheet_index}'
            )

            sample_filter_id = (
                f'sample-filter-{sheet_index}'
            )

            abundance_filter_id = (
                f'abundance-filter-{sheet_index}'
            )

            row_count_id = (
                f'row-count-{sheet_index}'
            )

            controls = []

            if pid_index is not None:

                controls.append(
                    '<div class="filter-control">'
                    '<label>PID filter</label>'
                    f'<select '
                    f'id="{pid_operator_id}" '
                    f'onchange="filterCurateTable('
                    f'{sheet_index}'
                    f')">'
                    '<option value="ge">≥</option>'
                    '<option value="gt">&gt;</option>'
                    '<option value="le">≤</option>'
                    '<option value="lt">&lt;</option>'
                    '<option value="eq">=</option>'
                    '</select>'
                    f'<input '
                    f'id="{pid_filter_id}" '
                    f'type="number" '
                    f'step="any" '
                    f'placeholder="e.g. 0.95" '
                    f'oninput="filterCurateTable('
                    f'{sheet_index}'
                    f')">'
                    '</div>'
                )

            if sample_indices:

                sample_options = [
                    '<option value="">'
                    'All samples'
                    '</option>'
                ]

                for idx in sample_indices:

                    sample_options.append(
                        f'<option value="{idx}">'
                        f'{html.escape(columns[idx])}'
                        f'</option>'
                    )

                controls.append(
                    '<div class="filter-control">'
                    '<label>Sample</label>'
                    f'<select '
                    f'id="{sample_filter_id}" '
                    f'onchange="filterCurateTable('
                    f'{sheet_index}'
                    f')">'
                    + ''.join(sample_options)
                    + '</select>'
                    '</div>'
                )

                controls.append(
                    '<div class="filter-control checkbox-control">'
                    f'<label>'
                    f'<input '
                    f'id="{abundance_filter_id}" '
                    f'type="checkbox" '
                    f'onchange="filterCurateTable('
                    f'{sheet_index}'
                    f')"> '
                    f'Abundance &gt; 0'
                    f'</label>'
                    '</div>'
                )

            controls.append(
                f'<div '
                f'class="filter-row-count" '
                f'id="{row_count_id}">'
                f'</div>'
            )

            head_html = ''.join(
                f'<th data-col-index="{idx}">'
                f'{html.escape(col)}'
                f'</th>'
                for idx, col
                in enumerate(columns)
            )

            body_rows = []

            for _, row in df.iterrows():

                cells = []

                for idx, col in enumerate(
                    df.columns
                ):

                    value = row[col]

                    display = (
                        ''
                        if value == ''
                        else str(value)
                    )

                    cells.append(
                        f'<td '
                        f'data-col-index="{idx}">'
                        f'{html.escape(display)}'
                        f'</td>'
                    )

                body_rows.append(
                    '<tr>'
                    + ''.join(cells)
                    + '</tr>'
                )

            table_html = (
                '<div class="curate-filter-bar">'
                + ''.join(controls)
                + '</div>'
                '<div class="table-wrap">'
                f'<table '
                f'id="{table_id}" '
                f'class="excel-table curate-filter-table" '
                f'data-pid-index="'
                f'{"" if pid_index is None else pid_index}'
                f'">'
                f'<thead><tr>{head_html}</tr></thead>'
                f'<tbody>{"".join(body_rows)}</tbody>'
                '</table>'
                '</div>'
            )

            sheet_blocks.append(
                '<div class="excel-sheet">'
                f'<h4>'
                f'Worksheet: '
                f'{html.escape(str(sheet_name))}'
                f'</h4>'
                f'{table_html}'
                '</div>'
            )

        body = ''.join(
            sheet_blocks
        )

    return (
        '<div class="table-card curate-table-card">'
        '<h3>'
        'Taxonomic assignments to validate'
        '</h3>'
        '<p class="table-help">'
        'Filter by PID, or select a sample and enable '
        '<strong>Abundance &gt; 0</strong> '
        'to show only taxa detected in that sample.'
        '</p>'
        f'{body}'
        f'<div class="file-name">{html.escape(str(path))}</div>'
        '</div>'
    )


# =========================================================
# Render supported files
# =========================================================

def render_file(
    path: Path,
    ordination_marker_size: float = 16
) -> str:

    suffix = (
        path.suffix.lower()
    )

    if suffix == '.svg':

        return svg_to_html(
            path
        )

    if (
        suffix == '.html'
        and path.name
        == 'krona_class_order_family.html'
    ):

        return krona_html_to_html(
            path
        )

    if suffix == '.html':

        return interactive_html_to_html(
            path,
            ordination_marker_size=
                ordination_marker_size
        )

    if suffix in {
        '.tsv',
        '.txt'
    }:

        return tsv_to_html(
            path
        )

    if suffix == '.xlsx':

        return xlsx_to_html(
            path
        )

    return (
        '<div class="file-card">'
        f'<h3>{html.escape(path.name)}</h3>'
        '<p>'
        'This file type is listed but not embedded.'
        '</p>'
        '</div>'
    )


# =========================================================
# Match output files
# =========================================================

def matching_files(
    directory: Path,
    patterns: list[str]
) -> list[Path]:

    found = []

    seen = set()

    for pattern in patterns:

        for path in sorted(
            directory.glob(pattern)
        ):

            if (
                not path.is_file()
                or path in seen
            ):

                continue

            if (
                path.suffix.lower() == '.svg'
                and path.with_suffix(
                    '.html'
                ).is_file()
            ):

                continue

            found.append(
                path
            )

            seen.add(
                path
            )

    return found


# =========================================================
# Build report section
# =========================================================

def section_html(
    section_id: str,
    title: str,
    description: str,
    files: list[Path],
    ordination_marker_size: float = 16
) -> str:

    if not files:

        content = (
            '<p class="empty">'
            'No output files were found for this analysis.'
            '</p>'
        )

    else:

        content = ''.join(
            render_file(
                path,
                ordination_marker_size=
                    ordination_marker_size
            )
            for path in files
        )

    return (
        f'<section '
        f'id="{html.escape(section_id)}">'
        f'<h2>{html.escape(title)}</h2>'
        f'<p class="section-description">'
        f'{html.escape(description)}'
        f'</p>'
        f'{content}'
        '</section>'
    )


# =========================================================
# Build report
# =========================================================

def build_report(
    results_dir: Path,
    output: Path,
    title: str,
    include_curate: bool = False,
    include_rf: bool = False,
    curate_xlsx: Path | None = None,
    rf_dir: Path | None = None,
    ordination_marker_size: float = 16
) -> None:

    plots_dir = (
        results_dir
        / "plots"
    )

    if not plots_dir.is_dir():

        raise FileNotFoundError(
            f'PIMBA plots directory not found: '
            f'{plots_dir}'
        )

    # Resolve optional components before creating report
    if include_curate:

        curate_xlsx = resolve_curate_xlsx(
            results_dir=results_dir,
            curate_xlsx=curate_xlsx
        )

    if include_rf:

        rf_dir = resolve_rf_dir(
            results_dir=results_dir,
            rf_dir=rf_dir
        )

    sections = []

    included_modules = [
        'PIMBA Plot'
    ]

    # =====================================================
    # Rarefaction
    # =====================================================

    sections.append(
        (
            'rarefaction',
            'Rarefaction',
            (
                'Sequencing-depth and '
                'observed-richness analysis.'
            ),
            matching_files(
                plots_dir,
                [
                    "rarefaction_curve.*"
                ]
            ),
        )
    )

    # =====================================================
    # Alpha diversity
    # =====================================================

    sections.append(
        (
            'alpha-diversity',
            'Alpha diversity',
            'Within-sample diversity results.',
            matching_files(
                plots_dir,
                [
                    "alpha_diversity_*"
                ]
            ),
        )
    )

    # =====================================================
    # Ordination
    # =====================================================

    sections.append(
        (
            'ordination',
            'Ordination',
            (
                'Community-level ordination analyses, '
                'including Bray-Curtis NMDS and '
                'UniFrac PCoA.'
            ),
            matching_files(
                plots_dir,
                [
                    "NMDS_*",
                    "PCoA_*",
                ]
            ),
        )
    )

    # =====================================================
    # Clustering
    # =====================================================

    sections.append(
        (
            'clustering',
            'Clustering',
            (
                'Community clustering and '
                'bootstrap-supported dendrogram results.'
            ),
            matching_files(
                plots_dir,
                [
                    "cluster_*"
                ]
            ),
        )
    )

    # =====================================================
    # Taxonomy
    # =====================================================

    sections.append(
        (
            'taxonomy',
            'Taxonomic composition',
            (
                'Taxonomic abundance plots and pivot '
                'tables from phylum through genus, '
                'plus the interactive Krona visualization.'
            ),
            matching_files(
                plots_dir,
                [
                    "phylum_barplots.*",
                    "phylum_pivot_table.tsv",
                    "class_barplots.*",
                    "class_pivot_table.tsv",
                    "order_barplots.*",
                    "order_pivot_table.tsv",
                    "family_barplots.*",
                    "family_pivot_table.tsv",
                    "genus_barplots.*",
                    "genus_pivot_table.tsv",
                    "krona_class_order_family.html",
                ]
            ),
        )
    )

    # =====================================================
    # Optional PIMBA-Curate
    # =====================================================

    if include_curate:

        sections.append(
            (
                'curate-validation',
                (
                    'PIMBA-Curate — '
                    'Taxonomic assignments to validate'
                ),
                (
                    'Taxonomic assignments flagged '
                    'by PIMBA-Curate for manual validation.'
                ),
                [
                    curate_xlsx
                ],
            )
        )

        included_modules.append(
            'PIMBA-Curate'
        )

    # =====================================================
    # Optional Random Forest
    # =====================================================

    if include_rf:

        sections.append(
            (
                'rf-importance',
                (
                    'Random Forest — '
                    'feature importance'
                ),
                (
                    'Top ASVs selected by the '
                    'Random Forest classifier and '
                    'their taxonomic annotation.'
                ),
                matching_files(
                    rf_dir,
                    [
                        "rf_top*_asvs_*.svg",
                        "rf_top_asvs_taxonomy_importance.tsv",
                    ]
                ),
            )
        )

        sections.append(
            (
                'rf-heatmap',
                (
                    'Random Forest — '
                    'ASV heatmap'
                ),
                (
                    'Abundance patterns for the '
                    'most informative ASVs.'
                ),
                matching_files(
                    rf_dir,
                    [
                        "rf_heatmap_*.svg"
                    ]
                ),
            )
        )

        sections.append(
            (
                'rf-pca',
                'Random Forest — PCA',
                (
                    'Interactive 3D PCA using '
                    'all ASVs and the selected '
                    'top ASVs.'
                ),
                matching_files(
                    rf_dir,
                    [
                        "rf_pca_*.html"
                    ]
                ),
            )
        )

        sections.append(
            (
                'rf-mds',
                (
                    'Random Forest — '
                    'Bray-Curtis MDS'
                ),
                (
                    'Interactive Bray-Curtis MDS '
                    'using all ASVs and the selected '
                    'top ASVs.'
                ),
                matching_files(
                    rf_dir,
                    [
                        "rf_mds_*.html"
                    ]
                ),
            )
        )

        included_modules.append(
            'Random Forest'
        )

    # =====================================================
    # Navigation
    # =====================================================

    nav_items = ''.join(
        (
            f'<a href="#{section_id}">'
            f'{html.escape(section_title)}'
            f'</a>'
        )
        for (
            section_id,
            section_title,
            _,
            _
        ) in sections
    )

    rendered_sections = ''.join(
        section_html(
            section_id,
            section_title,
            description,
            files,
            ordination_marker_size=
                ordination_marker_size
        )
        for (
            section_id,
            section_title,
            description,
            files
        ) in sections
    )

    report_subtitle = (
        'Results: '
        + ', '.join(
            included_modules
        )
    )

    # =====================================================
    # HTML document
    # =====================================================

    document = f'''<!DOCTYPE html>
<html lang="en">

<head>

<meta charset="utf-8">

<meta
    name="viewport"
    content="width=device-width, initial-scale=1"
>

<title>{html.escape(title)}</title>

<style>

    :root {{
        --bg: #f5f7fa;
        --panel: #ffffff;
        --text: #1f2937;
        --muted: #667085;
        --border: #d9e0e7;
        --accent: #176b56;
        --accent-dark: #0d4b3c;
    }}

    * {{
        box-sizing: border-box;
    }}

    html {{
        scroll-behavior: smooth;
    }}

    body {{
        margin: 0;
        font-family: Arial, Helvetica, sans-serif;
        background: var(--bg);
        color: var(--text);
        line-height: 1.5;
    }}

    header {{
        background:
            linear-gradient(
                135deg,
                var(--accent-dark),
                var(--accent)
            );
        color: white;
        padding: 34px 5vw;
    }}

    header h1 {{
        margin: 0 0 8px 0;
        font-size: 2rem;
    }}

    header p {{
        margin: 0;
        opacity: 0.9;
    }}

    nav {{
        position: sticky;
        top: 0;
        z-index: 20;
        display: flex;
        gap: 8px;
        overflow-x: auto;
        padding: 10px 5vw;
        background:
            rgba(
                255,
                255,
                255,
                0.96
            );
        border-bottom:
            1px solid var(--border);
        backdrop-filter: blur(8px);
    }}

    nav a {{
        flex: 0 0 auto;
        color: var(--accent-dark);
        text-decoration: none;
        font-weight: 600;
        font-size: 0.88rem;
        padding: 7px 11px;
        border-radius: 999px;
        background: #edf7f3;
    }}

    main {{
        width: min(1500px, 94vw);
        margin: 28px auto 60px auto;
    }}

    section {{
        background: var(--panel);
        border: 1px solid var(--border);
        border-radius: 14px;
        padding: 26px;
        margin-bottom: 28px;
        box-shadow:
            0 4px 16px
            rgba(
                16,
                24,
                40,
                0.04
            );
        scroll-margin-top: 70px;
    }}

    section h2 {{
        margin: 0 0 4px 0;
        color: var(--accent-dark);
        border-bottom:
            2px solid #e6f2ee;
        padding-bottom: 10px;
    }}

    .section-description {{
        color: var(--muted);
        margin: 10px 0 22px 0;
    }}

    .figure-card,
    .table-card,
    .file-card {{
        margin: 22px 0 30px 0;
        padding-top: 4px;
    }}

    .figure-card h3,
    .table-card h3,
    .file-card h3 {{
        margin-bottom: 12px;
    }}

    .excel-sheet {{
        margin: 18px 0 28px 0;
    }}

    .excel-sheet h4 {{
        margin: 0 0 10px 0;
        color: var(--accent-dark);
        font-size: 1rem;
    }}

    .svg-wrap {{
        width: 100%;
        overflow: auto;
        text-align: center;
        border: 1px solid var(--border);
        border-radius: 10px;
        background: white;
        padding: 12px;
    }}

    .svg-wrap svg {{
        max-width: 100%;
        height: auto;
    }}

    .interactive-wrap {{
        width: 100%;
        min-height: 650px;
        overflow: auto;
        border: 1px solid var(--border);
        border-radius: 10px;
        background: white;
        padding: 8px;
    }}

    .interactive-wrap .plotly-graph-div {{
        margin: 0 auto;
        max-width: 100%;
    }}

    .krona-wrap {{
        width: 100%;
        height: 820px;
        border: 1px solid var(--border);
        border-radius: 10px;
        overflow: hidden;
        background: white;
    }}

    .krona-frame {{
        width: 100%;
        height: 100%;
        border: 0;
    }}

    .table-wrap {{
        overflow: auto;
        max-height: 700px;
        border: 1px solid var(--border);
        border-radius: 10px;
    }}

    table {{
        border-collapse: collapse;
        width: 100%;
        font-size: 0.86rem;
        background: white;
    }}

    th,
    td {{
        padding: 8px 10px;
        border-bottom:
            1px solid #eaecf0;
        border-right:
            1px solid #f0f2f5;
        text-align: left;
        white-space: nowrap;
    }}

    th {{
        position: sticky;
        top: 0;
        background: #eff6f3;
        color: var(--accent-dark);
        z-index: 2;
    }}

    tbody tr:nth-child(even) {{
        background: #fafbfc;
    }}

    .curate-filter-bar {{
        display: flex;
        flex-wrap: wrap;
        align-items: end;
        gap: 12px;
        margin: 0 0 14px 0;
        padding: 12px;
        border: 1px solid var(--border);
        border-radius: 10px;
        background: #f8fbfa;
    }}

    .filter-control {{
        display: flex;
        align-items: center;
        gap: 7px;
    }}

    .filter-control label {{
        font-weight: 600;
        color: var(--accent-dark);
        white-space: nowrap;
    }}

    .filter-control input[type="number"],
    .filter-control select {{
        padding: 7px 9px;
        border: 1px solid var(--border);
        border-radius: 7px;
        background: white;
        color: var(--text);
        font: inherit;
    }}

    .checkbox-control {{
        padding-bottom: 4px;
    }}

    .filter-row-count {{
        margin-left: auto;
        padding: 7px 0;
        color: var(--muted);
        font-size: 0.86rem;
        white-space: nowrap;
    }}

    .table-help {{
        color: var(--muted);
        margin: -5px 0 14px 0;
    }}

    .file-name {{
        margin-top: 7px;
        color: var(--muted);
        font-family: monospace;
        font-size: 0.78rem;
    }}

    .empty {{
        color: var(--muted);
        font-style: italic;
    }}

    footer {{
        width: min(1500px, 94vw);
        margin: 0 auto 40px auto;
        color: var(--muted);
        font-size: 0.85rem;
        text-align: center;
    }}

    @media (max-width: 800px) {{

        header {{
            padding: 25px 4vw;
        }}

        main {{
            width: 96vw;
        }}

        section {{
            padding: 16px;
        }}

        .interactive-wrap {{
            min-height: 540px;
        }}

    }}

</style>

</head>

<body>

<header>

    <h1>
        {html.escape(title)}
    </h1>

    <p>
        {html.escape(report_subtitle)}
    </p>

</header>

<nav>

    {nav_items}

</nav>

<main>

    {rendered_sections}

</main>

<footer>

    Generated from
    {html.escape(str(results_dir))}

</footer>

<script>

function filterCurateTable(sheetIndex) {{

    const table =
        document.getElementById(
            `curate-table-${{sheetIndex}}`
        );

    if (!table) return;

    const pidInput =
        document.getElementById(
            `pid-filter-${{sheetIndex}}`
        );

    const pidOperator =
        document.getElementById(
            `pid-operator-${{sheetIndex}}`
        );

    const sampleSelect =
        document.getElementById(
            `sample-filter-${{sheetIndex}}`
        );

    const abundanceCheckbox =
        document.getElementById(
            `abundance-filter-${{sheetIndex}}`
        );

    const rowCount =
        document.getElementById(
            `row-count-${{sheetIndex}}`
        );

    const pidIndexRaw =
        table.dataset.pidIndex;

    const pidIndex =
        pidIndexRaw === ''
            ? null
            : Number(pidIndexRaw);

    const pidText =
        pidInput
            ? pidInput.value.trim()
            : '';

    const pidValue =
        pidText === ''
            ? null
            : Number(pidText);

    const operator =
        pidOperator
            ? pidOperator.value
            : 'ge';

    const sampleIndex =
        (
            sampleSelect
            && sampleSelect.value !== ''
        )
            ? Number(
                sampleSelect.value
            )
            : null;

    const requirePositiveAbundance =
        Boolean(
            abundanceCheckbox
            && abundanceCheckbox.checked
            && sampleIndex !== null
        );

    let visible = 0;

    const rows =
        table.querySelectorAll(
            'tbody tr'
        );

    rows.forEach(row => {{

        const cells =
            row.children;

        let keep = true;

        if (
            pidValue !== null
            && pidIndex !== null
            && cells[pidIndex]
        ) {{

            const value =
                Number(
                    cells[
                        pidIndex
                    ].textContent.trim()
                );

            if (
                !Number.isFinite(value)
            ) {{

                keep = false;

            }} else if (
                operator === 'ge'
            ) {{

                keep =
                    value >= pidValue;

            }} else if (
                operator === 'gt'
            ) {{

                keep =
                    value > pidValue;

            }} else if (
                operator === 'le'
            ) {{

                keep =
                    value <= pidValue;

            }} else if (
                operator === 'lt'
            ) {{

                keep =
                    value < pidValue;

            }} else if (
                operator === 'eq'
            ) {{

                keep =
                    value === pidValue;

            }}

        }}

        if (
            keep
            && requirePositiveAbundance
            && cells[sampleIndex]
        ) {{

            const abundance =
                Number(
                    cells[
                        sampleIndex
                    ].textContent.trim()
                );

            keep =
                Number.isFinite(abundance)
                && abundance > 0;

        }}

        row.style.display =
            keep
                ? ''
                : 'none';

        if (keep) {{

            visible += 1;

        }}

    }});

    if (rowCount) {{

        rowCount.textContent =
            `${{visible}} of ${{rows.length}} rows shown`;

    }}

}}

document.addEventListener(
    'DOMContentLoaded',
    () => {{

        document
            .querySelectorAll(
                '.curate-filter-table'
            )
            .forEach(table => {{

                const match =
                    table.id.match(
                        /curate-table-(\\d+)/
                    );

                if (match) {{

                    filterCurateTable(
                        Number(
                            match[1]
                        )
                    );

                }}

            }});

    }}
);

</script>

</body>

</html>
'''

    output.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    output.write_text(
        document,
        encoding='utf-8'
    )

    print(
        f'Report created: {output}'
    )

    print(
        'Included '
        f'{sum(len(files) for _, _, _, files in sections)} '
        'result file(s).'
    )


# =========================================================
# Main
# =========================================================

def main() -> None:

    args = parse_args()

    results_dir = (
        args.results_dir.resolve()
    )

    if not results_dir.is_dir():

        raise FileNotFoundError(
            f'Results directory not found: '
            f'{results_dir}'
        )

    output = (
        args.output
    )

    if output is None:

        output = (
            results_dir
            / "PIMBA_report.html"
        )

    else:

        output = (
            output.resolve()
        )

    curate_xlsx = (
        args.curate_xlsx.resolve()
        if args.curate_xlsx is not None
        else None
    )

    rf_dir = (
        args.rf_dir.resolve()
        if args.rf_dir is not None
        else None
    )

    build_report(
        results_dir=results_dir,
        output=output,
        title=args.title,
        include_curate=
            args.include_curate,
        include_rf=
            args.include_rf,
        curate_xlsx=
            curate_xlsx,
        rf_dir=
            rf_dir,
        ordination_marker_size=
            args.ordination_marker_size,
    )


if __name__ == '__main__':

    main()