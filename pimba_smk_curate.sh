#!/bin/bash

set -euo pipefail

# ============================================================================
# PIMBA 3.0 - r_curation
# ============================================================================

show_help() {
    echo " ___ ___ __  __ ___   _      ___ _   _ ___    _ _____ ___"
    echo "| _ \\_ _|  \\/  | _ ) /_\\    / __| | | | _ \\  /_\\_   _| __|"
    echo "|  _/| || |\\/| | _ \\/ _ \\  | (__| |_| |   / / _ \\| | | _|"
    echo "|_| |___|_|  |_|___/_/ \\_\\  \\___|\\___/|_|_\\/_/ \\_\\_| |___|"
    echo "                         v1.0"
    echo "Usage:"
    echo "  bash pimba_curate.sh -t|-threads <threads> -c|-config <config file> -d|-directory <working directory> [-u|-unlock]"
    echo ""
    echo "Options:"
    echo "  -h|-help       show help"
    echo "  -t|-threads    number of threads"
    echo "  -c|-config     config file"
    echo "  -d|-directory  working directory"
    echo "  -u|-unlock     unlock Snakemake working directory"
    echo ""
    echo "marker_gene examples:"
    echo "  marker_gene: 'COI-BOLD'"
    echo "  marker_gene: '16S-RDP'"
    echo "  marker_gene: '16S-GREENGENES'"
    echo "  marker_gene: '16S-SILVA'"
    echo "  marker_gene: 'ITS-FUNGI-UNITE'"
    echo "  marker_gene: 'NCBI'"
    echo "  marker_gene: '/path/to/custom_reference_folder'"
    echo ""
    exit 0
}

if [[ "${1:-}" == "-h" ]] || [[ "${1:-}" == "-help" ]]; then
    show_help
fi

ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -threads)
            ARGS+=("-t" "$2")
            shift 2
            ;;
        -config)
            ARGS+=("-c" "$2")
            shift 2
            ;;
        -directory)
            ARGS+=("-d" "$2")
            shift 2
            ;;
        -unlock)
            ARGS+=("-u")
            shift
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

set -- "${ARGS[@]}"

SETUNLOCK=""

while getopts "t:c:d:uh" opt; do
    case ${opt} in
        t ) THREADS=$OPTARG ;;
        c ) CONFIG=$OPTARG ;;
        d ) WORKDIR=$OPTARG ;;
        u ) SETUNLOCK="--unlock" ;;
        h ) show_help ;;
        \? ) echo "Invalid option: -$OPTARG" 1>&2; show_help ;;
    esac
done

if [ -z "${THREADS:-}" ] || [ -z "${CONFIG:-}" ] || [ -z "${WORKDIR:-}" ]; then
    echo "ERROR: Missing required arguments"
    show_help
fi

CONFIG=$(realpath "$CONFIG")
WORKDIR=$(realpath "$WORKDIR")

# --------------------------------
# Helper functions
# --------------------------------
clean_path() {
    echo "${1:-}" | tr -d '\r' | xargs
}

read_yaml_value() {
    local KEY="$1"
    local FILE="$2"

    grep -E "^${KEY}:" "$FILE" | \
        head -n 1 | \
        awk -F':' '{sub(/^[[:space:]]*/, "", $2); print $2}' | \
        tr -d '"' | \
        tr -d "'" | \
        tr -d '\r' | \
        xargs || true
}

require_directory() {
    local DIR="$1"
    local LABEL="$2"

    if [ -z "${DIR:-}" ]; then
        echo "ERROR: $LABEL is empty in config." >&2
        exit 1
    fi

    if [ ! -d "$DIR" ]; then
        echo "ERROR: $LABEL directory not found:" >&2
        echo "$DIR" >&2
        exit 1
    fi
}

detect_single_txt() {
    local DIR="$1"

    if [ ! -d "$DIR" ]; then
        echo "ERROR: Reference directory not found:" >&2
        echo "$DIR" >&2
        exit 1
    fi

    mapfile -t TXT_FILES < <(find "$DIR" -maxdepth 1 -type f -name "*.txt" | sort)

    if [ "${#TXT_FILES[@]}" -eq 0 ]; then
        echo "ERROR: No .txt taxonomy file found in reference directory:" >&2
        echo "$DIR" >&2
        exit 1
    fi

    if [ "${#TXT_FILES[@]}" -gt 1 ]; then
        echo "ERROR: More than one .txt taxonomy file found in reference directory:" >&2
        echo "$DIR" >&2
        echo "" >&2
        echo "Files found:" >&2
        printf '  %s\n' "${TXT_FILES[@]}" >&2
        echo "" >&2
        echo "Keep only one .txt reference taxonomy file in this folder." >&2
        exit 1
    fi

    echo "${TXT_FILES[0]}"
}

write_resolved_config() {
    local ORIGINAL_CONFIG="$1"
    local RESOLVED_CONFIG="$2"
    local REF_FILE="$3"

    grep -v '^ref_tax_file:' "$ORIGINAL_CONFIG" > "$RESOLVED_CONFIG"
    echo "ref_tax_file: '$REF_FILE'" >> "$RESOLVED_CONFIG"
}

# --------------------------------
# Read database directories from config
# These paths are user-defined in config.yaml.
# Each non-NCBI database directory must contain exactly one .txt taxonomy file.
# --------------------------------
BOLD=$(clean_path "$(read_yaml_value "COI-BOLD-DB" "$CONFIG")")
RDP=$(clean_path "$(read_yaml_value "16S-RDP-DB" "$CONFIG")")
SILVA=$(clean_path "$(read_yaml_value "16S-SILVA-DB" "$CONFIG")")
UNITE=$(clean_path "$(read_yaml_value "ITS-FUNGI-UNITE-DB" "$CONFIG")")
GREENGENES=$(clean_path "$(read_yaml_value "16S-GREENGENES-DB" "$CONFIG")")
NCBI_DB=$(clean_path "$(read_yaml_value "NCBI-DB" "$CONFIG")")

# Local NCBI taxizedb cache used by parsing_taxa.
# This is intentionally separate from NCBI-DB.
NCBI_TAXIZEDB=$(clean_path "$(read_yaml_value "ncbi_taxizedb" "$CONFIG")")

# --------------------------------
# Read marker_gene
# marker_gene controls database selection.
# Default databases use names such as COI-BOLD or 16S-RDP.
# Any non-default value is treated as a custom reference directory.
# --------------------------------
MARKER_GENE=$(clean_path "$(read_yaml_value "marker_gene" "$CONFIG")")

if [ -z "${MARKER_GENE:-}" ]; then
    echo "ERROR: marker_gene not found or empty in config."
    exit 1
fi

DB_FORMAT=""
DB_BIND=""
REF_FILE=""
CONFIG_ACTIVE="$CONFIG"

case "$MARKER_GENE" in

    "COI-BOLD")
        DB_FORMAT="COI-BOLD-DB"
        DB_DIR="$BOLD"
        require_directory "$DB_DIR" "COI-BOLD-DB"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"
        export PIMBA_BOLD_REF_FILE="$REF_FILE"
        ;;

    "16S-RDP")
        DB_FORMAT="16S-RDP-DB"
        DB_DIR="$RDP"
        require_directory "$DB_DIR" "16S-RDP-DB"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"
        export PIMBA_RDP_REF_FILE="$REF_FILE"
        ;;

    "16S-GREENGENES")
        DB_FORMAT="16S-GREENGENES-DB"
        DB_DIR="$GREENGENES"
        require_directory "$DB_DIR" "16S-GREENGENES-DB"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"
        export PIMBA_GREENGENES_REF_FILE="$REF_FILE"
        ;;

    "16S-SILVA")
        DB_FORMAT="16S-SILVA-DB"
        DB_DIR="$SILVA"
        require_directory "$DB_DIR" "16S-SILVA-DB"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"
        export PIMBA_SILVA_REF_FILE="$REF_FILE"
        ;;

    "ITS-FUNGI-UNITE")
        DB_FORMAT="ITS-FUNGI-UNITE-DB"
        DB_DIR="$UNITE"
        require_directory "$DB_DIR" "ITS-FUNGI-UNITE-DB"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"
        export PIMBA_UNITE_REF_FILE="$REF_FILE"
        ;;

    "NCBI")
        DB_FORMAT="NCBI-DB"

        if [ -n "${NCBI_DB:-}" ]; then
            require_directory "$NCBI_DB" "NCBI-DB"
            DB_BIND="-B $NCBI_DB:$NCBI_DB"
            export PIMBA_NCBI_DB_DIR="$NCBI_DB"
        else
            DB_BIND=""
            export PIMBA_NCBI_DB_DIR="NA"
        fi
        ;;

    *)
        DB_FORMAT="CUSTOM"
        DB_DIR=$(realpath "$MARKER_GENE")
        require_directory "$DB_DIR" "custom marker_gene"
        REF_FILE=$(detect_single_txt "$DB_DIR")
        DB_BIND="-B $DB_DIR:$DB_DIR"

        CONFIG_ACTIVE="$WORKDIR/.pimba_curate_config.resolved.yaml"
        write_resolved_config "$CONFIG" "$CONFIG_ACTIVE" "$REF_FILE"
        ;;
esac

export PIMBA_DB_FORMAT="$DB_FORMAT"

echo "Detected marker_gene: $MARKER_GENE"
echo "Internal database format: $DB_FORMAT"

if [ "$DB_FORMAT" != "NCBI-DB" ]; then
    echo "Reference directory: $DB_DIR"
    echo "Detected taxonomy file: $REF_FILE"
else
    echo "NCBI-DB path: ${NCBI_DB:-NA}"
fi

if [ "$DB_FORMAT" = "CUSTOM" ]; then
    echo "Resolved config file: $CONFIG_ACTIVE"
fi

# --------------------------------
# NCBI taxizedb cache bind
# Always mounted because parsing_taxa may validate species against NCBI.
# This uses ncbi_taxizedb, not NCBI-DB.
# --------------------------------
if [ -z "${NCBI_TAXIZEDB:-}" ]; then
    echo "ERROR: ncbi_taxizedb not found or empty in config."
    echo "This path is required because parsing_taxa validates species against the local NCBI taxizedb cache."
    exit 1
fi

require_directory "$NCBI_TAXIZEDB" "ncbi_taxizedb"

NCBI_BIND="-B $NCBI_TAXIZEDB:/ncbi_db"

# --------------------------------
# Run Snakemake
# --------------------------------
snakemake \
  --snakefile "workflow/Snakefile_curate" \
  --configfile "$CONFIG_ACTIVE" \
  --directory "$WORKDIR" \
  --use-singularity \
  --singularity-args "\
    -B $WORKDIR:$WORKDIR \
    -B $WORKDIR/workflow/scripts:/app \
    $DB_BIND \
    $NCBI_BIND" \
  --cores "$THREADS" \
  ${SETUNLOCK}