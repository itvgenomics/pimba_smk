configfile: "config.yaml"

import os

# ----------------------------------------------------------
# Prepare paths
# ----------------------------------------------------------
file_name_raw = config["outputprepare"]
strategy = config["strategy"]

current_path = os.getcwd()

input_dir = os.path.join(current_path, "results", "01-run", config["outputrun"])

db_format = os.environ["PIMBA_DB_FORMAT"]

db_label = (
    "BOLD" if db_format == "COI-BOLD-DB" else
    "RDP" if db_format == "16S-RDP-DB" else
    "SILVA" if db_format == "16S-SILVA-DB" else
    "NCBI" if db_format == "NCBI-DB" else
    "UNITE" if db_format == "ITS-FUNGI-UNITE-DB" else
    "CUSTOM"
)

output_dir = os.path.join(current_path, "results", "03-curate", config["outputrun"] + "_" + db_label + "_Curate")
if not os.path.exists(os.path.dirname(output_dir)):
    os.makedirs(os.path.dirname(output_dir))

rule r_curation:
    input:
        otu_asv_file = os.path.join(input_dir, file_name_raw + "_otu_table.txt" if strategy == "otu" else file_name_raw + "_asv_table.txt"),
        fasta_file = os.path.join(input_dir, file_name_raw + "_otus.fasta" if strategy == "otu" else file_name_raw + "_asvs.fasta"),
        hits_file = os.path.join(input_dir, file_name_raw + "_blast_ncbi.log")

    output:
        done = os.path.join(output_dir, "r_curation_" + ("multiSeq" if config["r_curation"]["mode"] == "multi" else "singleSeq") + ".done")

    log:
        os.path.join(output_dir, "logs", "r_curation_" + ("multiSeq" if config["r_curation"]["mode"] == "multi" else "singleSeq") + ".log")

    benchmark:
        os.path.join(output_dir, "benchmark", "r_curation_" + ("multiSeq" if config["r_curation"]["mode"] == "multi" else "singleSeq") + ".txt")

    singularity:
        "r_curation_v5.sif"

    params:
        bold_ref_file = os.environ.get("PIMBA_BOLD_REF_FILE", "NA"),
        rdp_ref_file = os.environ.get("PIMBA_RDP_REF_FILE", "NA"),
        silva_ref_file = os.environ.get("PIMBA_SILVA_REF_FILE", "NA"),
        unite_ref_file = os.environ.get("PIMBA_UNITE_REF_FILE", "NA"),
        ncbi_db_dir = os.environ.get("PIMBA_NCBI_DB_DIR", "NA"),
        ref_tax_file = ("NA" if str(config["r_curation"].get("ref_tax_file", "")).strip() == "" else config["r_curation"]["ref_tax_file"]),
        db_format = db_format,
        mode = config["r_curation"]["mode"]

    shell:
        r"""
        mkdir -p {output_dir}
        mkdir -p {output_dir}/logs
        mkdir -p {output_dir}/benchmark

        cd {output_dir}

        Rscript /app/main.R \
            {params.bold_ref_file} \
            {params.rdp_ref_file} \
            {params.silva_ref_file} \
            {params.unite_ref_file} \
            {input.otu_asv_file} \
            {input.fasta_file} \
            {input.hits_file} \
            {params.ref_tax_file} \
            {params.db_format} \
            {params.mode} \
            > {log} 2>&1

        touch {output.done}
        """