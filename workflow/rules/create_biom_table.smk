rule convert_otu_table_to_biom:
    input:
        table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt'),
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    output:
        biom = os.path.join(output_dir, file_name_raw + '_otu_table.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.biom')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'convert_otu_table_to_biom', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_convert_otu_table_to_biom.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container - convert: "
        grep -F -f <(awk 'NR > 1 {{print $1}}' {input.table}) {input.tax_assignments} > {input.tax_assignments}.temp
        mv {input.tax_assignments}.temp {input.tax_assignments}
        biom convert -i {input.table} -o {output.biom} --to-hdf5 --table-type="OTU table" &> {log}
        """

rule add_metadata_to_biom_table:
    input:
        biom = os.path.join(output_dir, file_name_raw + '_otu_table.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.biom'),
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    output:
        biom_with_metadata = os.path.join(output_dir, file_name_raw + '_otu_table_tax.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_tax.biom')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'add_metadata_to_biom_table', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_add_metadata_to_biom_table.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container"
        biom convert -i {input.biom} -o temp_table.tsv --to-tsv --header-key OTUID
        tail -n +2 temp_table.tsv | cut -f1 | sort > all_otus.txt
        cut -f1 {input.tax_assignments} | sort > assigned_otus.txt
        comm -23 all_otus.txt assigned_otus.txt | awk '{{print $1"\\tmissing\\tmissing"}}' > missing_assignments.txt
        cat {input.tax_assignments} missing_assignments.txt > tax_assignments_padded.txt
        biom add-metadata -i {input.biom} -o {output.biom_with_metadata} \
            --observation-metadata-fp tax_assignments_padded.txt \
            --observation-header OTUID,taxonomy,confidence \
            --sc-separated taxonomy --float-fields confidence &> {log}
        rm -f temp_table.tsv all_otus.txt assigned_otus.txt missing_assignments.txt tax_assignments_padded.txt
        """


rule summarize_biom_table:
    input:
        biom_with_metadata = os.path.join(output_dir, file_name_raw + '_otu_table_tax.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_tax.biom')
    output:
        summary = os.path.join(output_dir, file_name_raw + '_biom_table')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'summarize_biom_table', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_summarize_biom_table.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container - summarize-table: "
        biom summarize-table -i {input.biom_with_metadata} -o {output.summary} &> {log}
        """