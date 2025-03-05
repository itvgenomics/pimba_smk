rule convert_otu_table_to_biom:
    input:
        table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    output:
        biom = os.path.join(output_dir, file_name_raw + '_otu_table.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.biom')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'convert_otu_table_to_biom', file_name_raw + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container - convert: "
        biom convert -i {input.table} -o {output.biom} \
        --table-type="OTU table" --to-json &> {log}
        """

rule add_metadata_to_biom_table:
    input:
        biom = os.path.join(output_dir, file_name_raw + '_otu_table.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.biom'),
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    output:
        biom_with_metadata = os.path.join(output_dir, file_name_raw + '_otu_table_tax.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_tax.biom')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'add_metadata_to_biom_table', file_name_raw + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container - add-metadata: "
        biom add-metadata -i {input.biom} -o {output.biom_with_metadata} \
        --observation-metadata-fp {input.tax_assignments} --observation-header OTUID,taxonomy,confidence \
        --sc-separated taxonomy --float-fields confidence &> {log}
        """

rule summarize_biom_table:
    input:
        biom_with_metadata = os.path.join(output_dir, file_name_raw + '_otu_table_tax.biom') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_tax.biom')
    output:
        summary = os.path.join(output_dir, file_name_raw + '_biom_table')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'summarize_biom_table', file_name_raw + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_biom:v2.1.10"
    shell:
        """
        echo "Running the Biom-Format Container - summarize-table: "
        biom summarize-table -i {input.biom_with_metadata} -o {output.summary} &> {log}
        """