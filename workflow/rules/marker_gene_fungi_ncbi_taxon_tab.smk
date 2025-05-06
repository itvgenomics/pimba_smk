rule run_vsearch_usearch_global:
    input:
        raw_reads = os.path.join(current_path, 'results', '00-prepare', config['outputprepare'] + '.fasta'),
        filtered_fasta = os.path.join(output_dir, 'k__Fungi.fasta')
    output:
        uc_file = os.path.join(output_dir, file_name_raw + '_map_fungi.uc'),
        fungi_table = os.path.join(output_dir, file_name_raw + '_otu_table_fungi.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_fungi.txt')
    params:
        similarity = config['otu_similarity']
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'vsearch_usearch_global', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_vsearch_usearch_global.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.29.1"
    shell:
        """
        cd {output_dir}

        echo "Running the VSEARCH Container - --usearch_global: "
        vsearch --usearch_global {input.raw_reads} --db {input.filtered_fasta} --strand both \
        --id {params.similarity} --uc {output.uc_file} --otutabout {output.fungi_table} > {log} 2>&1
        """
rule run_create_fungi_taxon_table:
    input:
        blast_fungi_log = os.path.join(output_dir, file_name_raw + '_blast_fungi.log'),
        fungi_table = os.path.join(output_dir, file_name_raw + '_otu_table_fungi.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_fungi.txt')
    output:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
        txt_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    params:
        taxdump = config['taxdump'],
        strategy = strategy
    log:
        os.path.join(current_path, "results", "01-run", "logs", "create_taxon_table", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_create_taxon_table.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        mkdir -p {output_dir}/diversity_by_sample
        cd {output_dir}/diversity_by_sample
        echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
        python3.6 /qiimepipe/createTaxonTable_singleFile_smk.py {input.blast_fungi_log} {input.fungi_table} {params.taxdump}/rankedlineage.dmp {params.taxdump}/merged.dmp > {log} 2>&1
        cd ../
        mkdir -p {output_dir}/output
        mv *_otus_tax_assignments.txt {output.tax_assignments}
        cp {input.fungi_table} {output.txt_table}
        """