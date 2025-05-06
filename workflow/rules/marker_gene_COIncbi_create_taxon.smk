rule run_create_taxon_table:
    input:
        blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log'),
        txt_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    output:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    params:
        taxdump = config['taxdump']
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
        python3.6 /qiimepipe/createTaxonTable_singleFile_smk.py {input.blast_log} {input.txt_table} {params.taxdump}/rankedlineage.dmp {params.taxdump}/merged.dmp > {log} 2>&1
        cd ../
        mkdir -p {output_dir}/output
        mv *_otus_tax_assignments.txt {output.tax_assignments}
        """
