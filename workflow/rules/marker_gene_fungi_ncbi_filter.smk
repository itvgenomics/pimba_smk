rule run_qiimepipe_create_tax_assignment:
    input:
        blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log'),
        txt_table = os.path.join(output_dir, file_name_raw + '_pre_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_pre_asv_table.txt')
    output:
        tax_assignments = os.path.join(output_dir, 'diversity_by_sample_ncbi', 'ncbi_tax_assignments.txt')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_create_tax_assignment', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_create_tax_assignment.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        mkdir -p {output_dir}/diversity_by_sample_ncbi
        cd {output_dir}/diversity_by_sample_ncbi

        echo "Running the QiimePipe Container - create_otuTaxAssignment.py: "
        python3.6 /qiimepipe/create_otuTaxAssignment.py {input.blast_log} \
        {input.txt_table} {output.tax_assignments} > {log} 2>&1
        """
rule run_qiimepipe_filter:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta'),
        tax_assignments = os.path.join(output_dir, 'diversity_by_sample_ncbi', 'ncbi_tax_assignments.txt')
    output:
        filtered_fasta = os.path.join(output_dir, 'k__Fungi.fasta')
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_filter', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiimepipe_filter.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        cd {output_dir}
        echo "Running the QiimePipe Container - filterOTUs.py: "
        python3.6 /qiimepipe/filterOTUs.py {input.fasta_file} {input.tax_assignments} > {log} 2>&1
        """