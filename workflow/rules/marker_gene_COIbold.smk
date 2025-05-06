rule run_blast_assign_taxonomy:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
    output:
        log_file = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
    params:
        coibold_db = config['COI-BOLD-DB'],
        similarity_int_asg = similarity_int_asg,
        coverage_int = coverage_int,
        hits_per_subject = config['hits_per_subject'],
        e_value = config['e_value'],
        threads = config['num_threads']
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_blastn.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_blast:latest"
    shell:
        """
        echo "Running the BLAST Container - blastn: "
        blastn -query {input.fasta_file} -task megablast -db {params.coibold_db}/*.fasta \
        -perc_identity {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} \
        -max_hsps {params.hits_per_subject} -max_target_seqs {params.hits_per_subject} \
        -evalue {params.e_value} -parse_deflines -num_threads {params.threads} \
        -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.log_file}
        """
rule run_qiimepipe_create_taxon_table:
    input:
        blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log'),
        table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    output:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
        taxon_red_flagged = os.path.join(output_dir, 'output', file_name_raw + '_taxon_red_flagged.txt')
    params:
        coibold_db = config['COI-BOLD-DB']
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    log:
        os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_create_taxon_table', file_name_raw + '.log')
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_create_taxon_table.txt")
    shell:
        """
        mkdir -p {output_dir}/diversity_by_sample
        cd {output_dir}/diversity_by_sample

        echo "Running the QiimePipe Container - createTaxonTable_singleFile_flex.py: "
        python3.6 /qiimepipe/createTaxonTable_singleFile_flex.py {input.blast_log} \
        {input.table} {params.coibold_db}/*_tax.txt > {log} 2>&1

        cd {output_dir}
        mkdir -p output
        mv *_otus_tax_assignments.txt {output.tax_assignments}
        mv *_taxon_red_flagged.txt {output.taxon_red_flagged}
        """