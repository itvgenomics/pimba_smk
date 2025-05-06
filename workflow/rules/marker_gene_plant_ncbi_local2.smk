rule run_blastn_plants:
    input:
        filtered_fasta = os.path.join(output_dir, file_name_raw + '_otus_filtered.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_filtered.fasta')
    output:
        blast_plants_log = os.path.join(output_dir, file_name_raw + '_blast_plants.log')
    params:
        similarity_int_asg = similarity_int_asg,
        coverage_int = coverage_int,
        hits_per_subject = config['hits_per_subject'],
        evalue = config['e_value'],
        NCBI_DB = NCBI_DB,
        threads = config['num_threads'],
        db_type = config['db_type']
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_blastn_plants.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_blast:latest"
    shell:
        """
        echo "Running the BLAST Container - blastn: "
        cd {output_dir}
        blastn -query {input.filtered_fasta} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
            {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
            -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
            -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_plants_log}
        """