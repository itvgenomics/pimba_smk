rule run_blastn_fungi:
    input:
        filtered_fasta = os.path.join(output_dir, 'k__Fungi.fasta')
    output:
        blast_fungi_log = os.path.join(output_dir, file_name_raw + '_blast_fungi.log')
    params:
        similarity_int_asg = similarity_int_asg,
        coverage_int = coverage_int,
        hits_per_subject = config['hits_per_subject'],
        evalue = config['e_value']
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_blastn_fungi.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_blast:latest"
    shell:
        """
        echo "Running the BLAST Container - blastn: "
        cd {output_dir}
        blastn -query {input.filtered_fasta} -task megablast -db nt -remote -perc_identity \
            {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
            -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
            -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_fungi_log}
        """