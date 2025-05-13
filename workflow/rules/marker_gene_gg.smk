rule run_qiime_assign_taxonomy:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
    output:
        taxonomy = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    params:
        db_dir = config['16S-GREENGENES-DB'],
        similarity_int = similarity_int,
        similarity_assign = config['assign_similarity']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "qiime_assign_taxonomy", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiime_assign_taxonomy.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiime:latest"
    shell:
        """
        echo "Running the Qiime Container - assign_taxonomy.py: "
        assign_taxonomy.py -i {input.fasta_file} -o {output_dir}/output/ \
        -t {params.db_dir}/taxonomy/{params.similarity_int}_otu_taxonomy.txt \
        -r {params.db_dir}/rep_set/{params.similarity_int}_otus.fasta \
        --similarity={params.similarity_assign} > {log} 2>&1
        """
rule run_qiime_align_seqs:
    input:
        fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta'),
    output:
        rep_set_align = os.path.join(output_dir, 'rep_set_align', file_name_raw + '_otus_aligned.fasta') if strategy == "otu" else os.path.join(output_dir, 'rep_set_align', file_name_raw + '_asvs_aligned.fasta'),
    params:
        similarity_int = similarity_int,
        db_dir = config['16S-GREENGENES-DB']
    log:
        os.path.join(current_path, "results", "01-run", "logs", "qiime_align_seqs", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiime_align_seqs.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiime:latest"
    shell:
        """
        mkdir -p {output_dir}/rep_set_align
        echo "Running the Qiime Container - align_seqs.py: "
        align_seqs.py -i {input.fasta_file} -o {output_dir}/rep_set_align \
        -t {params.db_dir}/rep_set_aligned/{params.similarity_int}_otus.fasta > {log} 2>&1
        """
rule run_qiime_filter_alignment:
    input:
        rep_set_align = os.path.join(output_dir, 'rep_set_align', file_name_raw + '_otus_aligned.fasta') if strategy == "otu" else os.path.join(output_dir, 'rep_set_align', file_name_raw + '_asvs_aligned.fasta'),
    output:
        filtered_alignment = os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_otus_aligned_pfiltered.fasta') if strategy == "otu" else os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_asvs_aligned_pfiltered.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "qiime_filter_alignment", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiime_filter_alignment.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiime:latest"
    shell:
        """
        mkdir -p {output_dir}/filtered_alignment
        echo "Running the Qiime Container - filter_alignment.py: "
        filter_alignment.py -i {input.rep_set_align} -o {output_dir}/filtered_alignment > {log} 2>&1
        """
rule run_qiime_make_phylogeny:
    input:
        filtered_alignment = os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_otus_aligned_pfiltered.fasta') if strategy == "otu" else os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_asvs_aligned_pfiltered.fasta')
    output:
        phylogeny_tree = os.path.join(output_dir, 'rep_set.tre')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "qiime_make_phylogeny", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_qiime_make_phylogeny.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiime:latest"
    shell:
        """
        echo "Running the Qiime Container - make_phylogeny.py: "
        make_phylogeny.py -i {input.filtered_alignment} -o {output.phylogeny_tree} > {log} 2>&1
        """
rule run_create_abundance_file:
    input:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
        otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
    output:
        abundance_file_dir = directory(os.path.join(output_dir, 'diversity_by_sample'))
    log:
        os.path.join(current_path, "results", "01-run", "logs", "create_abundance_file", file_name_raw + ".log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "run_create_abundance_file.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        mkdir -p {output.abundance_file_dir}
        cd {output.abundance_file_dir}
        echo "Running the QiimePipe Container - createAbundanceFile.py: "
        python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
        """