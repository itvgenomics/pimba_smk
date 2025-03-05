if marker_gene == "16S-SILVA":
    rule run_qiime_assign_taxonomy:
        input:
            fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
        output:
            taxonomy = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
        params:
            db_dir = config['16S-SILVA-DB'],
            similarity_int = similarity_int,
            similarity_assign = config['assign_similarity']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_assign_taxonomy", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiime:latest"
        shell:
            """
            echo "Running the Qiime Container - assign_taxonomy.py: "
            assign_taxonomy.py -i {input.fasta_file} -o {output_dir}/output/ \
            -t {params.db_dir}/taxonomy/16S_only/{params.similarity_int}/taxonomy_7_levels.txt \
            -r {params.db_dir}/rep_set/rep_set_16S_only/{params.similarity_int}/silva_132_{params.similarity_int}_16S.fna \
            --similarity={params.similarity_assign} > {log} 2>&1
            """
    rule run_qiime_align_seqs:
        input:
            fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta'),
        output:
            rep_set_align = os.path.join(output_dir, 'rep_set_align', file_name_raw + '_otus_aligned.fasta') if strategy == "otu" else os.path.join(output_dir, 'rep_set_align', file_name_raw + '_asvs_aligned.fasta')
        params:
            similarity_int = similarity_int,
            db_dir = config['16S-SILVA-DB']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_align_seqs", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiime:latest"
        shell:
            """
            mkdir -p {output_dir}/rep_set_align
            echo "Running the Qiime Container - align_seqs.py: "
            align_seqs.py -i {input.fasta_file} -o {output_dir}/rep_set_align \
            -t {params.db_dir}/rep_set_aligned/{params.similarity_int}/{params.similarity_int}_alignment.fna > {log} 2>&1
            """
    rule run_qiime_filter_alignment:
        input:
            rep_set_align = os.path.join(output_dir, 'rep_set_align', file_name_raw + '_otus_aligned.fasta') if strategy == "otu" else os.path.join(output_dir, 'rep_set_align', file_name_raw + '_asvs_aligned.fasta')
        output:
            filtered_alignment = os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_otus_aligned_pfiltered.fasta') if strategy == "otu" else os.path.join(output_dir, 'filtered_alignment', file_name_raw + '_asvs_aligned_pfiltered.fasta')
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_filter_alignment", file_name_raw + ".log")
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
            reference_tree = os.path.join(output_dir, 'rep_set.tre')
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_make_phylogeny", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiime:latest"
        shell:
            """
            echo "Running the Qiime Container - make_phylogeny.py: "
            make_phylogeny.py -i {input.filtered_alignment} -o {output.reference_tree} > {log} 2>&1
            """
    rule run_create_abundance_file:
        input:
            tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
            otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
        output:
            abundance_file_dir = directory(os.path.join(output_dir, 'diversity_by_sample'))
        log:
            os.path.join(current_path, "results", "01-run", "logs", "create_abundance_file", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:latest"
        shell:
            """
            mkdir -p {output.abundance_file_dir}
            cd {output.abundance_file_dir}
            echo "Running the QiimePipe Container - createAbundanceFile.py: "
            python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
            """

elif marker_gene == "16S-GREENGENES":
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
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:latest"
        shell:
            """
            mkdir -p {output.abundance_file_dir}
            cd {output.abundance_file_dir}
            echo "Running the QiimePipe Container - createAbundanceFile.py: "
            python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
            """

elif marker_gene == "16S-RDP":
    rule run_qiime_assign_taxonomy:
        input:
            fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
        output:
            taxonomy = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
        params:
            db_dir = config['16S-RDP-DB'],
            similarity_assign = config['assign_similarity']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_assign_taxonomy", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiime:latest"
        shell:
            """
            echo "Running the Qiime Container - assign_taxonomy.py: "
            assign_taxonomy.py -i {input.fasta_file} -o {output_dir}/output/ \
            -t {params.db_dir}/*.txt \
            -r {params.db_dir}/*.fa \
            --similarity={params.similarity_assign} > {log} 2>&1
            """
    rule run_create_abundance_file:
        input:
            tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
            otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
        output:
            abundance_file_dir = directory(os.path.join(output_dir, 'diversity_by_sample'))
        log:
            os.path.join(current_path, "results", "01-run", "logs", "create_abundance_file", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:latest"
        shell:
            """
            mkdir -p {output.abundance_file_dir}
            cd {output.abundance_file_dir}
            echo "Running the QiimePipe Container - createAbundanceFile.py: "
            python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
            """

elif marker_gene == "16S-NCBI":
    if remote == 'yes':
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    else:
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
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

elif marker_gene == "COI-NCBI":
    if remote == 'yes':
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    else:
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
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

elif marker_gene == "COI-BOLD":
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

elif marker_gene == "ITS-PLANTS-NCBI":
    if remote == 'yes':
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    else:
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    rule run_qiimepipe_create_tax_assignment:
        input:
            blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log'),
            txt_table = os.path.join(output_dir, file_name_raw + '_pre_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_pre_asv_table.txt')
        output:
            tax_assignments = os.path.join(output_dir, 'diversity_by_sample_ncbi', 'ncbi_tax_assignments.txt')
        log:
            os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_create_tax_assignment', file_name_raw + '.log')
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
            filtered_fasta = os.path.join(output_dir, file_name_raw + '_otus_filtered.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_filtered.fasta')
        log:
            os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_filter', file_name_raw + '.log')
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:v2"
        shell:
            """
            cd {output_dir}
            echo "Running the QiimePipe Container - filterOTUs.py: "
            python3.6 /qiimepipe/filterOTUs.py {input.fasta_file} {input.tax_assignments} > {log} 2>&1
            """
    if remote == 'yes':
        rule run_blastn_plants:
            input:
                filtered_fasta = os.path.join(output_dir, file_name_raw + '_otus_filtered.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_filtered.fasta')
            output:
                blast_plants_log = os.path.join(output_dir, file_name_raw + '_blast_plants.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.filtered_fasta} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_plants_log}
                """
    else:
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
    rule run_vsearch_usearch_global:
        input:
            raw_reads = os.path.join(current_path, 'results', '00-prepare', config['outputprepare'] + '.fasta'),
            filtered_fasta = os.path.join(output_dir, file_name_raw + '_otus_filtered.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_filtered.fasta')
        output:
            uc_file = os.path.join(output_dir, file_name_raw + '_map_plants.uc'),
            plants_table = os.path.join(output_dir, file_name_raw + '_otu_table_plants.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_plants.txt')
        params:
            similarity = config['otu_similarity']
        log:
            os.path.join(current_path, 'results', '01-run', 'logs', 'vsearch_usearch_global', file_name_raw + '.log')
        singularity:
            "docker://itvdsbioinfo/pimba_vsearch:v2.29.1"
        shell:
            """
            cd {output_dir}

            echo "Running the VSEARCH Container - --usearch_global: "
            vsearch --usearch_global {input.raw_reads} --db {input.filtered_fasta} --strand both \
            --id {params.similarity} --uc {output.uc_file} --otutabout {output.plants_table} > {log} 2>&1
            """
    rule run_create_plants_taxon_table:
        input:
            blast_plants_log = os.path.join(output_dir, file_name_raw + '_blast_plants.log'),
            plants_table = os.path.join(output_dir, file_name_raw + '_otu_table_plants.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table_plants.txt'),
            filtered_fasta = os.path.join(output_dir, file_name_raw + '_otus_filtered.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_filtered.fasta')
        output:
            tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
            final_fasta = os.path.join(output_dir, file_name_raw + '_otus_plants.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs_plants.fasta'),
            txt_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
        params:
            taxdump = config['taxdump'],
            strategy = strategy
        log:
            os.path.join(current_path, "results", "01-run", "logs", "create_taxon_table", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:v2"
        shell:
            """
            mkdir -p {output_dir}/diversity_by_sample
            cd {output_dir}/diversity_by_sample
            echo "Running the QiimePipe Container - createTaxonTable_singleFile.py: "
            python3.6 /qiimepipe/createTaxonTable_singleFile_smk.py {input.blast_plants_log} {input.plants_table} {params.taxdump}/rankedlineage.dmp {params.taxdump}/merged.dmp > {log} 2>&1
            cd ../
            mkdir -p {output_dir}/output
            mv *_otus_tax_assignments.txt {output.tax_assignments}
            cp {input.filtered_fasta} {output.final_fasta}
            cp {input.plants_table} {output.txt_table}
            """

elif marker_gene == "ITS-FUNGI-UNITE":
    rule run_qiime_assign_taxonomy:
        input:
            fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
        output:
            taxonomy = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
        params:
            db_dir = config['ITS-FUNGI-UNITE-DB'],
            similarity_int = similarity_int,
            similarity_assign = config['assign_similarity']
        log:
            os.path.join(current_path, "results", "01-run", "logs", "qiime_assign_taxonomy", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiime:latest"
        shell:
            """
            echo "Running the Qiime Container - assign_taxonomy.py: "
            assign_taxonomy.py -i {input.fasta_file} -o {output_dir}/output/ \
            -t {params.db_dir}/*{params.similarity_int}*.txt \
            -r {params.db_dir}/*{params.similarity_int}*.fasta \
            --similarity={params.similarity_assign} > {log} 2>&1
            """
    rule run_create_abundance_file:
        input:
            tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt'),
            otu_table = os.path.join(output_dir, file_name_raw + '_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asv_table.txt')
        output:
            abundance_file_dir = directory(os.path.join(output_dir, 'diversity_by_sample'))
        log:
            os.path.join(current_path, "results", "01-run", "logs", "create_abundance_file", file_name_raw + ".log")
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:latest"
        shell:
            """
            mkdir -p {output.abundance_file_dir}
            cd {output.abundance_file_dir}
            echo "Running the QiimePipe Container - createAbundanceFile.py: "
            python /qiimepipe/createAbundanceFile.py {input.tax_assignments} {input.otu_table} > {log} 2>&1
            """

elif marker_gene == "ITS-FUNGI-NCBI":
    if remote == 'yes':
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    else:
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    rule run_qiimepipe_create_tax_assignment:
        input:
            blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log'),
            txt_table = os.path.join(output_dir, file_name_raw + '_pre_otu_table.txt') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_pre_asv_table.txt')
        output:
            tax_assignments = os.path.join(output_dir, 'diversity_by_sample_ncbi', 'ncbi_tax_assignments.txt')
        log:
            os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_create_tax_assignment', file_name_raw + '.log')
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
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:v2"
        shell:
            """
            cd {output_dir}
            echo "Running the QiimePipe Container - filterOTUs.py: "
            python3.6 /qiimepipe/filterOTUs.py {input.fasta_file} {input.tax_assignments} > {log} 2>&1
            """
    if remote == 'yes':
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
    else:
        rule run_blastn_fungi:
            input:
                filtered_fasta = os.path.join(output_dir, 'k__Fungi.fasta')
            output:
                blast_fungi_log = os.path.join(output_dir, file_name_raw + '_blast_fungi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.filtered_fasta} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_fungi_log}
                """
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

elif marker_gene == "ALL-NCBI":
    if remote == 'yes':
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db nt -remote -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
    else:
        rule run_blastn:
            input:
                fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
            output:
                blast_log = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
            params:
                similarity_int_asg = similarity_int_asg,
                coverage_int = coverage_int,
                hits_per_subject = config['hits_per_subject'],
                evalue = config['e_value'],
                NCBI_DB = NCBI_DB,
                threads = config['num_threads'],
                db_type = config['db_type']
            singularity:
                "docker://itvdsbioinfo/pimba_blast:latest"
            shell:
                """
                echo "Running the BLAST Container - blastn: "
                cd {output_dir}
                blastn -query {input.fasta_file} -task megablast -db {params.NCBI_DB}/{params.db_type} -num_threads {params.threads} -perc_identity \
                    {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} -max_hsps {params.hits_per_subject} \
                    -max_target_seqs {params.hits_per_subject} -evalue {params.evalue} \
                    -outfmt "6 qseqid sscinames sseqid staxids stitle pident qcovs evalue" > {output.blast_log}
                """
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

else:
    rule run_blast_assign_taxonomy:
        input:
            fasta_file = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
        output:
            log_file = os.path.join(output_dir, file_name_raw + '_blast_ncbi.log')
        params:
            db_dir = marker_gene,
            similarity_int_asg = similarity_int_asg,
            coverage_int = coverage_int,
            hits_per_subject = config['hits_per_subject'],
            e_value = config['e_value'],
            threads = config['num_threads']
        singularity:
            "docker://itvdsbioinfo/pimba_blast:latest"
        shell:
            """
            echo "Running the BLAST Container - blastn: "
            blastn -query {input.fasta_file} -task megablast -db {params.db_dir}/*.fasta \
            -perc_identity {params.similarity_int_asg} -qcov_hsp_perc {params.coverage_int} \
            -max_hsps {params.hits_per_subject} -max_target_seqs {params.hits_per_subject} \
            -evalue {params.e_value} -num_threads {params.threads} \
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
            db_dir = marker_gene
        singularity:
            "docker://itvdsbioinfo/pimba_qiimepipe:v2"
        log:
            os.path.join(current_path, 'results', '01-run', 'logs', 'qiimepipe_create_taxon_table', file_name_raw + '.log')
        shell:
            """
            mkdir -p {output_dir}/diversity_by_sample
            cd {output_dir}/diversity_by_sample

            echo "Running the QiimePipe Container - createTaxonTable_singleFile_flex.py: "
            python3.6 /qiimepipe/createTaxonTable_singleFile_flex.py {input.blast_log} \
            {input.table} {params.db_dir}/*_tax.txt > {log} 2>&1

            cd {output_dir}
            mkdir -p output
            mv *_otus_tax_assignments.txt {output.tax_assignments}
            mv *_taxon_red_flagged.txt {output.taxon_red_flagged}
            """