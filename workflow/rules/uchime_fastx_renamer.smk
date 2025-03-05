rule run_vsearch_uchime_denovo:
    input:
        uchime_input = lambda wildcards: os.path.join(output_dir, file_name_raw + '_otus1.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs1.fasta')
    output:
        chimeras = os.path.join(output_dir, file_name_raw + '_chim.fasta'),
        nonchimeras = os.path.join(output_dir, file_name_raw + '_noChim.fasta'),
        otusnochim = os.path.join(output_dir, file_name_raw + '_uchime_nochim')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "vsearch_uchime_denovo", file_name_raw + ".log")
    singularity:
        "docker://itvdsbioinfo/pimba_vsearch:v2.15.2"
    shell:
        """
        echo "Running the VSEARCH Container - --uchime_denovo:"
        vsearch --uchime_denovo {input.uchime_input} --fasta_width 0 --uchimeout {output.otusnochim} \
        --chimeras {output.chimeras} --nonchimeras {output.nonchimeras} > {log} 2>&1
        """

rule run_fastx_formatter:
    input:
        nonchimeras = os.path.join(output_dir, file_name_raw + '_noChim.fasta')
    output:
        formatted = os.path.join(output_dir, file_name_raw + '_noChim_formatted.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "fastx_formatter", file_name_raw + ".log")
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        echo "Creating and running a fastxtoolkit Container: "
        fasta_formatter -i {input.nonchimeras} -o {output.formatted} > {log} 2>&1
        """

rule run_bmp_renamer:
    input:
        formatted = os.path.join(output_dir, file_name_raw + '_noChim_formatted.fasta')
    output:
        renamed = os.path.join(output_dir, file_name_raw + '_otus.fasta') if strategy == "otu" else os.path.join(output_dir, file_name_raw + '_asvs.fasta')
    log:
        os.path.join(current_path, "results", "01-run", "logs", "bmp_renamer", file_name_raw + ".log")
    singularity:
        "docker://itvdsbioinfo/pimba_perl:v5"
    shell:
        """
        echo "Creating and running a Perl Container: "
        perl /data/bmp-otuName.pl -i {input.formatted} -o {output.renamed} > {log} 2>&1
        """