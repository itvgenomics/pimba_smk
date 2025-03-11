rule run_fastx_barcode_splitter_forbol:
    input:
        rawdata = raw_fastq,
        barcodes_3end_txt = barcodes_3end_txt
    output:
        forbol_file = expand(os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexforbol_{sample}.fastq"), sample=SAMPLES),
        unmatchedforbol = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexforbol_unmatched.fastq")
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'fastx_barcode_splitter_forbol', config["outputprepare"] + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        echo "Creating and running a fastxtoolkit Container: "
        mkdir -p {current_path}/results/00-prepare/prepare_output
        cat {input.rawdata} | fastx_barcode_splitter.pl --bcfile {input.barcodes_3end_txt} --prefix {current_path}/results/00-prepare/prepare_output/indexforbol_ --suffix .fastq --bol --exact > {log} 2>&1
        """

rule run_fastx_barcode_splitter_reveol:
    input:
        unmatchedforbol = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexforbol_unmatched.fastq"),
        barcodes_3end_rev = config["barcodes_3end_rev"]
    output:
        reveol_file = expand(os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexreveol_{sample}.fastq"), sample=SAMPLES),
        unmatchedreveol = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexreveol_unmatched.fastq")
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'fastx_barcode_splitter_reveol', config["outputprepare"] + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        cat {input.unmatchedforbol} | fastx_barcode_splitter.pl --bcfile {input.barcodes_3end_rev} --prefix {current_path}/results/00-prepare/prepare_output/indexreveol_ --suffix .fastq --eol --exact > {log} 2>&1
        """

rule run_fastx_reverse_complement_reveol:
    input:
        reveol_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexreveol_{sample}.fastq")
    output:
        creveol_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexcreveol_{sample}.fastq")
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        fastx_reverse_complement -i {input.reveol_file} > {output.creveol_file}
        rm {input.reveol_file}
        """

rule concatenate_barcodes:
    input:
        forbol_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexforbol_{sample}.fastq"),
        creveol_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "indexcreveol_{sample}.fastq")
    output:
        index_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "index_{sample}.fastq")
    shell:
        """
        cat {input.forbol_file} {input.creveol_file} > {output.index_file}
        rm {input.forbol_file} {input.creveol_file}
        """