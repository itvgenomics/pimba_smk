rule run_qiimepipe_strip:
    input:
        index_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "index_{sample}.fastq")
    output:
        clipped_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_clipped.fastq")
    params:
        reverse_adapter = config["reverse_adapter"],
        barcodes_3end_fasta = config["barcodes_3end_fasta"]
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        echo "Running the QiimePipe Container: "
        python3.6 /qiimepipe/fastq_strip_barcode_relabel2.py {input.index_file} {params.reverse_adapter} {params.barcodes_3end_fasta} Ex > {output.clipped_file}
        rm {input.index_file}
        """

rule run_prinseq:
    input:
        clipped_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_clipped.fastq")
    output:
        filtered_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_min50.fastq")
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'prinseq', config["outputprepare"] + '_{sample}' + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_prinseq:v0.20.4"
    shell:
        """
        echo "Creating and running a Prinseq Container: "
        prinseq-lite.pl -fastq {input.clipped_file} -min_len 50 -out_good {current_path}/results/00-prepare/prepare_output/{wildcards.sample}_min50 > {log} 2>&1
        rm {input.clipped_file}
        """

rule run_fastx_reverse_complement_filtered:
    input:
        filtered_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_min50.fastq")
    output:
        revclipped_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_revclipped.fastq")
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'fastx_reverse_complement_filtered', config["outputprepare"] + '_{sample}' + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        echo "Creating and running a fasxttoolkit Container: "
        fastx_reverse_complement -i {input.filtered_file} -o {output.revclipped_file} > {log} 2>&1
        rm {input.filtered_file} {current_path}/results/00-prepare/prepare_output/{wildcards.sample}_clipped_prinseq_bad_*.fastq
        """

rule run_qiimepipe_strip_revclipped:
    input:
        revclipped_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_revclipped.fastq")
    output:
        fastq_sample = os.path.join(current_path, "results", "00-prepare", "prepare_output", "samples_{sample}.fastq")
    params:
        forward_adapter = config["forward_adapter"],
        barcodes_5end_dir = config["barcodes_5end_dir"]
    singularity:
        "docker://itvdsbioinfo/pimba_qiimepipe:v2"
    shell:
        """
        echo "Running the QiimePipe Container: "
        python3.6 /qiimepipe/fastq_strip_barcode_relabel2.py {input.revclipped_file} {params.forward_adapter} {params.barcodes_5end_dir}/barcodes_{wildcards.sample}.fasta Ex > {output.fastq_sample}
        rm {input.revclipped_file}
        """