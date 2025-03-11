rule run_fastx_barcode_splitter:
    input:
        rawdata = raw_fastq,
        barcodes_5end_txt = barcodes_5end_txt
    output:
        forbol_file = expand(os.path.join(current_path, "results", "00-prepare", "prepare_output", prefix + "_{sample}.fastq"), sample=SAMPLES),
        unmatchedforbol = os.path.join(current_path, "results", "00-prepare", "prepare_output", prefix + "_unmatched.fastq")
    params:
        prefix = prefix
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'fastx_barcode_splitter', config["outputprepare"] + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_fastxtoolkit:v0.0.14"
    shell:
        """
        echo "Creating and running a fastxtoolkit Container: "
        mkdir -p {current_path}/results/00-prepare/prepare_output
        cat {input.rawdata} | fastx_barcode_splitter.pl --bcfile {input.barcodes_5end_txt} --prefix {current_path}/results/00-prepare/prepare_output/{params.prefix}_ --suffix .fastq --bol --exact > {log} 2>&1
        """

rule concatenate_samples:
    input:
        fastq_sample = expand(os.path.join(current_path, "results", "00-prepare", "prepare_output", prefix + "_{sample}.fastq"), sample=SAMPLES)
    output:
        fastq_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", config["outputprepare"] + ".fastq")
    params:
        prefix = prefix
    shell:
        """
        cat {current_path}/results/00-prepare/prepare_output/{params.prefix}_R* > {output.fastq_file}
        """

rule run_adapter_removal:
    input:
        fastq_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", config["outputprepare"] + ".fastq")
    output:
        discarded_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_min.discarded").format(output_name=config["outputprepare"]),
        truncated_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_min.truncated").format(output_name=config["outputprepare"]),
        settings_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_min.settings").format(output_name=config["outputprepare"])
    params:
        num_threads = config["num_threads"],
        minphred = config["minphred"],
        minlength = config["minlength"],
        output_name = config["outputprepare"]
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'adapter_removal', config["outputprepare"] + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_adapterremoval:v2.2.3"
    shell:
        """
        echo "Creating and running an AdapterRemoval Container: "
        AdapterRemoval --file1 {input.fastq_file} --threads {params.num_threads} \
        --minlength {params.minlength} --qualitymax 64 --basename {current_path}/results/00-prepare/prepare_output/{params.output_name}_min > {log} 2>&1
        """