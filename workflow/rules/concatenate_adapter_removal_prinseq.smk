rule concatenate_samples:
    input:
        fastq_sample = expand(os.path.join(current_path, "results", "00-prepare", "prepare_output", "samples_{sample}.fastq"), sample=SAMPLES)
    output:
        fastq_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", config["outputprepare"] + ".fastq")
    shell:
        """
        cat {current_path}/results/00-prepare/prepare_output/samples_R* > {output.fastq_file}
        """

rule run_adapter_removal:
    input:
        fastq_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", config["outputprepare"] + ".fastq")
    output:
        discarded_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_good.discarded").format(output_name=config["outputprepare"]),
        truncated_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_good.truncated").format(output_name=config["outputprepare"]),
        settings_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_good.settings").format(output_name=config["outputprepare"])
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
        AdapterRemoval --file1 {input.fastq_file} --threads {params.num_threads} --trimwindows 10 \
        --minquality {params.minphred} --minlength {params.minlength} --qualitymax 64 --basename {current_path}/results/00-prepare/prepare_output/{params.output_name}_good > {log} 2>&1
        """

rule run_prinseq_truncated:
    input:
        truncated_file = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{output_name}_good.truncated").format(output_name=config["outputprepare"])
    output:
        good_fasta = os.path.join(current_path, "results", "00-prepare", "{output_name}.fasta").format(output_name=config["outputprepare"])
    params:
        output_name = config["outputprepare"]
    log:
        os.path.join(current_path, 'results', '00-prepare', 'logs', 'prinseq_truncated', config["outputprepare"] + '.log')
    singularity:
        "docker://itvdsbioinfo/pimba_prinseq:v0.20.4"
    shell:
        """
        echo "Creating and running a Prinseq Container: "
        prinseq-lite.pl -fastq {input.truncated_file} -out_format 1 -out_good {current_path}/results/00-prepare/prepare_output/{params.output_name} > {log} 2>&1

        mv {current_path}/results/00-prepare/prepare_output/{params.output_name}.fasta {current_path}/results/00-prepare
        rm {current_path}/results/00-prepare/prepare_output/indexforbol_unmatched.fastq
        rm {current_path}/results/00-prepare/prepare_output/indexreveol_unmatched.fastq
        """