rule prepare_directories:
    input:
        full_path_raw = config["rawdatadir"]
    output:
        r1_files = [os.path.join(current_path, "results", "00-prepare", "prepare_output", f"R1/{sample}_R1{suffix}") for sample, suffix in sample_suffixes.items()],
        r2_files = [os.path.join(current_path, "results", "00-prepare", "prepare_output", f"R2/{sample}_R2{suffix}") for sample, suffix in sample_suffixes.items()]
    shell:
        """
        mkdir -p {current_path}/results/00-prepare/prepare_output/R1 {current_path}/results/00-prepare/prepare_output/R2
        cp {input.full_path_raw}/*_R1* {current_path}/results/00-prepare/prepare_output/R1/
        cp {input.full_path_raw}/*_R2* {current_path}/results/00-prepare/prepare_output/R2/
        """

rule run_adapter_removal:
    input:
        r1 = lambda wildcards: os.path.join(current_path, "results", "00-prepare", "prepare_output", f"R1/{wildcards.sample}_R1{sample_suffixes[wildcards.sample]}"),
        r2 = lambda wildcards: os.path.join(current_path, "results", "00-prepare", "prepare_output", f"R2/{wildcards.sample}_R2{sample_suffixes[wildcards.sample]}")
    output:
        os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.discarded"),
        os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair1.truncated"),
        os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair2.truncated"),
        os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.settings"),
        os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.singleton.truncated")
    params:
        adapters = config["adapters"],
        num_threads = config["num_threads"],
        min_phred = config["minphred"],
        min_length = config["minlength"]
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "adapter_removal", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_adapter_removal.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_adapterremoval:v2.2.3"
    shell:
        """
        cd {current_path}/results/00-prepare/prepare_output/
        AdapterRemoval --file1 {input.r1} --file2 {input.r2} --threads {params.num_threads} --mate-separator " " --adapter-list {params.adapters} \
        --trimwindows 10 --minquality {params.min_phred} --minlength {params.min_length} --qualitymax 64 --basename {wildcards.sample}_good --mm 5 \
        > {log} 2>&1
        rm {input.r1} {input.r2}
        """

rule run_pear:
    input:
        pair1 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair1.truncated"),
        pair2 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair2.truncated"),
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.assembled.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.discarded.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.unassembled.forward.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "pear", "{sample}", "{sample}.unassembled.reverse.fastq")
    params:
        num_threads = config["num_threads"]
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "pear", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_pear.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_pear:v0.9.10"
    shell:
        """
        mkdir -p {current_path}/results/00-prepare/assemblies/pear/{wildcards.sample}
        pear -j {params.num_threads} -f {input.pair1} -r {input.pair2} -o {current_path}/results/00-prepare/assemblies/pear/{wildcards.sample}/{wildcards.sample} \
        > {log} 2>&1
        for i in {current_path}/results/00-prepare/assemblies/pear/{wildcards.sample}/*.assembled.fastq; do sed -i "s/ /_/g" $i; done
        """