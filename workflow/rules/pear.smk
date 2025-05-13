rule run_pear:
    input:
        pair1 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair1.truncated"),
        pair2 = os.path.join(current_path, "results", "00-prepare", "prepare_output", "{sample}_good.pair2.truncated"),
    output:
        os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.assembled.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.discarded.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.unassembled.forward.fastq"),
        os.path.join(current_path, "results", "00-prepare", "assemblies", "{sample}", "{sample}.unassembled.reverse.fastq")
    params:
        num_threads = config["num_threads"],
        minoverlap = config["minoverlap"]
    log:
        os.path.join(current_path, "results", "00-prepare", "logs", "pear", "{sample}.log")
    benchmark:
        os.path.join(current_path, "results", "benchmark", "{sample}_prepare_pear.txt")
    singularity:
        "docker://itvdsbioinfo/pimba_pear:v0.9.10"
    shell:
        """
        mkdir -p {current_path}/results/00-prepare/assemblies/{wildcards.sample}
        pear -j {params.num_threads} -f {input.pair1} -v {params.minoverlap} -r {input.pair2} -o {current_path}/results/00-prepare/assemblies/{wildcards.sample}/{wildcards.sample} \
        > {log} 2>&1
        for i in {current_path}/results/00-prepare/assemblies/{wildcards.sample}/*.assembled.fastq; do sed -i "s/ /_/g" $i; done
        """