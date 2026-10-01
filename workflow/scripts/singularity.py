import os
import subprocess
import argparse


parser = argparse.ArgumentParser(description="Process some Docker images.")
parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to the configuration YAML file"
)

args = parser.parse_args()


def read_yaml(file_path):
    sif_dir = None

    with open(file_path, "r") as file:
        for line in file:
            # Strip whitespace and check for the sif_dir key
            line = line.strip()

            if line.startswith("sif_dir:"):
                # Extract the value after the colon and strip whitespace
                sif_dir = line.split(":", 1)[1].strip()
                break

    return sif_dir


# Docker images and corresponding Singularity image names
docker_images = {
    "metashot/itsx:1.1.2-1": "itsx.sif",
    "itvdsbioinfo/pimba_swarm:v3.1.0": "pimba_swarm.sif",
    "itvdsbioinfo/pimba_fastxtoolkit:v0.0.14": "pimba_fastxtoolkit.sif",
    "itvdsbioinfo/pimba_perl:v7": "pimba_perl.sif",
    "itvdsbioinfo/pimba_vsearch:v2.29.1": "pimba_vsearch.sif",
    "itvdsbioinfo/pimba_r:latest": "pimba_r.sif",
    "itvdsbioinfo/pimba_qiime:latest": "pimba_qiime.sif",
    "itvdsbioinfo/pimba_qiimepipe:v4": "pimba_qiimepipe.sif",
    "itvdsbioinfo/pimba_blast:latest": "pimba_blast.sif",
    "itvdsbioinfo/pimba_biom:v2.1.10": "pimba_biom.sif",
    "itvdsbioinfo/pimba_adapterremoval:v2.2.3": "pimba_adapterremoval.sif",
    "quay.io/biocontainers/cutadapt:5.2--py313hd978853_2": "pimba_cutadapt.sif",
    "itvdsbioinfo/pimba_pear:v0.9.10": "pimba_pear.sif",
    "itvdsbioinfo/pimba_python_plot:latest": "pimba_python_plot.sif",
    "itvdsbioinfo/pimba_phyloseq:v2": "pimba_phyloseq.sif",
    "itvdsbioinfo/pimba_prinseqpp:v1": "pimba_prinseqpp.sif",
    "itvdsbioinfo/r_curation:v6": "r_curation.sif",
    "itvdsbioinfo/pimba_mumu:v1.1.4": "pimba_mumu.sif",
    "itvdsbioinfo/pimba_krona:v1": "pimba_krona.sif"
}


sif_dir = read_yaml(args.config)
sif_dir = sif_dir.replace('"', "").replace("'", "")

os.makedirs(os.path.abspath(sif_dir), exist_ok=True)


for image, sif_name in docker_images.items():

    print(f"INFO: Fetching {image} to {sif_dir}.")

    output_file = os.path.join(
        os.path.abspath(sif_dir),
        sif_name
    )

    # Check if the SIF file already exists
    if os.path.exists(output_file):
        print(
            f"INFO: {output_file} already exists. "
            f"Skipping pull for {image}."
        )
        continue

    command = (
        f"singularity pull "
        f'"{output_file}" '
        f'"docker://{image}"'
    )

    try:
        subprocess.run(
            command,
            shell=True,
            check=True
        )

        print(
            f"INFO: Successfully pulled {image} "
            f"to {output_file}"
        )

    except subprocess.CalledProcessError as e:
        print(
            f"ERROR: Error pulling {image}: {e}"
        )