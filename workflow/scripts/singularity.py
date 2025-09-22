import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Process some Docker images.")
parser.add_argument(
    "--config", type=str, required=True, help="Path to the configuration YAML file"
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


docker_images = [
    "metashot/itsx:1.1.2-1",
    "itvdsbioinfo/pimba_swarm:v3.1.0",
    "itvdsbioinfo/pimba_fastxtoolkit:v0.0.14",
    "itvdsbioinfo/pimba_perl:v7",
    "itvdsbioinfo/pimba_vsearch:v2.29.1",
    "itvdsbioinfo/pimba_r:latest",
    "itvdsbioinfo/pimba_qiime:latest",
    "itvdsbioinfo/pimba_qiimepipe:v2",
    "itvdsbioinfo/pimba_blast:latest",
    "itvdsbioinfo/pimba_biom:v2.1.10",
    "itvdsbioinfo/pimba_adapterremoval:v2.2.3",
    "itvdsbioinfo/pimba_pear:v0.9.10",
    "itvdsbioinfo/pimba_python_plot:latest",
    "itvdsbioinfo/pimba_phyloseq:v2",
    "itvdsbioinfo/pimba_prinseqpp:v1",
    "itvdsbioinfo/r_curation:v2"
]

sif_dir = read_yaml(args.config)
sif_dir = sif_dir.replace('"', "").replace("'", "")

os.makedirs(os.path.abspath(sif_dir), exist_ok=True)

for image in docker_images:
    print(f"INFO: Fetching {image} to {sif_dir}.")

    # Construct the singularity pull command with output directory
    output_file = os.path.join(
        os.path.abspath(sif_dir), f"{image.split('/')[-1].split(':')[0]}.sif"
    )

    # Check if the SIF file already exists
    if os.path.exists(output_file):
        print(f"INFO: {output_file} already exists. Skipping pull for {image}.")
        continue

    command = f"singularity pull {output_file} docker://{image}"

    # Execute the command
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"INFO: Successfully pulled {image} to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Error pulling {image}: {e}")
