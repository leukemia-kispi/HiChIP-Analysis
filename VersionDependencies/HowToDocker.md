
Optional: for maximum reproducibility

Lock package builds
Inside the container, after testing:

conda list --explicit > macs2_legacy-lock.txt

This freezes every build, including channel hashes.
Later, you can rebuild from it using:

conda create --name macs2_legacy --file macs2_legacy-lock.txt


Tag the image

docker tag macs2_legacy_env macs2_legacy_env:v1
docker save macs2_legacy_env:v1 -o macs2_legacy_env_v1.tar



For Legacy

docker run -it --rm -v /home/ubuntu/HiChIPLegacy_Outputs:/outputs hichip_legacy

To include versions, use this instead:

dpkg-query -W -f='${Package}=${Version}\n' > system-packages.txt

If you ever want to restore those packages in another system (like in a Docker container):

xargs -a system-packages.txt apt-get install -y



⚙️ 2. Docker-outside-of-Docker (DoD) is the preferred method

For most bioinformatics tools that "call Docker images" (e.g., Nextflow, Snakemake, nf-core pipelines), you only need access to the host’s Docker socket, not a Docker installation inside the container.

You simply run your container with:

✅ The correct solution: Mount the host Docker socket

When you start your legacy container, run:

docker run -it --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v /home/ubuntu/HiChIPLegacy_Outputs:/outputs \
  hichip_legacy



Recommended alternative
  Instead of copying everything:

Start from ubuntu:18.04 (or 20.04).

Explicitly install:

Conda

Your legacy environments (exported .yml)

CLI tools missing from YAML (BWA, Samtools, Bowtie2, Bedtools, Picard, MACS2 2.2.6, etc.)


hej