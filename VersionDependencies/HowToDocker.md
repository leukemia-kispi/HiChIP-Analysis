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