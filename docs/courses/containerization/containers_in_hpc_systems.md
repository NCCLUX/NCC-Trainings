# Containers in HPC systems

There are multiple container standards and even more container engines and storage formats that support them. Singularity type containers [[1](#ref_Kurtzer_2017a)] are popular in HPC systems due to the simplicity with which they integrate to HPC environments [[2](#ref_Mosciatti_2020a)]. The main container engines supporting Singularity type containers are are [Apptainer](https://apptainer.org/docs/user/latest/) and [Singularity](https://docs.sylabs.io/guides/latest/user-guide/). Both these container engines are built around the Singularity Image File (SIF) format of containers, and are interchangeable.

??? info "Relation between Apptainer and Singularity"

    Singularity began as an open-source project in 2015, when a team of researchers at Lawrence Berkeley National Laboratory. The original project joined the [Linux Foundation](https://www.linuxfoundation.org/) and was renamed to Apptainer, where as the original Singularity code has been forked by the company Sylabs and distributed under various licenses and names.

??? info "Why Singularity type containers are favored in HPC systems"

    The initial versions of [Docker](https://docs.docker.com/) needed significant adjustments to work well in HPC systems due to some design choices. Docker is a client-server architecture composed of a daemon service (dockerd) and client applications (Docker clients). The client must run in the host system with root privileges and shares resources to the clients that run the containers.

    While this designs works well in cloud architectures where the containers often run on VMs, provide services that are designed to run for an indefinite amount of time, it also poses difficulties in HPC systems for 2 main reasons.

    - Deployment: for every Docker container, the container process is spawned as a child of a root owned Docker daemon; while Docker has introduced a rootless mode in 2020, a Docker daemon is still necessary to spawn containers.

    - Security: As the user is able to directly interact with and control the Docker daemon, it is theoretically possible to coerce the daemon process into granting the users escalated privileges; this is particularly problematic for root owned daemons.

    Singularity avoid these issues by launching containers directly as a user process, in the same user namespace with [minimal isolation](https://apptainer.org/docs/user/main/docker_and_oci.html#namespace-device-isolation).

    - There is no daemon process controlling the container, so there no need to start an extra daemon process. Even daemons requiring no root privileges add extra complexity, as they add extra layers of resource management (namespaces for use mappings and cgroups for resource management) on top of the cluster scheduler (Slurm in Meluxina). Singularity containers directly inherit their resource constraints from the systems scheduler and run on the same namespace as the user process that launches the container.

    - There is no need for elevated user rights reducing the risk of privilege escalation, and the user does not interact with processes that have elevated user rights, reducing the attack surface of the system.

    A lot of the drawbacks of the Docker architecture are resolved in [Podman](https://podman.io/docs), a daemonless tool container engine. However, the minimal isolation that Singularity offers by default is better suited for HPC systems. Container engines like Docker and Podman that target cloud application tend to isolate the container environment by default, and apply extra layers of resource management which are not transparent to the system scheduler and must be explicitly deactivated. In contrast, Singularity container engines like Singularity and Apptainer, mounts host directories (like your home, `/tmp`, and device directories) inside the container by default, granting direct access to your local files and GPU hardware without extra configuration.

    A subtle advantage of the low complexity approach favored by Singularity type container engines is that the flat storage format of SIF images is often more efficient in HPC systems. While overlayed image formats where layers can be shared between container instance to save memory are useful in cloud applications where multiple containers are running on a single system, it is rarely useful in HPC where containers are used to package and run an application from a single container. Furthermore, large centers are often using CernMV-FS repositories that load container over a network file system in a lazy manner, only loading the components that are actually used. While storing deduplicated images and loading images in a lazy manner is easy for flat SIF files, extra processing is required to achieve the same for overlayed  OCI Images, resulting in minor performance degradation [[2](#ref_Mosciatti_2020a)].

## Container image formats

The containers are effectively a directory structure into which a chroot jail root together with a set of configuration files is written. However storing and shipping directories is cumbersome. There are image formats that allow serializing the container directory structure into a single binary tarball for easy distributions. However, not all image format are monolithic.

Docker tarballs and [OCI Image](https://github.com/opencontainers/image-spec/blob/main/image-layout.md) of the [Open Container Initiative](https://opencontainers.org/) are using [ovelayFS](https://docs.kernel.org/filesystems/overlayfs.html) to build containers as a stack of layers. For instance, when building a OCI Image one layer may include the Debian base image with `libc`, a second level may install Python, and a third level may install a python `venv`. The layers are assembled dynamically during runtime. This saves space in container registries by storing a single copy of each shared layer, and runtime memory and computation in container engines by loading a single copy of the required shared libraries at runtime.

The [Singularity Image Format (SIF)](https://github.com/apptainer/sif) is the format used natively by Singularity type containers. Similar to Singularity type container engines, SIF emphasizes simplicity. The is a single overlay layer in SIF files. Singularity type container engines can work with OCI Images and Docker tarballs, but the have to convert them in a SIF compatible format first. There is an OCI compatible format for SIF containers, [OCI-SIF](https://docs.sylabs.io/guides/latest/user-guide/oci_runtime.html#oci-sif-images) that stores the OCI image with its overlays in a SIF file. Launching a container from an OCI-SIF format directly is [supported in Singularity](https://docs.sylabs.io/guides/latest/user-guide/oci_runtime.html#oci-mode-oci), but [not in Apptainer](https://github.com/sylabs/singularity/discussions/2948) yet. Apptainer users must convert their containers to a format compatible with conventional SIF file to run them.

## Using Apptiner in HPC systems

To run containers in HPC systems, you need to access the Apptainer or Singularity executables and and a container file. Apptainer is already installed on Meluxina, and is available as a module. To use it, load the Apptainer module.

```console
module load env/release/2025.1
module load Apptainer
```

There are multiple [container registries](#container-registries), collections container images, where you can find images for many popular applications. We now demonstrate how you can fetch, store, and run an existing container image with Apptainer.

### Example: Pulling a container into a SIF file

Registries are data collections that can be accessed with programs designed to fetch data from the registry. Apptainer provides the `pull` function that fetches containers from registries and store them locally into a SIF file. In this example we fetch a container and store it into a SIF file locally.

1. Create a directory to store your container images.

  ```console
  mkdir -p ${HOME}/containers
  ```

1. Fetch the [`rancher/cowsay`](https://hub.docker.com/r/rancher/cowsay) from Docker hub and store it locally.

  ```console
  apptainer pull ${HOME}/containers/cow.sif docker://rancher/cowsay:latest
  ```

??? question "Managing container images"

    === "Question"

        What operations did the pull command perform?

    === "Answer"

        Apptainer fetched the image for `rancher/cowsay`.

        - The argument `pull` instructs the Apptainer to fetch the container image and store it locally. We could have used the commands `run` or `exec` to run the image without storing it locally.

        - The argument `${HOME}/containers/cow.sif` is the location where the container image is stored. By default, the container is stored as a SIF image.

        - The argument `docker://rancher/cowsay:latest` indicates the source of the container images, in this case the [Docke hub](https://hub.docker.com/) registry. Since the images in Docker hub are not stored in SIF format, Apptainer will convert the image into a SIF file before executing any further steps.

        Artifacts of image conversions are stored in `${HOME}/.apptainer/cache/`.

### Example: Running containers

You can now launch the container stored locally. Singularity type container engines provide two launch commands,

- the `exec` command that takes as argument the executable to execute inside the container, and
- the `run` command that run a default executable defined during the container creation.

!!! question "Run the `cowsay` command"

    === "Question"

        The cowsay container is designed to run the cowsay command. Run the container with the `run` command and argument `Hellow, world!`, and with the `exec` command and calling `cowsay` explicitly.

    === "Answer"

        Executing the container with the run command produces the following result.

        ```console
        $ apptainer run ${HOME}/containers/cow.sif 'Hello, world!'
         _______________ 
        < Hello, world! >
         --------------- 
                \   ^__^
                 \  (oo)\_______
                    (__)\       )\/\
                        ||----w |
                        ||     ||
        ```

        To run the same command with `exec` we call the container as follows.

        ```console
        $ apptainer exec ${HOME}/containers/cow.sif cowsay 'Hello, world!'
         _______________ 
        < Hello, world! >
         --------------- 
                \   ^__^
                 \  (oo)\_______
                    (__)\       )\/\
                        ||----w |
                        ||     ||
        ```

### Example: Building a container from a definition file

In this example we build Singularity type containers with Apptainer from a Singularity definition files. The definition file (`.def`) a list of instruction that a container builder program like Apptainer can use to compose a container binary file.

The first container we build, is a container for [DIA-NN](https://github.com/vdemichev/DiaNN), an automated software suite for data-independent acquisition (DIA) proteomics data processing. The software is [supported for a specific version of Linux Mint](https://github.com/vdemichev/DiaNN#installation), and thus requires a container to run in Meluxina that uses a RHEL based system. The `.def` contains the following instructions.

!!! abstract "DIA-NN-2.5.1.def"

    ```singularity
    BootStrap: docker
    From: linuxmintd/mint21.2-amd64:latest

    %arguments
      VERSION=2.5.1
      RELEASE=2.0

    %post
      export DIA_NN='DIA-NN-{{VERSION}}-Academia-Linux.zip'
      if [ -f "/tmp/${DIA_NN}" ]; then
        rm "/tmp/${DIA_NN}"
      fi
      wget "https://github.com/vdemichev/DiaNN/releases/download/{{RELEASE}}/${DIA_NN}" --directory-prefix=/tmp
      unzip "/tmp/${DIA_NN}" -d /opt
      chmod ugo+x /opt/diann-{{VERSION}}/diann-linux
      rm "/tmp/${DIA_NN}"

    %environment
      export LC_ALL=C
      export PATH=/opt/diann-{{VERSION}}:${PATH}

    %runscript
      /opt/diann-{{VERSION}}/diann-linux

    %labels
      Container DIA-NN
      Version {{VERSION}}
      Author hpc-team@uni.lu

    %help
      A a universal software suite for data-independent acquisition (DIA) proteomics data processing. Conceived at the University of Cambridge, UK, in the laboratory of Kathryn Lilley (Cambridge Centre for Proteomics), DIA-NN opened a new chapter in proteomics, introducing a number of algorithms which enabled reliable, robust and quantitatively accurate large-scale experiments using high-throughput methods. DIA-NN is currently being further developed in the laboratory of Vadim Demichev at the Charité (University Medicine Berlin, Germany).

      Source: https://github.com/vdemichev/DiaNN
    ```

The definition file, is text file with a description on how to build the _DIA-NN_ container:

1. bootstrap the container by pulling from the Docker image `linuxmintd/mint21.2-amd64:latest`;
1. download and extract the DIA-NN binaries, and cleanup any artifacts leftover;
1. set some environment variables to access the binaries;
1. define a command to run when executing the container.

A more complete list of section names for the definition is available in the [Apptainer documentation](https://apptainer.org/docs/user/main/definition_files.html).

From this definition file, you can build the container with the following command.

```bash
apptainer build ${HOME}/containers/DIA-NN.sif DIA-NN.def
```

This command takes the definition file `DIA-NN.def` to create a new container image in the file `DIA-NN.sif`.

!!! question "Building and running containers from definition files"

    === "Question"

        Build the container from the definition file. Then, fetch the example input file [`Aug2022_Mouse Tissue Contaminants.fasta`](https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Sample-type%20specific%20contaminant%20FASTA/Aug2022_Mouse%20Tissue%20Contaminants.fasta), and run the following command in the container.

        ```bash
        diann-linux
            --lib "" \
            --out "/tmp/report.parquet" \
            --out-lib "/tmp/report-lib.parquet" \
            --fasta "Aug2022_Mouse Tissue Contaminants.fasta" \
            --threads 32 \
            --verbose 1 \
            --qvalue 0.01 \
            --matrices \
            --gen-spec-lib \
            --predictor \
            --reannotate \
            --fasta-search \
            --min-fr-mz 200 \
            --max-fr-mz 1800 \
            --met-excision \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --min-pr-mz 300 \
            --max-pr-mz 1800 \
            --min-pr-charge 1 \
            --max-pr-charge 4 \
            --cut K*,R* \
            --missed-cleavages 1 \
            --unimod4 \
            --rt-profiling
        ```

    === "Answer"

        From this definition file, you can build the container with the following command.

        ```bash
        apptainer build ${HOME}/containers/DIA-NN.sif DIA-NN.def
        ```

        Then, execute the container with the following command.

        ```bash
        apptainer exec ${HOME}/containers/DIA-NN.sif diann-linux
            --lib "" \
            --out "/tmp/report.parquet" \
            --out-lib "/tmp/report-lib.parquet" \
            --fasta "Aug2022_Mouse Tissue Contaminants.fasta" \
            --threads 32 \
            --verbose 1 \
            --qvalue 0.01 \
            --matrices \
            --gen-spec-lib \
            --predictor \
            --reannotate \
            --fasta-search \
            --min-fr-mz 200 \
            --max-fr-mz 1800 \
            --met-excision \
            --min-pep-len 7 \
            --max-pep-len 30 \
            --min-pr-mz 300 \
            --max-pr-mz 1800 \
            --min-pr-charge 1 \
            --max-pr-charge 4 \
            --cut K*,R* \
            --missed-cleavages 1 \
            --unimod4 \
            --rt-profiling
        ```

!!! question "Provide an interface for using the DIA-NN container"

    === "Question"

        The current version of `DIA-NN.def` creates a container where the `run` argument does not accept any arguments and thus it is not very useful. Modify the container definition file so that it access arguments for the DIA-NN executable.

    === "Answer"

        We modify the `runscript` section to pass the arguments to the `run` command through to the call to `diann-linux` executable. The resulting definition file is the following.

        ```singularity
        BootStrap: docker
        From: linuxmintd/mint21.2-amd64:latest

        %arguments
          VERSION=2.5.1
          RELEASE=2.0

        %post
          export DIA_NN='DIA-NN-{{VERSION}}-Academia-Linux.zip'
          if [ -f "/tmp/${DIA_NN}" ]; then
            rm "/tmp/${DIA_NN}"
          fi
          wget "https://github.com/vdemichev/DiaNN/releases/download/{{RELEASE}}/${DIA_NN}" --directory-prefix=/tmp
          unzip "/tmp/${DIA_NN}" -d /opt
          chmod ugo+x /opt/diann-{{VERSION}}/diann-linux
          rm "/tmp/${DIA_NN}"

        %environment
          export LC_ALL=C
          export PATH=/opt/diann-{{VERSION}}:${PATH}

        %runscript
          /opt/diann-{{VERSION}}/diann-linux "$@"

        %labels
          Container DIA-NN
          Version {{VERSION}}
          Author hpc-team@uni.lu

        %help
          A a universal software suite for data-independent acquisition (DIA) proteomics data processing. Conceived at the University of Cambridge, UK, in the laboratory of Kathryn Lilley (Cambridge Centre for Proteomics), DIA-NN opened a new chapter in proteomics, introducing a number of algorithms which enabled reliable, robust and quantitatively accurate large-scale experiments using high-throughput methods. DIA-NN is currently being further developed in the laboratory of Vadim Demichev at the Charité (University Medicine Berlin, Germany).

          Source: https://github.com/vdemichev/DiaNN
        ```

        Observe the extra "$@" argument to the call to the `diann-linux` executable in the `runscript` section.

### Example: Containers that require elevated right to be built

The second example is based on a definition file from the [documentation of Apptainer](https://apptainer.org/docs/user/latest/build_a_container.html#building-containers-from-apptainer-definition-files). The example creates a _lolcow_ container which is a combination of the [cowsay](https://en.wikipedia.org/wiki/Cowsay), [lolcat](https://github.com/busyloop/lolcat), and [fortune](https://en.wikipedia.org/wiki/Fortune_(Unix)) tools. The tools are installed from the repositories of Ubuntu using `apt-get`, and thus require elevated privileges to install in the container, either root access or Apptainer setup to run with `fakeroot`. Thus, you want be able to build this container in Meluxina, and you need to build it in a machine where you have the appropriate rights.

The definition file `lolcow.def` for this container is the following.

!!! abstract "lolcow.def"

    ```singularity
    Bootstrap: docker
    From: ubuntu:20.04
    
    %post
        apt-get -y update
        apt-get -y install cowsay lolcat fortunes
    
    %environment
        export LC_ALL=C
        export PATH=/usr/games:$PATH
    
    %runscript
        if [ "$#" -gt 0 ]; then
            cowsay "$@" | lolcat
        else
            fortune | cowsay | lolcat
        fi
    ```

The definition file describes the procedure to build the _lolcow_ container:

1. bootstrap the container by pulling from the Docker image `ubuntu:20.04`;
1. install the package `cowsay` and `lolcat`;
1. set some environment variables;
1. define a command to run when executing the container.

From this definition file, you can build the container with this command:

```bash
apptainer build lolcow.sif lolcow.def
```

This command takes the definition file `lolcow.def` to create a new container image in the file `lolcow.sif`.

??? info "The `runscript` section"

    The optional `runscript` section is used to define a script that is called when apptainer is called with the `run` argument. The contents of the `runscript` section are copied in a bash script `/singularity` in the container file system, and the `/singularity` script is executed when the `run` argument is used.

    In the `lolcow.sif` container for instance, the contents of the `/singularity` are the following.

    !!! abstract "/singularity"
        ```bash
        #!/bin/bash
    
        if [ "$#" -gt 0 ]; then
            cowsay "$@" | lolcat
        else
            fortune | cowsay | lolcat
        fi
        ```

    Furthermore, a call to the `run` command

    ```console
    apptainer run lolcow.sif <arguments>
    ```

    is equivalent to the more verbose

    ```console
    apptainer exec lolcow.sif /singularity <arguments>
    ```

    call to the `exec` command. Thus the `runscript` section provides the means of defining a interface for using the container.

### Example: Unpacking containers into sandboxes

Sandboxes are expansions of image files into a regular directory on the file system. The user may access and modify the sandbox form the host system or they can launch the sandbox as a container and modify it from withing. Sandboxes can be packages into SIF images with the container application.

!!! question "Modifying a container in a sandbox"

    === "Question"

        The [`rancher/cowsay`](https://hub.docker.com/r/rancher/cowsay) container hangs when no argument is given to the `run` command. Modify the `/singularity` script of the container to print the date (output of the `date` command) when no argument is given.

    === "Answer"

        Pull and unpack the container in the working directory.

        ```console
        apptainer build --sandbox cowsay/ docker://rancher/cowsay:latest
        ```

        The modify the container in the host filesystem.

        ```console
        cat > cowsay/singularity <<EOF
        #!/bin/bash

        if [ "$#" -gt 0 ]; then
          cowsay "$@"
        else
          cowsay `date`
        fi
        EOF
        ```

        Then, repackage the container.

        ```
        apptainer build ${HOME}/containers/cow.sif cowsay/
        ```

        Finally, cleanup your sandbox.

        ```console
        chmod -R u+rwX cowsay
        rm -r cowsay
        ```

!!! tip "Usage of inodes"

    Note that unpacking an image into a sandbox generates a lot of files. In cluster file systems with quota on the number of inodes, you can exhaust your quote fast. Remember to repackage and cleanup your sandbox after you have finish you modifications.

## Resources

### References

1. <span id="ref_Kurtzer_2017a">Kurtzer, Gregory M., Vanessa, Sochat, and Michael W., Bauer. "Singularity: Scientific containers for mobility of compute". _PLOS ONE_ 12, no.5 (2017): e0177459, doi: [10.1371/journal.pone.0177459](https://www.doi.org/10.1371/journal.pone.0177459).</span>

1. <span id="ref_Mosciatti_2020a">Mosciatti, S. and Blomer, J. and Ganis, G. and Popescu, R. "CernVM-FS Container Image Integration". _Journal of Physics: Conference Series_ 1525, no. 1 (2020): 012058, doi: [10.1088/1742-6596/1525/1/012058](https://www.doi.org/10.1088/1742-6596/1525/1/012058).</span>

### Popular container registries { #container-registries }

- [Docker Hub](https://hub.docker.com/)
- [NVIDIA NGC Catalog](https://catalog.ngc.nvidia.com/containers)
- [AMD Infinity Hub](https://www.amd.com/en/developer/resources/infinity-hub.html)
- [Sylabs Cloud Library](https://cloud.sylabs.io/library)
- [RED HAT Quay.io](https://quay.io/)
- [BioContainers](https://biocontainers.pro/)
