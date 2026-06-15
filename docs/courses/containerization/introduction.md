# The Linux environment

Environment management in shared computing resources, like HPC systems and workstations, involves a collection of techniques that allow user processes to run in isolation from other user and system processes. This isolation allows

- non-kernel processes to run on environments isolated from other processes, and
- multiple processes to share system resources without interfering with each other.

This isolation is particularly useful in HPC systems where multiple users install and use software simultaneously.

## The Linux process interface

Processes in a GNU/Linux system are categorized into

- user processes,
- operating system processes, and
- the system kernel.

A user process interacts with other process, through function calls ([fig. 1](#fig-1)) and interprocess communication methods, such as streams and message queues.

<figure markdown="span">
  ![Linux process types and interactions](resources/linux_structure.png){ align=left width="400" }

  <figcaption id="fig-1"><b>Fig. 1:</b> Linux process types and interactions.</figcaption>
</figure>

The Linux operating system provides mechanisms to modify the view that a user process has of the operating system and the processes that are running on it. Containers and related methods utilize these mechanisms to modify the view that a process has of the system.

## Containerization technologies

The term _container_ describes a wide and loose set of technologies that enable software applications to run in isolated user spaces. The 3 main technologies used in containers are the following.

- [_Chroot_](https://en.wikipedia.org/wiki/Chroot): a system operation that changes the apparent root directory for the current running process and its children.
- [_Namespaces_](https://en.wikipedia.org/wiki/Linux_namespaces): a Linux kernel feature that limits the resources that a set of processes can access; examples of resource include process IDs, user IDs, and interprocess communication mechanisms and groups.
- [_Control group_ (`cgroup`)](https://en.wikipedia.org/wiki/Cgroups): a Linux kernel feature that limits, accounts for, and isolates the resource usage for collections of processes; resource include CPU cores, memory, and disk IOPs.

<figure markdown="span">
  ![Overview of Linux containers](images/linux-containers.png){ align=left width="400" }

  <figcaption><b>Fig. 2:</b> Overview of technologies used in Linux containers.</figcaption>
</figure>

The complete set of technologies used in containers and virtual environments involves a few more features of the operating system, like

- `PATH`, an environment variable used to set the location where applications look for executables,
- `LD_LIBRARY_PATH`, an environment variable used by the linker to determine the locations where to look for shared libraries, and
- tools such as `fakeroot` and `chroot` used to build isolated environments.

In the following sessions we will look into the details of these techniques, and how they are used in practice. An overview the various methods and their relations is presented in ([fig. 3](#fig-3)).

<figure markdown="span">
  ![Containerization techniques](resources/linux_interface.png){ align=left width="800" }

  <figcaption id="fig-3"><b>Fig. 3:</b> Linux application interface and containerization techniques. Containerization methods control the interface available between the users process, the operating system, and the operating system kernel. </figcaption>
</figure>

## Further resources

### Extra reading

- [Apptainer and MPI applications](https://apptainer.org/docs/user/latest/mpi.html)
- [GPU Support (NVIDIA CUDA & AMD ROCm)](https://apptainer.org/docs/user/latest/gpu.html)
- [Portable MPI containerization with the Process Management Interface (PMI)](https://ciq.com/blog/a-new-approach-to-mpi-in-apptainer/)
- [Docker vs Apptainer](https://www.linkedin.com/pulse/docker-vs-apptainer-anup-khanal-vxvxf/)
- [Container image formats under the hood](https://snyk.io/blog/container-image-formats/)

### Terms

- _Environments_ modify the user space to provide access to different sets of executables and libraries without modifying globally installed tools. Environment management tools such as Environment Modules, Conda, and application specific tools such Python Virtual Environments (`venv`) rely on environment variables such as `PATH` and `LD_LIBRARY_PATH`.

- _Containers_ on top of user space modifications can also modify aspects of the operating system. Singularity containers for instance, can use a custom GNU C Library (`glibc`), a wrapper around the system calls of the Linux kernel. Containers rely on a variety of tools to control the environment, access to files, and access to resources such as cores and memory.

- Anything that requires customization of the kernel and drivers must use _virtualisation_, which is beyond the scope of most HPC applications.
