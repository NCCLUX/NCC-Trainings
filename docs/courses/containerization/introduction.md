# The Linux environment

Environment management in shared computing resources, like HPC systems and workstations, involves a collection of techniques that allow user processes to run in isolation from other user and system processes. This isolation allows

- non-kernel processes to setup software environments isolated from other processes, and
- multiple processes to share system resources without interfering with each other.

## The Linux process interface

Processes in a Linux system can be categorized into

- user processes,
- operating system processes, and
- the system kernel.

A user process can interact with other process, through function calls ([fig. 1](#fig-1)) and interprocess communication methods, such as streams and message queues.

<figure markdown="span">
  ![Linux process types and interactions](resources/linux_structure.png){ align=left width="400" }
  <figcaption id="fig-1"><b>Fig. 1:</b> Linux process types and interactions</figcaption>
</figure>

The Linux kernel and operating systems provides mechanism that modify the interface exposed to processes. Containers and related methods utilize these mechanisms to modify the view that a process has of the system ([fig. 2](#fig-2)).

<figure markdown="span">
  ![Containerization techniques](resources/linux_interface.png){ align=left width="600" }
  <figcaption>Fig. 2: Linux application interface and containerization techniques</figcaption>
</figure>

## Terms

_Environments_ modify the user space to provide access to different sets of executables and libraries without modifying globally installed tools. Environment management tools such as Environment Modules, Conda, and application specific tools such Python Virtual Environments (`venv`) rely on environment variables such as `PATH` and `LD_LIBRARY_PATH`.

_Containers_ on top of user space modifications can also modify aspects of the operating system. Singularity containers for instance, can use a custom GNU C Library (`glibc`), a wrapper around the system calls of the Linux kernel. Containers rely on a variety of tools to control the environment, access to files, and access to resources such as cores and memory.

Anything that requires customization of the kernel and drivers must use _virtualisation_, which is beyond the scope of most HPC applications.
