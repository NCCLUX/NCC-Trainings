# The Linux environment

Environment management in shared computing resources like HPC systems and workstations involves a collection of techniques that allow user processes to run in isolation from other user and system processes. The isolation allows

- multiple processes to share the same system without interfering with each other and with system processes with respect to resource usage, and
- each processes can configure its own software environment without interfering with the other user space and operating system processes.

The isolation methods can be categorized in methods using the  

<figure markdown="span">
  ![Container architecture](resources/linux_interface.png){ align=left width="400" }
  <figcaption>Fig. 1: The linux application interface</figcaption>
</figure>


## Terms

_Environments_ modify the user space to provide access to different sets of executables and libraries without modifying globally installed tools. Environment management tools such as Environment Modules, Conda, and application specific tools such Python Virtual Environments (`venv`) rely on environment variables such as `PATH` and `LD_LIBRARY_PATH`.

_Containers_ on top of user space modifications can also modify aspects of the operating system. Singularity containers for instance, can use a custom GNU C Library (`glibc`), a wrapper around the system calls of the Linux kernel. Containers rely on a variety of tools to control the environment, access to files, and access to resources such as cores and memory.

Anything that requires customization of the kernel and drivers must use _virtualisation_, which is beyond the scope of most HPC applications.
