# Basic principles of containers

Containers use modern kernel features to enable running applications in isolated environments. The techniques used by container allow much better isolation of the environment from the underlying system that what is offered by virtual environments. Operating systems like GNU/Linux usually separate processes 2 groups, the kernel and the user space, that run in different memory regions. The kernel is a privileged process that schedules processes, manages resources such as CPU and memory usage, and handles devices such as storage, network, and I/O devices. Processes request services from the kernel invoking _system calls_ (see [fig. 1](#fig-1)). For instance, the kernel provides functions for applications to write to storage devices. The user space contains user processes and operating system libraries and executables required for user programs to run. For instance the `glibc` package of Debian provides some version of the `libc` library as part of the operating system. The `libc` library  defines functions such as `malloc`, that provides an interface for user programs to access memory in a manner that hides the complexities of memory management by the kernel.

<figure markdown="span">
  ![Containerization techniques](resources/linux_interface.png){ align=left width="800" }

  <figcaption id="fig-1"><b>Fig. 1:</b> Linux application interface and containerization techniques. Containerization methods control the interface available between the users process, the operating system, and the operating system kernel. </figcaption>
</figure>

Typically, each container has its own `\bin` and `\lib` directories for instance that completely override the host system directories inside the container. Since operating system binaries and libraries reside in `\bin` and `\lib` respectively, containers can run any operating system whose kernel is compatible with the host kernel. Despite each container running its own operating system, all containers use the kernel of the host system, as depicted in [fig. 2](#fig-2).

<figure markdown="span">
  ![Container architecture](resources/container_architecture.png){ align=left width="800" }

  <figcaption id="fig-2">Fig. 2: A diagram of a container engine running 2 container, one container with 2 processes and one container with a single process.</figcaption>
</figure>

??? info "Differences between containers and virtual machines"

    Containers perform similar but distinct functions to virtual machines. In virtualisation systems, multiple kernels run on top of a hypervisor, a lightweight operating system presenting guest operating systems with a virtual operating platform. Each virtual machine requires its own kernel, and thus is significantly more resource intensive than a container (see [fig (a)](#fig-a1)). On the other hand, virtual machines can use different kernels on the same machine.

    Running multiple kernels on the same machine is useful for instance in applications that require different device drivers. Device drivers are part of the kernel. Hypervisors typically expose physical devices to virtual machines. If a machine has 2 GPUs for instance, it can host one virtual machine per GPU, and each machine can use a different device driver for the GPU.

    <figure markdown="span">
      ![Comparison of virtual machines and containers](resources/container_VM_comparison.png){ align=left width="800" }

      <figcaption id="fig-a1">Fig. (a): A diagram of a container engine running 2 container, and a virtualisation server running the same workload on 2 virtual machines. Note that the virtualisation server must run 2 kernel instances even if the 2 machines use the exact same kernel.</figcaption>
    </figure>

## Example: chroot jail environment

Before working on containers for HPC systems let's explore one of the first technologies that allowed the creation of isolated process environments, the `chroot` command. You should run this example in a machine where you have super user privileges (`sudo`), as the command required `sudo` access.

1. Create the root of the `chroot` file system.

  ```bash
  mkdir --parents ${HOME}/jail/{bin,lib,lib64,home/jailuser}
  ```

1. Copy the executables that should be available in the `chroot` system.

  ```bash
  for binary in /bin/bash /bin/sh /bin/ls /bin/cat; do
    cp ${binary} ${HOME}/jail/bin/
  done; unset binary
  cp $(realpath /bin/which) ${HOME}/jail/bin/which
  ```

1. Copy the linker.

  ```bash
  for linker in /lib64/ld-linux-x86-64.so*; do
    cp ${linker} ${HOME}/jail/lib64/
  done; unset linker
  ```

1. Copy the libraries required by the executables.

  ```bash
  while IFS="" read -r library; do
    cp ${library} ${HOME}/jail/lib/
  done < <(ldd /bin/bash /bin/sh /bin/ls /bin/cat $(/bin/which) | grep -E '=>' | awk 'BEGIN {FS="(=>)|( +)"} {print $4}' | sort | uniq); unset library
  ```

1. Create a text file to test the executables.

  ```bash
  echo 'Welcome to chroot jail!' > ${HOME}/jail/home/jailuser/hello.txt
  ```

1. Your isolated environment is now created in `${HOME}/jail`. To use the chroot jail environment, simply call `chroot` on the root of the jail.

  ```console
  sudo chroot ${HOME}/jail /bin/bash
  ```

Any program running in the isolated environment cannot name and thus cannot access paths outside the environment. Thus the term _chroot jail_ is often used to describe an isolated environment created with `chroot`.

!!! question "Using a chroot jail"

    === "Question"

        Can you list the content of the home directory of the chroot jail (`/home/jailuser`), and display the contents of `hello.txt` which is located in the home directory? What are the contents of root directory viewed from inside the jailed system? Can you use the `which` command to locate the `ls` executable inside the chroot jail, and to which file does this executable correspond in the host system?

    === "Answer"

        To  list the contents of the home directory of the chroot jail (`/home/jailuser`), and display the contents of `hello.txt` which is located in the home directory, change to the chroot jail, and use the `ls` and `cat` commands whose executables were copied into the jail during step 2.

            ```console
            $ sudo chroot ${HOME}/jail /bin/bash
            bash-5.1# ls /home/myusername/          
            hello.txt
            bash-5.1# cat /home/myusername/hello.txt 
            Welcome to chroot jail!
            ```

        The contents of the chroot jail root are displayed using the `ls` command.

            ```console
            bash-5.1# ls -lahF /
            total 0
            drwxr-xr-x 1 0 0  30 Jun 14 20:07 ./
            drwxr-xr-x 1 0 0  30 Jun 14 20:07 ../
            drwxr-xr-x 1 0 0  66 Jun 14 20:11 bin/
            drwxr-xr-x 1 0 0  16 Jun 14 20:07 home/
            drwxr-xr-x 1 0 0 126 Jun 14 20:07 lib/
            drwxr-xr-x 1 0 0  40 Jun 14 20:07 lib64/
            ```

        The path for `ls` reported from inside the chroot jail is the following.

        ```console
        bash-5.1# which ls
        /bin/ls
        ```

        This corresponds to the following path in the host system.

        ```text
        ${HOME}/jail/bin/ls
        ```

The chroot jail we created follows the filesystem hierarchy standard inside the `${HOME}/jail` path, and this allows programs to operate as if they were operating in a system with `${HOME}/jail` as root.

The `chroot` executable is still used to create containers. However, chroot jails require root permission which is not available in share systems. Modern containers use namespaces to control the access to resources, and with the introduction of user namespaces no special privileges are required for most operations in containers.

### Accessing host system resources

In a typical system you may want to access resources other than simple executables, like network and storage devices. Linux operating systems expose such resources in the root path using pseudo filesystems like _sysfs_ and _procfs_. In many Linux operating systems,

- `/proc` contains the process information pseudo-filesystem,
- `/sys` is a mount point for the sysfs filesystem that provides information about the kernel like `/proc`, but better structured,
- `/dev` contains special or device files that refer to physical devices, and
- `/run` is a directory that contains information about processes that describes the system state.

You can use `man hier` to see more of the resource you may need to access inside an isolated environment. You can mount any resource files using bindings:

```bash
for fs in dev dev/pts proc sys sys/firmware/efi/efivars run; do
  sudo mount --bind "/${fs}" "${system_root}/${fs}"
done; unset fs
```

However, this process involves a lot of manual work, and exposes system components to the jailed system. Containers automate this process and use systems such as namespace and cgroup to expose resources in a more secure manner.

## Differences between containers and virtual environments

Containers and virtual environments both modify the environment where applications run. However, they generally differ in the manner the modification is made and thus the aspects of the environment they can modify. Virtual environment often rely on features of the operating system, such as the `PATH` variable used by executables and the `LD_LIBRARY_PATH` used by the operating system dynamic linker. If your application needs a specific version of `glibc` for instance, then it is much safer to use a container that overrides the default linker (`/lib64/ld-linux-x86-64.so.2`) and `glibc` (`/lib64/libc.so.6`).
