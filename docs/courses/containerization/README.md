# 1. Theory
**Main goal**:
	* 1. Understand Cgroups, Chroot and Namespaces in Linux and their role
	* 2. Differences between containers and environments. Use cases
	* 3. Differentces between VM and containers. Containers performance
	* Explain differences between Docker and Singularity. Check documentation 
## Cgroups
* Links: [1. Red Hat](https://docs.redhat.com/en/documentation/red_hat_enterprise_linux/7/html/resource_management_guide/chap-introduction_to_control_groups) [2. kernel.org](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html)
* *The control groups, abbreviated as cgroups in this guide, are a Linux kernel feature that allows you to allocate resources — such as CPU time, system memory, network bandwidth, or combinations of these resources — among hierarchically ordered groups of processes running on a system*
* *cgroup is a mechanism to organize processes hierarchically and distribute system resources along the hierarchy in a controlled and configurable manner.*
* Using Cgroups we can assign hardware resources to particular OS processes
* **systemd** provides a unified hierarchy to the cgroup tree
	* We can check **systemd** processes using `systemd-cgtop`. Three different unit types are distinguished
		* **service**: process or group of processes that have been started by systemd through a unit configuration file. The process encansulapted can be started and stopped as one set. Named as `name.service`.
		* **scope**: externally created processes, started and stopped through the `fork()` function, and registered by **systemd** at runtime, User sessions, containers and virtual machines are treated as scopes. Named as `name.scope`
		* **slice**: group of hierarchically organized units. They organize a hierarchy in which scopes and services are placed, they do not contain processes. If the name a slice look as `parent-name.slice` is a subslice of the `parent.slice`
	* Services, scopes and slices are created manually by the system administrator or dinamically by programs
	* There are four slices created by default:
		* **-.slice**: root slice
		* **system.slice**: default place for system services
		* **user.slice**: default place for all user sessions
		* **machine.slice**: default place for virtual machines and linux containers
### Resource controllers in Linux Kernel
* They represent single resources
* By default, the Linux kernel provides a range of resource controllers that are mounted automatically by **systemd**
* List of currently mounted resource controllers: `/proc/cgroups`
* List of available resources:
	* `blkio`: implements the block io controller by stablishing various kinds if IO control policies. It sets limits on input/output access to and from block devices
	* `cpu`: uses CPU scheduler to provide access to the CPU
	* `cpuacct`: create automatic reports on CPU resources used by tasks in a cgroup. To accounf for the CPU usage of the group tasks
	Both `cpu` and `cpuacct` controllers are mounted together
	* `cpuset`: assign individual CPUs and memory nodes to tasks in a CPU (assuming a mutlicore system)
	* `devices`: allows or denies access to devices for tasks in a cgroup. Implemented to track and enforce open and mknod restrictions on device files
	* `freezer`: suspends or resumes tasks in a cgroup. Useful to barch job managemetn system
	* `memory`: sets limits on memory use by tasks in a cgroup and generates automatic reports on memory resources used by those tasks
	* `net_cls`: tag network packets with `classid`. This allows the Linux traffic controller to identify packets from a particular cgroup task
		* The `net_filter`: subsystem allows the Linux firewall to identify packets originating from a particular cgroup tasks
	* `perf_event`: enables monitoring cgroups with the perf tool
	* `hugetlb`: allows to use virtual memory pages of larger sizes and to enforce resource limits on these pages
## Chroot
* Links: [1, man7](https://www.man7.org/linux/man-pages/man2/chroot.2.html), [2, Wikipedia](https://en.wikipedia.org/wiki/Chroot), [3, ArchWiki](https://wiki.archlinux.org/title/Chroot)
*  *It is a shell command and a system call on Unix and Unix-like operating systems that changes the apparent root directory for the current running process and its children.*
*  A program that run in this modified environment cannot name files outside the directory tree
*  The modified environment is called a **chroot jail**
*  We need a Linux installatiion or installation media with the same ISA as the system being *chrooted* into
## Namespaces
* Links: [1. Linux handbook](https://linuxhandbook.com/namespaces/), [2. Wikipedia](https://en.wikipedia.org/wiki/Linux_namespaces)
* Emerged in the early 2000s to address the growing need for resource isolation in multi-user systems
* They provide isolation which, together with granular resource management (cgroups), allows for the creation of lightweight and fast deployable containers
### Example (Not clear lsns)
* We run the `unshare` command in Linux. It is used to create new namespaces
```bash
sudo unshare --uts /bin/bash
```
We can change the hostname within this new shell. This change will be isolated to this UTS namespace
```bash
hostname isolated-box
exec bash
```
Using `lsns` (list information about the currently accessible namespaces or a specific namespace), I assume that these are the processes that we're looking for:
```bash
4026531832 mnt       209    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531833 net       165    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531834 time      234    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531835 cgroup    234    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531836 pid       191    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531837 user      142    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531838 uts       234    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
4026531839 ipc       186    7896 jdelguerrero /usr/bin/dbus-broker-launch --scope user
```
Executing `lsns` on new namespace:
```bash
4026531832 mnt       509       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531833 net       492       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531834 time      566       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531835 cgroup    566       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531836 pid       523       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531837 user      470       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531838 uts       549       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
4026531839 ipc       518       1 root            /usr/lib/systemd/systemd --system --deserialize=148 rhgb
```
* Two critical behaviors to understand about namespaces:
	* **Inheritance**: a process creates a new child process using `fork()`, the child inherits a copy of all its parent's namespaces
	* **Lifecycle**: A namespace is automatically erased by the kernel when the last process terminates of leaves it
### Architecture
* *hierarchical*: resources of new contexts are related to the context of the already created system. Namespaces are isolated from each other, but they can be mapped in such a way that the main context can know that other contexts are running
* *non-hierarchical*: does not related resources in different contexts
### PID Namespace
* It isolates the process ID number space
* Processes in different PID namespaces can have the same PID
* For instance, PID 1 is the first process started by the kernel
* Containers should look like a fresh Linux system, so it needs its own PID 1. This can be solved by isolating PID numbering
* Example
```
sudo unshare --pid --fork --mount-proc /bin/bash

ps aux
USER         PID %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND
root           1  0.1  0.0 233996  6328 pts/5    S    16:35   0:00 /bin/bash
root          49  0.0  0.0 234392  4556 pts/5    R+   16:35   0:00 ps aux

```
 Only two processes are seen the bash shell (PID 1) and the ps command (PID 49)
 * Inside our namespace, PID is 1
 * External PID
```
ps aux | grep 'unshare.*pid'
root     1193574  0.0  0.0 246128 10808 pts/8    S+   16:35   0:00 sudo unshare --pid --fork --mount-proc /bin/bash
root     1193623  0.0  0.0 246128  3184 pts/5    Ss   16:35   0:00 sudo unshare --pid --fork --mount-proc /bin/bash
root     1193624  0.0  0.0 230368  2104 pts/5    S    16:35   0:00 unshare --pid --fork --mount-proc /bin/bash
```
### Network Namespace
* Provide network isolation, i.e. a completely isolated network stack: network interfaces, IP addresses, routing tables, socket listings, tracking tables and firewall rules
```
sudo unshare --net /bin/bash
ip -br link list
lo               DOWN           00:00:00:00:00:00 <LOOPBACK> 
```
### Mount Namespace
* One of the most important and oldest namespaces
* It shapes the view of the filesystem. For instance, you run a process inside a container, the paths you see under `/` are not the same ones that the host sees
* This namespace allows each container to have its own filesystem layout that appears self contained
```
sudo unshare --mount /bin/bash
mkdir /mnt/test
mount --bind /etc /mnt/test
findmnt
[...]

/mnt/test  /dev/mapper/fedoravol-osbtrfs[/root/etc] btrfs               rw,relatime,seclabel,compress=zstd:1,ssd,space_cache=v2,subvolid=261,subvol=/root
```
* If we run `findmnt` in a regular shell, we won't find any trace of the mount point 
* Containerization software create mount namespaces so that they can build custom paths from chains of real directories, overlay filesystems, etc, all without affecting the host
### User Namespace
* Deternines who you are inside the filesystem by isolating user and group IDs
* As a consequence, processes inside Linux can be run as root inside the namespace, although there are not root privileges on the host
* It stablishes a mapping between internal and external UIDs
* As in the PID case, the internal UID is the one the process sees inside the namespace, and the external UID is the one that the kernel associates with the process on the host
*  To create a new user namespace `unshare -r /bin/bash`
*  We can check user id by executing `id`
*  When a container declares that it is running as a root, it is referred to the internal root identity. For the host kernel, the process is seen as an unprivileged user
*  If a user namespace owns a mount namespace, the processes insie can mount filesystems without host privileges
### Cgroup Namespace
* Isolates views of control groups and provides a private view of the cgroup hierarchy of a process
* Therefore, containers expect to manage their own cgroups without interfering with the host or other containers
* Processes inside the namespace cannot observe or influence cgroups outside their namespace.
* Create cgroup namespace:
```
sudo unshare --cgroup /bin/bash
cat /proc/self/cgroup
0::/
```
* It does not create or destroy cgroups, it simply controls visibility
### IPC Namespace
* IPC: Inter-Process Communication
* Isolates System V IPC objects and POSIX message queues; all of them mechanisms for communication between  processes via shared memory segments, semaphore arrays and message queues
* *Without isolation*: all processes can see each other's IPC objects and potentially interfere with them
* *With IPC Namespace*: only processes within the same namespace can access these objects
* Create memory segment in the host
```
ipcmk -M 1024
ipcs -m
----- Shared Memory Segments --------
key        shmid      owner      perms      bytes      nattch     status      
0x1c170818 1212470    jdelguerre 644        1024       0                       
```
* Create new IPC namespace and check for memory segments:
```
sudo unshare --ipc /bin/bash
ipcs -m
----- Shared Memory Segments --------
key        shmid      owner      perms      bytes      nattch     status      

```
### Time Namespace
* More recent addition to the Kernel (5.6, 2020)
* It allows for per-namespace offsets to the system clock
* Not widely used
### Building a container
#### 1. Create a root filesystem
We use [Busybox](https://man.archlinux.org/man/busybox.1.en)
```bash
ROOTFS="$HOME/rootfs"
mkdir -p "$ROOTFS"/{bin,proc,sys,dev}
cp /usr/bin/busybox "$ROOTFS/bin/"
 for cmd in sh mount umount ls mkdir ps ping hostname; do
    ln -sf busybox "$ROOTFS/bin/$cmd"
done
```
 We're setting up a minimal filesystem with essential binaries
 #### 2. Start a new set of namespaces
 New shell where we are operating in fresh namespaces
```bash
 sudo unshare --mount --pid --uts --ipc --net --fork /bin/bash 
```
#### 3. Prepare the mount namespace
* The `--bind` parameter allows us  to mount a directory on another directory.
*  `-t` stands for filesystem type

```
 mount --bind "$ROOTFS/dev" "$ROOTFS/dev"
mount -t proc proc "$ROOTFS/proc"
 mount -t sysfs sys "$ROOTFS/sys"
```

* Set hostname `hostname isolated-box`
* Bring up loopback networking for basic networking inside the namespace
```
ip link set lo up
```

* Enter the new root filesystem:
```
ROOTFS=/home/jdelguerrero/rootfs
chroot "$ROOTFS" /bin/sh
```
* We have an isolated filesystem environment with `/proc`, `/sys`, and `/dev` properly mounted
#### 4. Start PID 1 inside the namespace
* Exit the isolated filesystem
* Execute:
```
exec /bin/bash
ps aux
USER         PID %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND
root           1  0.0  0.0 234008  6464 ?        S    10:43   0:00 /bin/bash
root         365  0.0  0.0 234392  4812 ?        R+   11:48   0:00 ps aux
```
### A toy example with Chroot
Follow the example presented in [containers tutorial](https://ncclux.github.io/NCC-Trainings/courses/containerization/containers/#container-technologies)
#### 1. Create the root of the `chroot` file system

```
mkdir --parents ${HOME}/jail/{bin,lib,lib64,home/myusername}
```
#### 2. Copy the executables that should be available in the `chroot`system:

```
for binary in /bin/bash /bin/ls /bin/cat; do
  cp ${binary} ${HOME}/jail/bin/
done; unset binary
```
* ELF is a format for storing binaries, libraries and core dumps on disks in Linux and Unix-based systems
* Structure:
![Screenshot_20260529_173224.png](:/0546a4b230034a82bddf32b9f9b4a76b)

* Check  [What Is an ELF File?](https://www.baeldung.com/linux/executable-and-linkable-format-file)
#### 3. Copy the linkers.
```
for linker in /lib64/ld-linux-x86-64.so*; do
  cp ${linker} ${HOME}/jail/lib64/
done; unset linker
```

* The programs `ld.so` and `ld-linux.so*` find and load the shared libreries needed by a program, prepare the program to run, and then run it
* Linux binaries require dynamic linking (unless the `-static` option was given during compilation)
* `ld.so` handles `a.out` binarires
* `ld-linux.so*` handles binaries that are in the ELF format
* ELF file for an executable: it always has a part called the Program Header table. It provides information that is necessary while running the executable
* Dynamic linker tasks:
	1. Locating shared libraries
	2. Loading shared libraries into memory
	3. Resolving program's undefined symbols using the symbols in shared libraries
	4. Assembling the program so that the process can call the functions in the shared libraries at run-time
* To list the dynamic linkers:

```
 readelf -l /usr/bin/ls | head -20
 ldd /usr/bin/ls
```
* Links [1, Baeldung](https://www.baeldung.com/linux/dynamic-linker), [2, Linux manual page](https://www.systutorials.com/linux-manual-page-8-ld-linux/)
#### 4. Copy the libraries required by the executables
```
while IFS="" read -r library; do
  cp ${library} ${HOME}/jail/lib/
done < <(ldd /bin/bash /bin/ls /bin/cat | grep -E '=>' | awk 'BEGIN {FS="(=>)|( +)"} {print $4}' | sort | uniq); unset library

```
#### 5. Create a text file to test the executables
```
echo 'Welcome to chroot jail!' > ${HOME}/jail/home/myusername/hello.txt
```

## Building a container
* References [1, Shipping Your Machine: Building a Container in 50 Lines of Code](https://dev.to/yechielk/shipping-your-machine-building-a-container-in-60-lines-of-code-part-1-14ma), [2, Build Your Own Docker with Linux Namespaces, cgroups, and chroot: Hands-on Guide ](https://akashrajpurohit.com/blog/build-your-own-docker-with-linux-namespaces-cgroups-and-chroot-handson-guide/)
* 
