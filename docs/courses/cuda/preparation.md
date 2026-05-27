## 1. Setting up the connection

Please read the [instructions](https://docs.lxp.lu/first-steps/quick_start/) on how to get access to the MeluXina machine.

Then read the instructions on how to connect.

- [Windows users](https://docs.lxp.lu/first-steps/connecting/)
- [Linux/Mac users](https://docs.lxp.lu/first-steps/connecting/)

## 2. Use your username to connect to MeluXina

For example, the below example shows the user of `u100490`.

```console
ssh u100490@login.lxp.lu -p 8822
```

If an alias `meluxina` is set up:

```console
ssh meluxina
```

## 3. Check your access to your home and project directories

Once you have logged in, you will be in a default home directory.

```console
[u100490@login02 ~]$ pwd
/home/users/u100490
```

After that, go to the project directory.

```console
[u100490@login02 ~]$ cd /project/home/p201350
[u100490@login02 p201350]$ pwd
/project/home/p201350
```

## 4. Create a working folders

Create your own working folder under the project directory

For example, here is the process for user `u100490`.

```console
[u100490@login02 p201350]$ mkdir ${USER}
```

Now move into the working directory. For example, for user `u100490`:

```console
[u100490@login02 p201350]$ cd u100490
```

## 6. Prepare the course material

Copy the folder which has examples and source files to your working directory. For example, the user home directory `u100490` executes:

```console
[u100490@login03 u100490]$ cp -r /project/home/p201350/CUDA .
[u100490@login03 u100490]$ cd CUDA/
[u100490@login03 CUDA]$ pwd
/project/home/p201350/u100490/CUDA
[u100490@login03 CUDA]$ ls -lthr
total 20K
-rw-r-----. 1 u100490 p201350   51 Mar 13 15:50 module.sh
drwxr-s---. 2 u100490 p201350 4.0K Mar 13 15:50 Vector-addition
drwxr-s---. 2 u100490 p201350 4.0K Mar 13 15:50 Unified-memory
...
...
```

## 7. Reserve a compute node

Until now, you are in the login node; now it is time to do the dry run test. Reserve the interactive node for running/testing CUDA applications.

```console
salloc -A p201350 --partition=gpu --qos default -N 1 -t 01:00:00
```

??? "check if your reservation is allocated"
    ```
    [u100490@login03 ~]$ salloc -A p201350 --partition=gpu --qos default -N 1 -t 01:00:00
    salloc: Pending job allocation 296848
    salloc: job 296848 queued and waiting for resources
    salloc: job 296848 has been allocated resources
    salloc: Granted job allocation 296848
    salloc: Waiting for resource configuration
    salloc: Nodes mel2131 are ready for job
    ```

You can also check if you got the interactive node for your computations. For example, for the user `u100490`:

```console
[u100490@mel2131 ~]$ squeue -u u100490
            JOBID PARTITION     NAME     USER    ACCOUNT    STATE       TIME   TIME_LIMIT  NODES NODELIST(REASON)
           304381       gpu interact  u100490    p201350  RUNNING       0:37     01:00:00      1 mel2131
```

## 8. Accessing the CUDA examples

Now we need to check that a simple CUDA application is working. Go to folder `Dry-run-test`.

```console
[u100490@login03 CUDA]$ cd Dry-run-test/
[u100490@login03 Dry-run-test]$ ls
Hello-world.cu  module.sh
```

## 9. Loading the compilers

We need to load the compiler to test the GPU CUDA codes. We need a Nvidia HPC SDK compiler for compiling and testing CUDA code.

```console
module load env/staging/2023.1
module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-11.7.0
export NVCC_APPEND_FLAGS='-allow-unsupported-compiler'
```

We also provide a script to simplify the process for you.

```console
source module.sh
```

??? "check if the module is loaded properly"
    ```console
    [u100490@mel2131 ~]$ module load env/staging/2023.1
    [u100490@mel2131 ~]$ module load OpenMPI/4.1.5-NVHPC-23.7-CUDA-11.7.0
    [u100490@mel2131 ~]$ export NVCC_APPEND_FLAGS='-allow-unsupported-compiler'
    [u100490@mel2131 ~]$ module list

    Currently Loaded Modules:
    1) env/release/2022.1           (S)   6) numactl/2.0.14-GCCcore-11.3.0  11) libpciaccess/0.16-GCCcore-11.3.0  16) GDRCopy/2.3-GCCcore-11.3.0                  21) knem/1.1.4.90-GCCcore-11.3.0
    2) lxp-tools/myquota/0.3.1      (S)   7) CUDA/11.7.0                    12) hwloc/2.7.1-GCCcore-11.3.0        17) UCX-CUDA/1.13.1-GCCcore-11.3.0-CUDA-11.7.0  22) OpenMPI/4.1.5-NVHPC-23.7-CUDA-11.7.0
    3) GCCcore/11.3.0                     8) NVHPC/23.7-CUDA-11.7.0         13) OpenSSL/1.1                       18) libfabric/1.15.1-GCCcore-11.3.0
    4) zlib/1.2.12-GCCcore-11.3.0         9) XZ/5.2.5-GCCcore-11.3.0        14) libevent/2.1.12-GCCcore-11.3.0    19) PMIx/4.2.2-GCCcore-11.3.0
    5) binutils/2.38-GCCcore-11.3.0      10) libxml2/2.9.13-GCCcore-11.3.0  15) UCX/1.13.1-GCCcore-11.3.0         20) xpmem/2.6.5-36-GCCcore-11.3.0

    Where:
        S:  Module is Sticky, requires --force to unload or purge
    ```

## 10. Compile and test a simple CUDA application

Please compile and test a simple CUDA application. For example, in `Dry-run-test`:

```console
# compilation
$ nvcc -arch=compute_70 Hello-world.cu -o Hello-World-GPU

# execution
$ ./Hello-World-GPU
Hello World from GPU!
Hello World from GPU!
Hello World from GPU!
Hello World from GPU!
```

## 11. Check that you can reserve a node

Similarly, for the hands-on session, we need reserve a node. For example:

```console
salloc -A p201350 --partition=gpu --qos default -N 1 -t 02:15:00
```

??? "check if your reservation is allocated"
    ```
    [u100490@login03 ~]$ salloc -A p201350 --partition=gpu --qos default -N 1 -t 02:15:00
    salloc: Pending job allocation 296848
    salloc: job 296848 queued and waiting for resources
    salloc: job 296848 has been allocated resources
    salloc: Granted job allocation 296848
    salloc: Waiting for resource configuration
    salloc: Nodes mel2131 are ready for job
    ```

## 12. Check that you are ready to access the examples

We will continue with our Hands-on exercise. For example, in the `Hello World` example, we do the following steps:

```console
[u100490@mel2063 CUDA]$ pwd
/project/home/p201350/u100490/CUDA
[u100490@mel2063 CUDA]$ ls
Dry-run-test  Matrix-multiplication  Profiling      Unified-memory
Hello-world   module.sh              Shared-memory  Vector-addition
[u100490@mel2063 CUDA]$ source module.sh
[u100490@mel2063 CUDA]$ cd Hello-world
# compilation
[u100490@mel2063 CUDA]$ nvcc -arch=compute_70 Hello-world.cu -o Hello-World-GPU

# execution
[u100490@mel2063 CUDA]$ ./Hello-World-GPU
Hello World from GPU
```
