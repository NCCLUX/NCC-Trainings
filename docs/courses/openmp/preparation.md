## 1. Setting up the connetion

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
[u100490@login02 ~]$ cd /project/home/p200898
[u100490@login02 p200898]$ pwd
/project/home/p200898
```

## 4. Create a working folders

Create your own working folder under the project directory

For example, here is the process for user `u100490`.

```console
[u100490@login02 p200898]$ mkdir ${USER}
```
  
Now move into the working directory. For example, for user `u100490`:

```console
[u100490@login02 p200898]$ cd u100490
```

## 6. Prepare the course material

Copy the folder which has examples and source files to your working directory. For example, the user home directory `u100490` executes:

```console
[u100490@login03 u100490]$ cp -r /project/home/p200898/OpenMP .
[u100490@login03 u100490]$ cd OpenMP/
[u100490@login03 OpenMP]$ pwd
/project/home/p200898/u100490/OpenMP
[u100490@login03 OpenMP]$ ls -lthr
drwxr-s---. 2 u100490 p200898 4.0K May 27 21:42 Data-Sharing-Attribute
drwxr-s---. 2 u100490 p200898 4.0K May 28 00:35 Parallel-Region
drwxr-s---. 2 u100490 p200898 4.0K May 30 18:26 Dry-run-test
...
...
```

## 7. Reserve a compute node

Until now, you are in the login node; now it is time to do the dry run test. Reserve the interactive node for running/testing OpenMP applications.

```console
salloc -A p200898 --res ncc-openmp --partition=cpu --qos default -N 1 -t 01:00:00
```

??? "check if your reservation is allocated"
    ```
    [u100490@login03 ~]$ salloc -A p200898 --res ncc-openmp --partition=cpu --qos default -N 1 -t 01:00:00
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
           304381       cpu interact  u100490    p200898  RUNNING       0:37     01:00:00      1 mel2131
```

## 8. Accessing the OpneMP examples

Now we need to check that a simple OpenMP application is working. Go to folder `Dry-run-test`.

```console
[u100490@login03 OpenMP]$ cd Dry-run-test/
[u100490@login03 Dry-run-test]$ ls
source.sh  Test.cc  Test.f90
```

## 9. Loading the compilers

We need to load the compiler to test our OpenMP codes. We will work with GNU compiler.

``` console
source module.sh
```

??? "check if the module is loaded properly"
    ```
    [u100490@mel2131 ~]$ module list

    currently Loaded Modules:
    1) env/release/2022.1                (S)  19) libpciaccess/0.16-GCCcore-11.3.0    37) jbigkit/2.1-GCCcore-11.3.0        55) VTune/2022.3.0                          73) NSS/3.79-GCCcore-11.3.0
    2) lxp-tools/myquota/0.3.1           (S)  20) X11/20220504-GCCcore-11.3.0         38) gzip/1.12-GCCcore-11.3.0          56) numactl/2.0.14-GCCcore-11.3.0           74) snappy/1.1.9-GCCcore-11.3.0
    3) GCCcore/11.3.0                         21) Arm-Forge/22.0.4-GCC-11.3.0         39) lz4/1.9.3-GCCcore-11.3.0          57) hwloc/2.7.1-GCCcore-11.3.0              75) JasPer/2.0.33-GCCcore-11.3.0
    4) zlib/1.2.12-GCCcore-11.3.0             22) libglvnd/1.4.0-GCCcore-11.3.0       40) zstd/1.5.2-GCCcore-11.3.0         58) OpenSSL/1.1                             76) nodejs/16.15.1-GCCcore-11.3.0
    5) binutils/2.38-GCCcore-11.3.0           23) AMD-uProf/3.6.449                   41) libdeflate/1.10-GCCcore-11.3.0    59) libevent/2.1.12-GCCcore-11.3.0          77) Qt5/5.15.5-GCCcore-11.3.0
    6) ncurses/6.3-GCCcore-11.3.0             24) Advisor/2022.1.0                    42) LibTIFF/4.3.0-GCCcore-11.3.0      60) UCX/1.13.1-GCCcore-11.3.0               78) CubeGUI/4.7-GCCcore-11.3.0
    Where:
        S:  Module is Sticky, requires --force to unload or purge
    ```

## 10. Compile and test a simple OpenMP application

Please compile and test your OpenMP application. For example, in `Dry-run-test`:

```console
# compilation (C/C++)
$ g++ Test.cc -fopenmp

# compilation (FORTRAN)
$ gfortran Test.f90 -fopenmp

# execution
$ ./a.out
Hello world from the master thread
Hello world from thread id Hello world from thread id Hello world from thread
id Hello world from thread id Hello world from thread id 4 from the team size of
1 from the team size of 20 from the team size of  from the team size of 555
```

## 11. Check that you can reserve a node

Similarly, for the hands-on session, we need reserve a node. For example:

```console
salloc -A p200898 --res ncc-openmp --partition=cpu --qos default -N 1 -t 02:30:00
```

??? "check if your reservation is allocated"
    ```console
    [u100490@login03 ~]$ salloc -A p200898 --res ncc-openmp --partition=cpu --qos default -N 1 -t 02:30:00
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
[u100490@mel2063 OpenMP]$ pwd
/project/home/p200898/u100490/OpenMP
[u100490@mel2063 OpenMP]$ ls
[u100490@mel2063 OpenMP]$ ls
drwxr-s---. 2 u100490 p200898 4.0K May 27 21:42 Data-Sharing-Attribute
drwxr-s---  2 u100490 p200898 4.0K May 28 00:35 Parallel-Region
drwxr-s---  2 u100490 p200898 4.0K May 28 23:45 Worksharing-Constructs-Schedule
drwxr-s---. 2 u100490 p200898 4.0K May 29 00:57 Worksharing-Constructs-Other
drwxr-s---. 2 u100490 p200898 4.0K May 29 18:07 Worksharing-Constructs-Loop
drwxr-s---. 2 u100490 p200898 4.0K May 30 18:25 SIMD-Others
drwxr-s---. 2 u100490 p200898 4.0K May 30 18:37 Dry-run-test
-rw-r-----  1 u100490 p200898  241 May 30 18:41 module.sh
[u100490@mel2063 OpenMP]$ source module.sh
```
