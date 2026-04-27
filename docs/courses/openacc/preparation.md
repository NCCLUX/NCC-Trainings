#### 1. [How to login to MeluXina machine](https://docs.lxp.lu/first-steps/quick_start/)
- 1.1 [Please take a look if you are using Windows](https://docs.lxp.lu/first-steps/connecting/)
- 1.2 [Please take a look if you are using Linux/Mac](https://docs.lxp.lu/first-steps/connecting/)

#### 2. Use your username to connect to MeluXina
- 2.1 For example, the below example shows the user of `u100490` 
  ```
  $ ssh u100490@login.lxp.lu -p 8822
  ### or
  $ ssh meluxina 
  ```
#### 3. Once you have logged in
- 3.1 Once you have logged in, you will be in a default home directory 
  ```
  [u100490@login02 ~]$ pwd
  /home/users/u100490
  ```
- 3.2 After that, go to the project directory.
  ```
  [u100490@login02 ~]$ cd /project/home/p201030
  [u100490@login02 p201030]$ pwd
  /project/home/p201030
  ```
  
#### 4. And please create your own working folder under the project directory
- 4.1 For example, here is the user with `u100490`:
  ```
  [u100490@login02 p201030]$ mkdir $USER
  ### or 
  [u100490@login02 p201030]$ mkdir u100490  
  ```
#### 5. Now it is time to move into your home directory
- 5.1 For example, with user home directory `u100490` 
  ```
  [u100490@login02 p201030]$cd u100490
  ```

#### 6. Now it is time to copy the folder which has examples and source files to your home directory
- 6.1 For example, with user home directory `u100490`
  ```
  [u100490@login03 u100490]$ cp -r /project/home/p201030/OpenACC .
  [u100490@login03 u100490]$ cd OpenACC/
  [u100490@login03 OpenACC]$ pwd
  /project/home/p201030/u100490/OpenACC
  [u100490@login03 OpenACC]$ ls -lthr
  total 20K
  -rw-r-----. 1 u100490 p201030   51 Mar 13 15:50 module.sh
  drwxr-s---. 2 u100490 p201030 4.0K Mar 13 15:50 Vector-addition
  drwxr-s---. 2 u100490 p201030 4.0K Mar 13 15:50 Unified-memory
  ...
  ...
  ```
#### 7. Until now, you are in the login node; now it is time to do the dry run test
- 7.1 Reserve the interactive node for running/testing OpenACC applications 
  ```
  $ salloc -A p201030 --res gpu-openacc-morning --partition=gpu --qos default -N 1 -t 01:00:00
  ```
  
    ??? "check if your reservation is allocated"
        ```
        [u100490@login03 ~]$ salloc -A p201030 --res gpu-openacc-morning --partition=gpu --qos default -N 1 -t 01:00:00
        salloc: Pending job allocation 296848
        salloc: job 296848 queued and waiting for resources
        salloc: job 296848 has been allocated resources
        salloc: Granted job allocation 296848
        salloc: Waiting for resource configuration
        salloc: Nodes mel2131 are ready for job
        ```
     
 - 7.2 You can also check if you got the interactive node for your computations, for example, here with the user `u100490`:
 ```
 [u100490@mel2131 ~]$ squeue -u u100490
             JOBID PARTITION     NAME     USER    ACCOUNT    STATE       TIME   TIME_LIMIT  NODES NODELIST(REASON)
            304381       gpu interact  u100490    p201030  RUNNING       0:37     01:00:00      1 mel2131
 ```

#### 8. Now we need to check the simple OpenACC application if that is going to work for you:
 - 8.1 Go to folder `Dry-run-test`
```
[u100490@login03 OpenACC]$ cd Dry-run-test/
[u100490@login03 Dry-run-test]$ ls 
Hello-world.cu  module.sh
```

#### 9. Finally, we need to load the compiler to test the GPU OpenACC codes
 - 9.1 We need a Nvidia HPC SDK compiler for compiling and testing OpenACC code
 ```
 $ module load OpenMPI/5.0.3-NVHPC-24.9-CUDA-12.6.0
 ### or
 $ source module.sh
 ```

    ??? "check if the module is loaded properly"
    	```
        [u100490@mel2129 OpenACC]$ module load OpenMPI/5.0.3-NVHPC-24.9-CUDA-12.6.0
        [u100490@mel2129 OpenACC]$ ml
        
        Currently Loaded Modules:
        1) env/release/2024.1            (S)   7) CUDA/12.6.0                         13) OpenSSL/3                                   19) PMIx/5.0.2-GCCcore-13.3.0
        2) lxp-tools/myquota/0.3.3       (S)   8) NVHPC/24.9-CUDA-12.6.0              14) libevent/2.1.12-GCCcore-13.3.0              20) PRRTE/3.0.5-GCCcore-13.3.0
        3) GCCcore/13.3.0                      9) XZ/5.4.5-GCCcore-13.3.0             15) UCX/1.16.0-GCCcore-13.3.0                   21) UCC/1.3.0-GCCcore-13.3.0
        4) zlib/1.3.1-GCCcore-13.3.0          10) libxml2/2.12.7-GCCcore-13.3.0       16) GDRCopy/2.4.1-GCCcore-13.3.0                22) NCCL/2.22.3-GCCcore-13.3.0-CUDA-12.6.0
        5) binutils/2.42-GCCcore-13.3.0       11) libpciaccess/0.18.1-GCCcore-13.3.0  17) UCX-CUDA/1.16.0-GCCcore-13.3.0-CUDA-12.6.0  23) UCC-CUDA/1.3.0-GCCcore-13.3.0-CUDA-12.6.0
        6) numactl/2.0.18-GCCcore-13.3.0      12) hwloc/2.10.0-GCCcore-13.3.0         18) libfabric/1.21.0-GCCcore-13.3.0             24) OpenMPI/5.0.3-NVHPC-24.9-CUDA-12.6.0
        Where:
        S:  Module is Sticky, requires --force to unload or purge
        ```


#### 10. Please compile and test your CUDA application 
 - 10.1 For example, Dry-run-test
 ```
 // compilation
 $ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world.c -o Hello-World-GPU

 // execution
 $ ./Hello-World-GPU

 // output
 $ Hello World from GPU!
   Hello World from GPU!
   Hello World from GPU!
   Hello World from GPU!
 ```

#### 11. Similarly, for the hands-on session, we need to do the node reservation:
 - 10.1 For example, reservation
 ```
 $ salloc -A p201030 --res gpu-openacc-afternoon --partition=gpu --qos default -N 1 -t 02:30:00
 ```

    ??? "check if your reservation is allocated"
        ```
        [u100490@login03 ~]$ salloc -A p201030 --res gpu-openacc-afternoon --partition=gpu --qos default -N 1 -t 02:30:00
        salloc: Pending job allocation 296848
        salloc: job 296848 queued and waiting for resources
        salloc: job 296848 has been allocated resources
        salloc: Granted job allocation 296848
        salloc: Waiting for resource configuration
        salloc: Nodes mel2131 are ready for job
        ```

#### 12. We will continue with our Hands-on exercise
 - 12.1 For example, in the `Hello World` example, we do the following steps:
 ```
 [u100490@mel2063 OpenACC]$ pwd
 /project/home/p201030/u100490/OpenACC
 [u100490@mel2063 OpenACC]$ ls
 [u100490@mel2063 OpenACC]$ ls
 Dry-run-test  Matrix-multiplication  Profiling      Unified-memory
 Hello-world   module.sh              Vector-addition
 [u100490@mel2063 OpenACC]$ source module.sh
 [u100490@mel2063 OpenACC]$ cd Hello-world
 // compilation
 [u100490@mel2063 OpenACC]$ nvc -fast -acc=gpu -gpu=cc80 -Minfo=accel Hello-world.c -o Hello-World-GPU

 // execution
 [u100490@mel2063 OpenACC]$ ./Hello-World-GPU

 // output
 [u100490@mel2063 OpenACC]$ Hello World from GPU
 ```
