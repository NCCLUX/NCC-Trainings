## 1. Login to MeluXina machine

You can find instructions on how to access the MeluXina machine in [MeluXina documentation](https://docs.lxp.lu/first-steps/quick_start/).

Details for:

- [Windows users](https://docs.lxp.lu/first-steps/connecting/), and
- [Linux/Mac users](https://docs.lxp.lu/first-steps/connecting/).


## 2. Use your username to connect to MeluXina

For exmaple the below example shows the user of `u100490`
```
$ ssh u100490@login.lxp.lu -p 8822
```

## 3. Ensure that you can access your home and project directories

### 3.1 Check your home directory

Once you have logged in, you will be in a default home directory. Check with the following command.
 
```
[u100490@login02 ~]$ pwd
/home/users/u100490
```

### 3.2 Check the bootcamp project directory

Change to the project directory for Nvidia Bootcamp activities.

```
[u100490@login02 ~]$ cd /project/home/p200117
[u100490@login02 p200117]$ pwd
/project/home/p200117
```
  
## 4. Create a personal working folder in the bootcamp project directory

Please create your own working folder under the project directory. For example, here it is the process for user with user name `u100490`.

```
[u100490@login02 p200117]$ mkdir "${USER}"
```

### 4.1 Copy the files required by the practicals in your working directory

Now copy `climate.simg` and `climate.sh` from project directory to your user directory. For example, here is the process for user `u100490`.

```
[u100490@login02 p200117]$ cp /project/home/p200117/climate.simg /project/home/p200117/u100490
[u100490@login02 p200117]$ cp /project/home/p200117/climate.sh /project/home/p200117/u100490
```

### 4.2 Copy the files required by the practicals in your home directory

Copy `cfd.simg` and `cfd.sh` from project directory to your user directory. For example, here is the process for user `u100490`.

```
[u100490@login02 p200117]$ cp /project/home/p200117/cfd.simg /project/home/p200117/u100490
[u100490@login02 p200117]$ cp /project/home/p200117/cfd.sh /project/home/p200117/u100490
```

### 4.3 Check that all required file are available

Go to your home directory and check if all the necessary files are there (`.simg` and `.sh`). For example, here is the process for user `u100490`.

```
[u100490@login02 p200117]$ cd u100490
[u100490@login02 u100490]$ pwd
[u100490@login02 u100490]$ /project/home/p200117/u100490
[u100490@login02 u100490]$ ls -lthr
total 15G
-rw-r-x---. 1 u100490 p200117  736 Feb  8 18:59 climate.sh
-rwxr-x---. 1 u100490 p200117 7.2G Feb  8 19:19 climate.simg
-rwxr-x---. 1 u100490 p200117 6.9G Feb  8 19:21 cfd.simg
-rw-r-x---. 1 u100490 p200117  723 Feb  8 19:21 cfd.sh
```

## 5. For the dry run (9th February from 11:30-12:30)

Follow the in instruction below step-by-step.

- Launch the service:

  ```
  [u100490@login02 u100490]$ salloc -A p200117 --res gpudev -q dev -N 1 -t 01:00:0
  [u100490@mel2123 u100490]$ mkdir -p $PROJECT/$USER/workspace-climate
  [u100490@mel2123 u100490]$ module load Singularity-CE/3.10.2-GCCcore-11.3.0
  
  [u100490@mel2123 u100490]$ singularity run --bind $PROJECT/$USER $PROJECT/$USER/climate.simg cp -rT /workspace $PROJECT/$USER/workspace-climate
  INFO:    Converting SIF file to temporary sandbox...
  INFO:    Cleaning up image...
  [u100490@mel2123 u100490]$ singularity run --nv --bind $PROJECT/$USER $PROJECT/$USER/climate.simg jupyter lab --notebook-dir=$PROJECT/$USER/workspace-climate/python/jupyter_notebook --port=8888 --ip=0.0.0.0 --no-browser --NotebookApp.token=""
  INFO:    Converting SIF file to temporary sandbox...
  WARNING: underlay of /usr/bin/nvidia-smi required more than 50 (452) bind mounts
  [W 10:10:32.723 LabApp] All authentication is disabled.  Anyone who can connect to this server will be able to run code.
  [I 10:10:33.043 LabApp] jupyter_tensorboard extension loaded.
  [I 10:10:33.047 LabApp] JupyterLab extension loaded from /usr/local/lib/python3.8/dist-packages/jupyterlab
  [I 10:10:33.047 LabApp] JupyterLab application directory is /usr/local/share/jupyter/lab
  [I 10:10:33.048 LabApp] [Jupytext Server Extension] NotebookApp.contents_manager_class is (a subclass of) jupytext.TextFileContentsManager already - OK
  [I 10:10:33.048 LabApp] Serving notebooks from local directory: /mnt/tier2/project/p200117/u100490/workspace-climate/python/jupyter_notebook
  [I 10:10:33.048 LabApp] Jupyter Notebook 6.2.0 is running at:
  [I 10:10:33.048 LabApp] http://hostname:8888/
  [I 10:10:33.049 LabApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
  ```

- Open a new terminal on your local computer, and again login to MeluXina to access the service.

- Use the `squeue` command to determine the node running the service. Nodes are listed in the `NODELIST` field.

- Connect to the service. For instance, if the service is running on node `mel2123` then use the command:

  ```
  ssh -L8080:mel2123:8888 u100490@login.lxp.lu -p 8822
  ```

- Keep those terminals open/alive (please do not close them).

- Copy and paste `localhost:8080` to your browser either to Chrome or Firefox.

  ```
  http://localhost:8080
  ```

You should now have access to the service.
  
## 6. For the afternoon session (9th and 10th February)

If have missed the dry run session, then please go through the steps from 1-4.

- Now it is time to edit your batch script (`climate.sh`) before launching your Jupyter notebook. Please follow the following steps.
  ```
  [u100490@login02 u100490]$ emacs(emacs -nw)/vim climate.sh
  #!/bin/bash -l
  #SBATCH --partition=gpu 
  #SBATCH --ntasks=1
  #SBATCH --nodes=1    
  ############  day one ##########
  #######SBATCH --time=02:00:00         ## use this option for day one
  #######SBATCH --res ai_bootcamp_day1   ## use this option for day one
  ################################
   
  ############  day two ##########
  #SBATCH --time=03:30:00         ## use this option for day two
  #SBATCH --res ai_bootcamp_day2  ## use this option for day two
  ################################
  #SBATCH -A p200117
  #SBATCH --qos default
  
  mkdir -p $PROJECT/$USER/workspace-climate
  module load Singularity-CE/3.10.2-GCCcore-11.3.0
  
  singularity run --bind $PROJECT/$USER $PROJECT/$USER/climate.simg cp -rT /workspace $PROJECT/$USER/workspace-climate
  singularity run --nv --bind $PROJECT/$USER $PROJECT/$USER/climate.simg jupyter lab --notebook-dir=$PROJECT/$USER/workspace-climate/python/jupyter_notebook --port=8888 --ip=0.0.0.0 --no-browser --NotebookApp.token=""
  ```

- Once you have modified your `climate.sh`, please launch your batch script as below.

  ```
  [u100490@login03 u100490]$ sbatch climate.sh
  Submitted batch job 276009
  [u100490@login03 u100490]$ squeue 
  JOBID PARTITION     NAME     USER    ACCOUNT    STATE       TIME   TIME_LIMIT  NODES NODELIST(REASON)
  276009       gpu climate.  u100490    p200117  RUNNING       0:16        20:00      1 mel2077
  ```

- You have now initiated your singularity container which will help you to open the Jupyter nootebook.

  ```
  [u100490@login03 u100490]$ ls -lthr
  total 7.2G
  -rwxr-x---. 1 u100490 p200117 7.2G Feb  3 14:53 climate.simg
  -rw-r-----. 1 u100490 p200117  613 Feb  3 17:06 climate.sh
  -rw-r-x---. 1 u100490 p200117  724 Feb  8 19:41 cfd.sh
  -rwxr-x---. 1 u100490 p200117 6.9G Feb  8 19:42 cfd.simg
  -rw-r--r--. 1 u100490 p200117 1.1K Feb  3 17:58 slurm-276009.out
  ```

- You can also check if everything is OK by executing the command below.

  ```
  [u100490@login03 u100490]$ head -30 slurm-276009.out 
  INFO:    Converting SIF file to temporary sandbox...
  INFO:    Cleaning up image...
  INFO:    Converting SIF file to temporary sandbox...
  WARNING: underlay of /usr/bin/nvidia-smi required more than 50 (452) bind mounts
  [W 17:58:37.489 LabApp] All authentication is disabled.  Anyone who can connect to this server will be able to run code.
  [I 17:58:37.807 LabApp] jupyter_tensorboard extension loaded.
  [I 17:58:37.811 LabApp] JupyterLab extension loaded from /usr/local/lib/python3.8/dist-packages/jupyterlab
  [I 17:58:37.811 LabApp] JupyterLab application directory is /usr/local/share/jupyter/lab
  [I 17:58:37.813 LabApp] [Jupytext Server Extension] NotebookApp.contents_manager_class is (a subclass of) jupytext.TextFileContentsManager already - OK
  [I 17:58:37.813 LabApp] Serving notebooks from local directory: /mnt/tier2/project/p200117/u100490/workspace-climate/python/jupyter_notebook
  [I 17:58:37.813 LabApp] Jupyter Notebook 6.2.0 is running at:
  [I 17:58:37.813 LabApp] http://hostname:8888/
  [I 17:58:37.813 LabApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
  ```

  If everything went OK your output should be similar.

- Open a new terminal on your local computer, and again login to MeluXina to access the service.

- Use the `squeue` command to determine the node running the service. Nodes are listed in the `NODELIST` field.

- Connect to the service. For instance, if the service is running on node `mel2077` then use the command:

  ```
  ssh -L8080:mel2077:8888 u100490@login.lxp.lu -p 8822
  ```

- Keep those terminals open/alive (please do not close them).

- Copy and paste `localhost:8080` to your browser either to Chrome or Firefox.

  ```
  http://localhost:8080
  ```
