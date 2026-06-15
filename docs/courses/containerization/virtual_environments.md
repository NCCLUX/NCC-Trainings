# Virtual environments

Environments typically rely on features of the operating system to modify the environment where processes run.

- The [`PATH` variable](https://pubs.opengroup.org/onlinepubs/9799919799/) is defined in the [POSIX standard](https://pubs.opengroup.org/onlinepubs/9799919799/basedefs/V1_chap01.html) of operating system interface and environment, as the sequence of path prefixes that functions and utilities apply in searching for an executable file.

- Similarly, the `LD_LIBRARY_PATH` is a list of paths where the dynamic linker (see [`ld.so`](https://man7.org/linux/man-pages/man8/ld.so.8.html)) searches for shared libraries that are linked dynamically to executables at runtime.

Environment systems typically install a partial filesystem hierarchy in some path, and prepend the `\bin` and `\lib` subpataths of the hierarchy to the `PATH` and `LD_LIBRARY_PATH` lists respectively. Thus, the programs and libraries available of the environment are added to the ones available in the `\bin` and `\lib` of the root directory.

## Filesystem hierarchy for isolated environments

The Linux operating systems rely on the [Filesystem Hierarchy Standard](https://refspecs.linuxfoundation.org/FHS_3.0/fhs-3.0.html), that places programs, drivers, and their components in specific paths. For instance,

- `\bin` contains commands that may be used by both the system administrator and by users,
- `\lib` contains those shared library images needed to run the commands in the root filesystem, and
- `\home` contains user home directories.

Every isolated environment and containerization system creates a file system hierarchy in some path of the file system that forms the _root_ of the isolated environment. You can use

- `man hier` to access a succinct description of the filesystem hierarchy of your system, and
- `man file-hierarchy` to access the systemd file system hierarchy requirements.

The latter is useful in container environment that need to support system services.

## Example: Python virtual environment (`venv`)

Start by loading the latest software release, and from the software release load the Python module to access the Python executable.

```console
module load  env/release/2025.1
module load Python
```

You can now use the python `venv` module to create and setup the virtual environment as a non-privileged user.

1. Create the virtual environment.

  ```console
  python -m venv ${HOME}/environments/env-numpy
  ```

1. Activate the virtual environment.

  ```console
  source ${HOME}/environments/env-numpy/bin/activate
  ```

1. Install the required packages. For instance, in this example we install `numpy` and `pyyaml`.

  ```console
  (env-numpy) $ pip install numpy pyyaml
  ```

After installing the packages the environment is ready for use. In future sessions, load the environment with the command

```console
source ${HOME}/environments/env-numpy/bin/activate
```

and unload the environment with the command:

```console
deactivate
```

To remove the environment, simply delete the environment directory.

```console
rm -r ${HOME}/environments/env-numpy
```

!!! question "Virtual environment setup and contents"

    === "Question"

        What are the contents of the virtual environment directory? Can you detect any of the path to the virtual environment directory in the `PATH` and `LD_LIBRARY_PATH` variables?

    === "Answer"

        The contents of the virtual environment directory follow the filesystem hierarchy structure.

        ```console
        $ ls -lahF $VIRTUAL_ENV
        total 28K
        drwxr-x---  5 u101064 hpcusers 4.0K Jun 14 18:12 ./
        drwx------ 15 u101064 hpcusers 4.0K Jun 14 18:12 ../
        drwxr-x---  2 u101064 hpcusers 4.0K Jun 14 18:12 bin/
        -rw-r-----  1 u101064 hpcusers   69 Jun 14 18:12 .gitignore
        drwxr-x---  3 u101064 hpcusers 4.0K Jun 14 18:12 include/
        drwxr-x---  3 u101064 hpcusers 4.0K Jun 14 18:12 lib/
        lrwxrwxrwx  1 u101064 hpcusers    3 Jun 14 18:12 lib64 -> lib/
        -rw-r-----  1 u101064 hpcusers  381 Jun 14 18:12 pyvenv.cfg
        ```
        The virtual environment sets the convenient variable `VIRTUAL_ENV` that points to the root of the virtual environment, `${HOME}/environments/env-numpy`. The `PATH` variable points to `$VIRTUAL_ENV/bin`.

        ```console
        $ echo $PATH | tr ':' '\n'
        /home/users/u101064/env-numpy/bin
        /apps/USE/easybuild/release/2025.1/software/Python/3.13.1-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/OpenSSL/3/bin
        /apps/USE/easybuild/release/2025.1/software/XZ/5.6.3-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/SQLite/3.47.2-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/Tcl/8.6.16-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/ncurses/6.5-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/bzip2/1.0.8-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/binutils/2.42-GCCcore-14.2.0/bin
        /apps/USE/easybuild/release/2025.1/software/GCCcore/14.2.0/bin
        /home/users/u101064/.local/bin
        /home/users/u101064/bin
        /usr/local/bin
        /usr/bin
        /usr/local/sbin
        /usr/sbin
        ```

        However, the `lib` directory of the environment does not appear in the `LD_LIBRARY_PATH`. The environment is a Python environment, and all external libraries are loaded by Python executable of the environment and not the system linker, so setting up the linker environment variable, `LD_LIBRARY_PATH`, is redundant. The Python executable in the environment,

        ```text
        /home/users/u101064/env-numpy/bin/python
        ```

        will look into the `pyvenv.cfg` at the root of the virtual environment to determine the location where Python will look for share libraries. The shared libraries of the environment are still stored in `$VIRTUAL_ENV/lib`. For instance the mathematical libraries used by NumPy are located in the following directory.

        ```console
        $ ls -lahF $VIRTUAL_ENV/lib/python3.13/site-packages/numpy.libs
        total 28M
        drwxr-x--- 2 u101064 hpcusers 4.0K Jun 14 18:12 ./
        drwxr-x--- 7 u101064 hpcusers 4.0K Jun 14 18:12 ../
        -rwxr-x--x 1 u101064 hpcusers 2.8M Jun 14 18:12 libgfortran-040039e1-0352e75f.so.5.0.0*
        -rwxr-x--x 1 u101064 hpcusers 246K Jun 14 18:12 libquadmath-96973f99-934c22de.so.0.0.0*
        -rwxr-x--x 1 u101064 hpcusers  25M Jun 14 18:12 libscipy_openblas64_-32a4b2a6.so*
        ```

        The get a list of paths where the Python interpreter is looking for shared libraries, you can use the following command.

        ```console
        python -c "import sys; from pathlib import PureWindowsPath; print('\n'.join(PureWindowsPath(p).as_posix() for p in sys.path))"
        ```
