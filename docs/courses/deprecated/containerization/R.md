# Managing packages in R

The R program has a built-in package manager. Assuming that you have access to an installation of R, the R system contains 2 utilities that allow the installation of packages in 3 different modes. First, there is the built-in package manager that can instal packages

- in system wide accessible locations for packages that should be available to all users (requires elevated privileges), or
- in user specific location for packages that are accessible to the current user only.

There are default locations where the built-in package manager searches for packages. The user specific locations take precedence over system locations. The package manager search path can be extended by the user to include bespoke locations. There is also the Packrat package manager which installs packages

- in project directories, with the packages being available in an environment isolated within the project directory.

The Packrat package manager is available as an R package. When creating an environment within a project directory, the environment is activated automatically when starting R in the project directory (but not in its subdirectories due to the implementation of Packrat).

In your local system you can install packages in any mode. In the HPC systems, you can only install packages in the user accessible location, so you are limited to user and project wide installations. Nevertheless, the HPC installation of R includes a number of commonly used packages, such as `dbplyr` and `tidyverse`. You should check if the package you require is installed and that the installed version provides the functionality you need before installing any packages locally. Remember, local package installations consume space and inodes against personal or project quota.

## Installing R packages locally and globally

Be default R installs packages system wide. When R detects that it does not have write access to the system directories it suggests installing packages for the current user only.

Start and interactive session and then load the R module and start R:
```bash
$ module load lang/R
$ R
```
You can list the directories where R is installing and looking for new packages using the function `.libPaths()`
```R
> .libPaths()
[1] "/mnt/irisgpfs/apps/resif/aion/2020b/epyc/software/R/4.0.5-foss-2020b/lib64/R/library"
```
If you haven't installed any libraries, only the system path appears in the path where R is looking for libraries. Now, try installing for instance the Packrat package globally with the `install.packages` command.
```R
> install.packages(c("packrat"))
Warning in install.packages(c("packrat")) :
  'lib = "/mnt/irisgpfs/apps/resif/aion/2020b/epyc/software/R/4.0.5-foss-2020b/lib64/R/library"' is not writable
Would you like to use a personal library instead? (yes/No/cancel) yes
Would you like to create a personal library
‘~/R/x86_64-pc-linux-gnu-library/4.0’
to install packages into? (yes/No/cancel) yes
--- Please select a CRAN mirror for use in this session ---
Secure CRAN mirrors
```
Select any mirror apart from `1: 0-Cloud [https]`; usually mirrors closer to your physical location will provide better bandwidth. After selecting a mirror the download and installation of the package proceeds automatically.

Note that after failing to install the package in the system directory, R creates an installation directory for the user in their home directory `~/R/x86_64-pc-linux-gnu-library/4.0` and installs the package for the user only. After the installation, you can check the path where R is looking for packages again.
```R
> .libPaths()
[1] "/mnt/irisgpfs/users/<user name>/R/x86_64-pc-linux-gnu-library/4.0"
[2] "/mnt/irisgpfs/apps/resif/aion/2020b/epyc/software/R/4.0.5-foss-2020b/lib64/R/library"
```
Now R will look for packages in the user directory first (`/mnt/irisgpfs` is another path for `/home` that appears in the `${HOME}` variable). Note by the naming convention that R uses when it creates the directory for installing user packages, you can have multiple minor versions of R installed, and their packages will not interfere with each other. For instance,

- R version 4.0.5 installs packages in `~/R/x86_64-pc-linux-gnu-library/4.0`, and
- R version 4.3.2 installs packages in `~/R/x86_64-pc-linux-gnu-library/4.3`.

Some useful commands for managing packages are,

- `installed.packages()` to list installed packages and various information regarding each package installation,
- `old.packages()` to list outdated packages,
- `update.packages()` to update installed packages, and
- `remove.packages(c("packrat"))` to remove packages.

To list the loaded packages, use the command
```R
search()
```
and to get a detailed description of the environment, use the command
```R
sessionInfo()
```
which provides information about the version of R, the OS, and loaded packages.

To load a library that has been installed use the command `library`. For instance,
```R
library(packrat)
```
where you cam notice that the use of quotes is optional and only a single can be loaded at a time. The `library` function causes an error when the loading of a package fails, so R provides the function `require` which returns the status of the package loading operation in a return variable, and is design for use inside R functions.

_Useful resources_

- [R Packages: A Beginner's Tutorial](https://www.datacamp.com/tutorial/r-packages-guide)
- [Efficient R programming: Efficient set-up](https://bookdown.org/csgillespie/efficientR/set-up.html)

## Configuring installation paths in R

So far we have only used the default installation paths of R. However, in a local installation where the user has rights to install in the system directories (e.g. in a Conda environment with R) the user installation directory is not created automatically. Open an R session in an interactive session in the HPC cluster or in your personal machine. To get the location where user packages are installed call
```R
> Sys.getenv("R_LIBS_USER")
[1] "/home/<user name>/R/x86_64-conda-linux-gnu-library/4.3"
```
which will print an environment variable, `R_LIBS_USER`, which is set by R and stores the default location for storing user packages. If you create the directory with
```bash
$ mkdir -p /home/<user name>/R/x86_64-conda-linux-gnu-library/4.3
```
then you can print the locations where R is searching for packages (after reloading R), and the default location should appear first in the list. For instance for o Conda installation of R using the Micromamba package manager, the paths printed are
```R
> .libPaths()
[1] "/home/<user name>/R/x86_64-conda-linux-gnu-library/4.3"
[2] "/home/<user name>/micromamba/envs/R/lib/R/library"
```
where R is installed in a Conda environment named `R` in the second entry of the search path.

There are now multiple locations where packages are stored. The location used by default is the first in the list. Thus, after creating the default location for user installed packages, packages are installed by default in user wide mode. For instance, installing the Packrat package,
```R
> install.packages(c("packrat"))
```
listing the user installation directory
```bash
$ ls "/home/<user name>/R/x86_64-conda-linux-gnu-library/4.3"
packrat
```
will show the directory with the installed Packrat package files. To install the package in a system wide installation, use the `lib`flag
```R
> install.packages(c("packrat"), lib="/home/<user name>/micromamba/envs/R/lib/R/library")
```
to specify the installation location. During loading, all directories in the path are searched consecutively until the package is located.

The package installation paths can also be used to maintain multiple independent environments in R. For instance, you can maintain a personal environment and project environment for your research group. Lets consider the case where you want the create an environment in a project directory. First, create a directory for the R environment
```bash
$ mkdir -p "${PROJECTHOME}<project name>/R-environment"
```
where the variable `PROJECTHOME` is defined in the UL HPC system environment to point to the home of the project directories (and includes a trailing slash '/'). To install a package in the project environment, call the installation function with the appropriate `lib` argument
```R
> install.packages( c("packrat"), lib=paste0( Sys.getenv("PROJECTHOME"), "<project name>/", "R-environment" ) )
```
and follow the typical instructions. To load the package, you now must also specify the location of the library,
```R
> library( packrat, lib.loc=paste0( Sys.getenv("PROJECTHOME"), "<project name>/", "R-environment" ) )
```
similar to the installation. Environment options can be used to extent the library paths and avoid having to specify the library path in each command.

A startup file mechanism is provided by R to set up user and project wide environment options. There are 2 kinds of file,

- `.Renviron` files used to set-up environment variables for R, and
- `.Rprofile` files used to run any R code during initialization.

Note that `.Renviron` files are simply a list of
```
key=value
```
assignment pairs which are read by R, not proper bash code (adding an `export` modifier is a syntax error). There are 2 locations where startup files appear,

- the home directory, `~/.Renviron` and `~/.Rprofile`, for user wide settings, and
- project directories for project wide settings.

The definitions in project `.Rprofile` files override the user wide definitions in `~/.Rprofile`. The definitions in `.Renviron` files supersede the definitions in `~/.Renviron`, that is if the project has an environment file, the user wide definitions are ignored. Note that R is designed to source setup files at the directory where R starts, and any setup files in parent or descendent directories are ignored.

Both the profile and environment startup files can setup a user wide environment. For instance, to use an environment setup in the project directories of the UL HPC systems add in the user wide environment setup file, `~/.Renviron`, the entry
```
R_LIBS=${PROJECTHOME}<project name>/R-environment
```
and then reload R. The new library path is
```R
> .libPaths()
[1] "/mnt/irisgpfs/projects/<project name>/R-environment"
[2] "/mnt/irisgpfs/users/<user name>/R/x86_64-pc-linux-gnu-library/4.0"
[3] "/mnt/irisgpfs/apps/resif/iris-rhel8/2020b/broadwell/software/R/4.0.5-foss-2020b/lib64/R/library"
```
assuming that all directories appearing in the path exist. Note that the setup file options precede any default options.

We can also use startup files to setup project wide libraries. For instance, assume that we are working on a project in a directory named `project` and the R packages are stored in a subdirectory `R-environment`. We use a project profile, to still be able to use any library paths defined in the user wide environment file. Add in a file `project/.Rprofile` the following definitions,
```R
project_path <- paste0( getwd(), "/R-environment" )
newpaths <- c( project_path, .libPaths() )
.libPaths( newpaths )
```
and then start R in the `project` directory. The new library path is
```R
> .libPaths()
[1] "/mnt/irisgpfs/users/<user name>/Documents/project/R-environment"
[2] "/mnt/irisgpfs/projects/<project name>/R-environment"
[3] "/mnt/irisgpfs/users/<user name>/R/x86_64-pc-linux-gnu-library/4.0"
[4] "/mnt/irisgpfs/apps/resif/iris-rhel8/2020b/broadwell/software/R/4.0.5-foss-2020b/lib64/R/library"
```
were the local project settings override the user and system wide settings. This is effectively a local project environment.

## Installing packages in R project directories with Packrat

The Packrat library is used to automate the creation and management of project based environments. Packrat also automates operations such as tracking the version of the packages installed in the environment with snapshots, and saving the snapshot information in a text file that can be version controlled. The R distribution available through the UL HPC modules has a fairly old version of Packrat, which nevertheless supports all the basic features. Packrat is a light package, so you can install a more modern version in a user wide mode or in some environment accessible to all the users of a UL HPC project.

To initialize the project, for instance in the directory `~/Documents/project`, use the commands:
```R
library(packrat)
packrat::init("~/Document/project")
```
The initialization command creates,
- a directory `~/Document/project/packrat` to store the packages, and
- a setup script `~/Document/project/.Rprofile` to initialize the project.
Therefore, start R within the project directory `~/Document/packrat`, to activate the project environment. After initializing the project or whenever you start R in the project directory, the `packrat` directory and its subdirectories will be the only ones appearing in the library paths:
```R
> .libPaths()
[1] "/mnt/irisgpfs/users/<user name>/Documents/project/packrat/lib/x86_64-pc-linux-gnu/4.0.5"
[2] "/mnt/irisgpfs/users/<user name>/Documents/project/packrat/lib-ext/x86_64-pc-linux-gnu/4.0.5"
[3] "/mnt/irisgpfs/users/<user name>/Documents/project/packrat/lib-R/x86_64-pc-linux-gnu/4.0.5"
```
Execute all package operations as usual. For instance, to install the `plyr` package, use the command:
```R
> install.packages(c("plyr"))
```
All packages are stored in the `packrat` subdirectory of the project.

Packrat stores the status of the project in the file `packrat/packrat.lock`. This file stores the precise package versions that were used to satisfy dependencies, including dependencies of dependencies, and should not be edited by hand. After any change in the installed packages run the command
```R
packrat::snapshot()
```
to update the file. You can use the command
```R
packrat::status()
```
to analyze the code in the project directory and get a report regarding the status of extraneous or missing packages. After running the `status` command, you can run
```R
packrat::clean()
```
to remove any unused packages. Finally, after restoring the `packrat/packrat.lock` file from a version control system, or if `status` detects a missing package, use the command
```R
packrat::restore()
```
to install any missing packages.

_Useful resources_

- [Official Packrat tutorial](https://rstudio.github.io/packrat/walkthrough.html)

## Issues with managing packages with the native R package managers

The native package manager of R is quite potent, and there are packages such as Packrat that further extend its capabilities. However, there are some drawbacks in installing packages with the native tools. Consider for instance installing the `hdf5r` package, a package used to read and write binary files, that is quite popular in HPC engineering applications. The installation mode is not important for our demonstration purposes, but assume that you are performing a user wide installation.
```R
> install.packages(c("hdf5r"))
```

During the installation, you can see that R is compiling the package components. This can be advantageous is the compilation process is tailored to optimize the build for the underlying system configuration. If you use the module available in the UL HPC systems, it is configured to use the main components of the FOSS tool chain (you can see that by calling `module list` after loading R), so the compiled packages are well optimized.

**N.B.** If you encounter any issues with missing packages load the whole FOSS tool chain module with the command,
```bash
module load toolchain/foss
```
as there are a few popular packages missing in the dependencies of R.

However, if you want to avoid compiling packages from source, which can be quite time consuming, you can use binary distributions of R. These include the distributions provided though native package managers in various Linux distributions, like APT and YUM, as well as Conda package managers like Mamba.



