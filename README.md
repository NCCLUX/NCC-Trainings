# NCC Trainings

This is the source of the [Supercomputing Luxembourg](https://supercomputing.lu/) trainings.

The courses can be provided to the [EVITA](https://www.eurohpc-ju.europa.eu/research-innovation/our-projects/evita_en) project. Please have a look on the [EVITA project template](https://github.com/ENCCS/evita-material-template/tree/main) and make sure your course format is compatible.

## Setup instructions

To setup the python environment you can use uv or any compatible package manager. The compilation is tested with Pthon 3.13.5. To setup the compilation environment with uv, run in the root directory of the repository:

```console
uv venv --python $(cat .python) ~/environments/ncc
source ~/environments/ncc/bin/activate
uv pip install --requirements requirements.txt
```

### Compiling the images in resource directories

There is code for LaTeX generated figures in some sections. The source code is located in directories named `resources`. The compilation of the code requires [TeX Live](https://tug.org/texlive/) and the [`convert`](https://linux.die.net/man/1/convert) tool. Compile the images manually by running the command

```console
make
```

inside the `resource` directories.

