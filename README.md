# pyGrappolo

## Steps
1. Clone the repository with `git clone https://github.com/neffrw/pyGrappolo.git --recurse-submodules`
    1. Run `pip install ./pyGrappolo`
    1. If you are running MacOS, ensure you are using a compiler with OpenMP support with `CC=<compiler> CXX=<compiler> pip install ./pyGrappolo`
    1. If you change the root folder from pyGrappolo to another name, change `./pyGrappolo` to that folder name instead.
1. In your python program, run `import grappolo`
    1. `grappolo.grappolo(args_dict)` (see test.py example)

## Troubleshooting
- If your build is failing on MacOS due to an `ld: unsupported tapi file type '!tapi-tbd' in yaml file` error within conda, your conda-supplied `ld` may be outdated compared to the system `ld`.
    - Please either run outside of the conda environment, or append the system's `ld` to the PATH environment variable after activating your conda environment with `export PATH=${PATH}:<path/to/system/ld>`.
    - _Hint: the system's `ld` location can be found by executing `which ld` outside of the conda environment._
