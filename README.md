# pyGrappolo

## Steps
1. Clone the repository with `git clone https://github.com/neffrw/pyGrappolo.git --recurse-submodules`
    1. Run `pip install ./pyGrappolo`
    1. If you are running MacOS, ensure you are using a compiler with OpenMP support with `CC=<compiler> CXX=<compiler> pip install ./pyGrappolo`
    1. If you change the root folder from pyGrappolo to another name, change `./pyGrappolo` to that folder name instead.
1. In your python program, run `import grappolo`
    1. `grappolo.grappolo(args_dict)` (see test.py example)
