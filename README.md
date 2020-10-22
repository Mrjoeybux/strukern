### Requirements
- Cmake
- ZLIB
- pybind11

### Compiling
- If you use a virtual python environment, make sure to activate it.
- Clone repo: `git clone --recurse-submodules -j4 git@github.com:Mrjoeybux/strukern.git`
- Let `<STRUKERN_ROOT>` be the path to the root directory of the strukern project.



##### Build Gedlib
- Navigate to `gedlib` directory i.e. `cd <STRUKERN_ROOT>/src/lib/gedlib`
- Exceute the following command `python install.py --lib gxl`

##### Build Dlib
- Navigate to `dlib` directory i.e. `cd <STRUKERN_ROOT>/src/lib/dlib`
- Execute the following commands

```
mkdir build
cd build
cmake ..
make
```

##### Build Strukern
- Navigate to `src` directory i.e. `cd <STRUKERN_ROOT>/src/`
- Exceute `python setup.py develop`