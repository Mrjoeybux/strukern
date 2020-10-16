- Clone repo
- Let `<STRUKERN_ROOT>` be the path to the root directory of the strukern project.

### Requirements
- Cmake
- ZLIB


#### Build Gedlib
- Navigate to `gedlib` directory i.e. `cd <STRUKERN_ROOT>/src/lib/gedlib`
- Exceute the following command `python install.py --lib gxl`

#### Build Dlib
- Navigate to `dlib` directory i.e. `cd <STRUKERN_ROOT>/src/lib/dlib`
- Execute the following commands
```
cd dlib/dlib
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=1 ..
make
sudo make install
```