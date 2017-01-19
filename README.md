# comphot
Comet photometry developed for the BAA comet section.

## Building
### Dependencies
#### cfitsio:
Download latest from http://heasarc.gsfc.nasa.gov/fitsio/
e.g. http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3410.tar.gz

Follow the included documentation for building on your platform, e.g. for Linux:
```
./configure --prefix=/usr/local
make shared
sudo make install
```

#### libgd:
Download latest from https://github.com/libgd/libgd/releases/
e.g. https://github.com/libgd/libgd/releases/download/gd-2.2.4/libgd-2.2.4.tar.gz

Follow the included documentation for building on your platform, e.g. for Linux:
```
./configure
make
sudo make install
```

### comphot
comphot is aimed at multiple platforms and therefore we use CMake to generate the appropriate build configuration. To get started follow the steps above to build and install the required dependencies, then fork and clone this repo:
````
git clone https://github.com/<username>/comphot.git
cd comphot
```
To build, create a clean directory, generate the build configuration and perform the build, e.g. for Linux:
```
mkdir build
cd build
cmake ../
make
```
Other platforms and other build systems are available using this more general approach:
```
cmake -G <generator> ../
cmake --build . --config <debug|release|...>
```
