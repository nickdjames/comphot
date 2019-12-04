# comphot
Comet photometry developed for the BAA comet section.

## Building
### Dependencies
#### Boost
Download latest from http://www.boost.org

Follow the included documentation for building on your platform, e.g. for Linux:
```
./bootstrap.sh
sudo ./b2 install
```

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
Comphot is aimed at multiple platforms and therefore we use CMake to generate the appropriate build configuration. To get started follow the steps above to build and install the required dependencies, then fork and clone this repo:
```
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
To build on other platforms or other build systems you can specify as required, e.g. for Windows using [Visual Studio Build Tools](https://visualstudio.microsoft.com/downloads/):
```
cmake -G "Visual Studio 14 2015" ../
cmake --build . --config Release
```

## Packaging
CMake includes CPack which can be used to package comphot for distribution, e.g. to produce a Windows installer:
```
cpack -G WIX -C Release
```
