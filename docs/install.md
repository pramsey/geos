---
layout: page
title: Install
nav_order: 4
permalink: /install/
---

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

## Download Source

| Release Date | Release | Download Link |
| ------- | ------ | ------ |
|  2020/12/10 | **3.9.0** | [geos-3.9.0.tar.bz2](http://download.osgeo.org/geos/geos-3.9.0.tar.bz2) |
|  2020/03/10 | **3.8.1** |  [geos-3.8.1.tar.bz2](http://download.osgeo.org/geos/geos-3.8.1.tar.bz2) |
|  2019/10/04 | **3.7.3** |  [geos-3.7.3.tar.bz2](http://download.osgeo.org/geos/geos-3.7.3.tar.bz2) |
|  2020/12/11 | **3.6.5** |  [geos-3.6.5.tar.bz2](http://download.osgeo.org/geos/geos-3.6.5.tar.bz2) |
|  2019/10/04 | **3.5.2** |  [geos-3.5.2.tar.bz2](http://download.osgeo.org/geos/geos-3.5.2.tar.bz2) |


## Install Packages

### Red Hat

There is a GEOS package in the EPEL (Extra Packages for Enterprise Linux) repository.

```bash
# Add the EPEL repository
yum -y install epel-release

# Install the GEOS runtime and development packages
rpm -Uvh geos geos-devel
```

### Ubuntu

The [Ubuntu GIS](https://wiki.ubuntu.com/UbuntuGIS) project maintains a collection of repositories with builds of common open source geospatial projects, including GEOS.

```bash
# Add the Ubuntu GIS PPA
sudo apt-get install python-software-properties
sudo add-apt-repository ppa:ubuntugis/ppa

# Install the packages

```

### Debian

The [Debian GIS](https://wiki.debian.org/DebianGis) project maintains [GEOS packages](https://tracker.debian.org/pkg/geos) and pushes them into the appropriate Debian respositories.

```bash
sudo apt-get install geos
```

### Amazon Linux

Amazon Linux is based on RH7, and can read from the EPEL repository. To enable using Amazon tools, use the `amazon-linux-extras` utility.

```bash
sudo yum install -y amazon-linux-extras
sudo amazon-linux-extras enable epel
sudo yum search geos
sudo yum install geos geos-devel
```

### Homebrew

For MacOS, GEOS can be installed using the [Homebrew](https://brew.sh/) package repository, which downloads source packages and builds them in place using a recipe to ensure all packages integrate with each other nicely.

First [install Homebrew](https://brew.sh/). Then:

```
brew install geos
```

### Macports

For MacOS, GEOS can be installed using the [MacPorts](https://www.macports.org/) package repository, which downloads source packages and builds them in place using a recipe to ensure all packages integrate with each other nicely.

First [install MacPorts](https://www.macports.org/install.php). Then:

```
port install geos
```


## Build From Source

### Build Requirements

* [​CMake](https://cmake.org/download/) 3.8 or later.
* C++11 compiler. We regularly test GCC, Clang and Microsoft Visual C++.
* [Doxygen](https://www.doxygen.nl/) to build the API documentation.

### Build

Builds with CMake are done "outside the tree" either in a build directory in the source tree or next to the tree.

```bash
# Unpack and setup build directory
tar xvfz geos-3.9.0.tar.bz2
cd geos-3.9.0
mkdir _build
cd _build
# Set up the build
cmake \
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    ..
# Run the build, test, install steps
make
ctest
make install
```

### Build Options

The GEOS build can be customized using build options.

| Option               | Default    | Note  |
| :------------------: | :--------: | :---: |
| CMAKE_BUILD_TYPE     | Release    | Use `Debug` to build with debug flags and optimizations off. Use `Release` for packaging and working installs. |
| CMAKE_INSTALL_PREFIX | /usr/local | Set to install root. Librarys end up in `./libs` headers in `./include` |
| BUILD_DOCUMENTATION  | ON         | Attempt to find `doxygen` executable and build API docs |
| BUILD_SHARED_LIBS    | ON         | Build dynamically linkable libraries. |
| DISABLE_GEOS_INLINE  | OFF        | Turn off inlining. This is bad for performance, only do this if you cannot build to pass tests on your platform with inlining on. |

### Test Options

It is possible to run ctest directly. This gives access to ctest command line options (see ctest --help for a listing).

```bash
ctest
ctest --verbose
```

A list of GEOS test suites is obtained by running `ctest --show-only`:

```
$ ctest --show-only
#
# Test project /home/dan/dev/libgeos/cmake-build-debug
#  Test #1: test_geos_unit
#  Test #2: test_xmltester
#  Test #3: test_bug234
#  Test #4: test_sweep_line_speed
```

A subset of test suites can be run using a regular expression (and in this case, running 4 jobs in parallel):

```
$ ctest --tests-regex test_ --parallel 4
```
