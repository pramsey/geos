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

```
# Add the EPEL repository
yum -y install epel-release

# Install the GEOS runtime and development packages
rpm -Uvh geos geos-devel
```

### Ubuntu

The [Ubuntu GIS](https://wiki.ubuntu.com/UbuntuGIS) project maintains a collection of repositories with builds of common open source geospatial projects, including GEOS.

```
# Add the Ubuntu GIS PPA
sudo apt-get install python-software-properties
sudo add-apt-repository ppa:ubuntugis/ppa

# Install the packages

```

### Debian

The [Debian GIS](https://wiki.debian.org/DebianGis) project maintains [GEOS packages](https://tracker.debian.org/pkg/geos) and pushes them into the appropriate Debian respositories.

```
sudo apt-get install geos
```

### Amazon Linux

Amazon Linux is based on RH7, and can read from the EPEL repository. To enable using Amazon tools, use the `amazon-linux-extras` utility.

```
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


## Build Source

### Build Requirements

* [​CMake](https://cmake.org/download/) 3.8 or later.
* C++11 compiler. We regularly test GCC, Clang and Microsoft Visual C++.
* [Doxygen](https://www.doxygen.nl/) to build the API documentation.

### CMake Configuration


