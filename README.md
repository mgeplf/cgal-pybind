# cgal-pybind

 [CGAL](http://cgal.org) Python binding with [pybind11](https://pybind11.readthedocs.io)

## Methods:
* Skeletonization/contraction
* Segmentation

## Installation
```bash
$> git clone https://github.com/CGAL/cgal.git $PATH_TO_GIT
$> git submodule init
$> git submodule update
$> export CGAL_DIR=$PATH_TO_GIT
$> pip install .
```

## Test
```bash
$> python example/test
```

## Requierements
* cmake > 3.0.9
* C++ compiler (with C++11 support)
* trimesh Python package
* CGAL header
