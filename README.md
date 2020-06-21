# cgal-pybind

Simple test of binding [cgal](http://cgal.org) with [pybind11](https://pybind11.readthedocs.io)

## Methods:
* Skeletonization/contraction
* Discrete_authalic_parameterizer_3

## Trying:
    $ sh build.sh
    $ cd examples
    $ python3 -mvenv venv
    $ venv/bin/pip install -U pip setuptools wheel
    $ venv/bin/pip install -r requirements.txt
    $ venv/bin/python test.py
