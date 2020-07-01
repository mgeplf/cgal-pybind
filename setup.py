import os
import sys
import subprocess

from setuptools import setup, Extension
from setuptools import find_packages
from setuptools.command.build_ext import build_ext

import imp
VERSION = imp.load_source("cgal_pybind.version", "src/cgal_pybind/version.py").VERSION


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            _ = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for auto-detection of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
        build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name='cgal_pybind',
    version='0.0.1',
    author='BBP',
    author_email='jean.jacquemier@epfl.ch',
    description='A Python binding for some CGAL classed',
    long_description='',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    ext_modules=[CMakeExtension('cgal_pybind._cgal_pybind')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[],
)

'''
packages = find_packages('src'),
# tell setuptools that all packages will be under the 'src' directory
# and nowhere else
package_dir = {'': 'src'},
# add an extension module named 'python_cpp_example' to the package
# 'python_cpp_example'
ext_modules = [CMakeExtension('python_cpp_example/python_cpp_example')],
'''