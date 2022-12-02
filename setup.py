from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages


with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='glassviewer',
    version='1.0',
    author='Sarath Menon',
    author_email='bggmail.com',
    description='Python library written in C++ for calculation of local atomic structural environment',
    long_description=readme,
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    # add an extension module named 'python_cpp_example' to the package
    headers=["src/glassviewer/atom.h", "src/glassviewer/system.h", "lib/voro++/voro++.hh","lib/wignerSymbols/include/wignerSymbols.h",'lib/fftw3/fftw3.h'],
    ext_modules=[
        Pybind11Extension(
            "glassviewer.catom",
            ["src/glassviewer/atom.cpp", "src/glassviewer/atom_binding.cpp"],
            language='c++',
            include_dirs=['lib/voro++']
        ),
        Pybind11Extension(
            "glassviewer.csystem",
            ["src/glassviewer/system.cpp", "src/glassviewer/system_binding.cpp", "src/glassviewer/atom.cpp", "lib/voro++/voro++.cc","lib/wignerSymbols/src/wignerSymbols-cpp.cpp"],
            language='c++',
            include_dirs=['lib/voro++','lib/wignerSymbols/include','lib/fftw3'],
            library_dirs=['lib/fftw3'],
            libraries=['libfftw3-3'],
        ),
    ],
    # add custom build_ext command
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    download_url = 'https://github.com/bgbaizh/GlassViewer',
    url = 'https://github.com/bgbaizh/GlassViewer',
    install_requires=['numpy', 'ase', 'plotly', 'ipywidgets'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    package_data={'glassviewer':['libfftw3-3.dll']}
)
