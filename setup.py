from setuptools import find_packages
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("forschungspraktikum.jiles_atherton.functions",
              ["forschungspraktikum/jiles_atherton/functions.pyx"],
              include_dirs=[
                  numpy.get_include(),
                  "forschungspraktikum/jiles_atherton/"
              ])
]

setup(
    name='forschungspraktikum',
    version='0.1',
    description='Quellcode für die Inhalte des Forschungspraktikums am LEMF, Uni Erlangen.',
    author='Adrian Schneider',
    author_email='post@adrian-schneider.de',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'algopy'],
    ext_modules=cythonize(extensions),
)