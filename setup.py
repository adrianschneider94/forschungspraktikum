from setuptools import setup

setup(
   name='forschungspraktikum',
   version='0.1',
   description='Quellcode für die Inhalte des Forschungspraktikums am LEMF, Uni Erlangen.',
   author='Adrian Schneider',
   author_email='post@adrian-schneider.de',
   packages=['forschungspraktikum'],
   install_requires=['numpy', 'scipy', 'algopy'],
)