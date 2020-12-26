from setuptools import setup

setup(
    name='Dummify',
    version='0.0.001',
    author='Aidan C',
    packages=['dummify'],
    scripts=['bin/tip4p_dummy.py'],
    url='http://pypi.python.org/pypi/PackageName/',
    description='An awesome package that does something',
    long_description=open('README.md').read(),
    install_requires=[
       "MDAnalysis >= 1.0.0",
       "numpy",
   ],
)
