from setuptools import setup, find_packages
import sys
import os

install_requires = ['pyfaidx','urllib3[secure]','six']
#if python is 2.6, requires argparse
if sys.version_info[0] == 2 and sys.version_info[1] < 7:
    install_requires.extend(['argparse'])


def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(
    name='TESQuIRE',
    provides='SQuIRE',
    version=get_version(open('squire/__init__.py').read()),
    author='Wan Rou Yang',
    author_email='wyang17@jhmi.edu',
    url='https://github.com/wyang17/SQuIRE',
    description='SQuIRE: Software for Quantifying Interspersed Repeat Expression',
    long_description=open('README.md').read(),
    license='GPL-3',
    packages=find_packages(),    
    include_package_data=True,
    package_data = {'squire': ['software/*']},
    install_requires=install_requires,
    entry_points={'console_scripts': ['squire = squire.cli:main']},
    classifiers=[
            "Development Status :: 4 - Beta",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Operating System :: Unix",
#            "Programming Language :: Python :: 3.4",
#            "Programming Language :: Python :: 3.3",
#            "Programming Language :: Python :: 3.2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 2.6",
#            "Programming Language :: Python :: Implementation :: PyPy",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
