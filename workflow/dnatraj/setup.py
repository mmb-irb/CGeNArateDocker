"""  Setup script. Used by easy_install and pip. """

import os
import sys
import re
from setuptools import dist, setup, find_packages

def read(*rnames):
    return open(os.path.join(os.path.dirname(__file__), *rnames)).read()

"""Discover the package version"""
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
VERSIONFILE = "dnatraj/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RunTimeError("Unable to find version string in {}.".format(VERSIONFILE))


"""Check Python version"""
if  sys.version_info[0:2] < (3, 4):
    raise RuntimeError('DNAtraj requires Python 3.4+')

setup_args = {
    'name':             "dnatraj",
    'version':          verstr,
    'description':      "Data I/O package for DNA trajectories.",
    'long_description_content_type':      "text/markdown",
    'long_description':      read('README.md'),
    'author':           "Charles Laughton",
    'author_email':     "charles.laughton@nottingham.ac.uk",
    'url':              "https:/?",
    'download_url':     "https://{}.tar.gz".format(verstr),
    'license':          "MIT license.",

    'classifiers': [
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix'
    ],

    'packages': find_packages(),

    'install_requires': [
                         'numpy',
                         'mdtraj',
                        ],

    'zip_safe': False,
}

setup(**setup_args)
