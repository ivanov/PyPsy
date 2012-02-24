"""pypsy version/release information"""

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
#_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 2 - Pre-Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               # XXX: I'm open to changing this, Dan
               "License :: OSI Approved :: BSD License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description = "pypsy: psychometric function fitting in Python"

# Note: this long_description is actually a copy/paste from the top-level
# README.txt, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = """
===================================================
 pypsy: psychometric function fitting in Python
===================================================

Website and mailing list
========================

We currently have neither.

Code
====

You can find our sources and single-click downloads:

* `Main repository`_ on Github.
* Download as a tar/zip file the `current trunk`_.
* Downloads of all `available releases`_.

.. _main repository: http://github.com/ivanov/pypsy
.. _Documentation: http://pypsy.pirsquared.org
.. _current trunk: http://github.com/ivanov/pypsy/archives/master
.. _available releases: http://github.com/ivanov/pypsy/downloads

Related Packages
================

* psignifit 3.0
    a rewrite of what of a popular Matlab package, now primarily Python
    focused

* palamedes
    a Matlab packages which accompanies Kingdom and Prins

License information
===================

pypsy is licensed under the terms of the new BSD license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2012, Dan Coates and Paul Ivanov
All rights reserved.
"""

NAME = "pypsy"
MAINTAINER = "Dan Coates and Paul Ivanov"
MAINTAINER_EMAIL = "{ dan.coates, pi } @  berkeley.edu"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/ivanov/PyPsy"
DOWNLOAD_URL = "http://github.com/ivanov/PyPsy/downloads"
LICENSE = "Simplified BSD"
AUTHOR = "Paul Ivanov"
AUTHOR_EMAIL = "{ dan.coates, pi } @  berkeley.edu"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['pypsy',
            'pypsy.tests',
            ]
PACKAGE_DATA = {"pypsy": ["LICENSE","tests/*.txt", "tests/*.npy"]}
REQUIRES = ["numpy", "matplotlib", "scipy"]
