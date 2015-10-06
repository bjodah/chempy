#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
from distutils.core import setup

pkg_name = "chempy"

CHEMPY_RELEASE_VERSION = os.environ.get('CHEMPY_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
if os.environ.get('CONDA_BUILD', '0') == '1':
    try:
        CHEMPY_RELEASE_VERSION = 'v' + open(
            '__conda_version__.txt', 'rt').readline().rstrip()
    except IOError:
        pass

release_py_path = os.path.join(pkg_name, '_release.py')

if (len(CHEMPY_RELEASE_VERSION) > 1 and
   CHEMPY_RELEASE_VERSION[0] == 'v'):
    TAGGED_RELEASE = True
    __version__ = CHEMPY_RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(open(release_py_path).read())

with open(pkg_name + '/__init__.py') as f:
    long_description = f.read().split('"""')[1]

submodules = [
    'chempy.kinetics',
    'chempy.util',
]

tests = [
    'chempy.tests',
    'chempy.kinetics.tests',
]

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.4',
]

setup_kwargs = {
    'name': pkg_name,
    'version': eval('__version__'),  # silence pyflakes
    'description': (
        'Package useful for (physical) chemistry'
    ),
    'long_description': long_description,
    'author': 'Björn Dahlgren',
    'author_email': 'bjodah@DELETEMEgmail.com',
    'license': 'BSD',
    'keywords': ("chemistry", "water properties"),
    'url': 'https://github.com/bjodah/' + pkg_name,
    'packages': [pkg_name] + submodules + tests,
    'classifiers': classifiers,
}

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set CHEMPY_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
