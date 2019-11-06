#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import sysconfig

with open('README.rst') as f:
    readme = f.read()

with open('HISTORY.rst') as f:
    history = f.read()

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pbpl-warp',
    version='0.1.0',
    description='Python package for working with WARP PIC code.',
    long_description=readme + '\n\n' + history,
    author='Brian Naranjo',
    author_email='brian.naranjo@gmail.com',
    url='https://github.com/bnara/pbpl-warp',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    license='MIT license',
    zip_safe=False,
    keywords='UCLA PBPL WARP particle cell beam relativistic',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    entry_points = {
        'console_scripts':
        ['pbpl-warp-phase-anim = pbpl.warp.phase_anim:main',
         'pbpl-warp-field-anim = pbpl.warp.field_anim:main']
    },
    console_scripts=[
        'bin/pbpl-warp-phase-anim',
        'bin/pbpl-warp-field-anim'
    ],
#    data_files=[('share', ['share/foo.ini'])],
    test_suite='tests',
    tests_require=test_requirements,
    namespace_packages=['pbpl']
)
