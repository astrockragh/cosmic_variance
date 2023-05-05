#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy>=1.17.0', 'scipy>=1.2.0', 'pandas>=1.0.0']

test_requirements = [ ]

setup(
    author="Christian Kragh Jespersen",
    author_email='ckragh@princeton.edu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Package to calculate cosmic variance in rectangular pencil-beam surveys",
    install_requires=requirements,
    license="MIT license",
    long_description = readme,
    long_description_content_type='text/x-rst',
    include_package_data = True,
    keywords = ['Cosmology', 'Galaxies', 'Statistics', 'Astrostatistics', 'Cosmic Variance'],
    name = 'cosmic_variance',
    packages=find_packages(include=['cosmic_variance', 'cosmic_variance.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/astrockragh/cosmic_variance',
    version='0.1.0',
    zip_safe=False,
)
