# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='replace_allele_counts',
    version='0.2.0',
    description='Package for replaceing allele counts in maf using a fill out produce with GetBaseCountMultiSample',
    long_description=readme,
    author='Ronak Shah',
    author_email='rons.shah@gmail.com',
    url='https://github.com/rhshah/remove_variants',
    license=license,
    install_requires=['nose==1.3.7', 'pandas', 'coloredlogs==5.2','codecov==2.0.5', 'coverage==4.3.4'],
    packages=find_packages(exclude=('tests', 'docs'))
)
