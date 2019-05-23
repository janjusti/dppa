# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="dppa",
    version="0.1.0",
    description="Deep Protein Polarity Analyser",
    license="MIT",
    author="Jan Justi",
    author_email='jan.ma.justi@gmail.com',
    packages=find_packages(),
    install_requires=[
        "tqdm==4.31.1",
        "anytree==2.6.0",
        "pandas==0.24.2",
        "numpy==1.16.3",
        "biopython==1.73",
        "openpyxl==2.6.2",
        "styleframe==2.0.3"
    ],
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
    ],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'rundppa = dppa.core:main'
        ]
    }
)