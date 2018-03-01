# coding: utf-8
from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='aimd',
    version='0.1',
    description='Python project for diffusion analysis from AIMD simulations',
    long_description=long_description,
    url='https://github.com/mogroupumd/aimd',
    author='Yifei Mo Research Group at UMD',
    author_email='hxfthu@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7'
    ],
    keywords=["diffusion", "AIMD"],
    packages=find_packages(),
    install_requires=['pymatgen>=2017.12.30','argparse','PrettyTable'],
    extras_require={},
    package_data={},
    data_files=[],
    entry_points={
        'console_scripts': [
            'analyze_aimd=aimd.script.analyze_aimd:main',
        ],
    },
    project_urls={
        'Bug Reports': 'https://github.com/mogroupumd/aimd/issues',
        'Source': 'https://github.com/mogroupumd/aimd',
    },
)