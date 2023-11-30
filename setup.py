import os
from setuptools import find_packages, setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='wfomc_tuttepoly',
    version='0.1',
    url='https://github.com/l2l7l9p/wfomc_tuttepoly',
    author='Qipeng Kuang',
    author_email='kuangqipeng@connect.hku.hk',
    license='MIT',

    packages=find_packages(include=['sampling_fo2', 'sampling_fo2.*']),
    install_requires=["symengine",
                      "sympy",
                      "lark",
                      "numpy",
                      "networkx",
                      "contexttimer",
                      "logzero",
                      "pandas",
                      "pysat",
                      "tqdm",
                      "dataclasses",
                      "PrettyPrintTree"],
    python_requires='>=3.8',
    description = ("Combining WFOMC and various polynomials such as Tutte Polynomial"),
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: MIT License",
    ],
    keywords="WFOMC Tutte-polynomial",
)
