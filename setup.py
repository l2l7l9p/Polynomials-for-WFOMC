import os
from setuptools import find_packages, setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='Polynomials-for-WFOMC',
    version='0.1',
    url='https://github.com/l2l7l9p/Polynomials-for-WFOMC',
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
    description = ("Computing polynomials for WFOMC for a C2 Sentence with cardinality constraints"),
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: MIT License",
    ],
    keywords="Weighted First-Order Model Counting, Axiom, Enumerative Combinatorics, Tutte Polynomial",
)
