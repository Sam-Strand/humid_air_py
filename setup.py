from setuptools import setup, find_packages

setup(
    name="humid_air_py",
    version="1.0.2",
    url="https://github.com/Sam-Strand/humid_air_py",
    author="Садовский М.К.",
    author_email="i@maxim-sadovskiy.ru",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "numba>=0.55.0"
    ],
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
    ],
)
