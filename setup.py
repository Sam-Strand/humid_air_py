from setuptools import setup, find_packages

setup(
    name="humid_air_py",
    version="1.0.2",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "numba>=0.55.0"
    ]
)
