""""The package to manage instances for my research problem."""

from setuptools import setup, find_packages

setup(
    name="instance",
    version="0.1.dev",
    author="Thiago Jobson",
    description=__doc__,
        packages = find_packages(
        where = 'src',
        include = ['instances',],
    #    exclude = ['additional',]
    ),
    package_dir = {"":"src"},
    include_package_data=False,
    install_requires=[
        "numpy>=1.20.1",
        "pytest>=6.2.2",
    ],
)