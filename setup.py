from setuptools import setup, find_packages

setup(
    name="mlplib",
    version="0.1.0",
    description="Materials simulation library with machine learning potentials",
    author="nkitamuraQC",
    packages=find_packages(),
    python_requires='>=3.9',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU LESSER GENERAL PUBLIC LICENSE",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
)