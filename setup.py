from setuptools import find_packages, setup

setup(
    name="laespec",
    version="0.1.0",
    description="Read and manipulate simulated LAE Lyman-alpha damping-wing spectra",
    packages=find_packages(exclude=["tests"]),
    package_data={"laespec": ["manifest.json"]},
    python_requires=">=3.7",
    install_requires=["numpy>=1.17"],
    extras_require={"plot": ["matplotlib"], "test": ["pytest"]},
    entry_points={"console_scripts": ["laespec=laespec.cli:main"]},
)
