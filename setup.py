from setuptools import setup, find_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="fits_metadata_extractor",  # Replace with your package name
    version="0.1.0",  # Update version as needed
    author="Bibal Sobeaux Pierre Gabriel",
    author_email="pierre.bibal-sobeaux@etu.unistra.fr",
    description="A package for processing, visualizing, and searching metadata from FITS files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pierregab/fits_metadata_extractor",  # Replace with your GitHub repository URL
    packages=find_packages(),  # Automatically find all packages in the directory
    install_requires=[
        "astropy>=5.0",
        "numpy>=1.21",
        "pandas>=1.3",
        "matplotlib>=3.4",
        "filelock>=3.0",
        "regions>=0.5",
        "shapely>=1.8",
        "mocpy>=0.10",
        "tqdm>=4.62",
        "reproject>=0.9",
        "adjustText>=0.8",
        "astroquery>=0.4"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires=">=3.7",  # Specify the minimum Python version
    entry_points={
        "console_scripts": [
            "fits-metadata-extractor=fits_metadata_extractor.cli:main",  # Optional CLI entry point
        ],
    },
    include_package_data=True,  # Include additional non-code files specified in MANIFEST.in
)
