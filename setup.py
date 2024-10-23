from setuptools import setup, find_packages

setup(
    name='fits_metadata_extractor',
    version='1.0.0',
    description='FITS Metadata Extraction and Homogenization Tool',
    author='Bibal Sobeaux Pierre Gabriel',
    author_email='your_email@example.com',  # Replace with your actual email
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'requests',
        'astropy',
        'astroquery',
        'mocpy',
        'filelock',
        # Add other dependencies as needed
    ],
    entry_points={
        'console_scripts': [
            'fits_metadata_extractor=fits_metadata_extractor.extractor:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  # Update if you use a different license
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
