from setuptools import setup, find_packages

setup(
    name='termol',
    version='0.1.4',
    packages=find_packages(),
    install_requires=[
        'numpy<2',
        'rdkit-pypi',
    ],
    entry_points={
        'console_scripts': [
            'termol=termol.cli:termol_cli',
            'termol-showcase=termol.cli:showcase_cli',
        ],
    },
    package_data={
        'termol': ['smiles_1000.csv'],
    },
    author='Nicholas Freitas',
    author_email='Nicholas.Freitas@ucsf.edu',
    description='A simple molecular renderer for the terminal using RDKit.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Nicholas-Freitas/TerMol',  # Replace with your GitHub URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
