from setuptools import setup, find_packages
from os import path
import glob

from io import open

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(

    name='chlamdb',
    version='3.0',
    description='ChlamDB',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/metagenlab/chlamdb',
    author='Trestan Pillonel',
    author_email='trestan.pillonel@gmail.com',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='genome orthology comparative genomics bacteria KEGG COG EC Pfam Interpro',
    packages=find_packages(exclude=['chlamdb', 'db_setup', 'tests']),
    python_requires='>=3.0',
    scripts=[script for script in glob.glob('db_setup/*')],
    project_urls={ 
        'Bug Reports': 'https://github.com/metagenlab/chlamdb/issues',
        'Source': 'https://github.com/metagenlab/chlamdb',
    },
)
