from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='EM2Tools',
    version='0.1.0',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/fplewniak/EM2Tools',
    license='CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21',
    author='Frédéric PLEWNIAK',
    author_email='f.plewniak@unistra.fr',
    description='Python tools for environmental genomics data analysis and manipulation.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.6', install_requires=['biopython', 'gffpandas', 'pandas']
)
