from setuptools import setup, find_packages

setup(
    name='MatPhylobi',
    version='0.2.0',
    author='Hania Kranas',
    author_email='matphylobi@gmail.com',
    url='https://github.com/hansiu/MatPhylobi',
    license='LICENSE.txt',
    packages=find_packages(exclude=['tests']),
    include_package_data=True,
    description='MatPhylobi is a command line tool for automatic construction of molecular data matrices for'
                'phylogenetic inference based on GenBank records.',
    long_description=open('README.md').read(),
    entry_points={'console_scripts': ['MatPhylobi=MatPhylobi.__main__:main']},
    install_requires=['biopython>=1.67'],
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.4.*, <4',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords=["bioinformatics", "phylogenetics", "genbank", "supermatrix"],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    zip_safe=True,
)
