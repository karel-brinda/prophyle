# see https://github.com/pypa/sampleproject
"""Installation script for ProPhyle.

To include binaries into the package, run make and set the system variable
PROPHYLE_PACKBIN to a non-zero value, e.g.,
PROPHYLE_PACKBIN=1 python3 setup.py install
"""

import glob
import os
import setuptools
import sys

if sys.version_info < (3, 4):
    sys.exit('Minimum supported Python version is 3.4')

bwa_dir = 'prophyle/prophyle_index/bwa'
if len(glob.glob(os.path.join(bwa_dir, "*.c"))) == 0 or len(glob.glob(os.path.join(bwa_dir, "*.h"))) == 0:
    sys.exit("BWA submodule is missing. Run 'make submodules' to download it.")

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Get the current version
exec(open("prophyle/version.py").read())

try:
    packbin = os.environ['PROPHYLE_PACKBIN']
except KeyError:
    packbin = False

if packbin:
    print("Adding executables and *.o files to the package", file=sys.stderr)

prophyle_files = [
    'Makefile',
    '*.py',
    'prophyle_assembler/*.cpp',
    'prophyle_assembler/*.h',
    'prophyle_assembler/Makefile',
    'prophyle_index/*.c',
    'prophyle_index/*.h',
    'prophyle_index/Makefile',
    'prophyle_index/bwa/*.c',
    'prophyle_index/bwa/*.h',
    'prophyle_index/bwa/Makefile',
    'prophyle_assignment/*.cpp',
    'prophyle_assignment/*.c',
    'prophyle_assignment/*.h',
    'prophyle_assignment/Makefile',
    'trees/*.nw',
] + (
    [
        'prophyle_index/prophyle_index',
        'prophyle_index/*.o',
        'prophyle_assembler/prophyle_assembler',
        'prophyle_assembler/*.o',
        'prophyle_assignment/prophyle_assignment',
        'prophyle_assignment/*.o',
        'prophyle_index/bwa/bwa',
        'prophyle_index/bwa/*.o',
    ] if packbin else []
)

setuptools.setup(
    name='prophyle',
    version=VERSION,
    description='ProPhyle metagenomic classifier',
    long_description=long_description,
    #
    url='https://github.com/prophyle/prophyle',
    download_url="https://github.com/prophyle/prophyle/releases",
    #
    author='Karel Brinda, Kamil Salikhov, Simone Pignotti, Gregory Kucherov',
    author_email=
    'kbrinda@hsph.harvard.edu, salikhov.kamil@gmail.com, pignottisimone@gmail.com, gregory.kucherov@univ-mlv.fr',
    license='MIT',
    #
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 3 :: Only',
        'Operating System :: Unix',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    #
    keywords='metagenomics classification NGS',
    #
    packages=["prophyle"],
    #
    package_data={'prophyle': prophyle_files},
    #
    entry_points={
        'console_scripts': [
            'prophyle = prophyle.prophyle:main',
            'prophyle_analyze.py = prophyle.prophyle_analyze:main',
            'prophyle_assignment.py = prophyle.prophyle_assignment:main',
            'prophyle_ncbi_tree.py = prophyle.prophyle_ncbi_tree:main',
            'prophyle_otu_table.py = prophyle.prophyle_otu_table:main',
            'prophyle_paired_end.py = prophyle.prophyle_paired_end:main',
            'prophyle_plot_tree.py = prophyle.prophyle_plot_tree:main',
            'prophyle_propagation_makefile.py = prophyle.prophyle_propagation_makefile:main',
            'prophyle_propagation_postprocessing.py = prophyle.prophyle_propagation_postprocessing:main',
            'prophyle_propagation_preprocessing.py = prophyle.prophyle_propagation_preprocessing:main',
            'prophyle_split_allseq.py = prophyle.prophyle_split_allseq:main',
        ],
    },
    #
    install_requires=[
        'ete3',
        'wheel',
        'bitarray',
        'psutil',
        'pysam',
        'scipy',
        'six',
    ],
    #
)
