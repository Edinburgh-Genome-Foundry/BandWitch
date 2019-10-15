import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('bandwitch/version.py').read()) # loads __version__

setup(
    name='bandwitch',
    version=__version__,
    author='Zulko',
    url='https://github.com/Edinburgh-Genome-Foundry/BandWitch',
    description='Enzyme selection for DNA verification and identification',
    long_description=open('pypi-readme.rst').read(),
    license='MIT',
    keywords="Restriction enzyme synthetic biology DNA band patterns",
    packages= find_packages(exclude='docs'),
    include_package_data=True,
    install_requires=('tqdm', 'biopython', 'scipy', 'matplotlib',
                      'bandwagon', 'dna_features_viewer', 'flametree',
                      'pandas')
)
