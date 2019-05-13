import setuptools
#from setuptools import find_packages, setup

with open('README.md') as fh:
    long_description = fh.read()

setuptools.setup(
    name='ncortho',
    version='1.0',
    author='Andreas Blaumeiser',
    author_email='',
    description='ncRNA orthology search',
    long_description= long_description,
    url='',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': ['ncortho = ncortho.ncortho:main'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: ',
        'Operating System :: Linux',
    ],
)



#setup, check if the required third-party tools are available

tools = ['blastn', 'infernal', 't_coffee', 'mafft', 'alifold']
