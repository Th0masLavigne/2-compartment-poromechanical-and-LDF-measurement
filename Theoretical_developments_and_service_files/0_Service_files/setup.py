from setuptools import setup, find_packages

setup(
    name='Porous_media_package',
    version='0.9.0',
    packages=find_packages(),
    description='Porous media mechanics within FEniCSx developed during Thomas Lavigne P.h.D',
    author='Thomas Lavigne',
    author_email='lavignethomas@hotmail.fr',
    install_requires=[],
    url='https://github.com/Th0masLavigne/2-compartment-poromechanical-and-LDF-measurement.git',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)