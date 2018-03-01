from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
exec(open('nanosv/version.py').read())

def readme():
    with open('README.rst') as f:
        return f.read()


setup( name='NanoSV',
       version=__version__,
       description='Structural variation detection tool for Oxford Nanopore data.',
       long_description=readme(),
       classifiers=[
           'Development Status :: 4 - Beta',
           'Intended Audience :: Science/Research',
           'License :: OSI Approved :: MIT License',
           'Programming Language :: Python :: 3',
           'Topic :: Scientific/Engineering :: Bio-Informatics'
       ],
       keywords='nanosv sv-caller nanopore',
       url='https://github.com/mroosmalen/nanosv',
       author='Mark van Roosmalen',
       author_email='m.vanroosmalen-2@umcutrecht.nl',
       license='MIT',
       python_requires='>=3',
       install_requires=['pysam','pyvcf','configparser'],
       packages=['nanosv','nanosv.utils','nanosv.classes'],
       package_data={'nanosv': ['config.ini','bedfiles/*.bed']},
       including_package_data=True,
       entry_points={
           'console_scripts': [
               'NanoSV=nanosv.NanoSV:main'
           ]
       }
)