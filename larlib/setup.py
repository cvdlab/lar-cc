from setuptools import setup

import sys,os,math

setup(name='larlib',
	  version='0.1',
	  description='Linear algebraic representation of topology and geometry as a complex of chain spaces',
	  url='https://github.com/cvdlab/lar-cc',
	  author='Alberto Paoluzzi',
	  author_email='paoluzzi@dia.uniroma3.it',
	  license='MIT',
	  packages=['larlib'],
	  zip_safe=False,
	  install_requires=[
						'scipy',
						#'pyplasm',
						'triangle',
						#'matplotlib',
						'Cython',
      ],
	  include_package_data=True)