from setuptools import setup

setup(name='codon_manipulator',
      version='0.1',
      description='Toolkit to manipulate synonymous codon usage in various ways',
      url='http://github.com/clauswilke/codon_manipulator',
      author='Claus O. Wilke',
      author_email='wilke@austin.utexas.edu',
      license='MIT',
      package_dir = {'codon_manipulator':'src'},
      install_requires=['Biopython'],
      packages=['codon_manipulator'],
      zip_safe=False
      )
