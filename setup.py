from setuptools import setup

setup(name='codon_tools',
      version='0.2',
      description='Toolkit to manipulate synonymous codon usage in various ways',
      url='http://github.com/clauswilke/codon_tools',
      author='Claus O. Wilke',
      author_email='wilke@austin.utexas.edu',
      license='MIT',
      test_suite='test',
      package_dir = {'codon_tools':'codon_tools'},
      install_requires=['Biopython'],
      packages=['codon_tools']
      )
