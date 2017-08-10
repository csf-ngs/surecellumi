from setuptools import setup

setup(name='surecellumi',
      version='0.1',
      description='SureCell WTA UMI extraction',
      url='http://github.com/csf/surecellumi',
      author='Ido M. Tamir',
      author_email='ido.tamir@vbcf.ac.at',
      license='GPLv3',
      packages=['surecellumi'],
      install_requires=[
          'cutadapt>=1.14',
          'Distance>=0.1.3'
      ],
      setup_requires=['pytest-runner'],
      tests_requires=['pytest'],
      include_package_data=True,
      scripts=['bin/mergeUmi'],
      zip_safe=False)
