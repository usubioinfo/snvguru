from setuptools import setup, find_packages

setup(
    name='VirnaSNV',
    version='0.0.1',
    license='MIT',
    author="David Guevara-Barrientos",
    author_email='david.guevara@usu.edu',
    packages=find_packages(),
    package_dir={"": "src"},
    url='http://bioinfo.usu.edu/virnasnv',
    install_requires=[
          'pandas==1.4.1', 'joblib', 'matplotlib', 'seaborn', 'biopython', 'xlrd==1.2.0'
      ],
    python_requires='>=3.8',
    include_package_data=True,
    package_data={"config": ["*.config"]}
)