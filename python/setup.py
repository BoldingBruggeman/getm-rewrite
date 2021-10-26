from setuptools import setup, find_packages

setup(
    name='pygetm',
    version='0.1.1',
    author='Bolding-Bruggeman ApS',
    author_email='jorn@bolding-bruggeman.com',
    license='GPL',
    packages=find_packages(include=['pygetm*']),
    package_data={'pygetm': ['*.so', '*.dll', '*.dylib', '*.pyd']},
    zip_safe=False
)


