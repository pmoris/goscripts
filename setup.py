from setuptools import setup

setup(
    name='goscripts',
    version='1.0',
    description='Tools for gene ontology enrichment analysis',
    url='https://github.com/pmoris/goscripts',
    author='Pieter Moris',
    author_email='pieter.moris@uantwerpen.be',
    license='MIT',
    packages=['goscripts'],
    install_requires=['pandas', 'numpy', 'statsmodels'],
    # entry_points={
    #     "console_scripts": [
    #         "go_enrichment_script = goscripts.go_enrichment_script:main",
    #     ]
    # },
    zip_safe=False)
