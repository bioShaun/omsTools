from setuptools import setup, find_packages


version = '0.1dev'

print '''
--------------------------------
Installing omsTools version {v}
--------------------------------
'''.format(v=version)


setup(
    name='omsTools',
    version=version,
    author='lx Gui',
    author_email='guilixuan@gmail.cn',
    keywords=['bioinformatics', 'NGS', 'RNAseq'],
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'envoy',
        'pandas',
        'HTSeq',
        'click'
    ],
    entry_points={
        'console_scripts': [
            "oms_lncRNA_classify=RNAseq.lncRNA.classify.lncRNA_classify:main"
        ]
    }

)

print'''
--------------------------------
omsTools installation complete!
--------------------------------
'''
