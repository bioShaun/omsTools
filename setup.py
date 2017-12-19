from setuptools import setup, find_packages


version = '0.1dev'

print '''
--------------------------------
Installing omsTools version {v}
--------------------------------
'''.format(v=version)


setup(
    name='omstools',
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
        'click',
        'distribute'
    ],
    entry_points={
        'console_scripts': [
            "oms_lncRNA_classify=omstools.RNAseq.lncRNA.classify.lncRNA_classify:main",
            "oms_gtf_split=omstools.general.gtf_split:main",
        ],
        'gtf': [
            "func_tr_from_gtf_lines=omstools.utils.transcript:transcripts_from_gtf_lines",
            "func_strand_int_to_str=omstools.utils.transcript:strand_int_to_str",
            "dict_GENCODE_CATEGORY_MAP=omstools.utils.transcript:GENCODE_CATEGORY_MAP",
            "func_parse_gtf=omstools.utils.transcript:parse_gtf",
            "class_GTFFeature=omstools.utils.gtf:GTFFeature",
        ]
    }

)

print'''
--------------------------------
omsTools installation complete!
--------------------------------
'''
