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
        'pandas==0.22.0',
        'numpy=1.13.3',
        'HTSeq',
        'click',
        'distribute',
        'tabulate',
        'biopython'
    ],
    scripts=['scripts/omstools'],
    entry_points={
        'console_scripts': [
            "oms_lncRNA_classify=omstools.RNAseq.lncRNA.classify.lncRNA_classify:main",
            "oms_gtf_split=omstools.general.gtf_split:main",
            "oms_transcript_feature=omstools.RNAseq.lncRNA.assembly.transcript_feature:main",
            "oms_gtf2bed=omstools.general.gtf2bed:main",
            "oms_kill=omstools.general.killjob:main",
            'oms_fa_split=omstools.general.fa_split_v2:main',
            'oms_merge=omstools.general.merge_files:main',
            'oms_rename=omstools.general.rename_batch_files:main'
        ],
        'gtf': [
            "func_tr_from_gtf_lines=omstools.utils.transcript:transcripts_from_gtf_lines",
            "func_strand_int_to_str=omstools.utils.transcript:strand_int_to_str",
            "func_parse_gtf=omstools.utils.transcript:parse_gtf",
            "func_to_formatted_gtf=omstools.utils.transcript:to_formatted_gtf",
            "func_sort_gtf=omstools.utils.gtf:sort_gtf",
            "func_get_gene_type=omstools.utils.gtf:get_gene_type",
            "func_get_tr_type=omstools.utils.gtf:get_tr_type",
            "func_merge_sort_gtf_files=omstools.utils.gtf:merge_sort_gtf_files",
            "class_GTFFeature=omstools.utils.gtf:GTFFeature",
            "error_GTFError=omstools.utils.gtf:GTFError",
            "dict_GENCODE_CATEGORY_MAP=omstools.utils.gtf:GENCODE_CATEGORY_MAP"
        ]
    }

)

print'''
--------------------------------
omsTools installation complete!
--------------------------------
'''
