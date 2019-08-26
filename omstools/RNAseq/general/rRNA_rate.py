import sys
import os
import re

if not len(sys.argv) == 4:
    print 'python ' + sys.argv[0] + ' sample_file mapping_dir summary'
    sys.exit(0)

sample_file = sys.argv[1]
mapping_dir = sys.argv[2]
summary = sys.argv[3]

sample_list = [each.strip().split()[0] for each in open(sample_file)]

summary_info = open(summary, 'w')
summary_info.write('Sample_ID\trRNA_rate\n')
for n, each_sample in enumerate(sample_list):
    each_sample_mapping_file = os.path.join(mapping_dir,
                                            '%s.log' % each_sample)
    with open(each_sample_mapping_file) as each_sample_mapping_file_info:
        for eachline in each_sample_mapping_file_info:
            if 'overall alignment rate' in eachline:
                mapping_rate = re.search(r'(.*)%', eachline).groups()[0]
                mapping_rate = float(mapping_rate) / 100
                summary_info.write('%s\t%s\n' % (each_sample, mapping_rate))
summary_info.close()
