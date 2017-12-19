echo start
date

oms_lncRNA_classify \
    -m Sus_scrofa.Sscrofa11.1.90.gtf \
    -l novel.lncrna.can.gtf

date
echo finished
