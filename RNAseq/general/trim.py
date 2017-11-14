import luigi
import envoy
import pandas as pd
import os

ADAPTER = '/public/software/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa'


class trim(luigi.Task):
    sample = luigi.Parameter()
    raw_dir = luigi.Parameter()
    out_dir = luigi.Parameter()

    def run(self):
        cmd = 'java -jar /public/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar \
               PE -phred33 -threads 6 \
               {t.raw_dir}/{t.sample}_1.clean.fq.gz \
               {t.raw_dir}/{t.sample}_2.clean.fq.gz \
               {t.out_dir}/{t.sample}_1.clean.fq.gz \
               {t.out_dir}/{t.sample}_1.clean.unpair.fq.gz \
               {t.out_dir}/{t.sample}_2.clean.fq.gz \
               {t.out_dir}/{t.sample}_2.clean.unpair.fq.gz \
               ILLUMINACLIP:{a}:2:30:10 \
               LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100'.format(
            t=self, a=ADAPTER)
        cmd_inf = envoy.run(cmd)
        with self.output().open('w') as run_inf:
            run_inf.write(cmd_inf.std_err)

    def output(self):
        return luigi.LocalTarget(
            '{t.raw_dir}/logs/{t.sample}.log'.format(t=self))


class trim_summary(luigi.Task):
    sample_inf = luigi.Parameter()
    raw_dir = luigi.Parameter()
    out_dir = luigi.Parameter()

    def requires(self):
        sample_df = pd.read_table(self.sample_inf, header=None, index_col=1)
        log_dir = os.path.join(self.raw_dir, 'logs')
        try:
            os.makedirs(log_dir)
        except OSError as e:
            print e
        return [
            trim(
                sample=each_sample, raw_dir=self.raw_dir, out_dir=self.out_dir)
            for each_sample in sample_df.index
        ]

    def run(self):
        pass


if __name__ == '__main__':
    luigi.run()
