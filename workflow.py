"""
Usage: workflow.py <working_dir> <prefix>
"""


from docopt import docopt
import logging
import os
import luigi
import sciluigi as sl

log = logging.getLogger('sciluigi-interface')


class MyWorkflow(sl.WorkflowTask):
    working_dir = sl.Parameter()
    prefix = sl.Parameter()

    def workflow(self):
        # Initialize tasks:
        input_task = self.new_task("input", InputTask,
                                   working_dir=self.working_dir,
                                   prefix=self.prefix)
        fastqc_task_before = self.new_task('fastQC before', FastqcTask,
                                           working_dir=self.working_dir,
                                           subdir="fastqc_before"
                                           )
        sickle_task = self.new_task('sickle', SickleTask,
                                    working_dir=self.working_dir)
        fastqc_task_after = self.new_task('fastQC after', FastqcTask,
                                          working_dir=self.working_dir,
                                          subdir="fastqc_after")

        # Connecting outputs to inputs:
        sickle_task.in_fastq1 = input_task.out_fastq1
        sickle_task.in_fastq2 = input_task.out_fastq2
        fastqc_task_before.in_fastq1 = input_task.out_fastq1
        fastqc_task_before.in_fastq2 = input_task.out_fastq2
        fastqc_task_after.in_fastq1 = sickle_task.out_fastq1
        fastqc_task_after.in_fastq2 = sickle_task.out_fastq2

        # Return the last task(s) in the workflow chain.
        return [fastqc_task_after, fastqc_task_before]


class InputTask(sl.ExternalTask):
    working_dir = sl.Parameter()
    prefix = sl.Parameter()

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.working_dir), str(self.prefix) + "_1.fq.gz"))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.working_dir), str(self.prefix) + "_2.fq.gz"))


class FastqcTask(sl.Task):
    working_dir = sl.Parameter()
    subdir = sl.Parameter()
    in_fastq1 = None
    in_fastq2 = None

    def out_dir(self):
        return sl.TargetInfo(self, os.path.join(str(self.working_dir), str(self.subdir)))

    def run(self):
        self.out_dir().target.fs.mkdir(self.out_dir().path)
        self.ex('fastqc --noextract -o {outdir} {in1} {in2}'.format(
            in1=self.in_fastq1().path,
            in2=self.in_fastq2().path,
            outdir=self.out_dir().path)
        )
        with self.out_flag().open("w") as fp:
            fp.write("Done")





class SickleTask(sl.Task):
    working_dir = sl.Parameter()
    in_fastq1 = None
    in_fastq2 = None

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.working_dir),
                                                "sickle_trimmed_{}".format(os.path.basename(self.in_fastq1().path))))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.working_dir),
                                                "sickle_trimmed_{}".format(os.path.basename(self.in_fastq2().path))))

    def run(self):
        self.ex("sickle pe -g -f {reads1} -r {reads2} -o {out1} -p {out2} -t sanger -s /dev/null".format(
            reads1=self.in_fastq1().path,
            reads2=self.in_fastq2().path,
            out1=self.out_fastq1().path,
            out2=self.out_fastq2().path))


if __name__ == '__main__':
    arguments = docopt(__doc__)
    log.debug(str(arguments))
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow,
              cmdline_args=[f'--working-dir={arguments["<working_dir>"]}',
                            f'--prefix={arguments["<prefix>"]}']
                 )
