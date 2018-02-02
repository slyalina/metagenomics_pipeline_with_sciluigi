"""
A basic metagenomics workflow using SciLuigi.

Usage: workflow.py <working_dir> <prefix>
"""


import hashlib
import logging
import os
import sys

import luigi
import sciluigi as sl
from docopt import docopt
from six import reraise as raise_
from datetime import datetime

log = logging.getLogger('sciluigi-interface')
traceback = sys.exc_info()[2]
BUF_SIZE = 128*1024


def get_md5(filepath):
    """
    Buffered calculation of standard md5sum on contents of provided file
    :param filepath: file to be opened and checksummed
    :return: md5 string
    """
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            md5.update(data)
    return md5.hexdigest()


class SingleSampleWorkflow(sl.WorkflowTask):
    working_dir = sl.Parameter()
    prefix = sl.Parameter()
    tool_log_subdir = sl.Parameter(default="tool_logs")
    contaminant_removal_method = sl.Parameter(default="bbsplit")

    def workflow(self):
        # Initialize tasks:
        input_task = self.new_task("input", Input,
                                   prefix=self.prefix)
        fastqc_task_before = self.new_task('fastQC before', Fastqc,
                                           subdir="fastqc_before"
                                           )
        sickle_task = self.new_task('sickle', TrimReads)
        fastqc_task_after = self.new_task('fastQC after', Fastqc,
                                          subdir="fastqc_after")

        contam_index = self.new_task("index", CreateIndexForContamRemoval)
        # Connecting outputs to inputs:
        sickle_task.in_fastq1 = input_task.out_fastq1
        sickle_task.in_fastq2 = input_task.out_fastq2
        fastqc_task_before.in_fastq1 = input_task.out_fastq1
        fastqc_task_before.in_fastq2 = input_task.out_fastq2
        fastqc_task_after.in_fastq1 = sickle_task.out_fastq1
        fastqc_task_after.in_fastq2 = sickle_task.out_fastq2

        # Return the last task(s) in the workflow chain.
        return [fastqc_task_after, fastqc_task_before, contam_index]


class Input(sl.ExternalTask):

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir), str(self.workflow_task.prefix) + "_1.fq.gz"))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir), str(self.workflow_task.prefix) + "_2.fq.gz"))


class Fastqc(sl.Task):
    subdir = sl.Parameter()
    in_fastq1 = None
    in_fastq2 = None

    def outfile(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir), str(self.subdir), "sentinel.txt"))

    def run(self):
        outdir=os.path.dirname(self.outfile().path)
        self.outfile().target.fs.mkdir(outdir)
        self.ex('fastqc --noextract -o {outdir} {in1} {in2}'.format(
            in1=self.in_fastq1().path,
            in2=self.in_fastq2().path,
            outdir=outdir)
        )
        with open(self.outfile().path,'w') as fp:
            fp.write('Completed task {task} at {dt}'.format(
                task=self.get_instance_name(),
                dt=str(datetime.now())
            ))


class TrimReads(sl.Task):
    in_fastq1 = None
    in_fastq2 = None
    method = sl.Parameter(default="sickle")

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                "sickle_trimmed_{}".format(os.path.basename(self.in_fastq1().path))))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                "sickle_trimmed_{}".format(os.path.basename(self.in_fastq2().path))))

    def run(self):
        self.ex("sickle pe -g -f {reads1} -r {reads2} -o {out1} -p {out2} -t sanger -s /dev/null".format(
            reads1=self.in_fastq1().path,
            reads2=self.in_fastq2().path,
            out1=self.out_fastq1().path,
            out2=self.out_fastq2().path))


class CreateIndexForContamRemoval(sl.Task):

    # TODO: Figure out how to track if ref has been indexed before

    genomes = luigi.ListParameter(default=["/Users/ipqb/Downloads/GCF_000195995.1_ASM19599v1_genomic.fna.gz",
                                           "/Users/ipqb/Downloads/GCF_000005845.2_ASM584v2_genomic.fna.gz"])
    ref_info_subdir = sl.Parameter(default= "ref_data")
    # genomes = luigi.ListParameter(default=[
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/cat_masked.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/dog_masked.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/mouse_masked.fa.gz"])

    def outfile(self):
        method = self.workflow_task.contaminant_removal_method
        outdir = os.path.join(self.workflow_task.working_dir, self.ref_info_subdir, method)
        if method == "bbsplit":
            ret_val = {}
            for genome in self.genomes:
                label = genome.split("/")[-1].split(".")[0] # Do this more robustly
                ret_val[str(genome)] = sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir), "ref", label))
            return ret_val

    def run(self):
        method = self.workflow_task.contaminant_removal_method
        os.mk
        if method == "bbsplit":
            for genome in self.genomes:
                label = genome.split("/")[-1].split(".")[0] # Do this more robustly
                self.ex("bbsplit.sh ref={genome} path={path} -Xmx16g".format(
                    genome=genome,
                    path=self.out_dirs()[genome].path
                ))
        elif method == "bowtie2":
            self.ex("cat seqs > tmp/combo.fa;"
                    "bowtie2-build tmp/combo.fa ref")


class ContaminantRemoval(sl.Task):

    # TODO: finish implementing dependency on index being made beforehand

    in_fastq1 = None
    in_fastq2 = None
    ref = None # Ref has to handle both bt2 and bbsplit index
    method = sl.Parameter()
    subdir = sl.Parameter(default="decontaminated")

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                self.subdir,
                                                os.path.basename(self.in_fastq1().path)))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                self.subdir,
                                                os.path.basename(self.in_fastq2().path)))

    def run(self):
        self.out_dir().target.fs.mkdir(self.out_dir().path)
        if self.method == "bowtie2":
            self.ex("bowtie2 -x {ref} -1 {reads1} -2 {reads2} \
                --un-conc-gz \
                {out}%.fq.gz > /dev/null 2>{logs}/{sample}/bowtie2_decontamination.log".format(
                reads1=self.in_fastq1().path,
                reads2=self.in_fastq2().path,
                out=os.path.dirname(self.out_fastq1().path),
                ref=self.ref,
                logs=os.path.join(self.workflow_task.working_dir,
                                  self.workflow_task.tool_log_subdir),
                sample=self.workflow_task.prefix

            ))
        elif self.method == "bbsplit":

            self.ex('''
            bbsplit.sh \
                in={reads1} \
                in2={reads2} \
                outu={out1} \
                outu2={out2} \
                statsfile={logs}/{sample}/bbsplit_decontamination.log \
                minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch \
                fast minhits=2 qtrim=rl trimq=10 untrim -Xmx140g \
                path={path}
                '''.format(
                reads1=self.in_fastq1().path,
                reads2=self.in_fastq2().path,
                out1=self.out_fastq1().path,
                out2=self.out_fastq2().path,
                path=os.path.basename(self.ref().path),
                logs=os.path.join(self.workflow_task.working_dir,self.workflow_task.tool_log_subdir),
                sample=self.workflow_task.prefix
            )
                    )
        else:
            raise_(ValueError, "Unimplemented host removal method chosen", traceback)


class ReSyncPairs(sl.Task):
    in_fastq1 = None
    in_fastq2 = None
    subdir = sl.Parameter(default="resynced")
    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                self.subdir,
                                                self.in_fastq1().task.instance_name,
                                                os.path.basename(self.in_fastq1().path)))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                self.subdir,
                                                self.in_fastq2().task.instance_name,
                                                os.path.basename(self.in_fastq2().path)))

    def run(self):
        self.ex('''
                repair.sh \
	            in1={reads1} \
	            in2={reads2} \
	            out1={out1} \
	            out2={out2}
        '''.format(
            reads1=self.in_fastq1().path,
            reads2=self.in_fastq2().path,
            out1=self.out_fastq1().path,
            out2=self.out_fastq2().path

        ))




class RunMetaphlan2(sl.Task):
    in_fastq1 = None
    in_fastq2 = None
    nproc = sl.Parameter(default=1)
    subdir = sl.Parameter(default="metaphlan2_out")

    def outfile(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                               self.subdir,
                                               self.workflow_task.prefix+"_metaphlan2_output.txt"))

    def run(self):
        self.ex('''
        metaphlan2.py \
        {reads1},{reads2} \
        --no_map -t rel_ab_w_read_stats \
        --sample_id {prefix} \
        --nproc {nproc} --input_type fastq -o {out}
        '''.format(
            reads1=self.in_fastq1().path,
            reads2=self.in_fastq2().path,
            out=self.outfile().path,
            nproc=self.nproc,
            prefix=self.workflow_task.prefix

        ))


class RunMIDAS(sl.Task):
    midas_subtask = sl.Parameter()
    in_fastq1 = None
    in_fastq2 = None
    nproc = sl.Parameter(default=1)
    subdir = sl.Parameter(default="midas_out")

    def outfile(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.working_dir),
                                                self.subdir,
                                                self.workflow_task.prefix, self.midas_subtask + "sentinel.txt"))

    def run(self):
        self.ex('''
        run_midas.py {subtask} \
        {outdir} \
        -1 {in1} \
        -2 {in2} \
        -t {nproc}
        '''.format(
            subtask=self.subtask,
            outdir=os.path.dirname(self.out_put().path),
            in1=self.in_fastq1().path,
            in2=self.in_fastq2().path,
            nproc=self.nproc
        ))
        with open(self.out_put().path, 'w') as fp:
            fp.write('Completed task {task} at {dt}'.format(
                task=self.get_instance_name(),
                dt=str(datetime.now())
            ))






if __name__ == '__main__':
    arguments = docopt(__doc__)
    log.debug(str(arguments))
    luigi.run(local_scheduler=True, main_task_cls=SingleSampleWorkflow,
              cmdline_args=['--scheduler-host=localhost',
                            f'--working-dir={arguments["<working_dir>"]}',
                            f'--prefix={arguments["<prefix>"]}']
              )
