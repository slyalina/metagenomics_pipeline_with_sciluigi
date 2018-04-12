# coding=utf-8
"""
A basic metagenomics workflow using SciLuigi.

Usage: workflow.py <working_dir> <prefix>
    [--midas_db=<midas_db_path>] [--ref_data_dir=<ref_dir>] [--genome_to_filter=<human.fa>]...


Options:
--midas_db=<midas_db_path>     Location of the MIDAS database
--ref_data_dir=<ref_dir>       Location of where precomputed indices for filtering may exist.
                               If they don't, they'll be created
--genome_to_filter=<human.fa>  Can be specified multiple times to filter against multiple fastas


"""

import hashlib
import json
import six
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
BUF_SIZE = 128 * 1024


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


def get_md5_of_unordered_fileset(list_of_filepaths):
    hashes = [get_md5(filepath) for filepath in list_of_filepaths]
    return hashlib.md5(";".join(sorted(hashes)).encode('utf-8')).hexdigest()


class AutoOutput(object):
    """
    This mixin is adapted from https://github.com/pyannote/pyannote-workflows/blob/master/pyannote_workflows/utils.py
    """

    def _output_from_hash(self):

        # working directory within which all automatic outputs will be stored
        workdir = self.workflow_task.workdir

        description = {}

        # add one {key: value} per in_xxxx method
        # key = 'in_xxxx'
        # value = F(in_xxxx().path)
        for attrname, attrval in six.iteritems(self.__dict__):
            if 'in_' == attrname[0:3]:
                path = attrval().path
                if path.startswith(workdir):
                    path = path[len(workdir):]
                description[attrname] = path

        # add one {key: value} per task parameter
        # key = parameter name
        # value = parameter value
        params = self.get_params()
        params = [name for name, _ in params]
        for param_name in params:
            # do not take 'instance_name' and 'workflow_task' into account
            if param_name in ['instance_name', 'workflow_task']:
                continue
            description[param_name] = getattr(self, param_name)

        # hash the resulting dictionary
        digest = hashlib.sha1(
            json.dumps(description, sort_keys=True).encode('utf-8')).hexdigest()

        # generate out_put path automatically
        template = '{workdir}/{workflow_name}/{instance_name}/{digest}'

        output_path = template.format(
            workdir=workdir,
            instance_name=self.instance_name,
            workflow_name=self.workflow_task.__class__.__name__,
            digest=digest)

        # save it into parent workflow task
        if not hasattr(self.workflow_task, 'auto_output'):
            self.workflow_task.auto_output = {}
        self.workflow_task.auto_output[self.instance_name] = output_path

        return output_path


class AutoTxt(AutoOutput):
    def out_put(self):
        # automagically get out_put path
        path = self._output_from_hash() + ".txt"
        target = sl.TargetInfo(self, path)
        target.target.makedirs()
        return target


class AutoSentinel(AutoOutput):
    def log_info(self, fp):
        json.dump({"Task": self.task_id,
                   "Task type": self.get_task_family(),
                   "Time": str(datetime.now())},
                  fp)

    def out_put(self):
        # automagically get out_put path
        path = self._output_from_hash() + ".json"
        target = sl.TargetInfo(self, path)
        target.target.makedirs()
        return target


class AutoPairedFastq(AutoOutput):
    def out_fastq1(self):
        path = self._output_from_hash()
        target = sl.TargetInfo(self, os.path.join(os.path.dirname(path),
                                                  os.path.basename(self.in_fastq1().path)))
        target.target.makedirs()
        return target

    def out_fastq2(self):
        path = self._output_from_hash()
        target = sl.TargetInfo(self, os.path.join(os.path.dirname(path),
                                                  os.path.basename(self.in_fastq2().path)))
        target.target.makedirs()
        return target


class SingleSampleWorkflow(sl.WorkflowTask):
    midas_db = sl.Parameter()
    workdir = sl.Parameter()
    prefix = sl.Parameter()
    tool_log_subdir = sl.Parameter(default="tool_logs")
    contaminant_removal_method = sl.Parameter(default="bbsplit")
    filter_genomes = luigi.ListParameter()
    ref_info_dir = sl.Parameter()
    ref_combo_hash = sl.Parameter()

    # genomes = luigi.ListParameter(default=[
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/cat_masked.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/dog_masked.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz",
    #             "/pollard/home/slyalina/work/reference/masked_genomes_for_filtering/mouse_masked.fa.gz"])


    def workflow(self):

        # Initialize tasks:
        input_task = self.new_task("input", Input,
                                   prefix=self.prefix)
        fastqc_task_before = self.new_task('fastQC_before', Fastqc)
        sickle_task = self.new_task('sickle', TrimReads)
        fastqc_task_after = self.new_task('fastQC_after', Fastqc)
        resync_task = self.new_task("resync", ReSyncPairs)
        metaphlan_task = self.new_task("metaphlan", RunMetaphlan2)
        midas_task = self.new_task("midas_species", RunMIDAS,
                                   midas_subtask="species")
        if len(self.filter_genomes) > 0:
            index_task = self.new_task("ref_index", CreateIndexForContamRemoval)
            decontam_task = self.new_task("decontam", ContaminantRemoval)

        # Connecting outputs to inputs:
        sickle_task.in_fastq1 = input_task.out_fastq1
        sickle_task.in_fastq2 = input_task.out_fastq2

        fastqc_task_before.in_fastq1 = input_task.out_fastq1
        fastqc_task_before.in_fastq2 = input_task.out_fastq2

        fastqc_task_after.in_fastq1 = sickle_task.out_fastq1
        fastqc_task_after.in_fastq2 = sickle_task.out_fastq2

        resync_task.in_fastq1 = sickle_task.out_fastq1
        resync_task.in_fastq2 = sickle_task.out_fastq2

        end_tasks = [fastqc_task_after, fastqc_task_before, metaphlan_task]

        if len(self.filter_genomes) > 0:
            decontam_task.ref = index_task.out_put
            decontam_task.in_fastq1 = resync_task.out_fastq1
            decontam_task.in_fastq2 = resync_task.out_fastq2

            metaphlan_task.in_fastq1 = decontam_task.out_fastq1
            metaphlan_task.in_fastq2 = decontam_task.out_fastq2

            midas_task.in_fastq1 = decontam_task.out_fastq1
            midas_task.in_fastq2 = decontam_task.out_fastq2

            end_tasks.append(index_task)
        else:
            metaphlan_task.in_fastq1 = resync_task.out_fastq1
            metaphlan_task.in_fastq2 = resync_task.out_fastq2

            midas_task.in_fastq1 = resync_task.out_fastq1
            midas_task.in_fastq2 = resync_task.out_fastq2

        # Return the last task(s) in the workflow chain.
        return end_tasks


class MultiSampleWorkflow(sl.WorkflowTask):
    midas_db = sl.Parameter()
    dataset_description = sl.Parameter()
    working_dir = sl.Parameter()
    contaminant_removal_method = sl.Parameter(default="bbsplit")
    filter_genomes = luigi.ListParameter()
    ref_info_dir = sl.Parameter()
    ref_combo_hash = sl.Parameter()

    def workflow(self):
        dataset_spec = json.load(self.dataset_description)
        tasks = []
        for sample in dataset_spec["samples"]:
            wf = self.new_task('SampleWorkflow', SingleSampleWorkflow,
                               working_dir=self.working_dir,
                               prefix=sample["prefix"])
            tasks.append(wf)
        return tasks


class Input(sl.ExternalTask):

    def out_fastq1(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.workdir),
                                                str(self.workflow_task.prefix) + "_1.fq.gz"))

    def out_fastq2(self):
        return sl.TargetInfo(self, os.path.join(str(self.workflow_task.workdir),
                                                str(self.workflow_task.prefix) + "_2.fq.gz"))


class Fastqc(sl.Task, AutoSentinel):
    in_fastq1 = None
    in_fastq2 = None

    def run(self):
        output = self.out_put()
        outdir = os.path.dirname(output.path)
        self.ex('fastqc --noextract -o "{outdir}" "{in1}" "{in2}"'.format(
            in1=self.in_fastq1().path,
            in2=self.in_fastq2().path,
            outdir=outdir)
        )
        with output.open("w") as fp:
            self.log_info(fp)


class TrimReads(sl.Task, AutoPairedFastq):
    in_fastq1 = None
    in_fastq2 = None
    method = sl.Parameter(default="sickle")

    def run(self):
        if self.method == "sickle":
            self.ex('sickle pe -g -f "{reads1}" -r "{reads2}" -o "{out1}" -p "{out2}" -t sanger -s /dev/null'.format(
                reads1=self.in_fastq1().path,
                reads2=self.in_fastq2().path,
                out1=self.out_fastq1().path,
                out2=self.out_fastq2().path))
        elif self.method == "fastq-mcf":
            pass  # TODO: Add fastq-mcf
        elif self.method == "bbduk":
            self.ex('bbduk.sh in="{reads1}" in2="{reads2}" '
                    'out="{out1}" out2="{out2}" ref=adapters,phix stats="{outstats}'.format(
                        reads1=self.in_fastq1().path,
                        reads2=self.in_fastq2().path,
                        out1=self.out_fastq1().path,
                        out2=self.out_fastq2().path,
                        outstats=os.path.join(os.path.dirname(self.out_fastq2().path), "stats.txt")))
        else:
            raise_(ValueError, "Unimplemented trimming method chosen", traceback)


class CreateIndexForContamRemoval(sl.Task, AutoSentinel):

    def run(self):
        method = self.workflow_task.contaminant_removal_method
        output = self.out_put()
        outdir = os.path.dirname(output.path)
        # Location of either existing ref data with subdirs per method. Must have trailing slash:
        ref_dir_for_wf = os.path.join(self.workflow_task.ref_info_dir, self.workflow_task.ref_combo_hash, method, "")

        # If the given ref dir has necessary indices, just symlink them to the local dir
        if not os.path.exists(ref_dir_for_wf):
            # TODO: check more robustly that index actually exists
            output.target.fs.mkdir(ref_dir_for_wf)
            # Else make them
            if method == "bbsplit":
                self.ex("bbsplit.sh ref={genome} path={path}".format(
                    genome=",".join(self.workflow_task.filter_genomes),
                    path=ref_dir_for_wf
                ))
            elif method == "bowtie2":
                self.ex("bowtie2-build {genome} {path}".format(
                    genome=",".join(self.workflow_task.filter_genomes),
                    path=ref_dir_for_wf)
                )
        output.target.fs.mkdir(outdir)
        os.symlink(os.path.abspath(ref_dir_for_wf), os.path.join(outdir, method))
        with output.open("w") as fp:
            self.log_info(fp)


class ContaminantRemoval(sl.Task, AutoPairedFastq, AutoTxt):
    in_fastq1 = None
    in_fastq2 = None
    ref = None  # Ref is the sentinel file that is in the dir created via symlink pointing at ref dirs for each method

    def run(self):
        method = self.workflow_task.contaminant_removal_method
        ref_dir = os.path.join(os.path.dirname(self.ref().path), method, "")  # Again have to ensure trailing slash
        outtxt = self.out_put()
        outdir = os.path.dirname(self.out_fastq1().path)

        if method == "bowtie2":
            self.ex("bowtie2 -x {ref} -1 {reads1} -2 {reads2} \
                    --un-conc-gz \
                    {out}/tmp%.fq.gz > /dev/null 2>{log};"
                    "mv {out}/tmp1.fq.gz {out1};"
                    "mv {out}/tmp1.fq.gz {out2}".format(
                        reads1=self.in_fastq1().path,
                        reads2=self.in_fastq2().path,
                        out=outdir,
                        ref=ref_dir,
                        out1=self.out_fastq1(),
                        out2=self.out_fastq2(),
                        log=outtxt.path,
                        sample=self.workflow_task.prefix
                    ))

        elif method == "bbsplit":

            self.ex('''
            bbsplit.sh \
                in={reads1} \
                in2={reads2} \
                outu={out1} \
                outu2={out2} \
                statsfile={log} \
                minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch \
                fast minhits=2 qtrim=rl trimq=10 untrim -Xmx140g \
                path={path}
                '''.format(
                reads1=self.in_fastq1().path,
                reads2=self.in_fastq2().path,
                out1=self.out_fastq1().path,
                out2=self.out_fastq2().path,
                path=ref_dir,
                log=outtxt.path,
            ))
        else:
            raise_(ValueError, "Unimplemented host removal method chosen", traceback)


class ReSyncPairs(sl.Task, AutoPairedFastq):
    in_fastq1 = None
    in_fastq2 = None

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


class RunMetaphlan2(sl.Task, AutoTxt):
    in_fastq1 = None
    in_fastq2 = None
    nproc = sl.Parameter(default=1)

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
            out=self.out_put().path,
            nproc=self.nproc,
            prefix=self.workflow_task.prefix

        ))


class RunMIDAS(sl.Task, AutoSentinel):
    midas_subtask = sl.Parameter()
    in_fastq1 = None
    in_fastq2 = None
    nproc = sl.Parameter(default=1)

    def run(self):
        output = self.out_put()
        outdir = os.path.dirname(output.path)
        self.ex('''
        run_midas.py {subtask} \
        {outdir} \
        -d \
        -1 {in1} \
        -2 {in2} \
        -t {nproc} --remove_temp
        '''.format(
            subtask=self.midas_subtask,
            outdir=outdir,
            in1=self.in_fastq1().path,
            in2=self.in_fastq2().path,
            nproc=self.nproc
        ))
        with output.open("w") as fp:
            self.log_info(fp)


if __name__ == '__main__':
    arguments = docopt(__doc__)
    if not arguments["--midas_db"]:
        arguments["--midas_db"] = os.getenv("MIDAS_DB")
    if not arguments["--midas_db"]:
        raise_(ValueError, "Need to provide location of MIDAS_DB (or set it as an env var)", traceback)

    if not arguments["--ref_data_dir"]:
        arguments["--ref_data_dir"] = os.path.join(arguments["<working_dir>"], "ref_data")

    if len(arguments["--genome_to_filter"]) > 0:
        ref_hash = get_md5_of_unordered_fileset(arguments["--genome_to_filter"])
    else:
        ref_hash = None
    log.info(str(arguments))

    luigi.run(local_scheduler=True, main_task_cls=SingleSampleWorkflow,
              cmdline_args=['--scheduler-host=localhost',
                            '--workdir={}'.format(arguments["<working_dir>"]),
                            '--prefix={}'.format(arguments["<prefix>"]),
                            '--midas-db={}'.format(arguments["--midas_db"]),
                            '--ref-info-dir={}'.format(arguments["--ref_data_dir"]),
                            '--ref-combo-hash={}'.format(ref_hash),
                            '--filter-genomes={}'.format(json.dumps(arguments["--genome_to_filter"]))]
              )
