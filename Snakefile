from os.path import realpath
from os.path import split as pathsplit
import subprocess
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
import sys

# Block annoying warnings
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

# META
__author__ = "Matt Lawlor"

# SETUP
shell.executable("/bin/bash")

def is_pe(name):
    fqs = MY_SAMPLES.get(name,None).get("fastq",None)
    return(len(fqs) == 2)

# DETERMINE REMOTE OR LOCAL RESOURCE
def determine_resource(path):
    if "gs://" in path:
         return GSRemoteProvider().remote(path.replace("gs://",""))
    elif "ftp://" in path:
         return FTPRemoteProvider().remote(path)
    elif "s3://" in path:
         return S3RemoteProvider().remote(path.replace("s3://",""))
    elif "http://" in path:
         return HTTPRemoteProvider().remote(path.replace("http://",""))
    elif "https://" in path:
         return HTTPRemoteProvider().remote(path.replace("https://",""))
    else:
        return path

PROJECT_NAME = config.get("project","rna")
MY_SAMPLES = config.get("samples",None)

GENCODE = config["txome"].get("gencode",False)
SALMON_BOOTSTRAPS = config["salmon_config"].get("bootstraps",50)

RUN_GALIGN = False
RUN_SALMON = True

target_files=[]

# we pretty much always want these
target_files.append("ref/fa/txfa")

if RUN_SALMON:
    target_files.append("salmon/{p}_salmon_dds.rds".format(p=PROJECT_NAME))

ruleorder: trim_se > trim_pe

rule target:
    input:
        target_files

include: "rules/trimming.smk"


rule get_txome_fasta:
    output:
        "ref/fa/txfa"
    params:
        fl=config["txome"].get("fasta",None)
    shell:
        "curl -J -L {params.fl} > {output[0]}"

rule get_gtf:
    output:
        "ref/tx.gtf"
    params:
        fl=config["genome"].get("gtf",None)
    shell:
        "curl -J -L {params.fl} > {output[0]}"

rule salmon_index:
    input:
        "ref/fa/txfa"
    output:
        directory("ref/idx/salmon/transcripts_index")
    params:
        k = 31,
        gencode = "--gencode" if GENCODE else ""
    conda:
        "envs/salmon.yaml"
    singularity:
        "docker://quay.io/biocontainers/salmon:0.14.2--ha0cc327_0"
    threads:
        8
    shell:
        "salmon index -t {input} -i {output} "
        "-k {params.k} {params.gencode} -p {threads}"


def get_proper_ended_salmon_call(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "-r {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "-1 {r1} -2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))



rule salmon_quant:
    input:
        fq  = lambda wc: get_fqs_for_aln(wc.s),
        idx = "ref/idx/salmon/transcripts_index",
    output:
        "salmon/{s}/quant.sf"
    params:
        bootstraps = "--numBootstraps {n}".format(n=SALMON_BOOTSTRAPS),
        read_arg = lambda wc: get_proper_ended_salmon_call(wc.s)
    threads:
        8
    conda:
        "envs/salmon.yaml"
    singularity:
        "docker://quay.io/biocontainers/salmon:0.14.2--ha0cc327_0"
    shell:
        "salmon quant --libType A -i {input.idx} "
        "{params.read_arg} "
        "-p {threads} {params.bootstraps} "
        "--seqBias --gcBias --posBias "
        "--validateMappings -o salmon/{wildcards.s}"

rule tx2gene_rds:
    input:
        gtf = rules.get_gtf.output
    output:
        "salmon/tx2gene.rds".format(p=PROJECT_NAME)
    conda:
        "envs/make_salmon_dds.yaml"
    script:
        "scripts/make_tx2gene.R"

rule salmon_dds:
    input:
        qsf = expand("salmon/{s}/quant.sf",s=MY_SAMPLES),
        gtf = rules.get_gtf.output,
        tx2g = rules.tx2gene_rds.output
    output:
        "salmon/{p}_salmon_dds.rds".format(p=PROJECT_NAME)
    params:
        names = [x for x in MY_SAMPLES],
        conds = {x:MY_SAMPLES[x]["condition"] for x in MY_SAMPLES},
        ctrl = config.get("control_level",None)
    conda:
        "envs/make_salmon_dds.yaml"
    script:
        "scripts/make_salmon_dds.R"
