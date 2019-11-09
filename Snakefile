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

include: "rules.smk"
