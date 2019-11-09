rule concat_fqs:
    input:
        lambda wc: [determine_resource(x) for x in config["samples"][wc.s]["fastq"][wc.end]]
    output:
        temp("fastq/{s}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

def get_fqs_for_trim(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.fq.gz"]
    else:
        return ["fastq/{s}_r1.fq.gz", "fastq/{s}_r2.fq.gz"]

def get_fqs_for_aln(x):
    if (len(config["samples"][x]["fastq"].keys()) == 1):
        return ["fastq/{s}_r1.trimmed.fq.gz"]
    else:
        return ["fastq/{s}_r1.trimmed.fq.gz", "fastq/{s}_r2.trimmed.fq.gz"]


def get_proper_ended_fastp_call(x):
    fqs = get_fqs_for_trim(x)
    if len(fqs) == 1:
        return "--in1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "--in1 {r1} --in2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))

def get_proper_ended_fastp_out(x):
    fqs = get_fqs_for_aln(x)
    if len(fqs) == 1:
        return "--out1 {r1}".format(r1=fqs[0].format(s=x))
    else:
        return "--out1 {r1} --out2 {r2}".format(r1=fqs[0].format(s=x), r2=fqs[1].format(s=x))

rule trim_se:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.s)
    output:
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"

rule trim_pe:
    input:
        fq = lambda wc: get_fqs_for_trim(wc.s)
    output:
        r1 = "fastq/{s}_r1.trimmed.fq.gz",
        r2 = "fastq/{s}_r2.trimmed.fq.gz",
        html = "fastq/{s}_fastp.html",
        json = "fastq/{s}_fastp.json"
    threads:
        2
    params:
        call_in = lambda wc: get_proper_ended_fastp_call(wc.s),
        call_out = lambda wc: get_proper_ended_fastp_out(wc.s)
    conda:
        "envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp {params.call_in} "
        "{params.call_out} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.s}_fastp"


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
