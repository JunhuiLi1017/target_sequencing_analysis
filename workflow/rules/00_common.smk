import pandas as pd
import glob
import os
import yaml
from snakemake.utils import validate

configfile: "config/config.yaml"
validate(config, schema = "../schemas/config.schema.yaml")

with open(config["resource"]) as f:
    resource = yaml.safe_load(f)
validate(resource, schema = "../schemas/resource.schema.yaml")

units = pd.read_table(config["units"], dtype = str).set_index(
    ["sample", "library", "flowlane"], drop = False
)

container_image = {}
if "container" in config:
    try:
        with open(config["container"]) as f:
            container_image = yaml.safe_load(f) or {}
    except Exception as e:
        print(f"Warning: Could not load container config from {config['container']}: {e}")
        container_image = {}

units['trim_front1'] = pd.to_numeric(units['trim_front1'], errors='coerce').fillna(0).astype(int)
units['trim_front2'] = pd.to_numeric(units['trim_front2'], errors='coerce').fillna(0).astype(int)
units['trim_tail1'] = pd.to_numeric(units['trim_tail1'], errors='coerce').fillna(0).astype(int)
units['trim_tail2'] = pd.to_numeric(units['trim_tail2'], errors='coerce').fillna(0).astype(int)
units['pcr_based'] = units['pcr_based'].map({'Yes': True, 'No': False, 'TRUE': True, 'FALSE': False, '1': True, '0': False, 'True': True, 'False': False}).astype(bool)
units['is_umi'] = units['is_umi'].map({'Yes': True, 'No': False, 'TRUE': True, 'FALSE': False, '1': True, '0': False, 'True': True, 'False': False}).astype(bool)


validate(units, schema = "../schemas/units.schema.yaml")

sample_num=units['sample'].unique().tolist()
Sample=units['sample']
Outpath = config['outpath']
intervals_dir=config["interval"]
Ref_version=config["ref_version"]
Caller=config["callers"]

# Get the list of all interval files
interval_files = glob.glob(os.path.join(intervals_dir, "*.intervals.list"))
CHROMOSOMES = [os.path.basename(f).replace(".intervals.list", "") for f in interval_files]
standard_order = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
    "chr20", "chr21", "chr22", "chrX", "chrY"
]

CHROMOSOMES = sorted(CHROMOSOMES, key=lambda x: standard_order.index(x) if x in standard_order else len(standard_order))

def is_pair_end(wildcards):
    try:
        fastqs = units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), ["fq1", "fq2"]].dropna()
        if len(fastqs) == 2:
            return True
        return False
    except KeyError:
        raise KeyError(f"No entry in units for sample={wildcards.sample}, library={wildcards.library}, flowlane={wildcards.flowlane}")

def get_fastq(wildcards):
    """Get fastq files for given sample, library, and flowlane."""
    fastqs = units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), ["fq1", "fq2"]].dropna()
    if is_pair_end:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_fastqc_input(wildcards):
    if is_pair_end:
        return [
            f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz",
            f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R2.fastq.gz"
        ]
    else:
        return [f"{Outpath}/01_multiqc/fastp/{wildcards.sample}_{wildcards.library}_{wildcards.flowlane}.R1.fastq.gz"]

def get_multiqc_input(wildcards):
    if is_pair_end:
        GROUP = ["R1", "R2"]
    else:
        GROUP = ["R1"]
    return expand(
            [
                f"{Outpath}/01_multiqc/fastqc/{{sample}}_{{library}}_{{flowlane}}.{{group}}_fastqc.html",
                f"{Outpath}/01_multiqc/fastqc/{{sample}}_{{library}}_{{flowlane}}.{{group}}_fastqc.zip"
            ],
            sample=[u.sample for u in units.itertuples()],
            library=[u.library for u in units.itertuples()],
            flowlane=[u.flowlane for u in units.itertuples()],
            group = GROUP
        )

def get_vcf_inputs(wildcards):
    callers = config["callers"]
    input_list = []
    for caller in callers:
        input_list.append(f"{Outpath}/03_variants/{caller}/04_pass/{wildcards.sample}.{Ref_version}.pass.vcf.gz")
    return input_list

def get_raw_bam(wildcards):
    sample_units = units.loc[units['sample'] == wildcards.sample]
    bam_files = [
        f"{Outpath}/02_Map/bwa/sort/{sample}/{row['sample']}_{row['library']}_{row['flowlane']}.sort.bam"
        for _, row in sample_units.iterrows()]
    return bam_files

def get_reads_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{library}_{flowlane}\tSM:{sample}\tPL:{platform}\tLB:{library}'".format(
        sample=wildcards.sample,
        library=wildcards.library,
        flowlane=wildcards.flowlane,
        platform=units.loc[(wildcards.sample, wildcards.library, wildcards.flowlane), "platform"]
    )

def get_venn_input(wildcards, input, output):
    callers = config["callers"]
    txt = f"{Outpath}/04_overlap/01_plot/{wildcards.sample}.all_callers.txt"
    lines = []
    for i, (caller, vcf) in enumerate(zip(callers, input)):
        short = caller.split("_", 1)[-1]
        # Extract variant id and prepend caller name
        lines.append(
            f"zcat {vcf} | grep -v '^#' | cut -f 1,2,4,5 | sed 's/\\t/_/g' | awk '{{print \"{short}\\t\"$1}}' >> {output.txt}"
        )
    return " && ".join(lines)


def get_scatter_cmd(wildcards, input, output):
    # input is now a list of VCF files
    callers = config["callers"]
    txt = output[0] if isinstance(output, (list, tuple)) else output.txt
    lines = []
    for i, (caller, vcf) in enumerate(zip(callers, input)):
        # Remove prefix for output variable name, e.g., "01_mutect2" -> "mutect2"
        short = caller.split("_", 1)[-1]
        if short == "pisces":
            lines.append(f"zcat {vcf} | cut -f 1,2,4,5,10 | sed 's/:/\t/g' | cut -f 1-4,8,9 | grep -v '^#' | awk '{{print \"pisces\\t\"$1\"_\"$2\"_\"$3\"_\"$4\"\t\"$5\"\t\"$6}}' >> {txt}")
        elif short == "mutect2":
            lines.append(f"zcat {vcf} | cut -f 1,2,4,5,10 | sed 's/:/\t/g' | cut -f 1-4,7,8 | grep -v '^#' | awk '{{print \"mutect2\\t\"$1\"_\"$2\"_\"$3\"_\"$4\"\t\"$6\"\t\"$5}}' >> {txt}")
        elif short == "replow":
            lines.append(f"zcat {vcf} | grep -v '^#' | sed 's/;/\t/g' | cut -f 1,2,4,5,8,12 | sed 's/DP=//' | sed 's/EBAF=//g' | awk '{{print \"replow\\t\"$1\"_\"$2\"_\"$3\"_\"$4\"\t\"$5\"\t\"$6}}' >> {txt}")
    return " && ".join(lines)


def get_dedup_files(wildcards):
    # Retrieve the UMI status from your sample sheet (assumed to be a pandas DataFrame)
    is_umi = units.loc[wildcards.sample, "is_umi"].any()
    
    # This logic forces Snakemake to pick the correct path
    if is_umi == True:
        return f"{wildcards.outpath}/02_map/03_rmdup/{wildcards.sample}/{wildcards.sample}.03_rmdup.umi.bam"
    else:
        return f"{wildcards.outpath}/02_map/03_rmdup/{wildcards.sample}/{wildcards.sample}.03_rmdup.gatk.bam"

def get_bam_for_stats(wildcards):
    # If we are looking at the sorted (not deduped) BAM, the path is standard
    if wildcards.bam_type == "02_sort":
        return f"{wildcards.outpath}/02_map/02_sort/{wildcards.sample}/{wildcards.sample}.02_sort.bam"
    
    # If we are looking at the deduped BAM, we must use the UMI/GATK logic
    if wildcards.bam_type == "03_rmdup":
        is_umi = units.loc[wildcards.sample, "is_umi"].any()
        if is_umi:
            # Match the output path from rule remove_dup_umi
            return f"{wildcards.outpath}/02_map/03_rmdup/{wildcards.sample}/{wildcards.sample}.03_rmdup.umi.bam"
        else:
            # Match the output path from rule remove_dup_gatk
            return f"{wildcards.outpath}/02_map/03_rmdup/{wildcards.sample}/{wildcards.sample}.03_rmdup.gatk.bam"