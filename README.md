# Target Sequence Analysis Pipeline

## 1. Introduction

The Target Sequence Analysis Pipeline is a comprehensive Snakemake-based workflow designed for analyzing targeted sequencing data. This pipeline provides end-to-end processing from raw FASTQ files to variant calling and analysis, with support for multiple variant callers and comprehensive quality control.

The pipeline is specifically designed for:
- **Targeted sequencing data** (exome, gene panels, custom targets)
- **Multiple variant callers** (Mutect2, Pisces, RePlow)
- **Flexible sample configuration** with per-sample parameters
- **Comprehensive quality control** and reporting
- **High-performance computing** environments

## 2. Pipeline Summary

The pipeline consists of the following main stages:
### stage 0: requirment for the fastq with umi information
- **merge_umi_to_header.py**: use the customed script to set umi sequence followed to the header of reads separated by ":"  usega:  python merge_umi_to_header.py -u umi.fastq.gz -i input.fastq.gz -o output.fastq.gz

for example:

    --------------------------------
    - first 2 reads of umi.fastq.gz
    --------------------------------

    @LH00160:653:2373YKLT4:8:1101:1380:1098 2:N:0:GATGTGTG+TAGCCATG
    ACCAAGGCC
    +
    9IIIIII9I
    @LH00160:653:2373YKLT4:8:1101:4309:1098 2:N:0:GATGTGTG+TAGCCATG
    CGGGTTGAG
    +
    IIIII9III

    --------------------------------
    - first 2 reads in input.fastq.gz
    --------------------------------

    @LH00160:653:2373YKLT4:8:1101:1380:1098 1:N:0:GATGTGTG+TAGCCATG

    @LH00160:653:2373YKLT4:8:1101:4309:1098 1:N:0:GATGTGTG+TAGCCATG

    --------------------------------
    - first 2 reads in output.fastq.gz
    --------------------------------

    @LH00160:653:2373YKLT4:8:1101:1380:1098:ACCAAGGCC 1:N:0:GATGTGTG+TAGCCATG

    @LH00160:653:2373YKLT4:8:1101:4309:1098:CGGGTTGAG 1:N:0:GATGTGTG+TAGCCATG


### Stage 1: Quality Control & Preprocessing
- **fastp(https://github.com/opengene/fastp)**
    A tool designed to provide ultrafast all-in-one preprocessing and quality control for FastQ data.

    ** Key functions:

    1. filter out bad reads (too low quality, too short, or too many N...)
    
    2. trim all reads in front and tail

    3. cut adapters. Adapter sequences can be automatically detected, which means you don't have to input the adapter sequences to trim them.


- **fastqc(https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**: Raw read quality assessment
- **MultiQC(https://docs.seqera.io/multiqc/#:~:text=MultiQC%20is%20a%20reporting%20tool,log%20files%20and%20console%20outputs.)**: Aggregated QC reports

### Stage 2: Read Mapping & Processing
- **BWA-MEM**: Read alignment to reference genome
- **SAMtools**: BAM file processing and indexing
- **umitools**: duplication reads remove for samples with umi is true/yes
- **picard**: duplication reads remove for samples with umi is not true/yes
- **GATK(https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR)**: Base Quality Score Recalibration (BQSR)
    Systematic bias can originate from library preparation, sequencing, manufacturing defects in the flowcell chips, sequencer variation, and sequencing chemistry. It results in over- or underestimation of quality scores. GATK recalibrates base quality scores by building an error model using known covariates from all base calls (BaseRecalibrator), then applying adjustment to the dataset based on the model (ApplyBQSR).
- **R script with bedtools output**: coverage, depth, duplication rate

### Stage 3: Variant Calling
- **Mutect2**: GATK's somatic variant caller
- **Pisces**: Illumina's variant caller
- **RePlow**: Alternative variant caller
- **Variant filtering**: Quality-based filtering

### Stage 4: Analysis & Reporting
- **Variant overlap analysis**: Venn diagrams and scatter plots
- **Annotation**: Variant annotation and prioritization
- **Final report**: Comprehensive HTML report

## 3. Tools and Versions

| Tool | Version | Purpose |
|------|---------|---------|
| **Snakemake** | 7.0+ | Workflow management |
| **BWA** | 0.7.17 | Read alignment |
| **SAMtools** | 1.15+ | BAM processing |
| **GATK** | 4.6.1.0 | Variant calling and BQSR |
| **FastQC** | 0.11.9 | Quality control |
| **FastP** | 0.23.2 | Adapter trimming |
| **MultiQC** | 1.14 | Report aggregation |
| **Pisces** | Latest | Variant calling |
| **RePlow** | Latest | Variant calling |
| **R** | 4.2+ | Statistical analysis and plotting |
| **Python** | 3.8+ | Scripting and data processing |

## 4. Quick Start

### Prerequisites
- Snakemake 7.0 or higher
- Singularity or Conda for environment management
- Sufficient storage space for your data
- Access to reference genome files

### Installation
```bash
# Clone the repository
git clone <repository-url>
cd target_sequence_analysis

# Install dependencies (if using conda)
conda env create -f environment.yml
conda activate target_seq_analysis
```

### Configuration
1. **Edit config file**: Modify `config/config.yaml` with your paths and parameters
2. **Prepare units file**: Create a TSV file with sample information (see example below)
3. **Set up reference files**: Ensure all reference files are accessible

### Example units.tsv format:
```tsv
sample	library	flowlane	fq1	fq2	platform	trim_front1	trim_front2	trim_tail1	trim_tail2	pcr_based is_umi
SAMPLE1	LIB1	LANE1	/path/to/SAMPLE1_R1.fastq.gz	/path/to/SAMPLE1_R2.fastq.gz	ILLUMINA	0	0	0	0	No True
SAMPLE2	LIB2	LANE2	/path/to/SAMPLE2_R1.fastq.gz	/path/to/SAMPLE2_R2.fastq.gz	ILLUMINA	5	5	0	0	Yes False
```

### Running the Pipeline
```bash
# Dry run to check the workflow
snakemake -s workflow/Snakefile.hg38 --configfile config/config.hg38.yaml -n

# Run the pipeline
snakemake -s workflow/Snakefile.hg38 --configfile config/config.hg38.yaml -j 4

# Run on cluster (example with SLURM)
snakemake -s workflow/Snakefile.hg38 --configfile config/config.hg38.yaml -j 99 --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M"
```

## 5. Features

### 1. Flexible Sample Configuration
- **PCR-based parameter specification**: Each sample can have different `-f` (trim_front) and `-t` (trim_tail) parameters
- **Per-sample customization**: All trimming parameters are configurable in the units file
- **Dynamic parameter application**: Parameters are automatically applied based on sample metadata

### 2. Automatic Read Type Detection
- **Single-end detection**: Automatically detects single-end reads when `fq2` is missing
- **Pair-end support**: Full support for paired-end sequencing data
- **Dynamic workflow adaptation**: Pipeline adapts based on detected read type

### 3. Resource Management
- **Three-tier resource allocation**: High, medium, and low resource levels
- **Rule-specific resource assignment**: Each rule can use appropriate resource levels
- **Memory and thread optimization**: Efficient resource utilization for different computational tasks

### 4. Comprehensive Reporting
- **MultiQC integration**: Aggregated quality control reports
- **Variant overlap analysis**: Interactive Venn diagrams and scatter plots
- **Final HTML report**: Complete pipeline summary with all results
- **Interactive visualizations**: Dynamic plots for variant comparison

## 6. Contribution and Support

### Contributing
We welcome contributions! Please follow these steps:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

### Reporting Issues
- Use the GitHub issue tracker
- Include error messages and log files
- Provide minimal reproducible examples
- Specify your environment (OS, versions, etc.)

### Getting Help
- Check the documentation in the `docs/` folder
- Review existing issues on GitHub
- Contact the maintainers for specific questions

---

**Note**: This pipeline is designed for research use. Please ensure compliance with your institution's data handling policies and obtain necessary permissions for data processing.

