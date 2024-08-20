RePlow
=======

**RePlow** is a Java based variant caller designed for detecting somatic single nucleotide variants (SNVs) from the replicated set of high-depth sequencing data. RePlow is highly specialized for the identification of somatic mutations with low variant allele frequency (VAF ~1%). RePlow accurately detects such low-level mutations based on the probabilistic model that jointly analyzes library-level replicates, regardless of the sequencing platform. The important features of RePlow are:

* On-the-fly estimation of the VAF distribution of background errors
* Improved low-level SNV calling based on the estimated error profiles and the observation concordance between replicates

Get the most recent version of RePlow here:
[Download](https://sourceforge.net/projects/replow/files/latest/download)
<br>

Prerequisite softwares
=======
**R packages** 

* [Rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html)
* [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/fitdistrplus.pdf)

<br>

Installation of RePlow
=======

RePlow was developed using JAVA JDK 8 64bit. To run RePlow, Java Runtime Environment (JRE) version 1.8.x or later is required. Download the most recent version of RePlow program in the code page and extract the gzipped archive.

    tar -zxvf RePlow-1.x.y.tar.gz

It contains the followings:

1. **RePlow.jar**
    : Executable JAR file of RePlow
+  **README.md**
    : Same as the wiki page
*  **lib/**
    : Directory containing other JAR libraries for RePlow
*  **testdata/**
    : Directory containing a small demo dataset for installation verification
*  **errorProfile.R**
    : R code to estimate error profiles in the variant detection process
*  **highMapQ.bed**
    : Predefined high mapQ regions for the estimation of background error rate. For computational efficiency, only shared regions between high mapQ and given target regions are used for the error rate estimation if the total size of target region > 1Mb. 
*  **RePlow.license**
    : Software licensing terms

Or you can directly build RePlow from the source code by yourself. To compile RePlow, Java JDK 8 and Maven 3.0+ are required. You'll get RePlow-1.x.y.jar in your folder if the source code is successfully compiled.

    # make a new source directory
    mkdir RePlow
    cd RePlow
    
    # get the source code
    git clone https://junhokim@git.code.sf.net/p/replow/code replow-code
    
    # build RePlow
    mvn package    
    
    # test running. You should correct the path for REF and RSCRIPT in run.sh
    cd testdata
    ./run.sh

<br>

Running RePlow
=======

To run RePlow, you can simply run the JAR executable file like below:

    java -jar RePlow.jar

The -h or -? options will bring the following usage of RePlow. If you see this, you are ready to run RePlow.

-------------------------------------

    RePlow: Replicate-based caller for the accurate detection of low-frequency SNV
    Version: 1.1.0
    Usage: java -jar RePlow.jar -r <reference_genome.fa> -b <replicate_sorted_indexed_BAM_list> -T <replicate_target_BED_file_list> -R <Rscript_path> [options]

Input options:

    -r, --reference_sequence <PATH>             Path for indexed reference genome FASTA [.fa/.fasta]
    -b, --bam_file_list <PATH>                  Comma-separated BAM file list
    -R, --Rscript <PATH>                        Path for Rscript (e.g. /usr/bin/Rscript)
    -T, --target_bed_list <PATH>                Comma-separated target BED file list

Analysis parameters for background error profiling:

    -q, --bQ_thres <INT>                        BaseQ threshold [15]
    -Q, --mQ_thres <INT>                        MapQ threshold [20]
    -N, --normal_BAM <PATH>                     BAM file of the matched-control sample
    -v, --germline_SNP <PATH>                   Germline SNP call (.vcf) (e.g. GATK SNP call)
    -V, --somatic_SNV <PATH>                    Comma-separated somatic SNV call list (e.g. MuTect SNV rawcall results)
    -f, --germlineBAFthres <FLOAT>              Threshold for maxBAF (germline-filtering for non-matched samples) [0.4]

Analysis parameters for variant detection:

    -n, --minBalleleCnt <INT>                   MinBalleleCnt [1]
    -p, --pval_for_normal <FLOAT>               P-value cutoff for matched normal (germline filter) [0.01]
    -e, --min_MOS_raction <FLOAT>               Required minimum fraction of MOS (MOS / total number of B alleles) to filter out sequencing errors [0.5]
    -d, --minDiseaseDepth <INT>                 Required minimum depth at the variant site from disease sample [50]
    -D, --minNormalDepth <INT>                  Required minimum depth at the variant site from matched-control sample [10]
    -y, --homopolymer_length <INT>              Length to determine homopolymeric site [4]
    -I, --minMisReadCnt <INT>                   Minimum number of B-allele to apply fraction-based filters [5]
    -c, --clustered_pos_thres <INT>             Threshold for the clustered position from the end of read [3]
    -C, --clusteredFractionThres <FLOAT>        Allowable B-allele fraction with clustered position [0.8]
    -S, --successiveFractionThres <FLOAT>       Allowable B-allele fraction with successive mutation [0.8]
    -M, --lowMQFraction <FLOAT>                 Allowable B-allele fraction with low mapping quality [1.0]
    -B, --lowBQFraction <FLOAT>                 Allowable B-allele fraction with low base-call quality [0.5]
    -g, --indelFraction <FLOAT>                 Allowable fraction of the reads with indels at the position [0.5]
    -G, --BalleleFractionWithIndel <FLOAT>      Allowable B-allele fraction with the indels within the indel offset [0.4]
    -i, --indel_offset <INT>                    Indel offset [3]
    -s, --NMThres <FLOAT>                       Allowable mismatch fraction compared to the read-length from a single read [0.05]
    -P, --NMFraction <FLOAT>                    Allowable fraction of the reads with excessive mismatches [0.3]
    -l, --BalleleCountWithIndel <INT>           Minimum number of B-allele with indels to determine proximity gap [3]
    -E, --mutFreeThres <FLOAT>                  Allowable fraction of mutation-free replicate [0.5]
    -F, --maxVAFdiff <FLOAT>                    Allowable VAF difference between replicates [0.05]
    -m, --mutationRate <FLOAT>                  Somatic mutation rate [3E-06]

Output options:

    -o, --output_directory <PATH>               Path for output directory
    -L, --label <STRING>                        Labeling the replicate set for output filename (label.snv.call)
    -t, --tempDir <PATH>                        Path for temporary directory
    -k, --keep_intermediate_files               Keep intermediate files

-------------------------------------------

<br>

Test running
=======

We have provided a small demo dataset for installation verification. You can find the data and do test run under testdata folder. Before running this, you should change the path for REF and RSCRIPT in run.sh file.

    cd testdata
    ./run.sh

You can also test the running command directly.

    java -jar ../RePlow.jar \
	-r Reference.fasta \
	-b rep1.bam,rep2.bam \
	-N normal.bam \
	-T test.bed \
	-R RSCRIPT_PATH \
	-o OUTPUT_DIR \
	-L LABEL_FOR_OUTPUT

You will get an output file (out.snv.call) under OUTPUT folder  when your running is successfully finished.


<br>

Resource requirements
=======

To analyze a 500X whole-exome sequencing dataset (~25 Gb), RePlow takes ~9 hours per sample for the entire process with a single core and ~6.5 Gb RAM usage. 
RePlow is not suitable for whole-genome datasets since it is designed to focus on detecting low-VAF mutations (~1% VAF), which requires deep enough depth for sequencing data, like 400X at least.

<br>

Inputs and options
=======
####Mandatory inputs:####
There are four **mandatory** inputs for RePlow.

Input                  | Option    | Description
---------------------- | --------- | -----------
**Reference sequence** | -r, --reference_sequence | FASTA formatted reference sequence file. *The reference must be BWA indexed.*
**List for replicate sequence data** | -b, --bam_file_list | Comma-separated list of BAM formatted alignment files for replicates. The BAM file must be (coordinate) sorted and indexed.
**Rscript executable file**   | -R, --Rscript | Rscript executable file path (absolute path).
**List for target BED files**   | -T, --target_bed_list | Comma-separated list of target BED files for each replicate. Just one BED file can also be given if target regions are identical for all replicates.

<br>

####Options:####
There are several options you can give to RePlow for more accurate running. Note that lenient filtration thresholds were applied as default settings compared with conventional methods, to apply our model for datasets regardless of the sequencing platform. Most false positives will be filtered by the calculated probability based on the estimated background error profiles.

**Analysis parameters for background error profiling:**

Option  | Default Value                       | Description
------- | ----------------------------------- | ------------------
-q, --bQ_thres | 15 | Minimum base call quality to consider a alternative base as a mismatch (not a sequencing error)
-Q, --mQ_thres | 20 |  Minimum mapping quality to consider a read for analysis
-N, --normal_BAM |  | BAM file of matched-control sample
-v, --germline_SNP |  | Germline SNP call (.vcf) to exclude genomic positions from the error profiling (e.g. SNP call from the GATK HaplotypeCaller)
-V, --somatic_SNV |  | Comma-separated somatic SNV call list to exclude repeatedly called positions from the error profiling (e.g. MuTect SNV calls for each replicate)
-f, --germlineBAFthres | 0.4 | VAF threshold for germline SNP filtering (for non-matched samples)

**Analysis parameters for variant detection:**

Option  | Default Value                       | Description
------- | ----------------------------------- | ------------------
-n, --minBalleleCnt | 1 | Required minimum B allele count to consider as a variant candidate
-p, --pval_for_normal | 0.01 | P-value cutoff to filter out non-somatic variants. If variant alleles are also observed at the corresponding site of matched control data and the probability that a background error will cause a larger VAF is less than this threshold, the candidate will be filtered out (germline SNP, contamination, mosaicism, etc).
-I, --minMisReadCnt | 5 |  Required minimum number of B-allele to apply fraction-based filters. 
-e, --min_MOS_fraction | 0.5 | Required minimum fraction of MOS (MOS / total number of B alleles) to filter out variant candidates with low base call quality (sequencing errors)
-d, --minDiseaseDepth | 50 | Required minimum depth to consider as a variant candidate from the disease sample
-D, --minNormalDepth | 10 | Required minimum depth at the correspoding site of the matched-control sample
-y, --homopolymer_length | 4 | Minimum length to consider as a homopolymeric region. If four or more identical reference bases are continuously connected from the corresponding site and the variant allele is the same as the first inconsistent base next to them, the candidate will be filtered out as a homopolymeric artifact.
-c, --clustered_pos_thres | 3 | Threshold for the clustered position from the end of read
-C, --clusteredFractionThres | 0.8 | Allowable B-allele fraction with clustered position (e.g. If more than 80% of B alleles are located within 3 bp from the end of read, the candidate will be filtered out)
-S, --successiveFractionThres | 0.8 | Allowable B-allele fraction with successive mutation (e.g. If more than 80% of B alleles are successively located beside another mismatch, the candidate will be filtered out)
-M, --lowMQFraction | 1.0 | Allowable B-allele fraction with low mapping quality
-B, --lowBQFraction | 0.5 | Allowable B-allele fraction with low base-call quality
-g, --indelFraction | 0.5 | Allowable fraction of the reads with indels at the corresponding position (e.g. If more than 50% of reads have indels at the corresponding site, the candidate will be filtered out)
-G, --BalleleFractionWithIndel | 0.4 | Allowable B-allele fraction with the indels within the indel offset (e.g. If more than 40% of the reads with B allele have indels within 3 bp from the corresponding site, the candidate will be filtered out)
-i, --indel_offset | 3 | Indel offset
-l, --BalleleCountWithIndel | 3 | Required minimum number of B-allele with indels to apply the --BalleleFractionWithIndel filter (-G)
-s, --NMThres | 0.05 | Allowable mismatch fraction compared to the read-length from a single read (e.g. If a given 100 bp-long read contains five of more mismatches, it is considered as a read with excessive mismatches.)
-P, --NMFraction | 0.3 | Allowable fraction of the reads with excessive mismatches
-E, --mutFreeThres | 0.5 | Allowable fraction of mutation-free replicate at the corresponding site
-F, --maxVAFdiff | 0.05 | Allowable VAF difference between replicates at the corresponding site
-m, --mutationRate | 3E-06 | Somatic mutation rate

**Output options:**

Option  | Default Value                       | Description
------- | ----------------------------------- | ------------------
-o, --output_directory | Current working directory | Path for output directory
-L, --label | | Labeling the replicate set for output filename (label.snv.call)


<br>

Output
=======
RePlow outputs the following files:

1. **ID.snv.call**
    : List of the variant candidates with the probability scores. All rejected candidates are also included with their filtered reasons. Final somatic candidates are tagged as 'PASS' in the Filter field. Labeled name (with -L) or the name of the first BAM file from the given input list will be used as an ID for the output files.
+  **intersection.bed (optional)**
    : If the target regions are different between given replicates, an intersection is calculated first and saved as a file (intersection.bed). Only the regions included in the intersection file is considered for the entire analysis.

<br>

Single running mode in RePlow
=======

Although RePlow is designed to simultaneously analyze multiple replicates, it can also be used for a single replicate. On-the-fly estimation of background errors provides great benefits to remove false positives even for the analysis of a single library (See the figure below). If you are trying to detect low-VAF mutations with deep-enough data and you don't have any replication for this, you can try our single running mode to eliminate false positive error calls. However, when the sequencing data contains a large number of background errors with VAFs similar to true mutations (e.g. ILA cases in the figure), the VAF cutoff becomes high enough to also remove true mutations due to the high estimated error rate, resulting in significant decrease of sensitivity. With replicates of those error-rich data, RePlow can recover the removed true mutations based on the VAF concordance between replicates therefore showed far higher sensitivity. In addition, RePlow with replicates provides more sophisticated filtration of false positives, considering the error profile concordance; false positives with exceptionally high VAFs (therefore beyond the cutoff) from the analysis of single library can be subsequently filtered out when multiple replicates are considered together, based on their concordance of VAF and error probability. So all modules including the ones work well for the single library have benefits of replication, thus we recommend to generate library replicates or just use single running mode to prioritize high confidence candidates for accurate detection of low-level mutations.

<p align="center">
<img src="https://sourceforge.net/p/replow/wiki/Home/attachment/Sup_fig14.png" width="600"  />
</p>
<br>
