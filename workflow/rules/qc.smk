rule trim:
    input:
        unpack(get_fastq)
    output:
        ["result/01_multiqc/fastp_trim/{sample}.R1.fastq.gz", "result/01_multiqc/fastp_trim/{sample}.R2.fastq.gz"] if paired_end else "result/01_multiqc/fastp_trim/{sample}.R1.fastq.gz"
    log:
        "logs/fastqc/{sample}.log"
    run:
        if paired_end:
            shell("fastp -f 10 -F 10 -i {input.r1} -I {input.r2} -o {output[0]} -O {output[1]} > {log} 2>&1")
        else:
            shell("fastp -f 10 -F 10 -i {input.r1} -o {output[0]} > {log} 2>&1")

rule fastqc:
    input:
        get_trim_fastq
    output:
        ["result/01_multiqc/fastqc/{sample}.R1_fastqc.html", "result/01_multiqc/fastqc/{sample}.R2_fastqc.html"] if paired_end else "result/01_multiqc/fastqc/{sample}.R1_fastqc.html",
        ["result/01_multiqc/fastqc/{sample}.R1_fastqc.zip", "result/01_multiqc/fastqc/{sample}.R2_fastqc.zip"] if paired_end else "result/01_multiqc/fastqc/{sample}.R1_fastqc.zip"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc -o result/01_multiqc/fastqc {input}"

rule samtools_stats:
    input:
        #get_bqsr_bam
        "result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam"
    output:
        "result/02_Map/bqsr_stat/{sample}.sort.rmdup.bqsr.stat"
    params:
        TargeRegion=config['TargeRegion']
    log:
         "logs/bwa/{sample}.bqsr.stat.log"
    shell:
        "samtools stats {input} -t {params.TargeRegion} > {output}"

rule multiqc:
    input:
        get_multiqc_input
    output:
        "result/01_multiqc/multiqc_report.html"
    conda:
        "../envs/multiqc_env.yaml"
    log:
        "logs/multiqc/multiqc.log"
    shell:
        """
        multiqc -o result/01_multiqc result/01_multiqc/fastqc result/02_Map/bqsr_stat/
        """