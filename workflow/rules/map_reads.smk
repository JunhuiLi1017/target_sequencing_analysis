rule map_reads:
    input:
        get_trim_fastq
    output:
        o1=temp("result/02_Map/bwa/{sample}.raw.bam")
    log:
        l1="logs/bwa/{sample}.bwa.log"
    threads:
        8
    params:
        mem="8G",
        rg=get_reads_group,
        ref=config['ref']
    run:
        if paired_end:
            shell("bwa mem -t {threads} -M {params.rg} {params.ref} {input[0]} {input[1]} | samtools view -b -o {output.o1} > {log.l1} 2>&1")
        else:
            shell("bwa mem -t {threads} -M {params.rg} {params.ref} {input[0]} | samtools view -b -o {output.o1} > {log.l1} 2>&1")

rule sort_bam:
    input:
        "result/02_Map/bwa/{sample}.raw.bam"
    output:
        o2="result/02_Map/bwa/{sample}.sort.bam",
        o3="result/02_Map/bwa/{sample}.sort.stat"
    log:
        l2="logs/bwa/{sample}.bwa.sort.log"
    threads:
        12
    params:
        mem="4000",
        tmpdir="02_Map/bwa",
        TargeRegion=config['TargeRegion']
    shell:
        """
        samtools sort -@ {threads} -m 8G -O bam -o {output.o2} {input}> {log.l2} 2>&1
        samtools stats {output.o2} -t {params.TargeRegion} > {output.o3}
        samtools index {output.o2}
        """  

rule remove_dup:
    input:
        #"result/02_Map/bwa/{sample}.rename.rx.sort.bam"
        "result/02_Map/bwa/{sample}.sort.bam"
    output:
        o1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    log:
        l1="logs/bwa/{sample}.removeduplicate.log",
        l2="logs/bwa/{sample}.removeduplicate.process.log"
    threads:
        8
    params:
        mem="4000",
        prefix="result/02_Map/dup/{sample}.dedup.umi",
        outdir="result/02_Map/dup"
        #/share/pkg/jdk/1.8.0_77/bin/java -XX:ParallelGCThreads=12 -Xmx14G -jar /share/pkg/picard/2.23.3/picard.jar UmiAwareMarkDuplicatesWithMateCigar --INPUT {input} --ASSUME_SORT_ORDER  coordinate --METRICS_FILE {output.o2} --OUTPUT {output.o1} --UMI_METRICS_FILE {output.o3} --REMOVE_DUPLICATES true > {log} 2>&1
    shell:
        """
        umi_tools dedup --stdin={input} --log={log.l1} --output-stats={params.prefix} --stdout={output.o1} --umi-separator=":" --paired > {log.l2} 2>&1
        samtools index {output.o1}
        """

rule BaseRecalibrator:
    input:
        "result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="result/02_Map/bqsr/{sample}.recal_data.table"
    log:
        "logs/bwa/{sample}.BQSR.log"
    threads:
        8
    params:
        mem="4000",
        ref=config['ref'],
        dpsnp138=config['dpsnp138'],
        known_indels=config['known_indels'],
        Mills_and_1000G=config['Mills_and_1000G'],
        outdir="result/02_Map/bqsr"
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar BaseRecalibrator -I {input} -O {output.o1} -R {params.ref} --known-sites {params.dpsnp138} --known-sites {params.known_indels} --known-sites {params.Mills_and_1000G} > {log} 2>&1
        """

rule ApplyBQSR:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam",
        i2="result/02_Map/bqsr/{sample}.recal_data.table"
    output:
        o1="result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam",
        o2="result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam.bai"
    log:
        "logs/bwa/{sample}.applyBQSR.log"
    threads:
        8
    params:
        mem="4000",
        ref=config['ref'],
        prefix="result/02_Map/bqsr/{sample}"
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar ApplyBQSR -I {input.i1} -O {output.o1} -R {params.ref} --bqsr-recal-file {input.i2} > {log} 2>&1
        samtools flagstat {output.o1} > {params.prefix}.txt
        samtools index {output.o1} {output.o2}
        """