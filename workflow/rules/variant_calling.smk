rule variants_Mutect2:
    input:
        i1="result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam"
    output:
        o1="result/03_Variants/Mutect2/{sample}.mutect2.vcf.gz"
    params:
        PON=config['PON'],
        af_only_gnomad=config['af_only_gnomad'],
        REFERENCE=config['ref'],
        sample="{sample}"
    log:
        "logs/variants/{sample}.mutect2.log"
    threads:
        8
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar Mutect2  -R {params.REFERENCE} -I {input.i1} --pon {params.PON} --tumor {params.sample} --germline-resource {params.af_only_gnomad} --interval-padding 100 -O {output.o1} > {log} 2>&1        
        """

rule FilterMutectCall:
    input:
        vcf="result/03_Variants/Mutect2/{sample}.mutect2.vcf.gz"
    output:
        "result/03_Variants/Mutect2/{sample}.mutect2.filter.vcf.gz"
    log:
        "logs/variants/{sample}.FilterMutectCalls.log"
    threads:
        8
    params:
        REFERENCE=config['ref']
    shell:
        '''
        java -Xms10g -XX:ParallelGCThreads={threads} -jar /home/junhui.li11-umw/anaconda3/envs/snakemake/share/gatk4-4.1.8.1-0/gatk-package-4.1.8.1-local.jar FilterMutectCalls -R {params.REFERENCE} -V {input.vcf} -O {output} > {log} 2>&1
        '''

rule correct_vcf:
    input:
        vcf="result/03_Variants/Mutect2/{sample}.mutect2.filter.vcf.gz"
    output:
        o1="result/03_Variants/Mutect2_correct/{sample}.mutect2.filter.vcf.gz"
    shell:
        """
        zcat {input.vcf} | sed 's/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=1/g' | bgzip > {output.o1} && tabix {output.o1}
        """

rule merge_all_vcf:
    input:
        expand(["result/03_Variants/Mutect2_correct/{u.sample}.mutect2.filter.vcf.gz"], u=units.itertuples())
    output:
        o1="result/03_Variants/01_merge_vcf/all.mutect2.filter.vcf.gz"
    params:
        meta_sample=config['meta_sample']
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        conda activate bcftools
        bcftools merge {input} -o {output.o1}
        conda deactivate
        """

rule variants_pysam:
    input:
        i1="result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam"
    output:
        o1="result/03_Variants/02_pysam_depth/{sample}.allel.depth.txt" 
    params:
        REFERENCE=config['ref'],
        variants_bed=config['FinalRegion']
    shell:
        """
        python /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/allel_coverage_v3.py {input.i1} {params.variants_bed} {output.o1}
        """

rule merge_all_depth:
    input:
        expand(["result/03_Variants/02_pysam_depth/{u.sample}.allel.depth.txt"], u=units.itertuples())
    output:
        o1="result/03_Variants/02_pysam_merge_depth/all.pysam.depth.txt"
    shell:
        """
        cat {input} > {output.o1}
        """

'''
rule variants_RePlow:
    input:
        i1="result/02_Map/bqsr/{sample}.sort.rmdup.bqsr.bam"
    output:
        o1="result/03_Variants/RePlow/{sample}.RePlow.vcf.gz"
    params:
        outdir="result/03_Variants/RePlow"
        af_only_gnomad=config['af_only_gnomad'],
        REFERENCE=config['ref'],
    log:
        "logs/variants/{sample}.RePlow.log"
    threads:
        8
    shell:
        """
        java -Xms10g -XX:ParallelGCThreads={threads} -jar ../../software/RePlow.jar  -r {params.REFERENCE} -b <replicate_sorted_indexed_BAM_list> -T <replicate_target_BED_file_list> -R ~/anaconda3/envs/snakemake/bin/Rscript --output_directory {param.outdir} --label  > {log} 2>&1        
        """
'''