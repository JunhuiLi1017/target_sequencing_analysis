rule variants_pisces:
    input:
        bam="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam",
        bai="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam.bai"
    output:
        vcf=temp("{outpath}/03_variants/02_pisces/01_raw/{sample}.vcf")
    log:
        "{outpath}/03_variants/logs/{sample}.variants_pisces.log"
    params:
        outdir="{outpath}/03_variants/02_pisces/01_raw/",
        ref=config['genomefolders']
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    container:
        container_image["pisces_5.2.10"]
    shell:
        """
        #set +u
        #module load pisces/5.2.10
        pisces -g {params.ref} -b {input.bam} -CallMNVs false -gVCF false -o {params.outdir} > {log} 2>&1
        #set -u
        """

rule variants_pisces_index:
    input:
        vcf="{outpath}/03_variants/02_pisces/01_raw/{sample}.vcf"
    output:
        gz="{outpath}/03_variants/02_pisces/02_filter/{sample}.{ref_version}.vcf.gz",
        tbi="{outpath}/03_variants/02_pisces/02_filter/{sample}.{ref_version}.vcf.gz.tbi"
    log:
        "{outpath}/03_variants/logs/{sample}.{ref_version}.variants_pisces_index.log"
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    shell:
        """
        bgzip -c {input.vcf} > {output.gz} && tabix {output.gz}
        """
