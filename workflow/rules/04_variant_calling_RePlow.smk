rule variants_RePlow:
    input:
        i1="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
    output:
        o1="{outpath}/03_variants/03_replow/01_raw/{sample}.{ref_version}.snv.call"
    log:
        "{outpath}/03_variants/logs/{sample}.{ref_version}.variants_RePlow.log"
    params:
        outdir="{outpath}/03_variants/03_replow/01_raw",
        TargeRegion=config['TargeRegion'],
        ref=config['reference'],
        label="{sample}.{ref_version}",
        RePlow=config['RePlow'],
        command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    shell:
        """
        java -Xms{params.command_mem}m \
        -XX:ParallelGCThreads={threads} \
        -jar {params.RePlow} \
        -r {params.ref} \
        -b {input.i1} \
        -T {params.TargeRegion} \
        -R ~/anaconda3/envs/snakemake/bin/Rscript \
        --output_directory {params.outdir} \
        --label {params.label} > {log} 2>&1
        """

rule convert_tsv_to_vcf:
    input:
        "{outpath}/03_variants/03_replow/01_raw/{sample}.{ref_version}.snv.call"
    output:
        vcf=temp("{outpath}/03_variants/03_replow/02_filter/{sample}.{ref_version}.vcf"),
        gz="{outpath}/03_variants/03_replow/02_filter/{sample}.{ref_version}.vcf.gz",
        tbi="{outpath}/03_variants/03_replow/02_filter/{sample}.{ref_version}.vcf.gz.tbi"
    log:
        "{outpath}/03_variants/logs/{sample}.{ref_version}.convert_tsv_to_vcf.log"  
    params:
        convert_tsv_to_vcf="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/convert_tsv_to_vcf.py"    
    threads:
        resource['resource']['medium']['threads']
    resources:
        mem_mb=resource['resource']['medium']['mem_mb']
    shell:
        """
        python {params.convert_tsv_to_vcf} {input} {output.vcf} > {log} 2>&1
        bgzip -c {output.vcf} > {output.gz} && tabix {output.gz}
        """