rule run_deepvariant:
	message:
		"The BAM file must be also sorted and indexed. Duplicate marking may be performed, in our analyses there is almost no difference in accuracy except at lower (<20x) coverages. Finally, we recommend that you do not perform BQSR."
	input:
		bam="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam",
		bai="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam.bai"
	output:
		vcf="{outpath}/03_variants/04_deepvariant/01_raw/{sample}.{individual_chr}.{ref_version}.output.vcf.gz",
		gvcf="{outpath}/03_variants/04_deepvariant/01_raw/{sample}.{individual_chr}.{ref_version}.output.gvcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/01_raw/{sample}.{individual_chr}.{ref_version}.output.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{individual_chr}.{ref_version}.deepvariant.log"
	params:
		model_type=config["deepvariant_model"],
		ref=config['reference'],
		interval_list=intervals_dir + "/{individual_chr}.intervals.list",
		intermediate_dir="{outpath}/03_variants/04_deepvariant/01_raw/{sample}.{individual_chr}.{ref_version}.intermediate_results_dir",
		sample="{sample}",
		chr="{individual_chr}"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/mosaic_variants_calling_largeinput/workflow/envs/deepvariant_1.8.0.sif"
	shell:
		"""
		run_deepvariant \
		  --model_type={params.model_type} \
		  --ref={params.ref} \
		  --reads={input.bam} \
		  --output_vcf={output.vcf} \
		  --output_gvcf={output.gvcf} \
		  --num_shards={threads} \
		  --logging_dir={log} \
		  --intermediate_results_dir={params.intermediate_dir} \
		  --regions={params.chr}
		"""

rule merge_variant:
	input:
		vcf=lambda wildcards: [f"{wildcards.outpath}/03_variants/04_deepvariant/01_raw/{wildcards.sample}.{chr}.{wildcards.ref_version}.output.vcf.gz" for chr in CHROMOSOMES],
		tbi=lambda wildcards: [f"{wildcards.outpath}/03_variants/04_deepvariant/01_raw/{wildcards.sample}.{chr}.{wildcards.ref_version}.output.vcf.gz.tbi" for chr in CHROMOSOMES]
	output:
		vcf="{outpath}/03_variants/04_deepvariant/02_filter/{sample}.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/04_deepvariant/02_filter/{sample}.{ref_version}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{sample}.{ref_version}.04_deepvariant.merge.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/envs/bcftools_v1.10.2.sif"
	shell:
		"""
		bcftools concat -a {input.vcf} | bcftools sort | bgzip > {output.vcf}
		tabix -p vcf {output.vcf}
		"""