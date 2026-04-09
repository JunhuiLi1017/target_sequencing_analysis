# Mutect2 variant calling rules
rule variants_mutect2:
	input:
		i1="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
	output:
		o1="{outpath}/03_variants/01_mutect2/01_raw/{sample}.{ref_version}.vcf.gz",
		o2="{outpath}/03_variants/01_mutect2/01_raw/{sample}.{ref_version}.vcf.gz.stats"
	log:
		log="{outpath}/03_variants/logs/{sample}.{ref_version}.mutect2.log"
	params:
		pon=config['pon'],
		af_only_gnomad=config['af_only_gnomad'],
		ref=config['reference'],
		sample="{sample}",
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["gatk_4.6.1.0"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		Mutect2 \
		-R {params.ref} \
		-I {input} \
		--pon {params.pon} \
		-tumor {params.sample} \
		--germline-resource {params.af_only_gnomad} \
		-O {output.o1} > {log.log} 2>&1
		"""

rule filter_mutectcalls:
	input:
		vcf="{outpath}/03_variants/01_mutect2/01_raw/{sample}.{ref_version}.vcf.gz",
		stat="{outpath}/03_variants/01_mutect2/01_raw/{sample}.{ref_version}.vcf.gz.stats"
	output:
		vcf="{outpath}/03_variants/01_mutect2/02_filter/{sample}.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/01_mutect2/02_filter/{sample}.{ref_version}.vcf.gz.tbi"
	log:
		log="{outpath}/03_variants/logs/{sample}.{ref_version}.FilterMutectCalls.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image["gatk_4.6.1.0"]
	shell:
		'''
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		FilterMutectCalls \
		-R {params.ref} \
		-V {input.vcf} \
		-O {output.vcf} > {log.log} 2>&1
		'''
