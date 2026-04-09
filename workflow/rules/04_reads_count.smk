rule samtools_mpileup:
	input:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam",
		bai="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam.bai"
	output:
		mpileup="{outpath}/02_map/06_reads_count/01_samtools_mpileup/{sample}.{ref_version}.txt"
	log:
		log="{outpath}/02_map/logs/{sample}.{ref_version}.06_reads_count.log"
	params:
		ref=config['reference'],
		bed=config['MosaicRegion'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["samtools_1.20"]
	shell:
		"""
		samtools mpileup -l {params.bed} -f {params.ref} {input.bam} > {output.mpileup} 2> {log}
		"""

rule reads_count:
	input:
		mpileup="{outpath}/02_map/06_reads_count/01_samtools_mpileup/{sample}.{ref_version}.txt"
	output:
		"{outpath}/02_map/06_reads_count/02_reads_count/{sample}.{ref_version}.reads_count.txt"
	log:
		"{outpath}/02_map/logs/{sample}.{ref_version}.reads_count.log"
	params:
		bed=config['MosaicRegion'],
		reads_count="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/reads_count_v1.3.py"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["terra_py_tools"]
	shell:
		"""
		python {params.reads_count} --input {input.mpileup} --bed {params.bed} --output {output}
		"""