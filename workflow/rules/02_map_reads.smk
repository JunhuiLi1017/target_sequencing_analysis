rule map_reads:
	input:
		get_fastqc_input
	output:
		raw_bam=temp("{outpath}/02_map/01_raw/{sample}/{sample}.{library}.{flowlane}.raw.bam")
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.{library}.{flowlane}.map_reads.log"
	params:
		rg=get_reads_group,
		ref=config['reference'],
		fq2_expr=lambda wildcards, input: f" {input[1]}" if is_pair_end else ""
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	#container:
	#	container_image["bwa_samtools_0.7.17"]
	#container:
	#	"http://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_6"
	shell:
		"""
		bwa mem -t {threads} \
		-K 500000000 \
		-M {params.rg} \
		{params.ref} \
		{input[0]} \
		{params.fq2_expr} | samtools view \
		-@ {threads} \
		 -b -o {output.raw_bam} > {log} 2>&1
		"""
		
## merge multiple bams for each single sample
rule merge_and_sort:
	input:
		lambda wildcards: [
			f"{wildcards.outpath}/02_map/01_raw/{wildcards.sample}/{wildcards.sample}.{u.library}.{u.flowlane}.raw.bam"
			for u in units[units["sample"] == wildcards.sample].itertuples()
		]
	output:
		bam=temp("{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam"),
		bai=temp("{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam.bai")
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.merge_and_sort.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	params:
		tmpdir="{outpath}/02_map/02_sort/{sample}/tmpdir_{sample}",
		# Conditionally handle single vs. multiple BAM files
		merged_bam=lambda wildcards, input: f"{outpath}/02_map/02_sort/{wildcards.sample}/{wildcards.sample}.merged.bam"
		if len(input) > 1 else input[0],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	container:
		container_image["sambamba_1.0.1"]
	#	"http://depot.galaxyproject.org/singularity/sambamba:1.0.1--he614052_4"
	shell:
		"""
		if [ {input} != {params.merged_bam} ]; then
			sambamba merge -t {threads} {params.merged_bam} {input} > {log} 2>&1
		fi

		mkdir -p {params.tmpdir}

		sambamba sort \
		-t {threads} \
		-m {params.command_mem}M \
		-o {output.bam} \
		--tmpdir {params.tmpdir} \
		{params.merged_bam} > {log} 2>&1
		sambamba flagstat {output.bam}
		sambamba index {output.bam}
		""" 

rule remove_dup_umi:
	input:
		bam="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam",
		bai="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam.bai"
	output:
		bam=temp("{outpath}/02_map/03_rmdup/{sample}/{sample}.03_rmdup.umi.bam")
	log:
		"{outpath}/02_map/02_map/logs/{sample}/{sample}.remove_dup.log"
	params:
		umi_separator=":",
		output_stats="{outpath}/02_map/03_rmdup/{sample}/{sample}.sort.rmdup"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["umi_tools_1.1.5"]
		#"umi_tools:1.1.6--py39hbcbf7aa_0"
	shell:
		"""
		#source ~/anaconda3/etc/profile.d/conda.sh; conda activate umi_tools116
		umi_tools dedup --stdin={input.bam} \
			--log={log} \
			--output-stats={params.output_stats} \
			--stdout={output.bam} \
			--umi-separator={params.umi_separator} \
			--paired >> {log} 2>&1
		#conda deactivate
		"""

rule remove_dup_gatk:
	input:
		bam="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam",
		bai="{outpath}/02_map/02_sort/{sample}/{sample}.02_sort.bam.bai"
	output:
		bam=temp("{outpath}/02_map/03_rmdup/{sample}/{sample}.03_rmdup.gatk.bam"),
		mtx="{outpath}/02_map/03_rmdup/{sample}/{sample}.sort.rmdup.matrix"
	log:
		"{outpath}/02_map/02_map/logs/{sample}/{sample}.remove_dup.log"
	params:
		dedup=lambda wildcards: (
			f"--REMOVE_DUPLICATES true" 
			if units[units['sample'] == wildcards.sample]['pcr_based'].any() 
			else "--REMOVE_DUPLICATES false --REMOVE_SEQUENCING_DUPLICATES true --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500"
		),
		tmpdir="{outpath}/02_map/03_rmdup/{sample}/tmpdir_{sample}",
		input_args=lambda wildcards, input: " ".join(f"--INPUT {bam}" for bam in input.bam),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 8000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
		#mem_mb=5000
	container:
		container_image["gatk_4.6.1.0"]
	shell:
		"""
		mkdir -p {params.tmpdir}
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		MarkDuplicates --INPUT {input.bam} \
		--METRICS_FILE {output.mtx} \
		--OUTPUT {output.bam} \
		{params.dedup} \
		--TMP_DIR {params.tmpdir} \
		--CREATE_INDEX true > {log} 2>&1
		"""

rule base_recalibrator:
	input:
		bam=get_dedup_files
	output:
		o1="{outpath}/02_map/04_bqsr/{sample}/{sample}.recal_data.table"
	log:
		"{outpath}/02_map/logs/{sample}/{sample}.base_recalibrator.log"
	params:
		gatk=config['gatk_current_using'],
		ref=config['reference'],
		dbsnp138=config['dbsnp138'],
		g1000_known_indels=config['g1000_known_indels'],
		mills_and_1000g=config['mills_and_1000g'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["gatk_4.6.1.0"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		BaseRecalibrator \
		-I {input.bam} \
		-O {output.o1} \
		-R {params.ref} \
		--known-sites {params.dbsnp138} \
		--known-sites {params.g1000_known_indels} \
		--known-sites {params.mills_and_1000g} > {log} 2>&1
		"""

rule apply_bqsr:
	input:
		bam=get_dedup_files,
		recal_table="{outpath}/02_map/04_bqsr/{sample}/{sample}.recal_data.table"
	output:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
	log:
		"{outpath}/02_map/logs/{sample}.apply_bqsr.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["gatk_4.6.1.0"]
	shell:
		"""
		gatk --java-options "-Xms{params.command_mem}m -XX:ParallelGCThreads={threads}" \
		ApplyBQSR \
		-I {input.bam} \
		-O {output.bam} \
		-R {params.ref} \
		--bqsr-recal-file {input.recal_table} > {log} 2>&1
		"""

rule apply_bqsr_index:
	input:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
	output:
		bai="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam.bai",
		stats="{outpath}/02_map/07_summary/stats/{sample}.samtools.stats.txt",
		idxstats="{outpath}/02_map/07_summary/idxstats/{sample}.samtools.idxstats.txt",
		flagstat="{outpath}/02_map/07_summary/flagstat/{sample}.samtools.flagstat.txt"
	log:
		"{outpath}/02_map/logs/{sample}.apply_bqsr.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 4000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		container_image["samtools_1.20"]
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools idxstats {input.bam} > {output.idxstats}
		samtools flagstat {input.bam} > {output.flagstat}
		samtools index {input.bam}
		"""

rule multiqc_stats_bqsr:
	input:
		stats=expand(
			["{outpath}/02_map/07_summary/stats/{u.sample}.samtools.stats.txt"],
			outpath=Outpath,
			u=units.itertuples()
		 )			
	output:
		stats_report="{outpath}/02_map/07_summary/report_stats/multiqc_report.html"
	log:
		"{outpath}/02_map/logs/multiqc_stats.log"
	params:
		out_multiqc_stats="{outpath}/02_map/07_summary/report_stats",
		in_stats="{outpath}/02_map/07_summary/stats"
	threads: 
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image['multiqc_1.22.3']
	shell:
		"""
		multiqc -o {params.out_multiqc_stats} {params.in_stats} --force
		"""

rule multiqc_idxstats_bqsr:
	input:
		idxstats=expand(
			["{outpath}/02_map/07_summary/idxstats/{u.sample}.samtools.idxstats.txt"],
			outpath=Outpath,
			u=units.itertuples()
		)			
	output:
		idxstats_report="{outpath}/02_map/07_summary/report_idxstats/multiqc_report.html"
	log:
		"{outpath}/02_map/logs/multiqc_idxstats.log"
	params:
		out_multiqc_idxstats="{outpath}/02_map/07_summary/report_idxstats",
		in_idxstats="{outpath}/02_map/07_summary/idxstats"
	threads: 
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image['multiqc_1.22.3']
	shell:
		"""
		multiqc -o {params.out_multiqc_idxstats} {params.in_idxstats} --force
		"""


rule multiqc_flagstat_bqsr:
	input:
		flagstat=expand(
			["{outpath}/02_map/07_summary/flagstat/{u.sample}.samtools.flagstat.txt"],
			outpath=Outpath,
			u=units.itertuples()
		)			
	output:
		flagstat_report="{outpath}/02_map/07_summary/report_flagstat/multiqc_report.html"
	log:
		"{outpath}/02_map/logs/multiqc_flagstat.log"
	params:
		out_multiqc_flagstat="{outpath}/02_map/07_summary/report_flagstat",
		in_flagstat="{outpath}/02_map/07_summary/flagstat"
	threads: 
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image['multiqc_1.22.3']
	shell:
		"""
		multiqc -o {params.out_multiqc_flagstat} {params.in_flagstat} --force
		"""