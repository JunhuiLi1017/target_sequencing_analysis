# haplotypecaller variant calling rules
rule variants_haplotypecaller:
	input:
		i1="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
	output:
		o1="{outpath}/03_variants_germline/04_haplotypecaller/01_raw/{sample}.{ref_version}.g.vcf.gz",
		o2="{outpath}/03_variants_germline/04_haplotypecaller/01_raw/{sample}.{ref_version}.g.vcf.gz.tbi"
	log:
		log="{outpath}/03_variants_germline/logs/{sample}.{ref_version}.haplotypecaller.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk4.6.1.0.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		HaplotypeCaller \
		-R {params.ref} \
		-I {input} \
		-ERC GVCF \
		-O {output.o1} > {log.log} 2>&1
		conda deactivate
		"""

rule sample_name_map:
	input:
		gvcfs=expand("{outpath}/03_variants_germline/04_haplotypecaller/01_raw/{sample}.{ref_version}.g.vcf.gz", outpath=Outpath, sample=Sample, ref_version=Ref_version),
		idx=expand("{outpath}/03_variants_germline/04_haplotypecaller/01_raw/{sample}.{ref_version}.g.vcf.gz.tbi", outpath=Outpath, sample=Sample, ref_version=Ref_version)
	output:
		"{outpath}/03_variants_germline/04_haplotypecaller/01_raw/sample.list.txt"
	log:
		"{outpath}/03_variants_germline/logs/genomics_db_import/sample.list.log"
	params:
		sample_list=lambda wildcards, input: "\n".join([f"{sample}\t{gvcf}" for sample, gvcf in zip(Sample, input.gvcfs)])
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		echo "{params.sample_list}" > {output}
		"""

rule genomics_db_import:
	input:
		"{outpath}/03_variants_germline/04_haplotypecaller/01_raw/sample.list.txt"
	output:
		db=directory("{outpath}/03_variants_germline/04_haplotypecaller/02_genomicsdb/{chrom}/{chrom}_gdb")
	log:
		"{outpath}/03_variants_germline/logs/{chrom}.genomics_db_import.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		interval_list=intervals_dir + "/{chrom}.intervals.list",
		tmp_dir="{outpath}/03_variants_germline/04_haplotypecaller/02_genomicsdb/temp_dir_{chrom}",
		gvcf_list=lambda wildcards, input: " ".join([f"-V {gvcf}" for gvcf in input]),
		intervals=lambda wildcards: f"{wildcards.chrom}",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	shell:
		"""
		mkdir -p {params.tmp_dir}
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		GenomicsDBImport \
		--sample-name-map {input} \
		--genomicsdb-workspace-path {output.db} \
		--tmp-dir {params.tmp_dir} \
		--batch-size 90 \
		-R {params.ref} \
		--max-num-intervals-to-import-in-parallel 25 \
		--intervals {params.intervals} > {log} 2>&1
		conda deactivate
		"""

# Rule to perform joint genotyping on GenomicsDB workspace
rule genotype_gvcfs:
	input:
		db="{outpath}/03_variants_germline/04_haplotypecaller/02_genomicsdb/{chrom}/{chrom}_gdb"
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/03_joint/{chrom}.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/03_joint/{chrom}.vcf.gz.tbi"
	log:
		"{outpath}/03_variants_germline/logs/{chrom}.genotype_gvcfs.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['very_high']['threads']
	resources:
		mem_mb=resource['resource']['very_high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		GenotypeGVCFs \
		-R {params.ref} \
		-V gendb://{input.db} \
		-O {output.vcf} \
		--include-non-variant-sites false > {log} 2>&1
		conda deactivate
		"""

# Rule to merge per-chromosome VCFs into a single cohort VCF
rule merge_vcfs:
	input:
		vcfs=expand("{outpath}/03_variants_germline/04_haplotypecaller/03_joint/{chrom}.vcf.gz", outpath=Outpath, chrom=CHROMOSOMES),
		idx=expand("{outpath}/03_variants_germline/04_haplotypecaller/03_joint/{chrom}.vcf.gz.tbi", outpath=Outpath, chrom=CHROMOSOMES)
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz.tbi"
	log:
		"{outpath}/03_variants_germline/logs/merge_vcfs.log"
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	params:
		gatk=config['gatk_current_using'],
		vcf_list=lambda wildcards, input: " ".join([f"-I {vcf}" for vcf in input.vcfs]),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		MergeVcfs \
		{params.vcf_list} \
		-O {output.vcf} > {log} 2>&1
		conda deactivate
		"""

# Rule to perform Variant Quality Score Recalibration (VQSR) for SNPs
rule variant_recalibrator_snp:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz.tbi"
	output:
		recal="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.snp.recal",
		tranches="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.snp.tranches",
		rscript="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.snp.plots.R"
	log:
		"{outpath}/03_variants_germline/04_haplotypecaller/logs/variant_recalibrator_snp.log"
	params:
		ref=config['reference'],
		hapmap=config["hapmap"],
		omni=config["omni"],
		g1000_known_indels=config["g1000_known_indels"],
		dbsnp=config["dbsnp138"],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		VariantRecalibrator \
			-R {params.ref} \
			-V {input.vcf} \
			--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} \
			--resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} \
			--resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.g1000_known_indels} \
			--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
			-mode SNP \
			-O {output.recal} \
			--tranches-file {output.tranches} \
			--rscript-file {output.rscript} \
			--dont-run-rscript \
			--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 > {log} 2>&1
			conda deactivate
		"""

# Rule to perform VQSR for INDELs
rule variant_recalibrator_indel:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz.tbi"
	output:
		recal="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.indel.recal",
		tranches="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.indel.tranches",
		rscript="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.indel.plots.R"
	log:
		"{outpath}/03_variants_germline/logs/variant_recalibrator_indel.log"
	params:
		ref=config['reference'],
		mills=config["mills_and_1000g"],
		dbsnp=config["dbsnp138"],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		VariantRecalibrator \
			-R {params.ref} \
			-V {input.vcf} \
			--resource:mills,known=false,training=true,truth=true,prior=12.0 {params.mills} \
			--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} \
			-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
			-mode INDEL \
			-O {output.recal} \
			--tranches-file {output.tranches} \
			--rscript-file {output.rscript} \
			--dont-run-rscript \
			--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
			--max-gaussians 4 > {log} 2>&1
			conda deactivate
		"""

# Rule to apply VQSR for SNPs
rule apply_vqsr_snp:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/04_merge_vcf/all.vcf.gz.tbi",
		recal="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.snp.recal",
		tranches="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.snp.tranches",
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.snp.recalibrated.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.snp.recalibrated.vcf.gz.tbi"
	log:
		"{outpath}/03_variants_germline/logs/apply_vqsr_snp/all.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		ApplyVQSR \
			-R {params.ref} \
			-V {input.vcf} \
			-O {output.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranches} \
			-mode SNP \
			--truth-sensitivity-filter-level 99.0 > {log} 2>&1
		conda deactivate
		"""

# Rule to apply VQSR for INDELs
rule apply_vqsr_indel:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.snp.recalibrated.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.snp.recalibrated.vcf.gz.tbi",
		recal="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.indel.recal",
		tranches="{outpath}/03_variants_germline/04_haplotypecaller/05_vqsr/all.indel.tranches"
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.recalibrated.vcf.gz",
		idx="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.recalibrated.vcf.gz.tbi"
	log:
		"{outpath}/03_variants_germline/logs/apply_vqsr_indel.log"
	params:
		ref=config['reference'],
		gatk=config['gatk_current_using'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/gatk.sif"
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate gatk4.6.1.0
		java -Xms{params.command_mem}m -XX:ParallelGCThreads={threads} \
		-jar {params.gatk} \
		ApplyVQSR \
			-R {params.ref} \
			-V {input.vcf} \
			-O {output.vcf} \
			--recal-file {input.recal} \
			--tranches-file {input.tranches} \
			-mode INDEL \
			--truth-sensitivity-filter-level 99.0 > {log} 2>&1
		conda deactivate
		"""

rule split_multiallelic_germline:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.recalibrated.vcf.gz",
		tbi="{outpath}/03_variants_germline/04_haplotypecaller/06_apply_vqsr/all.recalibrated.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/07_split_mulalle/all.recalibrated.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants_germline/04_haplotypecaller/07_split_mulalle/all.recalibrated.split_mulalle.vcf.gz.tbi"
	log:
		"{outpath}/03_variants_germline/logs/all.split_multiallelic.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	singularity:
		"../envs/bcftools.sif"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools norm -m - -o {output.vcf} -Oz {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		conda deactivate
		'''

rule pass_filter_germline:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/07_split_mulalle/all.recalibrated.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants_germline/04_haplotypecaller/07_split_mulalle/all.recalibrated.split_mulalle.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/08_pass/all.recalibrated.pass.vcf.gz",
		tbi="{outpath}/03_variants_germline/04_haplotypecaller/08_pass/all.recalibrated.pass.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/all.04_haplotypecaller.pass.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	singularity:
		"../envs/bcftools.sif"
	shell:
		'''
		source ~/anaconda3/etc/profile.d/conda.sh; conda activate bcftools
		bcftools filter -i 'FILTER="PASS"' -o {output.vcf} -Oz {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		conda deactivate
		'''

rule annotate_clinvar_gnomad_germline:
	input:
		vcf="{outpath}/03_variants_germline/04_haplotypecaller/08_pass/all.recalibrated.pass.vcf.gz",
	output:
		txt="{outpath}/03_variants_germline/04_haplotypecaller/09_annovar/all.pass.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants_germline/all.{ref_version}.clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants_germline/04_haplotypecaller/09_annovar/all.pass",
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	singularity:
		"../envs/perl.sif"
	shell:
		'''
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{params.ref_version} \
		-buildver {params.ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol refGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		'''

rule annotate_rcnv_gnomadlof_germline:
	input:
		tier_anno="{outpath}/03_variants_germline/04_haplotypecaller/09_annovar/all.pass.{ref_version}_multianno.txt"
	output:
		sub="{outpath}/03_variants_germline/04_haplotypecaller/10_score/all.pass.{ref_version}.exonic_splicing_multianno.txt",
		txt="{outpath}/03_variants_germline/04_haplotypecaller/10_score/all.pass.{ref_version}.rcnv_gnomadlof_multianno.txt"
	log:
		"{outpath}/03_variants_germline/logs/all.{ref_version}.rcnv_lof.log"
	params:
		ref_version=config['ref_version'],
		gnomad_LoF=config['gnomad_LoF'],
		rCNV_gene_score=config['rCNV_gene_score'],
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		'''
		cat <(awk '{{if($7=="exonic"){{print $0}}}}' {input.tier_anno} | grep -E 'nonsynonymous|stop') <(awk '{{if($7=="splicing"){{print $0}}}}' {input.tier_anno}) > {output.sub}
		cat <(paste <(head -n 1 {input.tier_anno}) <(head -n 1 {params.gnomad_LoF}) <(head -n 1 {params.rCNV_gene_score})) <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.gnomad_LoF} <(awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$8]){{print $0"\t"c[$8]}}else{{print $0"\tNA\tNA\tNA"}}}}' {params.rCNV_gene_score} {output.sub})) > {output.txt}
		'''