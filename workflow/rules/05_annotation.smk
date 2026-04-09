rule split_multiallelic:
	input:
		vcf="{outpath}/03_variants/{caller}/02_filter/{sample}.{ref_version}.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/02_filter/{sample}.{ref_version}.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/{caller}/03_split_mulalle/{sample}.{ref_version}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/03_split_mulalle/{sample}.{ref_version}.split_mulalle.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{caller}/{sample}.{ref_version}.split_multiallelic.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/envs/bcftools_v1.10.2.sif"
	shell:
		'''
		bcftools norm -m - -o {output.vcf} -Oz {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		'''

rule pass_filter:
	input:
		vcf="{outpath}/03_variants/{caller}/03_split_mulalle/{sample}.{ref_version}.split_mulalle.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/03_split_mulalle/{sample}.{ref_version}.split_mulalle.vcf.gz.tbi"
	output:
		vcf="{outpath}/03_variants/{caller}/04_pass/{sample}.{ref_version}.pass.vcf.gz",
		tbi="{outpath}/03_variants/{caller}/04_pass/{sample}.{ref_version}.pass.vcf.gz.tbi"
	log:
		"{outpath}/03_variants/logs/{caller}/{sample}.{ref_version}.pass.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/envs/bcftools_v1.10.2.sif"
	shell:
		'''
		bcftools filter -i 'FILTER="PASS"' -o {output.vcf} -Oz {input.vcf} > {log} 2>&1
		tabix -p vcf {output.vcf}
		'''

rule annotate_clinvar_gnomad:
	input:
		vcf="{outpath}/03_variants/{caller}/04_pass/{sample}.{ref_version}.pass.vcf.gz"
	output:
		txt="{outpath}/03_variants/{caller}/05_annovar/{sample}.{ref_version}_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{caller}/{sample}.{ref_version}.clinvar.log"
	params:
		ref_version=config['ref_version'],
		annovar_dir=config['annovar_dir'],
		outputanno="{outpath}/03_variants/{caller}/05_annovar/{sample}",
		protocol=lambda wildcards: (
			"refGene,dbnsfp42a,clinvar_20240917,gnomad211_genome,gnomad211_exome"
			if wildcards.ref_version == "hg19"
			else "refGene,dbnsfp42a,clinvar_20240917,gnomad41_genome,gnomad41_exome"
		),
		command_mem=lambda wildcards, resources, threads: (resources.mem_mb * threads - 2000)
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	shell:
		'''
		perl {params.annovar_dir}/table_annovar.pl \
		{input.vcf} \
		{params.annovar_dir}/humandb_{params.ref_version} \
		-buildver {params.ref_version} \
		-out {params.outputanno} \
		-remove \
		-protocol {params.protocol} \
		-operation g,f,f,f,f \
		-nastring . \
		-vcfinput
		'''

rule annotate_rcnv_gnomadlof:
	input:
		tier_anno="{outpath}/03_variants/{caller}/05_annovar/{sample}.{ref_version}_multianno.txt"
	output:
		sub="{outpath}/03_variants/{caller}/06_score/{sample}.{ref_version}.exonic_splicing_multianno.txt",
		txt="{outpath}/03_variants/{caller}/06_score/{sample}.{ref_version}.rcnv_gnomadlof_multianno.txt"
	log:
		"{outpath}/03_variants/logs/{caller}/{sample}.{ref_version}.rcnv_lof.log"
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

rule caller_merge_vcf:
	input:
		passvcf=lambda wildcards: expand(["{{outpath}}/03_variants/{{caller}}/04_pass/{sample}.{ref_version}.pass.vcf.gz"], sample=[u.sample for u in units.itertuples()], ref_version=Ref_version),
		rcnv_gnomadlof=lambda wildcards: expand(["{{outpath}}/03_variants/{{caller}}/06_score/{sample}.{ref_version}.rcnv_gnomadlof_multianno.txt"], sample=[u.sample for u in units.itertuples()], ref_version=Ref_version)
	output:
		merged_vcf="{outpath}/03_variants/{caller}/07_merge_vcf/all.{caller}.pass.vcf.gz"
	log:
		"{outpath}/03_variants/logs/{caller}/all.{caller}.merge.log"
	params:
		is_single = lambda wildcards, input: len(input.passvcf) == 1
	threads:
		resource['resource']['high']['threads']
	resources:
		mem_mb=resource['resource']['high']['mem_mb']
	container:
		"/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/envs/bcftools_v1.10.2.sif"
	shell:
		'''
		if [ "{params.is_single}" == "True" ]; then
			echo "Single file detected. Using cp logic." > {log}
			cp {input.passvcf} {output.merged_vcf} >> {log} 2>&1
		else
			echo "Multiple files detected. Using bcftools merge." > {log}
			bcftools merge {input.passvcf} -Oz -o {output.merged_vcf} >> {log} 2>&1
		fi
		tabix -f -p vcf {output.merged_vcf} >> {log} 2>&1
		'''
