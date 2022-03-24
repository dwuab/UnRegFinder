import os

# label for the population under investigation
pop_name = config["pop_name"]

# paths to softwares
bcftools = config["software"]["bcftools"]
plink = config["software"]["plink"]
plink2 = config["software"]["plink2"]
Rscript = config["software"]["Rscript"]
bedtools = config["software"]["bedtools"]

# various paths
work_dir=config["path"]["work_dir"]
qc_dir=os.path.join(work_dir,"1_qc/")
outlier_dir=os.path.join(work_dir,"2_outlier/")
pca_dir=os.path.join(work_dir,"3_pca/")
iteration_dir=os.path.join(work_dir,"regions_removed_pca/")

rule step1a_SNP_QC:
	input:
		vcf=os.path.join(config["path"]["input_dir"], config["input_vcf_fn"])
	output:
		bi_vcf=os.path.join(qc_dir, pop_name+".chr{chr}.biallelic.vcf.gz"),
		flag=touch(config["path"]["work_dir"]+"/flags/"+pop_name+".chr{chr}.QC.OK"),
		thin_vcf=os.path.join(qc_dir, pop_name+".chr{chr}.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf.gz")
	params:
		plink_prefix=qc_dir+pop_name+".chr{chr}.biallelic.maf_0.05.hwe_1e-8.thin2k",
		mem="5G",
		job_name="SNP_QC_chr{chr}"
	run:
		shell("if [ ! -d " + qc_dir + " ]; then mkdir -p " + qc_dir + "; fi")
		# if only a subset of samples should be included (config["sample_list"] is specified)
		if config["sample_list"] != "":
			sample_filter="-S "+config["sample_list"]
		else:
			sample_filter=""

		# extract bi-allelic SNPs
		shell("{bcftools} view " + sample_filter +
	    	" {input.vcf} -Ou "
	    	"| {bcftools} view -i'MAF>=0.05' -Ou "
	    	"| {bcftools} annotate -x INFO,^FORMAT/GT -Ou "
	    	"| {bcftools} norm -m +any -Ou "
	    	"| {bcftools} view -f PASS --max-alleles 2 -v snps "
	    	"--threads 10 -Oz > {output.bi_vcf}")

		# thin the SNPs
		shell("{plink2} "
	    	"--vcf {output.bi_vcf} "
        	"--hwe 1e-8 "
        	"--bp-space 2000 "
        	"--geno 0.05 "
        	"--mind 0.05 "
        	"--export vcf bgz id-paste=iid --out {params.plink_prefix}"
        )


rule step1b_combine:
	input:
		flags=expand(config["path"]["work_dir"]+"/flags/"+pop_name+".chr{chr}.QC.OK", chr=range(1,23)),
		thin_vcfs=expand(qc_dir+pop_name+".chr{chr}.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf.gz",chr=range(1,23))
	output:
		flag=touch(config["path"]["work_dir"]+"/flags/"+pop_name+".combine.OK"),
		combined_vcf=qc_dir+pop_name+".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf.gz"
	params:
		job_name="combine",
		mem="10G"
	run:
		# combine 22 autosomes
		shell("{bcftools} concat {input.thin_vcfs} "
			"-Ou | {bcftools} annotate --set-id +'%CHROM\\_%POS' --threads 40 -Oz > {output.combined_vcf}")


rule step1c_remove_PCA_outlier:
	input:
		flag=config["path"]["work_dir"]+"/flags/"+pop_name+".combine.OK",
		combined_vcf = qc_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf.gz"
	output:
		flag=touch(config["path"]["work_dir"]+"/flags/"+pop_name+".remove_outlier.OK"),
		bed=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.bed",
		fam=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.fam",
		vcf=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.vcf.gz"
	params:
		job_name="remove_PCA_outliers",
		mem="10G"
	run:
		shell("if [ ! -d " + outlier_dir + " ]; then mkdir -p " + outlier_dir + "; fi")
		geno_prefix = outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf2geno"
		output_prefix = outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed"
		shell(config["software"]["vcf2geno"] + " --inVcf {input.combined_vcf} "
			"--out {geno_prefix}")

		# remove PCA outliers; an outlier is such a sample that at least one of its
		# top 10 PCs is 6 SD away from the rest
		shell("cd {outlier_dir} && {workflow.basedir}/robustPCA.sh {geno_prefix} 10 6 " + 
			config["software"]["laser"] + " {Rscript} {workflow.basedir} {outlier_dir}/Results")

		shell("cd {outlier_dir}/Results && "
			"awk '{{print $1}}' removedSamples.txt | sort -n | uniq > {pop_name}.6sd.removed.ID")

		shell("{bcftools} view -S ^{outlier_dir}/Results/{pop_name}.6sd.removed.ID " +
			"{qc_dir}/{pop_name}.chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.vcf.gz " +
			"-Oz --threads 10 > {output.vcf}")

		shell("{plink2} --vcf {output.vcf} --make-bed --out {output_prefix}")


rule step2_PCA_GWAS:
	input:
		flag=config["path"]["work_dir"]+"/flags/"+pop_name+".remove_outlier.OK",
		bed=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.bed",
		fam=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.fam"
	output:
		flag=touch(config["path"]["work_dir"]+"/flags/"+pop_name+".PCA_GWAS.OK")
	params:
		job_name="PCA_GWAS",
		mem="10GB"
	run:
		shell("if [ ! -d " + pca_dir + " ]; then mkdir -p " + pca_dir + "; fi")
		for line in shell("cat {input.fam} | wc -l", iterable=True):
			sample_num=int(line.strip())

		prefix = outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed"
		prefix2 = pca_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed"
		shell("{plink} --bfile {prefix} --pca {sample_num} --out {prefix}.allpc")

		# rename the .eigenval file for next step
		shell("if [ ! -e {prefix}.allpc.eigenvalues ]; then ln -s {prefix}.allpc.eigenval {prefix}.allpc.eigenvalues; fi")

		# step-2. GWAS of significant single PC
		line_no=1
		pc_num=0
		for line in shell(Rscript + " " + workflow.basedir + "/tracy_widom.R {sample_num} {outlier_dir} " + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.allpc.eigenvalues", iterable=True):
			if (line_no==9):
				pc_num=int(line.strip())
			line_no=line_no+1

		shell("{plink} --bfile {prefix} --pca {pc_num} --out {prefix}.pc")

		shell("cp {prefix}.pc.eigen* {pca_dir}")
		shell("cp {prefix}.b* {pca_dir}")
		shell("cp {prefix}.fam {pca_dir}")

		shell("sed \"1i\\\\$(eval echo FID IID PC{{1..{pc_num}..1}})\" {prefix2}.pc.eigenvec > "
			"{prefix2}.pc.eigenvec.pheno")

		shell("{plink} --allow-no-sex --bfile {prefix2} --linear --pheno {prefix2}.pc.eigenvec.pheno --all-pheno --out {prefix2}")

		# output: 
		# EUR.chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.PC*.assoc.linear

		# step-3. adjust lambdaGC for each GWAS of PC 
		# combine association signals for each SNP across significant PCs using the Sidak’s method

		shell("{Rscript} {workflow.basedir}/sidak.R {pca_dir} {pop_name}.chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.PC {pc_num}")
		

rule step3_clumping:
	input:
		flag=config["path"]["work_dir"]+"/flags/"+pop_name+".PCA_GWAS.OK"
	output:
		flag=touch(config["path"]["work_dir"]+"/flags/iteration_1.OK")
	params:
		job_name="clumping",
		mem="10G"
	run:
		shell("if [ ! -d {iteration_dir} ]; then mkdir -p {iteration_dir}; fi")
		prefix = pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed"
		shell("{plink} --bfile {outlier_dir}/{prefix} "
			"--clump {pca_dir}/{prefix}.PCall.assoc.linear "
			"--clump-field P_concat "
			"--clump-p1 5e-8 "
			"--clump-p2 0.001 "
			"--clump-kb 10000 "
			"--clump-r2 0.3 "
			"--clump-allow-overlap "
			"--out {pca_dir}/{pop_name}")

		shell("{bcftools} query -f '%CHROM %POS %ID\n' {outlier_dir}/{prefix}.vcf.gz > {pca_dir}/{pop_name}.rsID.bp.txt")

		shell("{Rscript} {workflow.basedir}/format_clumped_results.R {pca_dir} {pop_name}")

		shell("mkdir -p {iteration_dir}/round1")
		shell(config["software"]["bedtools"] + " merge -i <(awk '{{OFS=\"\\t\"}} NR != 1 {{if ($5>=5 && $7!=\"\" && $8!=\"\") print $1, $7*100, $8*100}}' {pca_dir}/{pop_name}.clump.cleaned.txt |sort -k1n -k2n -k3n)|awk '{{OFS=\"\\t\"}} {{print $1,$2*10000, $3*10000}}' > {iteration_dir}/round1/{pop_name}.round1.unusual.regions.txt")
		
		shell("cd {iteration_dir}/round1 && ln -s {pca_dir}/{pop_name}.clump.cleaned.txt {pop_name}.round1.clump.cleaned.txt")

rule step4_iteration:
	input:
		flag=lambda wildcards: config["path"]["work_dir"]+"/flags/iteration_"+str(int(wildcards.i)-1)+".OK",
		fam=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.fam",
		vcf=outlier_dir + pop_name + ".chr1_22.biallelic.maf_0.05.hwe_1e-8.thin2k.6sdOutlier_removed.vcf.gz"
	output:
		flag=touch(config["path"]["work_dir"]+"/flags/iteration_{i}.OK")
	wildcard_constraints:
		i="[2-9]"
	params:
		job_name="iteration_{i}",
		mem="10G"
	run:
		i=int(wildcards.i)
		cur_it_dir=os.path.join(iteration_dir, "round"+str(i)+"/")
		last_it_dir=os.path.join(iteration_dir, "round"+str(i-1)+"/")
		shell("mkdir -p {cur_it_dir}")

		for line in shell("cat {input.fam} | wc -l", iterable=True):
                        sample_num=int(line.strip())

		shell("{bcftools} view {input.vcf} -T ^{last_it_dir}/{pop_name}.round"+str(i-1)+".unusual.regions.txt"
			" -Oz > {cur_it_dir}/{pop_name}.removed.round{i}.vcf.gz")

		shell("{plink} --vcf {cur_it_dir}/{pop_name}.removed.round{i}.vcf.gz "
			"--pca {sample_num} "
			"--out {cur_it_dir}/{pop_name}.removed.round{i}.allpc")

		shell("mv {cur_it_dir}/{pop_name}.removed.round{i}.allpc.eigenval {cur_it_dir}/{pop_name}.removed.round{i}.allpc.eigenvalues")

		line_no=1
                pc_num=0
                for line in shell(Rscript + " " + workflow.basedir + "/tracy_widom.R {sample_num} {cur_it_dir} {pop_name}.removed.round{i}.allpc.eigenvalues", iterable=True):
                        if (line_no==9):
                                pc_num=int(line.strip())
                        line_no=line_no+1

		shell("{plink} --vcf {cur_it_dir}/{pop_name}.removed.round{i}.vcf.gz "
                        "--pca {pc_num} "
                        "--out {cur_it_dir}/{pop_name}.removed.round{i}.pc")

		shell("sed \"1i\\\\$(eval echo FID IID PC{{1..{pc_num}..1}})\" " +
			"{cur_it_dir}/{pop_name}.removed.round{i}.pc.eigenvec > " +
			"{cur_it_dir}/{pop_name}.removed.round{i}.pc.eigenvec.pheno")

		shell("{plink} --vcf {input.vcf} --allow-no-sex "
			"--linear "
			"--pheno {cur_it_dir}/{pop_name}.removed.round{i}.pc.eigenvec.pheno "
			"--all-pheno "
			"--out {cur_it_dir}/{pop_name}.removed.round{i}")

		shell("Rscript {workflow.basedir}/sidak.R {cur_it_dir} {pop_name}.removed.round{i}.PC {pc_num}")

		prefix = pop_name + ".removed.round"+str(i)
                shell("{plink} --vcf {cur_it_dir}/{prefix}.vcf.gz "
                        "--clump {cur_it_dir}/{prefix}.PCall.assoc.linear "
                        "--clump-field P_concat "
                        "--clump-p1 5e-8 "
                        "--clump-p2 0.001 "
                        "--clump-kb 10000 "
                        "--clump-r2 0.3 "
                        "--clump-allow-overlap "
                        "--out {cur_it_dir}/{pop_name}")

                shell("{bcftools} query -f '%CHROM %POS %ID\n' {cur_it_dir}/{prefix}.vcf.gz > {cur_it_dir}/{pop_name}.rsID.bp.txt")

                shell("{Rscript} {workflow.basedir}/format_clumped_results.R {cur_it_dir} {pop_name}")

		shell("mv {cur_it_dir}/{pop_name}.clump.cleaned.txt {cur_it_dir}/tmp && cat {cur_it_dir}/tmp <(tail -n +2 {last_it_dir}/{pop_name}.round"+str(i-1)+".clump.cleaned.txt) > {cur_it_dir}/{pop_name}.round{i}.clump.cleaned.txt")

                shell("{bedtools} merge -i <(awk '{{OFS=\"\\t\"}} NR != 1 {{if ($5>=5 && $7!=\"\" && $8!=\"\") print $1, $7*100, $8*100}}' {cur_it_dir}/{pop_name}.round{i}.clump.cleaned.txt |sort -k1n -k2n -k3n)|awk '{{OFS=\"\\t\"}} {{print $1,$2*10000, $3*10000}}' > {cur_it_dir}/{pop_name}.round{i}.unusual.regions.txt")

                shell("mkdir -p {iteration_dir}/round"+str(i+1))

rule all:
	input:
		flag=config["path"]["work_dir"]+"/flags/iteration_"+str(config["max_iteration"])+".OK"
