# Code for finding unusal genomic regions 

### Authors: Qiuxuan Liu, Degang Wu, [Chaolong Wang](https://chaolongwang.github.io)

#### License: [GNU General Public License v3.0 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html)
---

## 1. Description

**UnRegFinder** is a pipeline to identify unusual genomic regions that distort population structure inference. Users can use this pipeline to 1) reproduce the unusual regions found and presented in our paper; 2) identify unsual genomic regions in their own data (whole genome sequencing data preferrably).


## 2. Citation for our pipeline 

Details of this pipeline can be found in our paper:  
* Qiuxuan Liu, Degang Wu, Chaolong Wang. Identification of genomic regions distorting population structure inference in diverse continental groups, under review


## 3. Dependencies
* snakemake (version >= 7.2.1)
* bcftools (version >= 1.14)
* plink2 (version: PLINK v2.00a3LM 64-bit Intel (27 Jul 2020))
* R (version >= 4.1.2)
** with LEA package installed
** with tidyverse package installed
* [LASER](https://sph.umich.edu/csg/software_pages/laser.html) (version = 2.04)
* bedtools (version >= 2.29.1)


## 4. Installation and configuration

To start, you can download our pipeline by the following command:

`git clone https://github.com/xxx/UnRegFinder.git` 

In subsequent sections, I will refer to the directory of the pipeline as `${pipeline_dir}`.

### 4.1 Configure the pipeline for your run

Use `${pipeline_dir}/config.json` as a template, create a configuration file for your run. Here is a short description of the sections needed specifying:
- `software`: specify the path to the depedent software binaries here, if they cannot be found in the default search paths.
- `path`: specify the `${work_dir}` (where the output is placed, preferrably a directory other than `${pipeline_dir}`) and `${input_dir}` (where the input vcf files are located)
- `input_fn_fn`: the VCF file name pattern. If your input vcf files are `prefix_chr1.vcf.gz`, `prefix_chr2.vcf.gz`, ..., write down `prefix_chr{chr}.vcf.gz` here.
- `pop_name`: the prefix of the output
- `sample_list`: if the input vcf files contain samples other than those you are going to use for subsequent analysis, specify the file that contains all sample names intended for subsequent analyses, one sample name per row. Otherwise, write down `""`.
- `max_iteration`: the number of iterations performed in the final step. Must be an integer >=2 and <=9.


## 5. Running the pipeline

### 5.1 Reproduce the unusual regions presented in our paper

1. Download high-coverage 1000G dataset from [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/) to a directory ${input_dir} using, for example, `wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{1..22}.recalibrated_variants.vcf.gz`.

2. Create a working directory, referred to as `${work_dir}`, and copy `find_unusual_regions.sh` and `config.json` from `${pipeline_dir}` to `${work_dir}`. 

3. Modify `find_unusual_regions.sh` as you see fit. Also, change `m` in `-j m` to limit the number of jobs running simultaneously. Most rules in our pipeline run sequentially except the first rule.

4. Modify `path/work_dir`, `path/input_dir`, `pop_name` and `sample_list` in `config.json`. In our paper, we only include samples that were present in the original 1KGP3 dataset. We provide the sample lists used under `${work_dir}/pop_list`.

5. Run `find_unusual_regions.sh`. This may take a while. The first rule could be IO-intensive so limit the number of jobs that can be run simultaneously.

6. After snakemake indicates the pipeline run has ended successfully, access the results in `${work_dir}/${pop_name}/regions_removed_pca/round${n}/${pop_name}.round${n}.unusual.regions.txt`.

Run `./find_unusual_regions.sh`.

### 5.2 Search for unusual regions in your own data

The procedure is similar to section 5.1, except that you need to specify your own input directory, the VCF filename pattern in `input_vcf_fn` of `config.json`.


## 6. Frequently encountered problems

### 6.1 



## 7. Questions
For further questions, please raise issues through github (recommended), or contact Degang Wu <dwuab@alumni.ust.hk> or Qiuxuan Liu <qiuxuan686@163.com>.
