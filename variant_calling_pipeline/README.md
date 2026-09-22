# Variant calling pipeline

A Snakemake pipeline for calling, filtering and annotating SNPs, then reporting 
the mutations for the variants. This is based on the variant calling pipeline in
the compomics Y2 course in BSc Biochemistry at Imperial College London.

All settings (sample fastq files, reference genome fasta file, GATK filter settings, 
snpEff annotation genome, gene regions corresponding to the gene location in Ensembl) 
can be specified in the `config.yaml` file.


## Usage

1. Edit `config.yaml` to use your samples' FASTQ (.fq) files, reference FASTA (.fa) file, the snpeff genome, the dbSNP files and `regions`. 
   - BRCA1 and TP53 samples,  are examples and previously used in the compomics course.

2. Set working directory to this folder by running `cd /"insert path"/variant_calling_pipeline` in bash.

3. Run the following in bash, after activating the uv venv or conda environment (see Requirements):

```
snakemake -n                            # dry run to see what will run
snakemake --cores 8                     # run pipeline with the uv venv
snakemake --cores 8 --use-conda         # or run pipeline with conda
```

   Add `--rerun-incomplete` if the pipeline was aborted previously.

4. Outputs are in the `results` folder.



## Requirements

- Snakemake >= 7
- bwa, samtools, GATK, snpEff, htslib, bcftools, Python 

Install them in either a uv venv or conda.






