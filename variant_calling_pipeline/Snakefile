import os
configfile: "config.yaml"


######## Settings from config.yaml #########

SAMPLES = list(config["samples"])

REF = config["reference"]
if REF.endswith(".gz"):
    raise ValueError(
        "config.yaml: `reference` must be an uncompressed FASTA "
        f"(e.g. genome.fa), not {REF}")
REF_FAI = f"{REF}.fai"
REF_DICT = os.path.splitext(REF)[0] + ".dict"  # GATK wants genome.dict next to genome.fa
BWA_INDEX = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")


REGIONS = config["regions"]

THREADS = config.get("threads") or {}

SNPEFF_MEM = config.get("snpeff_memory", "4g")
if config.get("snpeff_jar"):
    _jar = os.path.expanduser(config["snpeff_jar"])
    _java = os.path.expanduser(config.get("snpeff_java") or "java")
    SNPEFF_CMD = f'{_java} -Xmx{SNPEFF_MEM} -jar "{_jar}"'
else:
    SNPEFF_CMD = f"snpEff -Xmx{SNPEFF_MEM}"


FILTER_ARGS = " ".join(
    f'-filter-name "{name}" -filter "{expr}"'
    for name, expr in config["snp_filters"].items()
)


wildcard_constraints:
    sample="[^/]+",
    region="[^/]+",


rule all:
    input:
        expand(
            "results/{sample}/filtered_recalibrated_annotated_dbsnp_SNPs.vcf",
            sample=SAMPLES,
        ),
        expand(
            "results/{sample}/regions/{region}_variants.vcf",
            sample=SAMPLES,
            region=REGIONS,
        ),
        expand(
            "results/{sample}/regions/{region}_missense.tsv",
            sample=SAMPLES,
            region=REGIONS,
        ),


######## Indexing #########

rule bwa_index:
    input:
        REF,
    output:
        BWA_INDEX,
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa index "{input}""


rule samtools_faidx:
    input:
        REF,
    output:
        REF_FAI,
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools faidx "{input}""


rule samtools_dict:
    input:
        REF,
    output:
        REF_DICT,
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools dict "{input}" > "{output}""


######## Alignment #########

rule bwa_mem:
    input:
        ref=REF,
        idx=BWA_INDEX,
        r1=lambda wc: config["samples"][wc.sample]["r1"],
        r2=lambda wc: config["samples"][wc.sample]["r2"],
    output:
        "results/{sample}/aligned_reads.sam",
    threads: THREADS.get("bwa_mem", 4)
    conda:
        "envs/mapping.yaml"
    shell:
        r"bwa mem -t {threads} -R '@RG\tID:{wildcards.sample}\tLB:{wildcards.sample}\tPL:ILLUMINA\tPM:HISEQ\tSM:{wildcards.sample}' "{input.ref}" "{input.r1}" "{input.r2}" > "{output}""


######## Mark duplicates #########

rule mark_duplicates_spark:
    input:
        sam="results/{sample}/aligned_reads.sam",
    output:
        bam="results/{sample}/aligned_deduplicated_reads.bam",
        metrics="results/{sample}/dedup_metrics.txt",
    threads: THREADS.get("markduplicates", 4)
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk MarkDuplicatesSpark -I "{input.sam}" -M "{output.metrics}" -O "{output.bam}""


######## First Variant Calling (for BQSR) #########

rule haplotype_caller:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bam="results/{sample}/aligned_deduplicated_reads.bam",
    output:
        "results/{sample}/variants.vcf",
    threads: THREADS.get("haplotypecaller", 4)
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk HaplotypeCaller -R "{input.ref}" -I "{input.bam}" -O "{output}""


rule select_snps:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        vcf="results/{sample}/variants.vcf",
    output:
        "results/{sample}/SNPs.vcf",
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk SelectVariants -R "{input.ref}" -V "{input.vcf}" -select-type SNP -O "{output}""


rule filter_snps:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        vcf="results/{sample}/SNPs.vcf",
    output:
        "results/{sample}/filtered_SNPs.vcf",
    params:
        filters=FILTER_ARGS,
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk VariantFiltration -R "{input.ref}" -V "{input.vcf}" -O "{output}" {params.filters}"


rule select_bqsr_snps:
    input:
        vcf="results/{sample}/filtered_SNPs.vcf",
    output:
        "results/{sample}/BQSR_SNPs.vcf",
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk SelectVariants --exclude-filtered -V "{input.vcf}" -O "{output}""


######## BQSR #########

rule base_recalibrator:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bam="results/{sample}/aligned_deduplicated_reads.bam",
        known_sites="results/{sample}/BQSR_SNPs.vcf",
    output:
        "results/{sample}/recalibrated_bases.table",
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk BaseRecalibrator -R "{input.ref}" -I "{input.bam}" --known-sites "{input.known_sites}" -O "{output}""


rule apply_bqsr:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bam="results/{sample}/aligned_deduplicated_reads.bam",
        table="results/{sample}/recalibrated_bases.table",
    output:
        "results/{sample}/recalibrated_reads.bam",
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk ApplyBQSR -R "{input.ref}" -I "{input.bam}" -bqsr "{input.table}" -O "{output}""


######## Recalibrated variant calling #########

rule haplotype_caller_recal:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bam="results/{sample}/recalibrated_reads.bam",
    output:
        "results/{sample}/recalibrated_variants.vcf",
    threads: THREADS.get("haplotypecaller", 4)
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk HaplotypeCaller -R "{input.ref}" -I "{input.bam}" -O "{output}""


rule select_snps_recal:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        vcf="results/{sample}/recalibrated_variants.vcf",
    output:
        "results/{sample}/recalibrated_SNPs.vcf",
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk SelectVariants -R "{input.ref}" -V "{input.vcf}" -select-type SNP -O "{output}""


rule filter_snps_recal:
    input:
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        vcf="results/{sample}/recalibrated_SNPs.vcf",
    output:
        "results/{sample}/filtered_recalibrated_SNPs.vcf",
    params:
        filters=FILTER_ARGS,
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk VariantFiltration -R "{input.ref}" -V "{input.vcf}" -O "{output}" {params.filters} "


######## snpEff annotation #########

rule snpeff_annotate:
    input:
        vcf="results/{sample}/filtered_recalibrated_SNPs.vcf",
    output:
        vcf="results/{sample}/filtered_recalibrated_annotated_SNPs.vcf",
        summary="results/{sample}/snpEff_summary.html",
        genes="results/{sample}/snpEff_summary.genes.txt",
    params:
        snpeff=SNPEFF_CMD,
        genome=config["snpeff_genome"],
    conda:
        "envs/snpeff.yaml"
    shell:
        "{params.snpeff} -v -stats "{output.summary}" {params.genome} "{input.vcf}" > "{output.vcf}""


######## dbSNP annotation #########

rule bgzip_annotated:
    input:
        "results/{sample}/filtered_recalibrated_annotated_SNPs.vcf",
    output:
        "results/{sample}/filtered_recalibrated_annotated_SNPs.vcf.gz",
    conda:
        "envs/vcf_tools.yaml"
    shell:
        "bgzip -k "{input}""


rule tabix_index:
    input:
        "results/{sample}/filtered_recalibrated_annotated_SNPs.vcf.gz",
    output:
        "results/{sample}/filtered_recalibrated_annotated_SNPs.vcf.gz.tbi",
    conda:
        "envs/vcf_tools.yaml"
    shell:
        "tabix "{input}""


rule bcftools_annotate_dbsnp:
    input:
        vcf="results/{sample}/filtered_recalibrated_annotated_SNPs.vcf.gz",
        tbi="results/{sample}/filtered_recalibrated_annotated_SNPs.vcf.gz.tbi",
        dbsnp=config["dbsnp"],
    output:
        "results/{sample}/filtered_recalibrated_annotated_dbsnp_SNPs.vcf",
    conda:
        "envs/vcf_tools.yaml"
    shell:
        "bcftools annotate -c ID -a "{input.dbsnp}" "{input.vcf}" > "{output}""


######## SNPs in the Region  #########

rule parse_vcf:
    input:
        vcf="results/{sample}/filtered_recalibrated_annotated_dbsnp_SNPs.vcf",
        script="parsing scripts/parseVCF.py",
    output:
        "results/{sample}/regions/{region}_variants.vcf",
    params:
        chrom=lambda wc: REGIONS[wc.region]["chrom"],
        start=lambda wc: REGIONS[wc.region]["start"],
        end=lambda wc: REGIONS[wc.region]["end"],
    conda:
        "envs/vcf_tools.yaml"
    shell:
        "python "{input.script}" "{input.vcf}" {params.start} {params.end} --chrom {params.chrom} -o "{output}""


rule parse_vcf_transcripts:
    input:
        vcf="results/{sample}/filtered_recalibrated_annotated_dbsnp_SNPs.vcf",
        script="parsing scripts/parseVCFTranscripts.py",
    output:
        "results/{sample}/regions/{region}_missense.tsv",
    params:
        chrom=lambda wc: REGIONS[wc.region]["chrom"],
        start=lambda wc: REGIONS[wc.region]["start"],
        end=lambda wc: REGIONS[wc.region]["end"],
    conda:
        "envs/vcf_tools.yaml"
    shell:
        "python "{input.script}" "{input.vcf}" {params.start} {params.end} --chrom {params.chrom} > "{output}""
