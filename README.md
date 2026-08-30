# Intra-Host Genomic Variation Analysis of Dengue Virus Serotype 2: Mild vs. Severe Infection

## Objective
This project investigates intra-host genomic diversity and mutation burden in Dengue Virus Serotype 2 (DENV-2) infections. By processing raw high-throughput sequencing data across clinical cohorts, the workflow compares viral variant counts, Transition/Transversion (Ti/Tv) ratios, and Allele Frequency (AF) spectrum distributions between **Mild (Primary)** and **Severe (Secondary)** infection cohorts.

---

## Datasets
The analysis compares two distinct cohorts of DENV-2 infected patients ($N=6$ total samples) sourced from SRA:

* **Group A: Mild / Primary Dengue Fever (DF)**
  * `DRR067444` (D2-2A-Hue-534-4MV)
  * `DRR067448` (D2-1A-Hue-519-3U)
  * `DRR067449` (D2-2A-Hue-534-4U)
* **Group B: Severe / Secondary Dengue (DHF)**
  * `DRR067462` (D2-2B-Hue-534-4U)
  * `DRR067463` (D2-3B-Hue-101-7U)
  * `DRR067465` (D2-5B-Hue-106-9U)
* **Reference Genome:** Dengue virus serotype 2 complete genome (`NC_001474.2`).

---

## Workflow Architecture

### 1. Preprocessing & Alignment (Linux / Terminal)
1. **Quality Control:** Run initial read inspection via `fastqc`.
2. **Trimming:** Filter low-quality bases (`-q 20`) and short reads (`--length_required 50`) using `fastp`.
3. **Alignment:** Index reference genome (`NC_001474.2`) and align trimmed FASTQ reads using `bwa mem`.
4. **BAM Processing:** Convert SAM to BAM, sort, and index via `samtools`.
5. **Variant Calling:** Call intra-host single nucleotide variants and indels using `freebayes` (`--min-alternate-fraction 0.05 --min-coverage 10 -p 1`).

```bash
# Quality trimming loop
for sample in DRR067444 DRR067448 DRR067449 DRR067462 DRR067463 DRR067465; do
  fastp -i rawfastq/${sample}.fastq.gz \
        -o trimmed_result/${sample}.trimmed.fastq.gz \
        -h qc/${sample}_fastp.html \
        -j qc/${sample}_fastp.json \
        -q 20 --length_required 50
done

# BWA Alignment, Indexing and sorting 
bwa index reference/dengue_ref.fasta
samtools faidx reference/dengue_ref.fasta

for sample in DRR067444 DRR067448 DRR067449 DRR067462 DRR067463 DRR067465; do
  bwa mem reference/dengue_ref.fasta trimmed_result/${sample}.trimmed.fastq.gz > aligned/${sample}.sam
  samtools view -bS aligned/${sample}.sam \vert{} samtools sort -o aligned/${sample}.sorted.bam
  samtools index aligned/${sample}.sorted.bam
  rm aligned/${sample}.sam
done

# Variant calling with Freebayes
for sample in DRR067444 DRR067448 DRR067449 DRR067462 DRR067463 DRR067465; do
  freebayes -f reference/dengue_ref.fasta \
            -p 1 \
            --min-alternate-fraction 0.05 \
            --min-coverage 10 \
            aligned/${sample}.sorted.bam > variants/${sample}.vcf
done
```

## 2. Downstream Analysis & Visualization (R / RStudio)
The resulting VCF files are parsed and analyzed in `R` using `vcfR`, `tidyverse`, and `ggplot2`. 
Quality Filtering: Variants filtered for QUAL >= 20. 
Transition/Transversion (Ti/Tv) Ratio: Evaluates substitution bias across cohorts. 
Allele Frequency (AF) Spectrum: Classifies variants into Major ($\text{AF} \ge 0.50$) and Minor ($\text{AF} < 0.50$) subpopulations. Statistical Significance: Evaluates variant count differences using Welch's Two-Sample $t$-test. 

The full, annotated R script is available in this repository.

```r
# Execute the R script file provided.
source("variants_analysis.R")
```

## Results and Data summaries

1. Total Raw variant counts per sample vs quality filtered variant counts.

| Sample ID | Cohort | Total Variant count | Quality filtered variant count |
| --- | ---| --- | --- |
| DRR067444 | Mild | 181 | 149 |
| DRR067448 | Mild | 210 | 128 |
| DRR067449 | Mild | 201 | 135 |
| DRR067462 | Severe | 510 | 258 |
| DRR067463 | Severe | 509 | 240 |
| DRR067465 | Severe | 472 | 260 |

2. Cohort Mutation Burden summary

| Cohort | Total samples | Total Variants | Avg. Variants per sample |
| --- | --- | --- | --- |
| **Mild** | 3 | 412 | 137.3 |
| **Severe** | 3 | 758 | 253.7 |

3. Transition vs Transversion (Ti/Tv) Ratio


|Cohort | Indel/Other| Transition| Transversion| Ti_Tv_Ratio|
|:------|-----------:|----------:|------------:|-----------:|
|Mild   |          65|        305|           42|        7.26|
|Severe |         167|        526|           65|        8.09|

4. Allele Frequency (AF) spectrum analysis.

|Cohort |variant_class              | total_count| avg_per_sample|
|:------|:--------------------------|-----------:|--------------:|
|Mild   |Major Variant (AF >= 0.50) |         412|          137.3 |
|Mild   |Minor Variant (AF < 0.05)  |         0  |          0   |
|Severe |Major Variant (AF >= 0.50) |         748|          249.3 |
|Severe |Minor Variant (AF < 0.05)  |         10 |          3.3 | 


## 3. Visualization
The plot below shows the distribution of variant counts between the two clinical groups. 

![Variant Burden Boxplot](variants_per_sample_type.png)

## Biological Conclusion:

* Significant Variant Increase: Severe Dengue samples showed an average of ~253 variants, whereas Mild samples averaged ~137 variants. 
* Statistical Validation: A Welch two-sample t- test confirmed that higher mutation burden in severe cases is statistically significant (p< 0.05)

The elevated variant count in severe Dengue cases reflects expanded **intra-host viral quasispecies diversity**. Higher rates of viral replication during severe infection lead to an accumulation of subclonal mutations, which may contribute to immune evasion and increased pathogenesis.


