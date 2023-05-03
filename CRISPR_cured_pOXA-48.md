# Setup, read trimming and genome assembly

Create a project directory from which run all commands and to store outputs:
```sh
mkdir CRISPR_cured_pOXA-48
cd CRISPR_cured_pOXA-48
```

**Data retrieval from BioProjects:**

* From BioProject **PRJNA626430**, download raw Illumina and Nanopore/PacBio reads to the corresponding directories, `../Reads_WT_pOXA-48/raw/` or `../Reads_long`, as indicated in the table below. Also, download the indicated draft assemblies to `../Assemblies_WT_pOXA-48/`. These are all wild-type (WT) strains.
* From BioProject **PRJNA838107**, download the closed assemblies and Illumina reads of C288, K153 and K163 to `../Closed_sequences/` and `../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/`, respectively.
* From BioProject **PRJNA962867**, download raw Illumina reads to `./reads/raw/`. These are WT strains that could be cured, and a subsample of nine cured strains (indicated with the suffix c1 or c2).

BioSample accessions are provided in Table S1 of the manuscript.


|BioProject|Read sequence code|Sample name|Destination|
|:----|:----|:----|:----|
|PRJNA626430|WTCHG_354241_281104|CF12|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|ccs_CF12_lbc70|CF12|../Reads_long/|
|PRJNA626430|WTCHG_354241_282116|CF13|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|ccs_CF13_lbc62|CF13|../Reads_long/|
|PRJNA626430|WTCHG_423853_206166|H53|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|H53_nanopore|H53|../Reads_long/|
|PRJNA626430|WTCHG_423853_208190|J57|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|J57_nanopore|J57|../Reads_long/|
|PRJNA626430|WTCHG_423853_218120|K147|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|K147_nanopore|K147|../Reads_long/|
|PRJNA626430|WTCHG_423853_261149|C325|../Reads_WT_pOXA-48/raw/|
|PRJNA626430|C325_nanopore|C325|../Reads_long/|
|PRJNA626430|-|C609|../Assemblies_WT_pOXA-48/|
|PRJNA626430|-|C642|../Assemblies_WT_pOXA-48/|
|PRJNA626430|-|C662|../Assemblies_WT_pOXA-48/|
|PRJNA626430|-|R10|../Assemblies_WT_pOXA-48/|
|PRJNA838107|C288WT_S139|C288|../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/|
|PRJNA838107|-|C288|../Closed_sequences/|
|PRJNA838107|K153WT_S103|K153|../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/|
|PRJNA838107|-|K153|../Closed_sequences/|
|PRJNA838107|K163WT_S99|K163|../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/|
|PRJNA838107|-|K163|../Closed_sequences/|
|PRJNA962867|2933_1_S37|AJ_K25|./reads/raw/|
|PRJNA962867|2933_5_S41|AJ_K127|./reads/raw/|
|PRJNA962867|2933_7_S43|AJ_K293|./reads/raw/|
|PRJNA962867|2933_8_S44|AJ_K308|./reads/raw/|
|PRJNA962867|2933_11_S47|AJ_K88|./reads/raw/|
|PRJNA962867|2933_13_S49|AJ_N23|./reads/raw/|
|PRJNA962867|2933_17_S53|AJ_K219|./reads/raw/|
|PRJNA962867|2933_18_S54|AJ_K244|./reads/raw/|
|PRJNA962867|2933_19_S55|AJ_C310|./reads/raw/|
|PRJNA962867|2933_21_S57|AJ_S8|./reads/raw/|
|PRJNA962867|2933_22_S58|AJ_J61|./reads/raw/|
|PRJNA962867|2933_24_S60|AJ_K198|./reads/raw/|
|PRJNA962867|2933_25_S61|AJ_K318|./reads/raw/|
|PRJNA962867|2933_26_S62|AJ_C163|./reads/raw/|
|PRJNA962867|2933_27_S1|AJ_C164|./reads/raw/|
|PRJNA962867|2933_28_S2|AJ_C527|./reads/raw/|
|PRJNA962867|2933_29_S3|AJ_Nh27|./reads/raw/|
|PRJNA962867|2933_30_S4|AJ_L38|./reads/raw/|
|PRJNA962867|2933_31_S5|AJ_Z29|./reads/raw/|
|PRJNA962867|2933_48_S22|AJ_C646|./reads/raw/|
|PRJNA962867|2933_49_S23|AJ_C728|./reads/raw/|
|PRJNA962867|2933_50_S24|AJ_N46|./reads/raw/|
|PRJNA962867|WTCHG_678479_73745350|CF13c1|./reads/raw/|
|PRJNA962867|WTCHG_678479_73755351|J57c1|./reads/raw/|
|PRJNA962867|WTCHG_678479_73765352|H53c1|./reads/raw/|
|PRJNA962867|WTCHG_678479_73775353|C325c1|./reads/raw/|
|PRJNA962867|WTCHG_678479_73785354|CF12c1|./reads/raw/|
|PRJNA962867|WTCHG_678479_73795355|K147c1|./reads/raw/|
|PRJNA962867|C288c2_Clon_2_S141|C288c2|./reads/raw/|
|PRJNA962867|K153c2_S104|K153c2|./reads/raw/|
|PRJNA962867|K163c1_S100|K163c1|./reads/raw/|


All Illumina reads were trimmed with Trim Galore v0.6.4 (Cutadapt v2.8) using the following commands. Note that the base name is the preffix WTCHG, SQCST or MiGS, followed by the strain name (e.g. WTCHG_CF12 or SQCST_AJ_K25):

```sh
# WT strains (PRJNA626430)
trim_galore --quality 20 --length 50 --fastqc --basename WTCHG_<strain_name> --output_dir ../Reads_WT_pOXA-48/trimmed --paired ../Reads_WT_pOXA-48/raw/<fq1> ../Reads_WT_pOXA-48/raw/<fq2>
# WT strains (PRJNA838107)
trim_galore --quality 20 --length 50 --fastqc --basename MiGS_<strain_name> --output_dir ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/ --paired ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/<fq1> ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_raw/<fq2>
# WT strains (PRJNA962867)
trim_galore --quality 20 --length 50 --fastqc --basename SQCST_<strain_name> --output_dir ./reads/trimmed --paired ./reads/raw/<fq1> ./reads/raw/<fq2>
# Cured strains (PRJNA962867: C325c1, CF12c1, CF13c1, H53c1, J57c1, K147c1)
trim_galore --quality 20 --length 50 --fastqc --basename WTCHG_<strain_name> --output_dir ./reads/trimmed --paired ./reads/raw/<fq1> ./reads/raw/<fq2>
# Cured strains (PRJNA962867: C288c2, K153c2, K163c1)
trim_galore --quality 20 --length 50 --fastqc --basename MiGS_<strain_name> --output_dir ./reads/trimmed --paired ./reads/raw/<fq1> ./reads/raw/<fq2>
# FastQC reports are generated for read quality analysis
```

The genomes of the WT and cured strains (BioProject PRJNA962867) were *de novo* assembled with SPAdes v3.15.2. Assembly quality was assessed with QUAST v5.0.2:

```sh
# Genome assembly
for fq1 in ./reads/trimmed/*val_1.fq.gz
do
	fq2=${fq1%%val_1.fq.gz}"val_2.fq.gz"
	strain=${fq1:22}
	strain=${strain::-12}
	spades.py --isolate --cov-cutoff auto -o ./assemblies_SPAdes/$strain -1 $fq1 -2 $fq2
done
# QUAST
for contigs in ./assemblies_SPAdes/*/contigs.fasta
do
	dir=${contigs::-14}
	quast.py $contigs -o $dir/QUAST/
done
```

Hybrid assemblies of WT strains were generated from trimmed Illumina reads and raw Nanopore/PacBio reads with Unicycler v0.4.9:

```sh
# CF12
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_2.fq.gz -l ../Reads_long/ccs_CF12_lbc70.fastq -o ./assemblies_unicycler/CF12
# CF13
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_2.fq.gz -l ../Reads_long/ccs_CF13_lbc62.fastq -o ./assemblies_unicycler/CF12
# C325
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_2.fq.gz -l ../Reads_long/C325_nanopore.fastq.gz -o ./assemblies_unicycler/C325
# H53
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_2.fq.gz -l ../Reads_long/H53_nanopore.fastq.gz -o ./assemblies_unicycler/H53
# J57
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_2.fq.gz -l ../Reads_long/J57_nanopore.fastq.gz -o ./assemblies_unicycler/J57
# K147
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_1.fq.gz -2 ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_2.fq.gz -l ../Reads_long/K147_nanopore.fastq.gz -o ./assemblies_unicycler/K147
```

Completeness was confirmed by visualizing the assemblies in Bandage v0.8.1. The chromosome of CF12 was fragmented. In case short-read processing could be affecting assembly quality, Illumina reads of CF12 and CF13 —since both were sequenced with PacBio— were alternatively trimmed with Trimmomatic v0.39, using the maximum information quality filtering algorithm (MAXINFO), which balances the benefits of keeping longer reads against the costs of retaining low-quality bases. The targeted minimum read length was set to 50 and the strictness to 0.8 (MAXINFO:50:0.8). Then, hybrid assemblies were obtained with Unicycler v0.4.9 as before:

```sh
# Read trimming
java -jar trimmomatic-0.39.jar PE -phred33 ../Reads_WT_pOXA-48/raw/WTCHG_354241_281104_1.fastq.gz ../Reads_WT_pOXA-48/raw/WTCHG_354241_281104_2.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_1_paired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_1_unpaired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_2_paired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_2_unpaired.fastq.gz ILLUMINACLIP:/home/ltorcel/Software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 MAXINFO:50:0.8
java -jar trimmomatic-0.39.jar PE -phred33 ../Reads_WT_pOXA-48/raw/WTCHG_354241_282116_1.fastq.gz ../Reads_WT_pOXA-48/raw/WTCHG_354241_282116_2.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_1_paired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_1_unpaired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_2_paired.fastq.gz ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_2_unpaired.fastq.gz ILLUMINACLIP:/home/ltorcel/Software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 MAXINFO:50:0.8
# Read quality analysis
fastqc ../Reads_WT_pOXA-48/trimmomatic_maxinfo/* -o ../Reads_WT_pOXA-48/trimmomatic_maxinfo/
# Hybrid assemblies
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_1_paired.fastq.gz -2 ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF12_2_paired.fastq.gz -l ../Reads_long/ccs_CF12_lbc70.fastq -o ./assemblies_unicycler/CF12_maxinfo
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_1_paired.fastq.gz -2 ../Reads_WT_pOXA-48/trimmomatic_maxinfo/WTCHG_CF13_2_paired.fastq.gz -l ../Reads_long/ccs_CF13_lbc62.fastq -o ./assemblies_unicycler/CF13_maxinfo
```

This way, the chromosomes and plasmids of both strains were closed. However, the assembly of strain CF13 WT was missing a 6 Kb ColRNAI plasmid that was present in the CF13 cured strain (NODE_30). The presence of the plasmid in CF13 WT was confirmed by PCR (forward primer: 5'-TCGTGATCGGCAAGGATACG-3'; reverse primer: 5'-CGATGAACGCAAGAAGCTCG-3'). Therefore, this plasmid was probably lost during DNA extraction for sequencing. Since CF12 and CF13 are clonal (isolated at different time points), the sequence of the ColRNAI plasmid from CF12 was pasted to the FASTA file of CF13. It was later confirmed that this plasmid is identical between CF12 and CF13 (see variant calling results). The FASTA sequences of the six complete assemblies were copied to `../Closed_sequences`.

The correct assembly of plasmid pOXA-48 was assessed by inspecting plasmid length: in agreement with DelaFuente *et al.* 2022 (https://doi.org/10.1038/s41559-022-01908-7), C325 and K147 carry the full-length plasmid (65499 bp) and the rest of strains carry a plasmid variant lacking the *ltrA* gene (14883-16369 deletion, 63589 bp). Additionally, the assemblies were aligned to the pOXA-48_K8 reference plasmid (MT441554) with minimap2 v2.21, confirming the *ltrA* deletion and the two SNPs in H53 (22416-C:T, 28222-A:G).

```sh
mkdir ./alignments_minimap2/
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/C325.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_C325.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_C325.bam
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/CF12.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_CF12.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_CF12.bam
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/CF13.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_CF13.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_CF13.bam
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/H53.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_H53.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_H53.bam
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/J57.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_J57.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_J57.bam
minimap2 -ax asm5 ../Closed_sequences/plasmids/pOXA-48_K8.fasta ../Closed_sequences/K147.fasta | samtools sort -o ./alignments_minimap2/ref_pOXA-48_K8_map_K147.bam
samtools index ./alignments_minimap2/ref_pOXA-48_K8_map_K147.bam
```



# Phylogenetic tree

A mash distance phylogenetic tree of cured strains was constructed with mashtree v1.2.0, with a bootstrap of 100:

```sh
mkdir ./phylogeny_mash/

# Rename contig files of the WT strains
for dir in ./assemblies_SPAdes/AJ*
do
	name=${dir:20}
	mv $dir/contigs.fasta $dir"/"$name"_contigs.fasta"
done

# Copy genome sequences of cured strains
cp ./assemblies_SPAdes/AJ_C646/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_C728/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_N46/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K127/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K293/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K308/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K88/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_N23/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K25/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K219/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K244/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_C310/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_S8/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_J61/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K198/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_K318/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_C163/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_C164/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_C527/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_Nh27/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_L38/*contigs.fasta ./phylogeny_mash/
cp ./assemblies_SPAdes/AJ_Z29/*contigs.fasta ./phylogeny_mash/
cp ../Assemblies_WT_pOXA-48/C642_contigs.fasta ./phylogeny_mash/
cp ../Assemblies_WT_pOXA-48/C609_contigs.fasta ./phylogeny_mash/
cp ../Assemblies_WT_pOXA-48/C662_contigs.fasta ./phylogeny_mash/
cp ../Assemblies_WT_pOXA-48/R10_contigs.fasta ./phylogeny_mash/
cp ../Closed_sequences/C325.fasta ./phylogeny_mash/
cp ../Closed_sequences/CF12.fasta ./phylogeny_mash/
cp ../Closed_sequences/CF13.fasta ./phylogeny_mash/
cp ../Closed_sequences/H53.fasta ./phylogeny_mash/
cp ../Closed_sequences/J57.fasta ./phylogeny_mash/
cp ../Closed_sequences/K147.fasta ./phylogeny_mash/
cp ../Closed_sequences/C288*.fna ./phylogeny_mash/C288.fasta
cp ../Closed_sequences/K153*.fna ./phylogeny_mash/K153.fasta
cp ../Closed_sequences/K163*.fna ./phylogeny_mash/K163.fasta

# Mash tree
mashtree_bootstrap.pl --numcpus 18 --reps 100 ./phylogeny_mash/*fasta > ./phylogeny_mash/mashtree_btrp.dnd
```



# Plasmid content of WT and cured strains

ABRicate v1.0.1 was used to find plasmid replicons in the genomes of the WT and cured strains. As expected, all cured strains lacked only the IncL pOXA-48 replicon (Table S3).

```sh
mkdir ./plasmids_abricate/
abricate --db plasmidfinder ../Closed_sequences/C325.fasta > ./plasmids_abricate/C325.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/C325c1/*contigs.fasta > ./plasmids_abricate/C325c1.tsv
abricate --db plasmidfinder ../Closed_sequences/CF12.fasta > ./plasmids_abricate/CF12.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/CF12c1/*contigs.fasta > ./plasmids_abricate/CF12c1.tsv
abricate --db plasmidfinder ../Closed_sequences/CF13.fasta > ./plasmids_abricate/CF13.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/CF13c1/*contigs.fasta > ./plasmids_abricate/CF13c1.tsv
abricate --db plasmidfinder ../Closed_sequences/H53.fasta > ./plasmids_abricate/H53.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/H53c1/*contigs.fasta > ./plasmids_abricate/H53c1.tsv
abricate --db plasmidfinder ../Closed_sequences/J57.fasta > ./plasmids_abricate/J57.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/J57c1/*contigs.fasta > ./plasmids_abricate/J57c1.tsv
abricate --db plasmidfinder ../Closed_sequences/K147.fasta > ./plasmids_abricate/K147.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/K147c1/*contigs.fasta > ./plasmids_abricate/K147c1.tsv
abricate --db plasmidfinder ../Closed_sequences/C288WT.fna > ./plasmids_abricate/C288.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/C288c2/*contigs.fasta > ./plasmids_abricate/C288c2.tsv
abricate --db plasmidfinder ../Closed_sequences/K153WT-30A7.fna > ./plasmids_abricate/K153.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/K153c2/*contigs.fasta > ./plasmids_abricate/K153c2.tsv
abricate --db plasmidfinder ../Closed_sequences/K163WT-28F1.fna > ./plasmids_abricate/K163.tsv
abricate --db plasmidfinder ./assemblies_SPAdes/K163c1/*contigs.fasta > ./plasmids_abricate/K163c1.tsv

abricate --summary ./plasmids_abricate/* > ./plasmids_abricate/summary.tsv
```




# Identifying mutations in cured strains

The closed genomes of the WT strains were annotated with PGAP v2021-07-01.build5508. The annotation GBK and GFF files were also copied to `../Closed_sequences`. The draft genomes of the cured strains were annotated with Prokka v1.14.6:

```sh
prokka --outdir ./assemblies_SPAdes/annotations_prokka/C325c1/ ./assemblies_SPAdes/C325c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/CF12c1/ ./assemblies_SPAdes/CF12c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/CF13c1/ ./assemblies_SPAdes/CF13c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/H53c1/ ./assemblies_SPAdes/H53c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/J57c1/ ./assemblies_SPAdes/J57c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/K147c1/ ./assemblies_SPAdes/K147c1/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/C288c2/ ./assemblies_SPAdes/C288c2/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/K153c2/ ./assemblies_SPAdes/K153c2/*contigs.fasta
prokka --outdir ./assemblies_SPAdes/annotations_prokka/K163c1/ ./assemblies_SPAdes/K163c1/*contigs.fasta
```

Variant calling was performed with breseq v0.35.6 or v0.35.7. First, to discard false positive calls due to misassemblies, the Illumina reads of the WT strains were mapped to their respective closed genomes. Then, the Illumina reads of the cured strains were mapped to their corresponding WT reference to identify mutations in cured strains. To further confirm or discard confusing mutations, the Illumina reads of the WT strains were mapped to the draft genomes of the cured strains.

```sh
# WT to WT
breseq -o ./variants_breseq/ref_C325_map_C325 -r ../Closed_sequences/C325.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF12_map_CF12 -r ../Closed_sequences/CF12.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF13_map_CF13 -r ../Closed_sequences/CF13.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_2.fq.gz
breseq -o ./variants_breseq/ref_H53_map_H53 -r ../Closed_sequences/H53.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_2.fq.gz
breseq -o ./variants_breseq/ref_J57_map_J57 -r ../Closed_sequences/J57.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_2.fq.gz
breseq -o ./variants_breseq/ref_K147_map_K147 -r ../Closed_sequences/K147.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_2.fq.gz
breseq -o ./variants_breseq/ref_C288_map_C288 -r ../Closed_sequences/C288WT.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_C288WT_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_C288WT_val_2.fq.gz
breseq -o ./variants_breseq/ref_K153_map_K153 -r ../Closed_sequences/K153WT-30A7.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_K153_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_K153_val_2.fq.gz
breseq -o ./variants_breseq/ref_K153_map_K153 -r ../Closed_sequences/K163WT-28F1.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_K163_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/reads_trimmed/MiGS_K163_val_2.fq.gz

# Cured to WT
breseq -o ./variants_breseq/ref_C325_map_C325c1 -r ../Closed_sequences/C325.gbk ./reads/trimmed/WTCHG_C325c1_val_1.fq.gz ./reads/trimmed/WTCHG_C325c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF12_map_CF12c1 -r ../Closed_sequences/CF12.gbk ./reads/trimmed/WTCHG_CF12c1_val_1.fq.gz ./reads/trimmed/WTCHG_CF12c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF13_map_CF13c1 -r ../Closed_sequences/CF13.gbk ./reads/trimmed/WTCHG_CF13c1_val_1.fq.gz ./reads/trimmed/WTCHG_CF13c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_H53_map_H53c1 -r ../Closed_sequences/H53.gbk ./reads/trimmed/WTCHG_H53c1_val_1.fq.gz ./reads/trimmed/WTCHG_H53c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_J57_map_J57c1 -r ../Closed_sequences/J57.gbk ./reads/trimmed/WTCHG_J57c1_val_1.fq.gz ./reads/trimmed/WTCHG_J57c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_K147_map_K147c1 -r ../Closed_sequences/K147.gbk ./reads/trimmed/WTCHG_K147c1_val_1.fq.gz ./reads/trimmed/WTCHG_K147c1_val_2.fq.gz
breseq -o ./variants_breseq/ref_C288_map_C288c2 -r ../Closed_sequences/C288WT.gbk ./reads/trimmed/MiGS_C288c2_Clon_2_val_1.fq.gz ./reads/trimmed/MiGS_C288c2_Clon_2_val_2.fq.gz
breseq -o ./variants_breseq/ref_K153_map_K153c2 -r ../Closed_sequences/K153WT-30A7.gbk ./reads/trimmed/MiGS_K153c2_val_1.fq.gz ./reads/trimmed/MiGS_K153c2_val_2.fq.gz
breseq -o ./variants_breseq/ref_K163_map_K163c1 -r ../Closed_sequences/K163WT-28F1.gbk ./reads/trimmed/MiGS_K163c1_val_1.fq.gz ./reads/trimmed/MiGS_K163c1_val_2.fq.gz

# WT to cured
breseq -o ./variants_breseq/ref_C325c1_map_C325 -r ./assemblies_SPAdes/annotations_prokka/C325c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_C325_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF12c1_map_CF12 -r ./assemblies_SPAdes/annotations_prokka/CF12c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_CF12_val_2.fq.gz
breseq -o ./variants_breseq/ref_CF13c1_map_CF13 -r ./assemblies_SPAdes/annotations_prokka/CF13c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_CF13_val_2.fq.gz
breseq -o ./variants_breseq/ref_H53c1_map_H53 -r ./assemblies_SPAdes/annotations_prokka/H53c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_H53_val_2.fq.gz
breseq -o ./variants_breseq/ref_J57c1_map_J57 -r ./assemblies_SPAdes/annotations_prokka/J57c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_J57_val_2.fq.gz
breseq -o ./variants_breseq/ref_K147c1_map_K147 -r ./assemblies_SPAdes/annotations_prokka/K147c1/*.gbk ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_1.fq.gz ../Reads_WT_pOXA-48/trimmed/WTCHG_K147_val_2.fq.gz
breseq -o ./variants_breseq/ref_C288c2_map_C288 -r ./assemblies_SPAdes/annotations_prokka/C288c2/*.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_C288WT_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_C288WT_val_2.fq.gz
breseq -o ./variants_breseq/ref_K153c2_map_K153 -r ./assemblies_SPAdes/annotations_prokka/K153c2/*.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_K153_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_K153_val_2.fq.gz
breseq -o ./variants_breseq/ref_K163c1_map_K163 -r ./assemblies_SPAdes/annotations_prokka/K163c1/*.gbk ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_K163_val_1.fq.gz ../Within_patient_evolution/within_patient_evol_cured/reads_trimmed/MiGS_K163_val_2.fq.gz
```

Mutations (SNPs, deletions and genomic rearrangements) were analyzed and summarized in Table S4. Additionally, breseq results further confirmed that no other plasmid was eliminated during pOXA-48 curing (missing coverage evidence, see Table S4).



# Confirming mutations are not due to off-target CRISPR activity

To determine whether the identified variants could be due to off-target cuts of the CRISPR-Cas9 system, the nucleotide sequences of sgOXA48 and sgPEMK were aligned to a 40 bp subregion enclosing the mutations using EMBOSS Needle v6.6.0.0:

```sh
mkdir ./analysis_mutations/

needle ./analysis_mutations/sgOXA-48_f.fasta ./analysis_mutations/seqs_enclosing_mutations.fasta -outfile ./analysis_mutations/needle_sgOXA-48_f.txt -gapopen 10.0 -gapextend 0.5
needle ./analysis_mutations/sgOXA-48_r.fasta ./analysis_mutations/seqs_enclosing_mutations.fasta -outfile ./analysis_mutations/needle_sgOXA-48_r.txt -gapopen 10.0 -gapextend 0.5
needle ./analysis_mutations/sgPemK_f.fasta ./analysis_mutations/seqs_enclosing_mutations.fasta -outfile ./analysis_mutations/needle_sgPemK_f.txt -gapopen 10.0 -gapextend 0.5
needle ./analysis_mutations/sgPemK_r.fasta ./analysis_mutations/seqs_enclosing_mutations.fasta -outfile ./analysis_mutations/needle_sgPemK_r.txt -gapopen 10.0 -gapextend 0.5
```
