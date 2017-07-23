# Contents of folders

- ```figures/```: figures
-```tables/```: tables
- ```supplementary_tables/```: supplementary tables
- ```scripts/```: additional scripts used for analysis
- ```supplementary_material/```: Output files of the analyses

--------------------------------------

# Analyses presented in the BlobTools manuscript

# 0 Program versions
- ```BlobTools```: v1.0
- ```ART```: v2.5.8
- ```bbmap shuffle.sh```: v37.02
- ```CLC assembler```: v5.0.0.142510-160525-215357
- ```BWA mem```: v0.7.15-r1140
- ```BLASTn```: v2.6.0+
- ```fastaqual_select.pl```: v1.0 from [GitHub](https://github.com/sujaikumar/assemblage/archive/v1.0.tar.gz)
- ```Diamond```: v0.9.5
- ```BUSCO```: v2.0.1

--------------------------------------

# 1 Preparing simulated data

## 1.1 generate reads

```
art_illumina -ss HS25 --id CELEG --fcov 50 -l 150 -m 500 -s 10 -i CELEG.fna -o CELEG -M
art_illumina -ss HS25 --id ECOLI --fcov 25 -l 150 -m 500 -s 10 -i ECOLI.fna -o ECOLI -M
art_illumina -ss HS25 --id HSAP19 --fcov 10 -l 150 -m 500 -s 10 -i HS19.fna -o HSAP19 -M
art_illumina -ss HS25 --id HSMT --fcov 250 -l 150 -m 500 -s 10 -i HSMT.fna -o HSMT -M
art_illumina -ss HS25 --id PAERU --fcov 100 -l 150 -m 500 -s 10 -i PAERU.fna -o PAERU -M
art_illumina -ss HS25 --id CELEG --fcov 25 -l 150 -m 500 -s 10 -i CELEG.fna -o CELEG.25. -M
```

## 1.2 concatenate into libraries

```
cat CELEG1.fq ECOLI1.fq HSAP191.fq HSMT1.fq > blobtools.dataset_A.1.fq
cat CELEG2.fq ECOLI2.fq HSAP192.fq HSMT2.fq > blobtools.dataset_A.2.fq
cat PAERU1.fq CELEG.25.1.fq > blobtools.dataset_B.1.fq
cat PAERU2.fq CELEG.25.2.fq > blobtools.dataset_B.2.fq
```

## 1.3 shuffle datasets

```
shuffle.sh in=blobtools.dataset_A.1.fq in2=blobtools.dataset_A.2.fq out=blobtools.dataset_A.1.shuffled.fq out2=blobtools.dataset_A.2.shuffled.fq
shuffle.sh in=blobtools.dataset_B.1.fq in2=blobtools.dataset_B.2.fq out=blobtools.dataset_B.1.shuffled.fq out2=blobtools.dataset_B.2.shuffled.fq
```

## 1.4 concatenate into one library (for mapping purposes)

```
cat blobtools.dataset_A.1.shuffled.fq blobtools.dataset_B.1.shuffled.fq > blobtools.dataset_both.1.shuffled.fq
cat blobtools.dataset_A.2.shuffled.fq blobtools.dataset_B.2.shuffled.fq > blobtools.dataset_both.2.shuffled.fq
```
--------------------------------------

# 2 Simulated read datasets by taxon

## 2.1 Assembly of simulated reads by taxonomic group

### 2.1.1 CLC

```
clc_assembler -o assembly.CELEG-SIM.fasta -p fb ss 300 700 -q -i CELEG.25.1.fq CELEG.25.2.fq -p fb ss 300 700 -q -i CELEG1.fq CELEG2.fq
clc_assembler -o assembly.ECOLI-SIM.fasta -p fb ss 300 700 -q -i ECOLI1.fq ECOLI2.fq
clc_assembler -o assembly.HSAPI-SIM.fasta -p fb ss 300 700 -q -i HSAP191.fq HSAP192.fq -p fb ss 300 700 -q -i HSMT1.fq HSMT2.fq
clc_assembler -o assembly.PAERU-SIM.fasta -p fb ss 300 700 -q -i PAERU1.fq PAERU2.fq
```

### 2.1.2 Rename sequences in assemblies

```
perl -i -pe "s/^>/>CELEG./g" assembly.CELEG-SIM.fasta
perl -i -pe "s/^>/>ECOLI./g" assembly.ECOLI-SIM.fasta
perl -i -pe "s/^>/>HSAPI./g" assembly.HSAPI-SIM.fasta
perl -i -pe "s/^>/>PAERU./g" assembly.PAERU-SIM.fasta
```

> ```supplementary_data/1_assembly_sim/assembly.CELEG-SIM.fasta```

> ```supplementary_data/1_assembly_sim/assembly.ECOLI-SIM.fasta```

> ```supplementary_data/1_assembly_sim/assembly.HSAPI-SIM.fasta```

> ```supplementary_data/1_assembly_sim/assembly.PAERU-SIM.fasta```

### 2.1.3 Concatenate into one file (for mapping purposes)

```
cat assembly.*.fasta > assembly.sim.all.fasta
```

## 2.2 Map read libraries

### 2.2.1 BWA

```
bwa index assembly.sim.all.fasta
bwa mem assembly.sim.all.fasta blobtools.dataset_both.1.shuffled.fq blobtools.dataset_both.2.shuffled.fq | samtools view -b - > blobtools.dataset_both.vs.assembly.sim.all.bam
```

### 2.2.2 Generate read counts by taxon for each sequence

```
samtools view -F 2304 blobtools.dataset_both.vs.assembly.sim.all.bam | cut -f1,3 | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | sed 's/HS19/HSAPI/g' | sed 's/HSMT/HSAPI/g' | sed 's/ENA|AE004091|AE004091/PAERU/g' | perl -lane 'if ($F[0] eq "*"){ print $F[0]."\t".(split /\./, $F[1])[0] }else{ print $F[0]."\t".(split /\./, $F[1])[0]}' | sort -Vk1 | uniq -c > blobtools.dataset_both.vs.assembly.sim.all.bam.read_count_by_reference.txt
```

> ```supplementary_data/1_assembly_sim/blobtools.dataset_both.vs.assembly.sim.all.bam.read_count_by_reference.txt```

### 2.2.3 Infer true taxonomy based on read counts using the script generate_table_based_on_read_counts_by_sequence.py (CHECK)

```
scripts/generate_table_based_on_read_counts_by_sequence.py -i blobtools.dataset_both.vs.assembly.sim.all.bam.read_count_by_reference.txt > assembly.sim.all.table_based_on_read_counts.txt
```

> ```supplementary_data/1_assembly_sim/assembly.sim.all.table_based_on_read_counts.txt```

--------------------------------------

# 3 Simulated read libraries for BlobTools

## 3.1 CLC assembly of both simulated read libraries

```
clc_assembler -o blobtools.assembly.A_B.fasta -p fb ss 300 700 -q -i blobtools.dataset_A.1.shuffled.fq blobtools.dataset_A.2.shuffled.fq -p fb ss 300 700 -q -i blobtools.dataset_B.1.shuffled.fq blobtools.dataset_B.2.shuffled.fq
```

> ```supplementary_data/2_simulated_libraries/blobtools.assembly.A_B.fasta```

## 3.2 Mapping of read libraries

### 3.2.1 BWA
- generates a BAM file for each read library mapped against the assembly of both simulated libraries

```
bwa index blobtools.assembly.A_B.fasta
bwa mem blobtools.assembly.A_B.fasta blobtools.dataset_A.1.shuffled.fq blobtools.dataset_A.2.shuffled.fq | samtools view -bS - > blobtools.dataset_A.vs.blobtools.assembly.A_B.bam
bwa mem blobtools.assembly.A_B.fasta blobtools.dataset_B.1.shuffled.fq blobtools.dataset_B.2.shuffled.fq | samtools view -bS - > blobtools.dataset_B.vs.blobtools.assembly.A_B.bam
```

### 3.2.2 Convert BAM to COV format using BlobTools ```map2cov```
- generates files containing coverage information in COV format which are used in construction of BlobDBs

```
parallel -j 2 'blobtools map2cov -i blobtools.assembly.A_B.fasta -b {}' ::: *.blobtools.assembly.A_B.bam
```

> ```supplementary_data/2_simulated_libraries/blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.cov```

> ```supplementary_data/2_simulated_libraries/blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.cov```

--------------------------------------

# 4 Assessment of efficiency of BlobTools taxonomic annotation based on similarity search results

## 4.1 Extracting base/read coverage information from BAM files

### 4.1.1 Generate read counts by taxon for each sequence
- used for assessing performance of taxonomic annotation of blobtools based on different similarity search results

#### 4.1.1.1 Generate list of sequence IDs by taxon from which reads originated, for each read library

```
samtools view -F 2304 blobtools.dataset_A.vs.blobtools.assembly.A_B.bam | cut -f1,3 | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | sed 's/HS19/HSAPI/g' | sed 's/HSMT/HSAPI/g' | sed 's/ENA|AE004091|AE004091/PAERU/g' | perl -lane 'if ($F[0] eq "*"){ print $F[0]."\t".(split /\./, $F[1])[0] }else{ print $F[0]."\t".(split /\./, $F[1])[0]}'  > blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.read_ids_by_contig_id.txt
samtools view -F 2304 blobtools.dataset_B.vs.blobtools.assembly.A_B.bam | cut -f1,3 | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | sed 's/HS19/HSAPI/g' | sed 's/HSMT/HSAPI/g' | sed 's/ENA|AE004091|AE004091/PAERU/g' | perl -lane 'if ($F[0] eq "*"){ print $F[0]."\t".(split /\./, $F[1])[0] }else{ print $F[0]."\t".(split /\./, $F[1])[0]}' > blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.read_ids_by_contig_id.txt
```

#### 4.1.1.2 join both files from previous step and get read counts

- generates file containing : counts, sequence_id ("*" for unmapped), and true origin of reads that mapped

```
cat *read_ids_by_contig_id.txt | sort -Vk1 | uniq -c > blobtools.dataset_A_B.vs.blobtools.assembly.A_B.bam.read_count_by_reference.txt
```

> ```supplementary_data/3_assessment_taxonomic_annotation/blobtools.dataset_A_B.vs.blobtools.assembly.A_B.bam.read_count_by_reference.txt```

## 4.2 Infer true taxonomy based on read counts using the script infer_true_taxonomy_based_on_read_mapping.py

```
generate_table_based_on_read_counts_by_sequence.py -i blobtools.dataset_A_B.vs.blobtools.assembly.A_B.bam.read_count_by_reference.txt > blobtools.assembly.A_B.true_taxonomy_by_contig.txt
```

> ```supplementary_data/3_assessment_taxonomic_annotation/blobtools.assembly.A_B.true_taxonomy_by_contig.txt```


## 4.3 Generating similarity search results
- ```MTS1``` : ```[-]-max_target_seqs 1``` (Diamond blastx, blastn)
- ```MTS10``` : ```[-]-max_target_seqs 10``` (Diamond blastx, blastn)
- ```HSP1``` : ```[-]-max_hsps 1``` (Diamond blastx, blastn)
- ```CUL10``` : ```-culling_limit 10``` (blastn)
- ```no-mask``` : search is performed against database without removing any sequences
- ```mask``` : search is performed against database removing sequences
 - ```supplementary_data/3_assessment_taxonomic_annotation/taxids_to_exlude.txt``` : file containing NCBI subtree TaxIDs for the following NCBI taxids: 9604 (Hominidae), 561 (Escherichia), 6239 (Caenorhabditis elegans), 286 (Pseudomonas), 28384 (other sequences). Subtree TaxIDs were retrieved through NCBI taxonomy web interface. Sequences belonging to these NCBI subtree TaxIDs were removed to generate masked UniProt Reference Proteomes Diamond database (uniprot_ref_proteomes.masked.diamond-v0.9.5) as specified in [MISC-section]

### 4.3.1 BLASTn searches

#### 4.3.1.1 ```no-mask```, ```CUL10```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.cul10.out -culling_limit 10
```

#### 4.3.1.2 ```no-mask```, ```MTS1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.mts1.out -max_target_seqs 1
```

#### 4.3.1.3 ```no-mask```, ```MTS10```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.mts10.out -max_target_seqs 10
```

#### 4.3.1.4 ```no-mask```, ```CUL10```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.cul10.max_hsp_1.out -culling_limit 10 -max_hsps 1
```

#### 4.3.1.5 ```no-mask```, ```MTS1```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.mts1.max_hsp_1.out -max_target_seqs 1 -max_hsps 1
```

#### 4.3.1.6 ```no-mask```, ```MTS10```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.no_mask.mts10.max_hsp_1.out -max_target_seqs 10 -max_hsps 1
```

#### 4.3.1.7 ```mask```, ```CUL10```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.cul10.out -culling_limit 10 -negative_gilist gis_to_exclude.txt
```

#### 4.3.1.8 ```mask```, ```MTS1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.mts1.out -max_target_seqs 1 -negative_gilist gis_to_exclude.txt
```

#### 4.3.1.9 ```mask```, ```MTS10```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.mts10.out -max_target_seqs 10 -negative_gilist gis_to_exclude.txt
```

#### 4.3.1.10 ```mask```, ```CUL10```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.cul10.max_hsp_1.out -culling_limit 10 -negative_gilist gis_to_exclude.txt -max_hsps 1
```

#### 4.3.1.11 ```mask```, ```MTS1```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.mts1.max_hsp_1.out -max_target_seqs 1 -negative_gilist gis_to_exclude.txt -max_hsps 1
```

#### 4.3.1.12 ```mask```, ```MTS10```, ```HSP1```

```
blastn -query blobtools.assembly.A_B.fasta -db ncbi.2017-06-13/nt -evalue 1e-25 -outfmt '6 qseqid staxids bitscore std' -out A_B.vs.nt.mask.mts10.max_hsp_1.out -max_target_seqs 10 -negative_gilist gis_to_exclude.txt -max_hsps 1
```

### 4.3.2 Diamond blastx searches

#### 4.3.2.1 ```no-mask```, ```MTS1```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.no_mask.mts1.out --max-target-seqs 1
```

#### 4.3.2.2 ```no-mask```, ```MTS10```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.no_mask.mts10.out --max-target-seqs 10
```

#### 4.3.2.3 ```no-mask```, ```MTS10```, ```HSP1```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.no_mask.mts10.max_hsp_1.out --max-target-seqs 10 --max-hsps 1
```

#### 4.3.2.4 ```mask```, ```MTS1```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.masked.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.mask.mts1.out --max-target-seqs 1
```

#### 4.3.2.5 ```mask```, ```MTS10```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.masked.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.mask.mts10.out --max-target-seqs 10
```

#### 4.3.2.6 ```mask```, ```MTS10```, ```HSP1```

```
diamond blastx --query 3_assembly/blobtools.assembly.A_B.fasta --db Reference_Proteomes_2017_07/uniprot_ref_proteomes.masked.diamond-v0.9.5.dmnd --sensitive --evalue 1e-25 --outfmt 6 --out A_B.vs.refprot.mask.mts10.max_hsp_1.out --max-target-seqs 10 --max-hsps 1
```

### 4.3.3 Add TaxIDs to Diamond blastx searches

#### 4.3.3.1 BlobTools taxify

```
parallel -j1 'blobtools taxify -f {} -m Reference_Proteomes_2017_07/uniprot_ref_proteomes.taxids -s 0 -t 2' ::: A_B.vs.refprot*.out
```

#### 4.3.3.2 Remove un-taxified Diamond blastx searches

```
rm A_B.vs.refprot.no_mask.mts1.out
rm A_B.vs.refprot.no_mask.mts10.out
rm A_B.vs.refprot.no_mask.mts10.max_hsp_1.out
rm A_B.vs.refprot.mask.mts1.out
rm A_B.vs.refprot.mask.mts10.out
rm A_B.vs.refprot.mask.mts10.max_hsp_1.out
```

> ```supplementary_data/3_assessment_taxonomic_annotation/search_results.tar.gz```

## 4.4 Generate BlobDBs using BlobTools ```create```

### 4.4.1 Using a single similarity search results

```
parallel -j1 'blobtools create -i blobtools.assembly.A_B.fasta -c blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.cov -c blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.cov -t {} -o {/.}' ::: *.out
```

### 4.4.2 Using two similarity search results

#### 4.4.2.1 ```no-mask```

```
parallel -j1 'blobtools create -i blobtools.assembly.A_B.fasta -c blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.cov -c blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.cov -x bestsumorder -t {1} -t {2} -o {1/.}.AND.{2/.}  ::: A_B.vs.nt.no_mask.* ::: A_B.vs.refprot.no_mask.*
```

#### 4.4.2.2 ```mask```

```
parallel -j1 'blobtools create -i blobtools.assembly.A_B.fasta -c blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.cov -c blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.cov -x bestsumorder -t {1} -t {2} -o {1/.}.AND.{2/.}  ::: A_B.vs.nt.mask.* ::: A_B.vs.refprot.mask.*
```

## 4.5 Generate tabular views of BlobDBs using BlobTools ```view```

### 4.5.1 Generate tabular views of single similarity search result BlobDBs

```
parallel -j1 'blobtools view -i {} -r superkingdom -r phylum -r order --hits' ::: `ls | grep "json" | grep -v 'AND'
```

### 4.5.2 Generate tabular views of two similarity search result BlobDBs, using taxrule 'bestsumorder'

```
parallel -j1 'blobtools view -i {} -x bestsumorder -r superkingdom -r phylum -r order --hits' ::: *AND*.json
```

> ```supplementary_data/3_assessment_taxonomic_annotation/tabular_views.tar.gz```

# 4.6 Evaluate results of taxonomic annotation of BlobTools
- ```supplementary_data/3_assessment_taxonomic_annotation/order_of_tables.txt```: contains filenames of tabular views of BlobDBs paired with search parameters

```
generate_taxonomy_tables.py -t blobtools.assembly.A_B.true_taxonomy_by_contig.txt -d 3_assessment_taxonomic_annotation/ -b order_of_tables.txt --taxrank order
```

> ```supplementary_data/3_assessment_taxonomic_annotation/blobtools_table_analysis/```

# 5 Visualising assembly of simulated read libraries using BlobTools
- BlobDB used is : ```supplementary_data/2_simulated_libraries/1_prefilter/A_B.vs.nt.mask.mts1.max_hsp_1.AND.A_B.vs.refprot.mask.mts1.taxified.blobDB.json```

## 5.1 Generating Blobplots and Readcovplots using 'BlobTools plot'

```
blobtools plot -i A_B.vs.nt.mask.mts1.max_hsp_1.AND.A_B.vs.refprot.mask.mts1.taxified.blobDB.json -x bestsumorder -r order --format png -o blobplot_png/```
```

## 5.2 Generating a Covplot using 'BlobTools covplot'

```
blobtools covplot -i A_B.vs.nt.mask.mts1.max_hsp_1.AND.A_B.vs.refprot.mask.mts1.taxified.blobDB.json --lib cov0 -c blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.cov --xlabel 'Library A' --ylabel 'Library B' --max 1e4 -x bestsumorder -r order
```

# 6 Filter sequence IDs in assembly based on tabular view of BlobDB using awk (CHECK)

## 6.1 "Rhabditida" sequences

```awk '$5>1 && $6>1' blobtools.A_B.blobDB.bestsumorder.table.txt | cut -f1 > rhadbditida.contig_ids.txt```

## 6.2 "Primates" sequences

```awk '$5>0.1 && $6<0.5 && $8!="Bacteria" && $12 !="Nematoda"' blobtools.A_B.blobDB.bestsumorder.table.txt | cut -f1 > primates.contig_ids.txt```

## 6.3 "Enterobacterales" sequences

```awk '$5>20 && $5<40 && $6<0.5 && $8!="Eukaryota"' blobtools.A_B.blobDB.bestsumorder.table.txt | cut -f1 > enterobacterales.contig_ids.txt```

## 6.4 "Pseudomonadales" sequences

```awk '$5<0.2 && $6>0.5' blobtools.A_B.blobDB.bestsumorder.table.txt | cut -f1 > pseudomonadales.contig_ids.txt```

# 7 Filter reads based on lists of sequence IDs using 'BlobTools bamfilter'

## 7.1 "Rhabditida" reads

### 7.1.1 "Rhabditida" reads in library A

```blobtools bamfilter -b blobtools.dataset_A.vs.blobtools.assembly.A_B.bam -i rhadbditida.contig_ids.txt -o rhadbditida_A```

### 7.1.2 "Rhabditida" reads in library B

```blobtools bamfilter -b blobtools.dataset_B.vs.blobtools.assembly.A_B.bam -i rhadbditida.contig_ids.txt -o rhadbditida_B```

## 7.2 "Primates" reads

```blobtools bamfilter -b blobtools.dataset_A.vs.blobtools.assembly.A_B.bam -i primates.contig_ids.txt -o primates```

## 7.3 "Enterobacterales" reads

```blobtools bamfilter -b blobtools.dataset_B.vs.blobtools.assembly.A_B.bam -i enterobacterales.contig_ids.txt -o enterobacterales```

## 7.4 "Pseudomondales" reads

```blobtools bamfilter -b blobtools.dataset_B.vs.blobtools.assembly.A_B.bam -i pseudomonadales.contig_ids.txt -o pseudomonadales```

# 8. Assembly of filtered reads by taxonomic group

## 8.1 CLC
- using only using read pairs where both reads mapped to sequences in lists ('InIn')

### 8.1.1 "Rhabditida" reads

```clc_assembler -o blobtools.assembly.rhadbditida-BT.fasta -p fb ss 300 700 -q rhadbditida_A.blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.InIn.fq -p fb ss 300 700 -q rhadbditida_B.blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.InIn.fq```

### 8.1.2 "Primates" reads

```clc_assembler -o blobtools.assembly.primates-BT.fasta -p fb ss 300 700 -q primates.blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.InIn.fq```

### 8.1.3 "Enterobacterales" reads

```clc_assembler -o blobtools.assembly.enterobacterales-BT.fasta -p fb ss 300 700 -q enterobacterales.blobtools.dataset_A.vs.blobtools.assembly.A_B.bam.InIn.fq```

### 8.1.4 "Pseudomonadales" reads

```clc_assembler -o blobtools.assembly.pseudomonadales-BT.fasta -p fb ss 300 700 -q pseudomonadales.blobtools.dataset_B.vs.blobtools.assembly.A_B.bam.InIn.fq```

## 8.2 Rename sequences in assemblies

```
perl -i -pe "s/^>/>rhabditida./g" blobtools.assembly.rhadbditida-BT.fasta
perl -i -pe "s/^>/>primates./g" blobtools.assembly.primates-BT.fasta
perl -i -pe "s/^>/>enterobacterales./g" blobtools.assembly.enterobacterales-BT.fasta
perl -i -pe "s/^>/>pseudomonadales./g" blobtools.assembly.pseudomonadales-BT.fasta
```

> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.assembly.rhadbditida-BT.fasta```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.assembly.primates-BT.fasta```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.assembly.enterobacterales-BT.fasta```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.assembly.pseudomonadales-BT.fasta```

## 8.3 Concatenate into one file (for mapping purposes)

```cat blobtools.assembly.*.fasta > blobtools.assembly.final.fasta```

# 9 Evaluation of 'cleaned' assemblies

## 9.1 Taxonomic evaluation based on mapping of original reads

### 9.1.1 BWA
```
bwa index blobtools.assembly.all.fasta
bwa mem blobtools.assembly.all.fasta blobtools.dataset_both.1.shuffled.fq blobtools.dataset_both.2.shuffled.fq | samtools view -b - > blobtools.dataset_both.vs.blobtools.assembly.all.bam
```
### 9.1.2 Generate read counts by taxon for each sequence
```
samtools view -F 2304 blobtools.dataset_both.vs.blobtools.assembly.all.bam | cut -f1,3 | awk ' { t = $1; $1 = $2; $2 = t; print; } ' | sed 's/HS19/HSAPI/g' | sed 's/HSMT/HSAPI/g' | sed 's/ENA|AE004091|AE004091/PAERU/g' | perl -lane 'if ($F[0] eq "*"){ print $F[0]."\t".(split /\./, $F[1])[0] }else{ print $F[0]."\t".(split /\./, $F[1])[0]}' | sort -Vk1 | uniq -c > blobtools.dataset_both.vs.blobtools.assembly.all.bam.read_ids_by_contig_id.txt
```
### 9.1.3 Generate table of taxonomic annotation based on read counts
```
generate_table_based_on_read_counts_by_sequence.py -i blobtools.dataset_both.vs.blobtools.assembly.all.bam.read_ids_by_contig_id.txt  > blobtools.assembly.all.table_based_on_read_counts.txt
```

> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.assembly.all.table_based_on_read_counts.txt```

### 9.1.4 Generate hits files based on table of taxonomic annotation based on read counts

```grep '^CELEG' blobtools.assembly.all.table_based_on_read_counts.txt | cut -f1,6 | perl -lane 'if($F[1] eq "CELEG"){print $F[0]."\t6239\t100"}elsif($F[1] eq "HSAPI"){print $F[0]."\t9606\t100"}elsif($F[1] eq "PAERU"){print $F[0]."\t287\t100"}elsif($F[1] eq "ECOLI"){print $F[0]."\t562\t100"}' > blobtools.assembly.rhadbditida.hits.txt```

```grep '^HSAPI' blobtools.assembly.all.table_based_on_read_counts.txt | cut -f1,6 | perl -lane 'if($F[1] eq "CELEG"){print $F[0]."\t6239\t100"}elsif($F[1] eq "HSAPI"){print $F[0]."\t9606\t100"}elsif($F[1] eq "PAERU"){print $F[0]."\t287\t100"}elsif($F[1] eq "ECOLI"){print $F[0]."\t562\t100"}' > blobtools.assembly.primates.hits.txt```

```grep '^ECOLI' blobtools.assembly.all.table_based_on_read_counts.txt | cut -f1,6 | perl -lane 'if($F[1] eq "CELEG"){print $F[0]."\t6239\t100"}elsif($F[1] eq "HSAPI"){print $F[0]."\t9606\t100"}elsif($F[1] eq "PAERU"){print $F[0]."\t287\t100"}elsif($F[1] eq "ECOLI"){print $F[0]."\t562\t100"}' > blobtools.assembly.enterobacterales.hits.txt```

```grep '^PAERU' blobtools.assembly.all.table_based_on_read_counts.txt | cut -f1,6 | perl -lane 'if($F[1] eq "CELEG"){print $F[0]."\t6239\t100"}elsif($F[1] eq "HSAPI"){print $F[0]."\t9606\t100"}elsif($F[1] eq "PAERU"){print $F[0]."\t287\t100"}elsif($F[1] eq "ECOLI"){print $F[0]."\t562\t100"}' > blobtools.assembly.pseudomonadales.hits.txt```

## 9.2 Generate individual blobplots for each 'cleaned' assembly using taxonomic annotation based on read counts

### 9.2.1 Convert BAM to COV format using 'BlobTools map2cov'

```blobtools map2cov -i blobtools.assembly.all.fasta -b blobtools.dataset_both.vs.blobtools.assembly.all.bam```

### 9.2.2 Subset COV file by taxonomic group

```
grep -Pv 'enterobacterales|primates|pseudomonadales' blobtools.dataset_both.vs.blobtools.assembly.all.bam.cov > blobtools.dataset_both.vs.blobtools.assembly.all.bam.rhabditida.cov
grep -Pv 'primates|pseudomonadales|rhabditida' blobtools.dataset_both.vs.blobtools.assembly.all.bam.cov > blobtools.dataset_both.vs.blobtools.assembly.all.bam.enterobacterales.cov
grep -Pv 'enterobacterales|pseudomonadales|rhabditida' blobtools.dataset_both.vs.blobtools.assembly.all.bam.cov > blobtools.dataset_both.vs.blobtools.assembly.all.bam.primates.cov
grep -Pv 'enterobacterales|primates|rhabditida' blobtools.dataset_both.vs.blobtools.assembly.all.bam.cov > blobtools.dataset_both.vs.blobtools.assembly.all.bam.pseudomonadales.cov
```

> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.dataset_both.vs.blobtools.assembly.all.bam.rhabditida.cov```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.dataset_both.vs.blobtools.assembly.all.bam.enterobacterales.cov```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.dataset_both.vs.blobtools.assembly.all.bam.primates.cov```
> ```supplementary_data/2_simulated_libraries/2_postfilter/blobtools.dataset_both.vs.blobtools.assembly.all.bam.pseudomonadales.cov```

## 9.3 Create BlobDBs by taxonomic group
- using assembly
- using COV file
- using taxonomic annotation based on read mapping

```
blobtools create -i blobtools.assembly.rhadbditida-BT.fasta -c blobtools.dataset_both.vs.blobtools.assembly.all.bam.CELEG.cov -t blobtools.assembly.rhadbditida.hits.txt -o blobtools.assembly.rhadbditida
```
```
blobtools create -i blobtools.assembly.primates-BT.fasta -c blobtools.dataset_both.vs.blobtools.assembly.all.bam.HSAPI.cov -t blobtools.assembly.primates.hits.txt -o blobtools.assembly.primates
```
```
blobtools create -i blobtools.assembly.enterobacterales-BT.fasta -c blobtools.dataset_both.vs.blobtools.assembly.all.bam.ECOLI.cov -t blobtools.assembly.enterobacterales.hits.txt -o blobtools.assembly.enterobacterales
```
```
blobtools create -i blobtools.assembly.pseudomonadales-BT.fasta -c blobtools.dataset_both.vs.blobtools.assembly.all.bam.PAERU.cov -t blobtools.assembly.pseudomonadales.hits.txt -o blobtools.assembly.pseudomonadales
```
## 9.4 Make BlobPlots for each BlobDB using defined colours
```
blobtools plot -i blobtools.assembly.rhadbditida.blobDB.json -o blobplots_png/ -r order --colours blobtools_colours.txt ; \
blobtools plot -i blobtools.assembly.primates.blobDB.json -o blobplots_png/ -r order --colours blobtools_colours.txt ; \
blobtools plot -i blobtools.assembly.enterobacterales.blobDB.json -o blobplots_png/ -r order --colours blobtools_colours.txt ; \
blobtools plot -i blobtools.assembly.pseudomonadales.blobDB.json -o blobplots_png/ -r order --colours blobtools_colours.txt
```

## 9.5 BUSCO analysis

### 9.5.1 BUSCO analysis of assemblies of filtered reads by taxonomic group

```
BUSCO.py -i blobtools.assembly.rhadbditida-BT.fasta -o blobtools.assembly.rhadbditida -m genome -l nematoda_odb9/ ; \
BUSCO.py -i blobtools.assembly.primates-BT.fasta -o blobtools.assembly.primates -m genome -l mammalia_odb9/ ; \
BUSCO.py -i blobtools.assembly.enterobacterales-BT.fasta -o blobtools.assembly.enterobacterales -m genome -l enterobacteriales_odb9/ ; \
BUSCO.py -i blobtools.assembly.pseudomonadales-BT.fasta -o blobtools.assembly.pseudomonadales -m genome -l gammaproteobacteria_odb9/
```

### 9.5.2 BUSCO analysis of assembly of original reads by taxonomic group

```
BUSCO.py -i CELEG.fasta -o CELEG_SIM -m genome -l nematoda_odb9/ ; \
BUSCO.py -i ECOLI.fasta -o ECOLI_SIM -m genome -l enterobacteriales_odb9/ ; \
BUSCO.py -i PAERU.fasta -o PAERU_SIM -m genome -l gammaproteobacteria_odb9 ; \
BUSCO.py -i HSAPI.fasta -o HSAPI_SIM -m genome -l mammalia_odb9/
```

### 9.5.3 BUSCO analysis of reference assemblies

```
BUSCO.py -i assembly.CELEG-SIM.fasta -o CELEG_REF -m genome -l nematoda_odb9/ ; \
BUSCO.py -i assembly.ECOLI-SIM.fasta -o ECOLI_REF -m genome -l enterobacteriales_odb9/ ; \
BUSCO.py -i assembly.PAERU-SIM.fasta -o PAERU_REF -m genome -l gammaproteobacteria_odb9 ; \
BUSCO.py -i assembly.HSAPI-SIM.fasta -o HSAPI_REF -m genome -l mammalia_odb9/
```

-----------------------------------------------------------------------------------------------------------------------

# X MISC

## X.1 preparation of Diamond databases

### X.1.1 Download UniProt Reference Proteomes

```wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2017_07.tar.gz```

### X.1.2 Unpack protein FASTAs for each kingdom

```parallel -j8 'gunzip {}' ::: `ls | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional'```

### X.1.3 Concatenate all protein sequences into ```uniprot_ref_proteomes.fasta```

```cat */*.fasta > uniprot_ref_proteomes.fasta```

### X.1.4 Change sequence IDs

```cat uniprot_ref_proteomes.fasta | sed -r 's/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d"|" > temp; mv temp uniprot_ref_proteomes.fasta```

### X.1.5 make "no-mask" database

```diamond makedb --in uniprot_ref_proteomes.fasta -d uniprot_ref_proteomes.diamond-v0.9.5```

### X.1.6 "mask" database

#### X.1.6.1 Subset mapping IDs to only contain TaxID entries

```cat */.idmapping | grep "NCBI_TaxID" > uniprot_ref_proteomes.taxids```

#### X.1.6.2 Get sequence IDs to exclude

```colgrep -f uniprot_ref_proteomes.taxids -i taxids_to_exlude.txt -c 3 | cut -f1 > sequence_ids_to_exclude.txt```

#### X.1.6.3 Exclude sequences based on list to exclude

```fastaqual_select.pl -f uniprot_ref_proteomes.fasta -e sequence_ids_to_exclude.txt > uniprot_ref_proteomes.masked.fasta```

#### X.1.6.4 make "mask" database

```diamond makedb --in uniprot_ref_proteomes.masked.fasta -d uniprot_ref_proteomes.masked.diamond-v0.9.5```
