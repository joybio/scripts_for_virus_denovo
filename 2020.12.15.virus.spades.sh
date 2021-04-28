#!/bin/bash
#coding:utf-8
#set -o pipefail
#set -e有一个例外情况，就是不适用于管道命令。只要一个子命令失败，整个管道命令就失败，脚本就会终止执行。

#####################################################################################
# Junbo Yang (yang_junbo_hi@126.com)
# Last modified: 7.15 12:34:47 2020
# Change log : #
# Change log : #
# Usage: sh script.sh path_to_sample

#####################################################################################
# Locations require to change (based to different computers)

# assign reference genome of the species
THREADS="20"
#REF="/media/ruizhi/database/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"
#HIAST_INDEX="/media/ruizhi/database/homo_sapiens/hisat2/GRCh38"
#BWA_INDEX="/media/ruizhi/database/homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa"
Bowtie_index="/media/ruizhi/database/ovis_canadensis/ovis"
#####################################################################################


#####################################################################################
#####################################################################################

# change illumina data format

#####################################################################################
mkdir sfastqc_results
mkdir -p genome_assembly/1.QC
# quality control
ls *_R1.fq.gz | while read id; do(mkdir -p sfastqc_results/$(basename $id '.fq.gz');fastqc -o sfastqc_results/$(basename $id '.fq.gz') $id);done
ls *_R2.fq.gz | while read id; do(mkdir -p sfastqc_results/$(basename $id '.fq.gz');fastqc -o sfastqc_results/$(basename $id '.fq.gz') $id);done
multiqc sfastqc_results/ -o genome_assembly/1.QC
#####################################################################################
#host genome
ls *_R1.fq.gz | while read id; do(i=${id%_*};bowtie2 -p 24 -x ${Bowtie_index} -1 $id -2 ${i}_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam);done
#ls */*_R1.fq.gz | while read id; do(i=${id%_*};bwa mem -t 2 -R '@RG\tID:lib_lane\tPL:Illumina\tLB:lib\tSM:L' -M -k 30 ${BWA_INDEX} $id ${i}_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam);done
#ls *_R1.fq.gz | while read id; do(hisat2 -p ${THREADS} --pen-noncansplice 1000000 -x ${HIAST_INDEX} -1 $id -2 $(basename $id '_R1.fq.gz')_R2.fq.gz -S $(basename $id '_R1.fq.gz').sam); done
#extract unmapped reads
ls *.sam | while read id;do(samtools view -bS -@ ${THREADS} -T ${REF} -h -u -f 12 -F 256 $id > $(basename $id 'sam')unmapped.bam);done
#sort by positon
ls *.unmapped.bam | while read id;do(samtools sort -@ ${THREADS} -n $id $(basename $id 'bam')sorted);done
#bam2fq
ls *.unmapped.sorted.bam | while read id;do(bamToFastq -i $id -fq $(basename $id '.sorted.bam').R1.fq -fq2 $(basename $id '.sorted.bam').R2.fq);done
#assembly
ls *R1.fq | while read id;do(mkdir -p genome_assembly/2.Assembly/$(basename $id ".unmapped.R1.fq");spades.py -k 21,33,55,77 --only-assembler --pe1-1 $id --pe1-2 $(basename $id "R1.fq")R2.fq -o genome_assembly/2.Assembly/$(basename $id ".unmapped.R1.fq"));done
#filter scaffolds (gt 500)
ls *R1.fq.gz | while read id;do(python /media/ruizhi/software/min_500.py -i genome_assembly/2.Assembly/$(basename $id "_R1.fq.gz")/scaffolds.fasta -o genome_assembly/2.Assembly/$(basename $id "_R1.fq.gz")/$(basename $id '.unmapped.R1.fq')contigs.500.fa);done
#ls *.unmapped.R1.fq | while read id;do(mkdir -p genome_assembly/2.Assembly/;megahit -1 $id -2 $(basename $id '.unmapped.R1.fq').unmapped.R2.fq -o genome_assembly/2.Assembly/$(basename $id '.unmapped.R1.fq') --out-prefix $(basename $id '.unmapped.R1.fq'));done

#filter scaffolds (gt 500)
#ls *R1.fq | while read id;do(python /media/ruizhi/software/min_500.py -i genome_assembly/2.Assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.fa -o genome_assembly/2.Assembly/$(basename $id ".unmapped.R1.fq")/$(basename $id '.unmapped.R1.fq').contigs.500.fa);done

ls *R1.fq | while read id;do(mkdir -p genome_assembly/3.Evaluation//$(basename $id '.unmapped.R1.fq');quast.py -o genome_assembly/3.Evaluation/$(basename $id '.unmapped.R1.fq')  -t 20 genome_assembly/2.Assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.500.fa);done

#orf prediction
#ls *R1.fq.gz | while read id;do(mkdir -p genome_assembly/orf_prediction/$(basename $id "_R1.fq.gz");prodigal -p meta -a genome_assembly/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").protein_seq.fasta -d genome_assembly/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").nucleotide_seq.fasta -o genome_assembly/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").gbk -s genome_assembly/orf_prediction/$(basename $id "_R1.fq.gz")/$(basename $id "_R1.fq.gz").poteintial.stat -i genome_assembly/spades/$(basename $id "_R1.fq.gz")/scaffolds.500.fasta);done
ls *R1.fq | while read id;do(mkdir -p genome_assembly/5.Annotation/;prokka --cpus 10 genome_assembly/2.Assembly/$(basename $id '.unmapped.R1.fq')/$(basename $id '.unmapped.R1.fq').contigs.500.fa --outdir genome_assembly/5.Annotation/$(basename $id '.unmapped.R1.fq')/ --prefix $(basename $id ".unmapped.R1.fq") --metagenome --kingdom Viruses);done
#non_redundunt_gene
mkdir -p genome_assembly/4.Non_refunfunt_geneset/
cat genome_assembly/5.Annotation/*/*.faa > genome_assembly/4.Non_refunfunt_geneset/total.redundunt.faa.fa
cat genome_assembly/5.Annotation/*/*.ffn > genome_assembly/4.Non_refunfunt_geneset/total.redundunt.ffn.fa
#cd-hit -M 8000 -i genome_assembly/non_redundunt/total.redundunt.faa -o genome_assembly/non_redundunt/non_redundunt.faa.fa
cd-hit -M 8000 -T 30 -i genome_assembly/4.Non_refunfunt_geneset/total.redundunt.faa.fa -o genome_assembly/4.Non_refunfunt_geneset/non_redundunt.nt.fa -c 0.9 -aS 0.9 -d 0
python /media/ruizhi/software/seq_format.py -i genome_assembly/4.Non_refunfunt_geneset/total.redundunt.ffn.fa -o genome_assembly/4.Non_refunfunt_geneset/total.redundunt.format.ffn.fa
python /media/ruizhi/software/seq_format.py -i genome_assembly/4.Non_refunfunt_geneset/non_redundunt.nt.fa -o genome_assembly/4.Non_refunfunt_geneset/non_redundunt.nt.format.fa

python /media/ruizhi/software/ffn.non.py -p genome_assembly/4.Non_refunfunt_geneset/non_redundunt.nt.format.fa -f genome_assembly/4.Non_refunfunt_geneset/total.redundunt.format.ffn.fa -o genome_assembly/4.Non_refunfunt_geneset/non_redundunt.fnn.format.fa
#############################################################################################################

#############################################################################################################
mkdir -p genome_assembly/6.Function_annotation/eggNOG
split -l 200000 -a 3 -d genome_assembly/4.Non_refunfunt_geneset/non_redundunt.nt.format.fa genome_assembly/6.Function_annotation/eggNOG/non_redundunt.nt.format.chunk_
ls genome_assembly/6.Function_annotation/eggNOG/*.chunk_0* | while read id;do(/media/ruizhi/software/eggnog-mapper/emapper.py -m diamond --no_annot --no_file_comments --cpu 10 -i $id -o $id);done
cat genome_assembly/6.Function_annotation/eggNOG/*.chunk_*.emapper.seed_orthologs > genome_assembly/6.Function_annotation/eggNOG/input_file.emapper.seed_orthologs
/media/ruizhi/software/eggnog-mapper/emapper.py --annotate_hits_table genome_assembly/6.Function_annotation/eggNOG/input_file.emapper.seed_orthologs --no_file_comments -o genome_assembly/6.Function_annotation/eggNOG/output_file --cpu 10
echo "query_name\tseed eggNOG ortholog\tseed ortholog evalue\tseed ortholog score\tPredicted taxonomic group\tPredicted protein name\tGene Ontology terms \tEC number\tKEGG_ko\tKEGG_Pathway\t KEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy \tBiGG Reaction\ttax_scope: eggNOG taxonomic level used for annotation\teggNOG OGs \tbestOG (deprecated, use smallest from eggnog OGs)\tCOG Functional Category\teggNOG free text description" > genome_assembly/6.Function_annotation/eggNOG/head.xls
cat genome_assembly/6.Function_annotation/eggNOG/head.xls genome_assembly/6.Function_annotation/eggNOG/output_file.emapper.annotations > genome_assembly/6.Function_annotation/eggNOG/output_file.emapper.annotations.xls




















