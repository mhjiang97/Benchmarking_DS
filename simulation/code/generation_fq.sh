## --------------------- INDEX --------------------- ##
##### generate STAR index using the primary assembly dna fasta file and the gtf file incorporating of only protein coding genes #####
## STAR 2.7.5a ##
conda activate sim
STAR --runThreadN 40 \
--runMode genomeGenerate \
--genomeDir ~/doc/reference/star_2.7.5a_pa \
--genomeFastaFiles ~/doc/reference/fa/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile ~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf
##### prepare RSEM index and only RSEM index #####
## RSEM 1.3.3 ##
rsem-prepare-reference \
--gtf ~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf \
~/doc/reference/fa/GRCh38.primary_assembly.genome.fa \
~/doc/reference/rsem_pa/index
##### prepare RSEM index with BOWTIE2 index #####
rsem-prepare-reference \
--bowtie2 \
--bowtie2-path ~/miniconda3/envs/assw_3.8/bin/ \
--gtf ~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf \
~/doc/reference/fa/GRCh38.primary_assembly.genome.fa \
~/doc/reference/rsem_pa_bowtie2/index
##### prepare Salmon index #####
cd ~/doc/reference/salmon_1.4.0_pa/
grep "^>" ~/doc/reference/fa/GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat ~/doc/reference/rsem_pa_bowtie2/index.idx.fa ~/doc/reference/fa/GRCh38.primary_assembly.genome.fa > gentrome.fa
gzip gentrome.fa
salmon index -t gentrome.fa.gz -d decoys.txt -p 40 -i index


##### download fastq files from sra #####
fastq-dump -Q 33 --split-e SRR493366


id=SRR493366
mkdir -p ~/projects/AS/data/real/rsem/pa/bowtie2/${id}/
rsem-calculate-expression \
-p 20 \
--paired-end --bowtie2 --seed 123 \
--bowtie2-path ~/miniconda3/envs/assw_3.8/bin/ \
~/projects/AS/data/real/${id}_1.fastq ~/projects/AS/data/real/${id}_2.fastq \
~/doc/reference/rsem_pa_bowtie2/index \
~/projects/AS/data/real/rsem/pa/bowtie2/${id}/${id} \
1> ~/projects/AS/data/real/rsem/pa/bowtie2/${id}/_log 2>&1


##### turn to R scripts for generating simulation files #####
Rscript simulation.R


##### simulate reads #####
RSEM=~/bin/RSEM/bin/

nproc=0
for nrep in {1..10}
do

if [ ! -d ../profiles/sim_${nrep}/fq/ ]
then
mkdir -p ../profiles/sim_${nrep}/fq/
fi

for i in {1..8}
do

nproc=$((${nproc}+1))
if [ ${nproc} -gt 5 ]
then
wait
nproc=0
fi

${RSEM}/rsem-simulate-reads \
~/doc/reference/rsem_pa_bowtie2/index \
../files/SRR493366.highQ.model \
../profiles/sim_${nrep}/sample${i}.txt \
0.05 40000000 \
../profiles/sim_${nrep}/fq/sample${i} \
--seed 0 && \
gzip ../profiles/sim_${nrep}/fq/sample${i}_1.fq \
../profiles/sim_${nrep}/fq/sample${i}_2.fq && \
mv ../profiles/sim_${nrep}/fq/sample${i}_1.fq.gz ../profiles/sim_${nrep}/fq/sample${i}_R1.fq.gz && \
mv ../profiles/sim_${nrep}/fq/sample${i}_2.fq.gz ../profiles/sim_${nrep}/fq/sample${i}_R2.fq.gz &

done
done




