RSEM=~/bin/RSEM/bin/
nrep=1


if [ ! -d ../../profiles/single_cell/20/sim_${nrep}/fq/ ]
then
    mkdir -p ../../profiles/single_cell/20/sim_${nrep}/fq/
fi

nproc=0
for i in {1..40}
do

    ${RSEM}/rsem-simulate-reads \
    ~/doc/reference/rsem_pa_bowtie2/index \
    ../../files/SRR493366.highQ.model \
    ../../profiles/single_cell/20/sim_${nrep}/dropout/sample${i}.txt \
    0.05 5000000 \
    ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i} \
    --seed 0 && \
    gzip ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_1.fq \
    ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_2.fq && \
    mv ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_1.fq.gz ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_R1.fq.gz && \
    mv ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_2.fq.gz ../../profiles/single_cell/20/sim_${nrep}/fq/sample${i}_R2.fq.gz &

    nproc=$((${nproc}+1))
    if [ ${nproc} -ge 5 ]
    then
        wait
        nproc=0
    fi

done




