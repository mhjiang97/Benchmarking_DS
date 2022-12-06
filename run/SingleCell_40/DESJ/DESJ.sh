Rscript merge.sj.r --sj_dir ~/projects/AS/data/sim/rsem/revision2/sim_1/star \
--outdir ~/DESJoutput/ --cpu 10 \
--min_cell 1 --min_read 5

perl Junction.ann.pl ~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf \
~/DESJoutput/Alljunction.filter.list.xls ~/projects/AS/data/sim/rsem/revision2/sim_1/star ~/DESJoutput/

Rscript DESJ-detetion.r --junction_matrix=~/DESJoutput/merge.count.txt \
--col_ann=~/DESJ/cellanno.txt --min_pct=0.2 \
--cellinfo=~/DESJ/cellinfo.txt --row_ann=~/DESJoutput/Alljunction.filter.list.ann.onegene.xls \
--clus=A,B --outdir=~/DESJoutput/ --cpu=20