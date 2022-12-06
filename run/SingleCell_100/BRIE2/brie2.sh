briekit-event -a ~/doc/reference/gtf/gencode.v32.annotation.primary_assembly.protein_coding2.gtf -o ~/brie/briekit/AS_events

brie-count -a ~/brie/SE.gff3.gz -S ~/brie/cell_list.txt -o ~/brie/brieSE -p 20

brie-quant -i ~/brieSE/brie_count.h5ad -c ~/brie/SampleCondition.txt \
-o ~/brie/brieSE/brie_quant_cellA.h5ad \
--interceptMode gene --LRTindex 0 --batchSize 800000

brie-count -a ~/brie/MXE.gff3.gz -S ~/brie/cell_list.txt -o ~/brie/brieMXE -p 20

brie-quant -i ~/brieMXE/brie_count.h5ad -c ~/brie/SampleCondition.txt \
-o ~/brie/brieMXE/brie_quant_cellA.h5ad \
--interceptMode gene --LRTindex 0 --batchSize 800000