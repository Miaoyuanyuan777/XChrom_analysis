#/bin/bash
python XChrom_train.py --input_folder ./train_data/ --cell_embedding_ad ../data/e18/processed_data/rna_pc32.h5ad --epochs 1000 --trackscore True --celltype 'pc32_leiden'