#/bin/bash

python scbasset_preprocess.py --ad_file ../m_brain/0_data/processed_data/atac_ad.h5ad --input_fasta /picb/bigdata/project/miaoyuanyuan/mm10.fa
python scbasset_train.py --input_folder processed/ --epochs 1000