#!/usr/bin/env python
import sys
sys.path.append('.')
import os
import h5py
import anndata
import configargparse
import numpy as np
import pandas as pd
import scanpy as sc
from utils import split_train_test_val, make_h5_sparse, make_bed_seqs_from_df, dna_1hot_2vec  
from scipy import sparse


def make_parser():
    parser = configargparse.ArgParser(
        description="Preprocess anndata to generate inputs for scBasset.")
    parser.add_argument('--ad_file', type=str,
                       help='Input scATAC anndata. .var must have chr, start, end columns. anndata.X must be in csr format.')
    parser.add_argument('--input_fasta', type=str,
                       help='Genome fasta file.')
    parser.add_argument('--out_path', type=str, default='./processed',
                       help='Output path. Default to ./processed/')

    return parser


def main():
    parser = make_parser()  ## 创建了一个命令行参数解析器，并配置了解析器的选项和参数。
    args = parser.parse_args()  ## 用于解析实际的命令行参数，并将解析结果存储在args变量中。
    
    input_ad = args.ad_file
    input_fasta = args.input_fasta
    output_path = args.out_path
    
    ad = anndata.read_h5ad(input_ad)

    # sample cells
    data_path = '.'
    os.makedirs(output_path, exist_ok=True)  ## 创建目录output_path，如果该目录不存在的话。True如果目录已经存在，不会引发异常，程序会继续执行。如果设置为False（默认值），并且目录已经存在，那么会引发 FileExistsError 异常。
    seq_len = 1344

    # save anndata
    ad.write('%s/ad.h5ad'%output_path)
    print('successful writing h5ad file.')  ## 将ad对象的内容保存到新文件中，而不会覆盖输入文件 input_ad

    # save peak bed file
    ad.var.loc[:,['chr','start','end']].to_csv('%s/peaks.bed'%output_path, sep='\t', header=False, index=False)
    print('successful writing bed file.')

    # save train, test, val splits
    train_ids, test_ids, val_ids = split_train_test_val(np.arange(ad.shape[1]))
    f = h5py.File('%s/splits.h5'%output_path, "w")
    f.create_dataset("train_ids", data=train_ids)
    f.create_dataset("test_ids", data=test_ids)
    f.create_dataset("val_ids", data=val_ids)
    f.close()
    print('successful writing split file.')

    # save labels (ad.X)
    m = ad.X.tocoo().transpose().tocsr()
    ## tocoo()：这是一个用于稀疏矩阵（sparse matrix）的方法，它将 ad.X 转换为 COO 格式（坐标格式）。COO 格式是一种表示稀疏矩阵的方式，其中只存储非零元素的坐标和值
    ## transpose()：这是一个用于矩阵的方法，它将 COO 格式的稀疏矩阵转置，即将行和列互换。
    ## tocsr()：这是另一个用于稀疏矩阵的方法，它将 COO 格式的稀疏矩阵转换为 CSR 格式（压缩稀疏行格式）。CSR 格式是一种用于高效存储和操作稀疏矩阵的格式。
    m_train = m[train_ids,:]
    m_val = m[val_ids,:]
    m_test = m[test_ids,:]
    sparse.save_npz('%s/m_train.npz'%output_path, m_train, compressed=False)
    sparse.save_npz('%s/m_val.npz'%output_path, m_val, compressed=False)
    sparse.save_npz('%s/m_test.npz'%output_path, m_test, compressed=False)
    print('successful writing sparse m.')
    ## SciPy库的函数，用于保存稀疏矩阵，保存为NumPy的稀疏矩阵文件 (*.npz格式)，compressed=False以确保不进行压缩

    # save sequence h5 file
    ad_train = ad[:,train_ids]
    ad_test = ad[:,test_ids]
    ad_val = ad[:,val_ids]
    make_h5_sparse(ad, '%s/all_seqs.h5'%output_path, input_fasta)
    make_h5_sparse(ad_train, '%s/train_seqs.h5'%output_path, input_fasta)
    make_h5_sparse(ad_test, '%s/test_seqs.h5'%output_path, input_fasta)
    make_h5_sparse(ad_val, '%s/val_seqs.h5'%output_path, input_fasta)

if __name__ == "__main__":
    main()
## 通过将 main() 函数封装在 if __name__ == "__main__": 条件下，您可以确保只有在直接运行脚本时才执行 main() 函数，而在被导入到其他脚本时不会执行。