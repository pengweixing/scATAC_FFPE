
## 1. install miniconda3.9

## 2. replace python3.7 in cellranger with python3.9
```shell
rm ./cellranger-atac-2.0.0/external/anaconda/bin/python
rm ./cellranger-atac-2.0.0/external/anaconda/bin/python3
ln -s miniconda3.9/python ./cellranger-atac-2.0.0/external/anaconda/bin/python
ln -s miniconda3.9/python ./cellranger-atac-2.0.0/external/anaconda/bin/python3
```
## 2. install the following packages under miniconda3.9
```
    lz4
    pyfaidx, pysam
    bedtools
    pybedtools
    biopython
    scipy
    matplotlib
    hdf5
    samtools
    pyBigWig
    sklearn
    h5py
    tables
    MOODS
    statsmodels
    mkl-service
    mkl
    scikit-learn
    tsne
    umap
```
And then update all packages

```conda update --all```

## 3. modify the 10x scripts
For the script of 
`./cellranger-atac-2.0.0/lib/python/cellranger/analysis/bhtsne.py`

Do

```Annotate: Line 78-82```



For the script of 
`./cellranger-atac-2.0.0/mro/atac/stages/processing/cell_calling/remove_barcode_multiplets/__init__.py`

Do

```
Add new lines to Line:
Line261:     part_b_count = 100
Line262:     part_a_count = 100

Line85:     part_a_count = 100
Line86:     part_b_count = 100

Line147,152,163,175,180,190,: change it to  "self_signal = part_a_linkage_matrix[1,1,1,1]"

annotate: Line309-311, Line314-316, Line320-322, Line325-327
```

Purpose: The matrix generated in this step is too large to 
load into memory. So we limited it to 100*100

