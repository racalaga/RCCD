# RNA CITE Concatenate Decoder
## introduction 
RNA CITE seq Concatenate Decoder aggreagte RNA-seq and CITE-seq sequance to reproduce CITE seq sequance
## installlation 
RCCD is python scrpits
```
git clone https://github.com/wheaton5/souporcell.git
```
It dependendent follw python packages:
```
-scanpy
-tensorflow
```
if use conda env, install package use conda:
```
conda install -c conda-forge scanpy
conda install -c conda-forge tensorflow
```

## To run the script

RCCD script use two csv file **config.csv** and **match.csv**

**config.csv**
```
anndata_path,rnacite.h5ad ## input file path
marker_path,match.csv ## match.csv file path
threshold_quality,30 ## train set
rae_epochs,10 ## RNA autoencoder epochs
pae_epochs,10 ## protein autoencoder epochs
gen_epochs,10 ## RCCD epochs
save_path,produce.h5ad ## output file path
```
**match.csv**
it composed RNA-id and protein-id match 
```
CD86,CD86-1
CD274,CD274-1
TNFRSF14,CD270
PVR,CD155
NECTIN2,CD112
CD47,CD47-1
CD48,CD48-1
CD40,CD40-1
CD40LG,CD154
CD52,CD52-1
CD3D,CD3
CD8A,CD8
NCAM1,CD56
CD19,CD19-1
CD33,CD33-1
ITGAX,CD11c
HLA-A,HLA-ABC
...
```
**Run scrpits**
```
python RCCD.py config.csv
```

## input and outputs

input file can be any type **anndata** data <https://anndata.readthedocs.io/en/latest/#>
```
AnnData object with n_obs × n_vars = 104236 × 36776
    var: 'gene_ids', 'feature_types', 'genome'
```

output file is saved **save_path** in **config.csv**
```
AnnData object with n_obs × n_vars = 104236 × 28003
    var: 'gene_ids', 'feature_types', 'genome', 'n_counts'
    uns: 'marker_quality_after', 'marker_quality_before', 'original_protein', 'produced_protein', 'protein_marker', 'rna', 'rna_marker'

RCCD save results in anndata.uns :

- marker_quality_before, marker_quality_after : cell by quality metric 
- original_protein : orignal protein counts
- produced_protein : Decoded protein counts
- rna : orignal rna counts
- protein_marker, rna_marker: markers pair definded by match.csv

```
![rae_train](./img/rae_train.png)
![pae_train](./img/pae_train.png)
![gen_train](./img/gen_train.png)
![cos_sim](./img/cos_sim.png)
![euclid_distance](./img/euclid_distance.png)