## This script contains the implementation of the scOnto-cpa pipeline in mouse
## embryonic and fibroblasts cells stimulated with interferon

## Dependencies
import sys
import os
import scanpy as sc
sys.path.append('/workspace/sconto-vae')
from sconto_vae.module.ontobj import *
from sconto_vae.module.utils import *
from sconto_vae.model.sconto_vae import *

## settings
path2project = '/workspace/sconto-cpa-validation/'
path2figures = path2project + '/figures'


## Loading data
scRNA_datapath = "/workspace/sconto-cpa-validation/data/IFNb_MEF_ESC_all_conditions_top10000genes.h5ad"
adata = sc.read_h5ad(scRNA_datapath)

path2figures = path2project + '/figures/01_exploratory_analysis'
if not os.path.exists(path2figures):
    os.makedirs(path2figures)



## plot cell line
with plt.rc_context():  
     sc.pl.umap(adata, color='cell_type', show=False)
     plt.savefig(path2figures + "/umap_PBMC_MEF_ESC_cell_type.pdf", bbox_inches="tight")
     plt.clf()



with plt.rc_context():  
     sc.pl.umap(adata, color='stimulation_time', show=False)
     plt.savefig(path2figures + "/umap_PBMC_MEF_ESC_stimulation.pdf", bbox_inches="tight")
     plt.clf()





ontobj = Ontobj()
ontobj.load('/workspace/sconto_vae/data/GO/GO.ontobj')
adata = sc.read_h5ad('/workspace/sconto_vae/data/SCZ_Puvogel.h5ad')



