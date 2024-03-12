## This script contains the implementation of the scOnto-cpa pipeline in mouse
## embryonic and fibroblasts cells stimulated with interferon

## Dependencies
import sys
import os
import scanpy as sc
import torch.nn as nn
import neptune
sys.path.append('/workspace/sconto-vae')
from sconto_vae.module.ontobj import *
from sconto_vae.module.utils import *
from sconto_vae.model.sconto_vae import *
from sconto_vae.model.scontovae_cpa import *
from sconto_vae.module.autotune import *
from ray import tune
from sconto_vae.module import autotune as at
from sconto_vae.model import scontovae_cpa as sv
import importlib
importlib.reload(at)
importlib.reload(sv)


## settings
path2project = '/workspace/sconto-cpa-testing/'
path2figures = path2project + '/figures'


## Loading data
scRNA_datapath = "/workspace/sconto-cpa-testing/data/IFNb_MEF_ESC_all_conditions_top10000genes.h5ad"
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



## subsetting the dataset
adata_subset = adata[adata.obs['stimulation_time'].isin(['0h', '1h'])].copy()


## Setting up 
adata_subset = setup_anndata_ontovae(adata_subset, 
                              ontobj,
                              top_thresh = 10000,
                              bottom_thresh = 10,
                              cpa_keys=['stimulation_time', 'cell_type'])




#train_adata, test_adata = split_adata(adata_subset)


ontobj = Ontobj()
ontobj.initialize_dag(gene_annot='/workspace/sconto-cpa-testing/data/grn/CollecTRItable_IFNb_MEF_ESC_overlapGenes_minRegsize10.tsv')

# trim the ontology
ontobj.trim_dag(top_thresh=10000, 
                bottom_thresh=10)


# create masks for decoder initialization
ontobj.create_masks(top_thresh=10000,
                 bottom_thresh=10)




model = OntoVAEcpa(adata_subset,
                   root_layer_latent=False,
                   use_batch_norm_dec=True,
                   use_activation_dec=True,
                   activation_fn_dec=nn.ReLU,
                   hidden_layers_enc=3,
                   hidden_layers_class=2,
                   neurons_per_class_layer=64,
                  )



run = neptune.init_run(
    project="scontocpa/scontocpa-testing",
    api_token="eyJhcGlfYWRkcmVzcyI6Imh0dHBzOi8vYXBwLm5lcHR1bmUuYWkiLCJhcGlfdXJsIjoiaHR0cHM6Ly9hcHAubmVwdHVuZS5haSIsImFwaV9rZXkiOiJmZmE2NjY0NS00YjUxLTQ4ZTktYmUzMi05YmUxMzAzY2UzNmMifQ==",
)  # your credentials
runid = str(vars(run)['_sys_id'])
modelpath = '/workspace/sconto-cpa-testing/data/models/' + runid
if not os.path.isdir(modelpath):
    os.mkdir(modelpath)



model.train_model(modelpath,  
                  train_size=0.85,
                  #mixup_lambda=1,
                  #kl_coeff=1e-2,
                  #lr_vae=1e-4,
                  #lr_adv=1e-3,
                  adv_coeff=100,
                  pen_coeff=2,
                  epochs=300,
                  adv_step=1,
                  run=run)


#
#path2model = '/workspace/sconto-cpa-testing/data/models/SCON-1'
#model = OntoVAEcpa.load(test_adata, path2model)
model = OntoVAEcpa.load(test_adata, modelpath)


embed = model.to_latent(test_adata)



if not os.path.exists(path2figures + '/' + runid):
    os.makedirs(path2figures + '/' + runid)


views = list(embed.keys())
for view in views:
    embedding = embed[view]
    covars = ['cell_type', 'stimulation_time']
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    # make scatterplot
    for c in range(len(covars)):
        # create color dict
        covar_categs = adata_subset.obs[covars[c]].unique().tolist()
        palette = sns.color_palette(cc.glasbey, n_colors=len(covar_categs))
        color_dict = dict(zip(covar_categs, palette))
        # make scatter plot
        sns.scatterplot(x=embedding[:,0],
                        y=embedding[:,1], 
                        hue=test_adata.obs[covars[c]],
                        palette=color_dict,
                        legend='full',
                        s=15,
                        rasterized=True,
                        ax=ax.flatten()[c])
    plt.tight_layout()
    plt.savefig(path2figures + '/' + runid + '/umap_' + view + '.pdf')
    run["images/UMAP/" + view].upload(fig)
    plt.close()





act = model.get_pathway_activities(test_adata)



####################################################################
## Model tuner
# create the tuner and show a summary
tuner = ModelTuner(OntoVAEcpa)
tuner.info()

# create the search space for as many parameters as desired from the tuner's summary
search_space = {
                #"mixup_lambda": tune.choice([0.0,2.0]), 
                #"kl_coeff": tune.loguniform(1e-4, 1e-2), 
                "lr_vae": tune.loguniform(1e-4, 1e-2) #,
                #"lr_adv": tune.loguniform(1e-4, 1e-2)
                }



cpa_keys = ['stimulation_time', 'cell_type']
# run the fit function of the tuner
results = tuner.fit(adata_subset, 
                    ontobj, 
                    search_space, 
                    epochs = 300, 
                    cpa_keys = cpa_keys, 
                    num_samples = 2, 
                    resources = {'gpu': 1})
# show the best hyperparameter settings
tuner.output(results[0], results[1])
