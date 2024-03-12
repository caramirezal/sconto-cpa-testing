## Based in Daria's code

# import packages
import sys
import os
import scanpy as sc
sys.path.append('/workspace')
from sconto_vae.module.ontobj import *
from sconto_vae.module.utils import *
from sconto_vae.model.scontovae_cpa import *
import neptune
from sklearn.decomposition import PCA
import umap
import pandas as pd

## settings
path2project = '/workspace/sconto-cpa-testing/'
path2figures = path2project + '/figures'


ontobj = Ontobj(description='CollecTRI_mouse_GRN')
#ontobj.initialize_dag(gene_annot='/workspace/sconto_vae/data/CollecTRItable_IFNb_MEF_ESC_overlapGenes_minRegsize10.tsv')
#ontobj.initialize_dag(gene_annot='/workspace/sconto-cpa-testing/data/grn/CollecTRItable_IFNb_MEF_ESC_overlapGenes_noFilter.tsv')
ontobj.initialize_dag(gene_annot='/workspace/sconto-cpa-testing/data/grn/CollecTRItable_IFNb_MEF_ESC_overlapGenes_minRegsize10.tsv')


ontobj.trim_dag(top_thresh=1000, 
                bottom_thresh=1)
ontobj.create_masks(top_thresh=1000,
                    bottom_thresh=1)
ontobj.save('/workspace/sconto-cpa-testing/data/CollecTRI_mouse_GRN_filt.ontobj')
ontobj = Ontobj()
#ontobj.load('/workspace/sconto_vae/data/CollecTRI_mouse_GRN_filt.ontobj')
ontobj.load('/workspace/sconto-cpa-testing/data/CollecTRI_mouse_GRN_filt.ontobj')


# load adata
adata = sc.read_h5ad('/workspace/sconto-cpa-testing/data/IFNb_MEF_ESC_all_conditions_top10000genes.h5ad')
adata.var_names = adata.var_names.astype('str')
adata.var_names_make_unique()
adata.var_names = [s.upper() for s in adata.var_names]
cpa_keys = ['stimulation_time', 'cell_type']

# setup the anndata
adata = setup_anndata_ontovae(adata, 
                              ontobj,
                              top_thresh=1000,
                              bottom_thresh=1,
                              cpa_keys=cpa_keys) 
# split into train-val and test


#train_adata, test_adata = split_adata(adata)
## Using the same train and test data from Anna
train = sc.read_h5ad('/workspace/sconto-cpa-testing/data/IFNb_MEF_ESC_train_preprocessed_top10000genes.h5ad')
test = sc.read_h5ad('/workspace/sconto-cpa-testing/data/IFNb_MEF_ESC_test_preprocessed_top10000genes.h5ad')
train_adata = adata[ adata.obs_names.isin(train.obs_names)]
test_adata = adata[ adata.obs_names.isin(test.obs_names)]
train_adata
test_adata


# create model
model = OntoVAEcpa(train_adata,
                   root_layer_latent=False,
                   activation_fn_dec=nn.ReLU,
                   neuronnum=1,
                   #use_batch_norm_dec=True,
                   #use_activation_dec=True,
                   #hidden_layers_enc=3,
                   #hidden_layers_class=2,
                   #neurons_per_class_layer=64
                  )



# init neptune run
run = neptune.init_run(
    project="scontocpa/scontocpa-testing",
    api_token="eyJhcGlfYWRkcmVzcyI6Imh0dHBzOi8vYXBwLm5lcHR1bmUuYWkiLCJhcGlfdXJsIjoiaHR0cHM6Ly9hcHAubmVwdHVuZS5haSIsImFwaV9rZXkiOiJmZmE2NjY0NS00YjUxLTQ4ZTktYmUzMi05YmUxMzAzY2UzNmMifQ=="
) 
#runid = 'SCON-15'
runid = str(vars(run)['_sys_id'])
modelpath = '/workspace/sconto-cpa-testing/data/models/' + runid
if not os.path.isdir(modelpath):
    os.mkdir(modelpath)

path2figures = path2project + "/figures" + '/' + runid
if not os.path.exists(path2figures + '/' + runid):
    os.makedirs(path2figures + '/' + runid)




'''
# train the model
model.train_model(modelpath,  
                   pos_weights = False,
                   #train_size=0.85,
                   #mixup_lambda=1,
                   #kl_coeff=0.001,
                   #lr_vae=1e-4,
                   #lr_adv=1e-3,
                   #adv_coeff=1000,
                   #pen_coeff=2,
                   #epochs=300,
                   #adv_step=1,
                   run=run)  
'''


model = OntoVAEcpa.load(test_adata, modelpath)
# compute different embeddings
embed = model.to_latent(test_adata)
views = list(embed.keys())
plt.rcParams.update({'font.size': 30})
for view in views:
    # run PCA and UMAP
    pca = PCA(n_components=30)
    pca_res = pca.fit_transform(embed[view])
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(pca_res)
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    # make scatterplot
    for c in range(len(cpa_keys)):
        # create color dict
        covar_categs = adata.obs[cpa_keys[c]].unique().tolist()
        palette = sns.color_palette(cc.glasbey, n_colors=len(covar_categs))
        color_dict = dict(zip(covar_categs, palette))
        # make scatter plot
        sns.scatterplot(x=embedding[:,0],
                        y=embedding[:,1], 
                        hue=test_adata.obs[cpa_keys[c]],
                        palette=color_dict,
                        legend='full',
                        s=15,
                        rasterized=True,
                        ax=ax.flatten()[c])
    plt.tight_layout()
    plt.savefig(path2figures + '/' + runid + '/umap_' + view + '.pdf')
    run["images/UMAP/" + view].upload(fig)
    plt.close()

run.stop()


tfa_test = model.get_pathway_activities(test_adata)
tfa_train = model.get_pathway_activities(train_adata)
tfa = model.get_pathway_activities(adata)
anns = adata.uns['_ontovae']['annot']



def save_results_to_tsv(path = '.',
                        tag = 'train_',
                        tfa_result = None):
    for st in tfa_result:
        print('Printing: '  + st)
        df = pd.DataFrame(tfa_result[st])
        df.to_csv(path + '/' + tag + st + '.tsv.gz', sep='\t')
        print('Printing to: ' + path + '/' + tag + st + '.tsv.gz')


## Saving results
save_results_to_tsv(path=modelpath, tag='train_', tfa_result=tfa_train)
save_results_to_tsv(path=modelpath, tag='test_', tfa_result=tfa_test)
save_results_to_tsv(path=modelpath, tag='all_', tfa_result=tfa)
anns.to_csv(modelpath + '/tfa_anns.tsv.gz', sep='\t')
test_adata.obs.to_csv(modelpath + '/test_cell_anns.tsv.gz', sep='\t')
train_adata.obs.to_csv(modelpath + '/train_cell_anns.tsv.gz', sep='\t')
adata.obs.to_csv(modelpath + '/all_cell_anns.tsv.gz', sep='\t')





##-----------------------------------------------------------------------------
## 
tfa_lin_layer_False_test = model.get_pathway_activities(test_adata,lin_layer=False)
save_results_to_tsv(path=modelpath, 
                    tag='test_lin_layer_False_', 
                    tfa_result=tfa_lin_layer_False_test)