
## Dependencies
library(ComplexHeatmap)
library(dplyr)
library(viridis)
library(ggplot2)
library(readr)


set.seed(333)


path2project <- '/home/bq_cramirez/sconto-cpa-testing'

## Source settings
source(paste0(path2project, "/scripts/R/settings.R"))


path2figures <- paste0(path2project, '/figures/03_comparing_cobra_vs_vega')
if ( ! dir.exists(path2figures)) {
    dir.create(path2figures)
}

## Reading VEGA data
vega_latent <- read_tsv(
    paste0(path2project, 
           '/data/vega/latent_space.tsv.gz')
)
dim(vega_latent)
vega_latent[1:5, 1:5]
vega_latent <- vega_latent[,2:ncol(vega_latent)]
## cell annotations
cell_anns <- read_tsv(
    paste0(path2project, 
           '/data/vega/cell_anns.tsv.gz')
)
dim(cell_anns)
cell_anns[1:5, 1:5]
## TF annotations
tf_anns <- read_tsv(
    paste0(path2project, 
           '/data/vega/tf_anns.tsv.gz')
)
dim(tf_anns)
tf_anns[1:5, 1:2]


vega_latent <- as.matrix(vega_latent)
colnames(vega_latent) <- tf_anns$`0`
rownames(vega_latent) <- cell_anns$`...1`


selected_tfs_vega <- apply(vega_latent, 2, var) %>%
                            sort(decreasing=TRUE) %>%
                            head(50) %>%
                            names()
selected_tfs_vega


column_anns_vega <- columnAnnotation(
    cell_type=cell_anns$cell_type,
    stimulation_time=cell_anns$stimulation_time,
    col = list(cell_type=cell_type_colors,
               stimulation_time=stimulation_time_colors),
    show_legend = FALSE)



#max_value <- sapply(tfa_mean.df.list, 
#       function(df) select(df, -stimulation_time, -cell_type) %>% max ) %>% max
#min_value <- 0
#col_fun = circlize::colorRamp2(c(min_value, (max_value/2), max_value), viridis(3))

pdf(paste0(path2figures, '/heatmap_vega.pdf'),
    height=12, width=10)
Heatmap(
    t(vega_latent[ , selected_tfs_vega]),
    top_annotation = column_anns_vega,
    col=viridis(30),
    column_title = 'VEGA',
    show_column_names = FALSE
                      )
dev.off()






##-------------------------------------------------------------------------------
## Loading COBRA results


test_lin_layer = FALSE


if ( test_lin_layer == FALSE ) {
    pattern = 'test_lin_layer_False_z_'
} else {
    pattern = 'test_z_'
}


## Three conditions
path2cobra_results <- paste0(path2project, '/data/models/SCON-27/')

tfa.files <- list.files(path2cobra_results,
                        pattern = pattern)
tfa.files <- tfa.files[!grepl('umap', tfa.files)]
tfa.files.full_names <- paste0(path2cobra_results, tfa.files)

if ( test_lin_layer == FALSE ) {
    names <- gsub('test_lin_layer_False_z_|.tsv.gz', '', tfa.files)
    names <- names[!grepl('umap', names)]
} else {
    names <- gsub('test_z_|.tsv.gz', '', tfa.files)
    names <- names[!grepl('umap', names)]
}



## Reading TF annotations
anns_file_name <- list.files(path2cobra_results, pattern = 'tfa_anns.tsv.gz')
anns <- read_tsv(paste0(path2cobra_results, anns_file_name))

## Reading cell annotations
cell_anns_file_name <- list.files(path2cobra_results, pattern = 'test_cell_anns.tsv.gz')
cell_anns <- read_tsv(paste0(path2cobra_results, cell_anns_file_name))


## Reading results from COBRA
tfa.df.list <- lapply(tfa.files.full_names, read_tsv)
tfa.df.list  <- lapply(tfa.df.list, function(df) select(df, -`...1`))     ## Removing artificially added extra column

## Adding cell and tf names
tfa.df.list  <- lapply(tfa.df.list,
    function(df){
        df <- as.data.frame(df)
        rownames(df) <- cell_anns$`...1`
        colnames(df) <- anns$Name
        df$"cell_type" <- cell_anns$"cell_type"
        df$"stimulation_time" <- cell_anns$"stimulation_time"
        return(df)
    }
)
names(tfa.df.list) <- names



tfa_stimulation_cobra <- tfa.df.list$"stimulation_time"
mtx <- select(tfa_stimulation_cobra, -cell_type, -stimulation_time)
mtx <- as.matrix(mtx)


selected_tfs_cobra <- apply(mtx, 2, var) %>%
                            sort(decreasing=TRUE) %>%
                            head(50) %>%
                            names()
selected_tfs_cobra


column_anns_cobra <- columnAnnotation(
    cell_type=cell_anns$cell_type,
    stimulation_time=cell_anns$stimulation_time,
    col = list(cell_type=cell_type_colors,
               stimulation_time=stimulation_time_colors),
    show_annotation_name = FALSE)


tfs <- selected_tfs_vega[selected_tfs_vega %in% colnames(mtx)]
pdf(paste0(path2figures, '/heatmap_cobra.pdf'),
    height=12, width=10)
Heatmap(
        t(scale(mtx[, tfs])),
        top_annotation = column_anns_cobra,
        column_title = "COBRA",
        col=viridis(30),
       show_column_names = FALSE
)
dev.off()





##--------------------------------------------------------------
## 

selected_tfs_vega <- apply(vega_latent, 2, var) %>%
                            sort(decreasing=TRUE) %>%
                            head(25) %>%
                            names()
selected_tfs_vega

selected_tfs_cobra <- apply(mtx, 2, var) %>%
                            sort(decreasing=TRUE) %>%
                            head(25) %>%
                            names()
selected_tfs_cobra


selected_join_tfs <- c(selected_tfs_vega, 
                       selected_tfs_cobra) 
selected_join_tfs <- selected_join_tfs[selected_join_tfs!="UNANNOTATED_0"]
selected_join_tfs <- unique(selected_join_tfs)

heatmap_cobra <- Heatmap(
        t(mtx[, selected_join_tfs]),
        top_annotation = column_anns_cobra,
        column_title = "COBRA",
    col=viridis(30),
    show_column_names = FALSE
)

heatmap_vega <- Heatmap(
    t(scale(vega_latent[ , selected_join_tfs])),
    top_annotation = column_anns_vega,
    col=viridis(30),
    column_title = 'VEGA',
    show_column_names = FALSE
                      )

heatmap <- heatmap_cobra + heatmap_vega

pdf(paste0(path2figures, '/heatmaps_tfs_cobra_vs_vega.pdf'),
    height=12, width=10)
heatmap_cobra + heatmap_vega
dev.off()




##-----------------------------------------------------------------
mtx.list <- lapply(names, 
     function(name){
        df <- tfa.df.list[name][[1]]
        mtx <- select(df, -cell_type, -stimulation_time)
       as.matrix(mtx)
     }
)
names(mtx.list) <- names 

heatmap_cobra.list <- lapply(names,
      function(name){
           mtx <- mtx.list[name][[1]]
           Heatmap(
               t(mtx[, selected_join_tfs]),
               top_annotation = column_anns_cobra,
               column_title = paste("COBRA - ", name),
               col=viridis(30),
               show_column_names = FALSE
        )
      }
)
names(heatmap_cobra.list) <- names



pdf(paste0(path2figures, '/heatmaps_tfs_cobra_vs_vega.pdf'),
    height=12, width=10)
heatmap_cobra.list$total +
    heatmap_cobra.list$stimulation_time + 
    heatmap_cobra.list$cell_type +
    heatmap_cobra.list$basal + heatmap_vega
dev.off()