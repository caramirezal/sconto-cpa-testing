## Visualization of the TFA 


## Loading libraries
## Dependencies
library(dplyr)
library(readr)
library(viridis)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ComplexHeatmap)
library(umap)


test_lin_layer = FALSE

set.seed(333)


path2project <- '/home/bq_cramirez/sconto-cpa-testing'

## Source settings
source(paste0(path2project, "/scripts/R/settings.R"))


#path2figures <- paste0(path2project, '/figures/01_tfa_visualization')
#path2figures <- paste0(path2project, '/figures/01_tfa_visualization_two_conditions')
path2figures <- paste0(path2project, '/figures/01_tfa_visualization_three_conditions')
if ( ! dir.exists(path2figures)){
    dir.create(path2figures)
}

if ( test_lin_layer == FALSE ) {
    pattern = 'test_lin_layer_False_z_'
} else {
    pattern = 'test_z_'
}


## Two conditions
path2cobra_results <- paste0(path2project, '/data/models/SCON-20/')

## Three conditions
#path2cobra_results <- paste0(path2project, '/data/models/SCON-15/')

tfa.files <- list.files(path2cobra_results,
                        pattern = pattern)
tfa.files.full_names <- paste0(path2cobra_results, tfa.files)

if ( test_lin_layer == FALSE ) {
    names <- gsub('test_lin_layer_False_z_|.tsv.gz', '', tfa.files)
} else {
    names <- gsub('test_z_|.tsv.gz', '', tfa.files)
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




#--------------------------------------------------------------------------------
## @Heatmap

## Mean TFA across conditions and cell lines
tfa_mean.df.list <- lapply(
    tfa.df.list,
    function(df){
        df %>%
            group_by(cell_type, stimulation_time) %>%
            summarise_all(mean) %>%
            ungroup
    }
)


view <- "basal"


max_value <- sapply(tfa_mean.df.list, 
       function(df) select(df, -stimulation_time, -cell_type) %>% max ) %>% max
min_value <- 0
col_fun = circlize::colorRamp2(c(min_value, (max_value/2), max_value), viridis(3))


plot_heatmap <- function(view=view){ 
        df <- tfa_mean.df.list[view][[1]]
        mtx <- select(df, -cell_type, -stimulation_time) %>%
                    t() %>%
                    as.matrix()
        anns_heatmap <- select(df, cell_type, stimulation_time)
        head(anns_heatmap)

        column_anns <- columnAnnotation(cell_type=anns_heatmap$cell_type,
                                        stimulation_time=anns_heatmap$stimulation_time,
                                        col = list(cell_type=cell_type_colors,
                                            stimulation_time=stimulation_time_colors))


        ## Selecting TFs by variability
        tfa_var <- apply(mtx, 1, var)
        n_top <- 50
        tfa_top <- sort(tfa_var, decreasing=TRUE) %>% 
                  head(n_top)   %>%
                  names()

        path2tfs <- paste0(path2project, 
                           '/analysis/top_var_tfs_',
                           'three_conditions')
        if ( ! dir.exists(path2tfs) ) {
                dir.create(path2tfs)
        }
        writeLines(tfa_top,
                   paste0(path2tfs, 
                         '/top_var_tfs_',
                         view,
                         '.txt'))

        h1 <- Heatmap(
                mtx,
                top_annotation = column_anns,
                col=col_fun,
                column_title = view
                )


        h2 <- Heatmap(
                      mtx[tfa_top, ],
                      top_annotation = column_anns,
                      col=col_fun,
                      column_title = view,
                      )

        return(list(heatmap_all_tfs=h1,
                    heatmap_top_tfs=h2))
}

names(names) <- names
heatmaps_list <- lapply(names, plot_heatmap)



pdf(paste0(path2figures, '/heatmap_mean_tfa_lin_layer_False_z.pdf'),
    height=65, width=20)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmaps_list[names[1]][[1]][[1]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heatmaps_list[names[2]][[1]][[1]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(heatmaps_list[names[3]][[1]][[1]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(heatmaps_list[names[4]][[1]][[1]], newpage = FALSE)
upViewport()

dev.off()





pdf(paste0(path2figures, '/heatmap_mean_tfa_top_50_lin_layer_False_z.pdf'),
    height=10, width=20)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 4)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmaps_list[names[1]][[1]][[2]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heatmaps_list[names[2]][[1]][[2]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(heatmaps_list[names[3]][[1]][[2]], newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(heatmaps_list[names[4]][[1]][[2]], newpage = FALSE)
upViewport()

dev.off()




##------------------------------------------------------------------------
## UMAP projections
umap.list <- lapply(tfa.df.list, 
                    function(df) {
                        mtx <- select(df, -cell_type, -stimulation_time) %>%
                                      as.matrix 
                        umap(mtx, n_neighbors = 300, min_dist = 0.1) 
                                  }
                        )

umap.df.list <- lapply(names, 
                       function(name){
                            res <- umap.list[name][[1]]
                            df <- res$layout %>% as.data.frame()
                            colnames(df) <- paste0('UMAP', 1:2)
                            df <- cbind(df, 
                                        select(tfa.df.list[name][[1]],
                                               -`cell_type`,
                                               -`stimulation_time`), 
                                               cell_anns)
                            return(df)

})

lapply(umap.df.list, colnames)

view <- 'basal'
plot_view <- function(view){
      umap.cell_type <- umap.df.list[view][[1]] %>%
        ggplot(aes(x=UMAP1, y=UMAP2,
                   colour=cell_type)) +
            geom_point(size=0.1) +
            theme_classic() +
            scale_colour_manual(values=cell_type_colors) +
            labs(colour='',
                 title=view)

      umap.stimulation_time <- umap.df.list[view][[1]] %>%
        ggplot(aes(x=UMAP1, y=UMAP2,
                   colour=stimulation_time)) +
            geom_point(size=0.1) +
            theme_classic() +
            scale_colour_manual(values=stimulation_time_colors) +
            labs(colour='')

       return(list(cell_type=umap.cell_type, 
                   stimulation_time=umap.stimulation_time))
}

umap_plot_list <- lapply(names, plot_view)

view <- 'basal'
pdf(paste0(path2figures, '/umap_views.pdf'),
    width=12, height=6)
grid.arrange(umap_plot_list$basal$cell_type,
             umap_plot_list$basal$stimulation_time,
             umap_plot_list$cell_type$cell_type,
             umap_plot_list$cell_type$stimulation_time,
             umap_plot_list$stimulation_time$cell_type,
             umap_plot_list$stimulation_time$stimulation_time,
             umap_plot_list$total$cell_type,
             umap_plot_list$total$stimulation_time,
             ncol=4
)
dev.off()

#tfa_stimulation_time <- umap.tfa_mean.df.list$"stimulation_time"
#selected_tfs <- tfa_stimulation_time %>%
#                    select(AHR:ZNF91) %>%
#                        apply(2, var) %>%
#                        sort(decreasing=TRUE) %>%
#                        head(50) %>%
#                        names()




#pdf(paste0(path2figures, '/umap_tfa_top_tfs_lin_layer_False_z.pdf'),
#     width=20, height=5)
#umap_list <- lapply(selected_tfs,
#       function(tf){
#           umap.df.list$"stimulation_time" %>%
#           ggplot(aes(x=UMAP1, y=UMAP2,
#                    colour=!!sym(tf))) +
#             geom_point(size=0.1) +
#            theme_void() +
#             scale_colour_viridis() +
#            labs(colour='',
#                 title=tf) +
#                 theme(legend.position="none")
#})
#gridExtra::grid.arrange(grobs=umap_list, ncol=10)
#dev.off()



