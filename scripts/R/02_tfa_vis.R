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




path2project <- '/home/bq_cramirez/sconto-cpa-testing'

## Source settings
source(paste0(path2project, "/scripts/R/settings.R"))


path2figures <- paste0(path2project, '/figures/01_tfa_visualization')
if ( ! dir.exists(path2figures)){
    dir.create(path2figures)
}

path2cobra_results <- paste0(path2project, '/data/models/SCON-15/')
tfa.files <- list.files(path2cobra_results,
                        pattern = 'test_z')
tfa.files.full_names <- paste0(path2cobra_results, tfa.files)


names <- gsub('test_z_|.tsv.gz', '', tfa.files)


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

        h1 <- Heatmap(
                mtx,
                top_annotation = column_anns,
                col=viridis(100, option="magma"),
                column_title = view
                )


        h2 <- Heatmap(
                      mtx[tfa_top, ],
                      top_annotation = column_anns,
                      col=viridis(100, option="magma"),
                      column_title = view
                      )

        return(list(heatmap_all_tfs=h1,
                    heatmap_top_tfs=h2))
}


heatmaps_list <- lapply(names, plot_heatmap)


pdf(paste0(path2figures, '/heatmap_mean_tfa.pdf'),
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





pdf(paste0(path2figures, '/heatmap_mean_tfa_top_50.pdf'),
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



