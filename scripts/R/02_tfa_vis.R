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

path2project <- '/home/bq_cramirez/sconto-cpa-testing'

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
        return(df)
    }
)
names(tfa.df.list) <- names


names(names) <- names
cell_anns <- as.data.frame(cell_anns)
rownames(cell_anns) <- cell_anns$`...1`
seurat_list <- lapply(names,
    function(name){
           CreateSeuratObject(counts=t(tfa.df.list[name][[1]]),
                              assay='TFA',
                              project=name,
                              min.cells=0,
                              min.features=0,
                              meta.data = cell_anns)
    }
)
names(seurat_list)


seurat_list


##---------------------------------------------------------------------------

## Differential TFA in MEFs 1 vs 0 hrs
test<-0.01
#test<-.5
dtfa_list.mef1vs0 <- lapply(
    seurat_list,
    function(seurat){
        s <- subset(seurat,
                    cell_type == 'MEF' &
                     stimulation_time %in% c('0h', '1h'))
        Idents(s) <- s$stimulation_time            
        FindMarkers(object=s,
                    slot='counts',
                    groups = 'stimulation_time',
                    ident.1 = '1h',
                    logfc.threshold = test,
                    min.pct = 0,                    
                    min.diff.pct = -Inf)
    }
)

add_annotations <- function(dtfa_list){
    lapply(dtfa_list,
       function(df){
             df$"regulon_size" <- plyr::mapvalues(x=rownames(df),
                                                  from=anns$Name,
                                                  to=anns$genes)
            df$"gene" <- rownames(df)
            return(df)
    })
}
dtfa_list.mef1vs0 <- add_annotations(dtfa_list.mef1vs0)


thresholds.mef1vs0 <- list(c(-0.075, 0.065), 
                   c(-0.2, 0.2),
                   c(-1, 0.5),
                   c(-1, 1.0))
limits.mef1vs0 <- list(c(-0.15, 0.15),
               c(-0.5, 1.5),
               c(-4, 2),
               c(-5.5, 4.5))
names(thresholds.mef1vs0) <- names
names(limits.mef1vs0) <- names

create_vulcano <- function(dtfa_list,
                thresholds=thresholds,
                limits=limits){
    lapply(names,
    function(name){
        df <- dtfa_list[name][[1]] %>%
            filter(p_val>0) %>%
            mutate(highlight=ifelse(avg_log2FC < thresholds[name][[1]][1] |
                                      thresholds[name][[1]][2] < avg_log2FC,
                                         TRUE, FALSE )) %>%
            mutate(label=ifelse(highlight==TRUE, gene, "")) %>%
            mutate(regulon_size=as.numeric(regulon_size)) %>%
            arrange(desc(regulon_size))                           
        df %>%    ggplot(aes(x=avg_log2FC,
                        y=-log10(p_val),
                        colour=regulon_size,
                        label=label)) +
                        geom_point() +
                        geom_text_repel(colour="red", max.overlaps=1000) +
                        theme_bw() +
                        labs(title=name) +
                        scale_colour_viridis() +
                        xlim(limits[name][[1]])
    })                
}

vulcano.list.mef1vs0 <- create_vulcano(dtfa_list.mef1vs0, 
                                        limits=limits.mef1vs0,
                                        thresholds=thresholds.mef1vs0)

pdf(file=paste0(path2figures, '/vulcano_plots_MEFS_1vs0.pdf'),
    height=3.5, width=18)
grid.arrange(grobs=vulcano.list.mef1vs0, ncol=4)
dev.off()




##---------------------------------------------------------------------------
## Differential TFA in MEFs 6 vs 0 hrs
dtfa_list.mef6vs0 <- lapply(
    seurat_list,
    function(seurat){
        s <- subset(seurat,
                    cell_type == 'MEF' &
                     stimulation_time %in% c('0h', '6h'))
        Idents(s) <- s$stimulation_time            
        FindMarkers(object=s,
                    slot='counts',
                    groups = 'stimulation_time',
                    ident.1 = '6h',
                    logfc.threshold = test,
                    min.pct = 0,                    
                    min.diff.pct = -Inf)
    }
)

## MEFs
thresholds.mef6vs0 <- list(c(-0.25, 0.25), 
                   c(-1.5, 1.5),
                   c(-0.45, 0.45),
                   c(-0.8, 0.5))
limits.mef6vs0 <- list(c(-0.5, 0.5),
               c(-3, 3),
               c(-1, 2),
               c(-3, 1))
names(thresholds.mef6vs0) <- names
names(limits.mef6vs0) <- names
dtfa_list.mef6vs0 <- add_annotations(dtfa_list.mef6vs0)
vulcano.list.mef6vs0 <- create_vulcano(dtfa_list.mef6vs0, 
                                        limits=limits.mef6vs0,
                                        thresholds=thresholds.mef6vs0)

pdf(file=paste0(path2figures, '/vulcano_plots_MEFS_6vs0.pdf'),
    height=3.5, width=18)
grid.arrange(grobs=vulcano.list.mef6vs0, ncol=4)
dev.off()




##---------------------------------------------------------------------------

## Differential TFA in ESCs 1 vs 0 hrs
dtfa_list.esc1vs0 <- lapply(
    seurat_list,
    function(seurat){
        s <- subset(seurat,
                    cell_type == 'ESC' &
                     stimulation_time %in% c('0h', '1h'))
        Idents(s) <- s$stimulation_time            
        FindMarkers(object=s,
                    slot='counts',
                    groups = 'stimulation_time',
                    ident.1 = '1h',
                    logfc.threshold = test,
                    min.pct = 0,                    
                    min.diff.pct = -Inf)
    }
)



dtfa_list.esc1vs0 <- add_annotations(dtfa_list.esc1vs0)

thresholds.esc1vs0 <- list(c(-0.05, 0.05), 
                   c(-0.2, 0.2),
                   c(-1, 1),
                   c(-1.5, 2.5))
limits.esc1vs0 <- list(c(-0.6, 0.6),
               c(-1, 1.5),
               c(-4, 2.5),
               c(-6, 5))
names(thresholds.esc1vs0) <- names
names(limits.esc1vs0) <- names

vulcano.list.esc1vs0 <- create_vulcano(dtfa_list.esc1vs0, 
                                        limits=limits.esc1vs0,
                                        thresholds=thresholds.esc1vs0)


pdf(file=paste0(path2figures, '/vulcano_plots_ESC_1vs0.pdf'),
    height=3.5, width=18)
grid.arrange(grobs=vulcano.list.esc1vs0, ncol=4)
dev.off()



##---------------------------------------------------------------------------

## Differential TFA in ESCs 6 vs 0 hrs
dtfa_list.esc6vs0 <- lapply(
    seurat_list,
    function(seurat){
        s <- subset(seurat,
                    cell_type == 'ESC' &
                     stimulation_time %in% c('0h', '6h'))
        Idents(s) <- s$stimulation_time            
        FindMarkers(object=s,
                    slot='counts',
                    groups = 'stimulation_time',
                    ident.1 = '6h',
                    logfc.threshold = test,
                    min.pct = 0,                    
                    min.diff.pct = -Inf)
    }
)




## ESCs
dtfa_list.esc6vs0 <- add_annotations(dtfa_list.esc6vs0)

thresholds.esc6vs0 <- list(c(-0.27, 0.27), 
                   c(-0.75, 0.5),
                   c(-0.9, 0.9),
                   c(-0.85, 0.85))
limits.esc6vs0 <- list(c(-0.6, 0.6),
               c(-2, 2),
               c(-2, 2),
               c(-1, 1))
names(thresholds.esc6vs0) <- names
names(limits.esc6vs0) <- names

vulcano.list.esc6vs0 <- create_vulcano(dtfa_list.esc6vs0, 
                                        limits=limits.esc6vs0,
                                        thresholds=thresholds.esc6vs0)


pdf(file=paste0(path2figures, '/vulcano_plots_ESC_6vs0.pdf'),
    height=3.5, width=18)
grid.arrange(grobs=vulcano.list.esc6vs0, ncol=4)
dev.off()


pdf(file=paste0(path2figures, '/vulcano_plots.pdf'),
    height=7, width=18)
grid.arrange(vulcano.list.mef1vs0$basal,
                   vulcano.list.mef1vs0$cell_type,
                   vulcano.list.mef1vs0$stimulation_time,
                   vulcano.list.mef1vs0$total, 
                   vulcano.list.mef6vs0$basal,
                   vulcano.list.mef6vs0$cell_type,
                   vulcano.list.mef6vs0$stimulation_time,
                   vulcano.list.mef6vs0$total,
                   ncol=4)
grid.arrange(vulcano.list.esc1vs0$basal,
                   vulcano.list.esc1vs0$cell_type,
                   vulcano.list.esc1vs0$stimulation_time,
                   vulcano.list.esc1vs0$total, 
                   vulcano.list.esc6vs0$basal,
                   vulcano.list.esc6vs0$cell_type,
                   vulcano.list.esc6vs0$stimulation_time,
                   vulcano.list.esc6vs0$total,
                   ncol=4)
dev.off()