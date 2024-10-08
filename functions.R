library(ggplot2)
library(tidyr)
library(dplyr)
library(knitr)
library(pandoc)
library(stringr)
library(readxl)
library(plotly)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(DT)
library(ggplotify)
library(tools)


ALL_DATASET <- "All"

RANKS <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species Adj", "Strain")
RANKS_SP <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# METRICS <- c("F1", "Sensitivity", "Precision", "TP", "TN", "FP", "FN", "FFP", "PearsonCorrelationIntersect",
#              "PearsonCorrelationUnion", "BrayCurtisIntersect", "BrayCurtisUnion", "Bray-Curtis Dissimilarity", 
#              "L2Intersect", "L2Union", "ShannonDiversity", "ShannonDiversityGold", "ShannonDiversityDiff")
METRICS <- c("F1", "Sensitivity", "Precision", "TP", "TN", "FP", "FN", 
             "PearsonCorrelation-TP", "PearsonCorrelation","SpearmanCorrelation-TP", "SpearmanCorrelation", 
             "Bray-Curtis similarity", "BrayCurtisIntersect", "BrayCurtis", "Bray-Curtis Dissimilarity", 
             "L2-TP", "L2", "ShannonDiversity", "ShannonDiversityGold", "ShannonDiversityDiff",
             "ShannonDiversityTP", "ShannonDiversityGoldTP", "ShannonDiversityDiffTP", "RichnessTP", "RichnessGoldTP", "RichnessTPDiff")

pairs <- c( "#a6cee3", "#1f78b4",
            "#fb9a99", "#e31a1c",
            "#b2df8a", "#33a02c",
            "#fdbf6f", "#ff7f00",
            "#cab2d6", "#6a3d9a",
            "#d0cd6c", "#a39047",
            "#ce81d1", "#ac31b0",
            "#b8d06c", "#90a347")

color_palette_paired <- function(pair) {
  
}

mean_sd_label <- function(ls, digits=3) {
  paste(round(mean(ls), digits), "±", round(sd(ls), digits), sep='')
  #paste(round(mean(ls), digits), "±", sd(ls), sep='')
}

harmonic_mean_list <- function(...) {
  arguments <- list(...)
  print(arguments)
  # arguments <- as.list(substitute(list(...)))[-1L]
  denominator <- sum(unlist(lapply(arguments, function(x) { 1/x })))
  result <- length(arguments) / denominator
  return(result)
}

harmonic_mean_vector <- function(vec) {
  arguments <- vec
  # arguments <- as.list(substitute(list(...)))[-1L]
  denominator <- sum(unlist(lapply(arguments, function(x) { 1/x })))
  result <- length(arguments) / denominator
  return(result)
}

read_table <- function(path) {
  extension <- file_ext(path)
  if (extension == "xlsx") {
    table <- read_xlsx(path)
  } else {
    table <- read.csv(path, header=TRUE, sep='\t', row.names=NULL)
  }
  return(table)
}

load_data <- function(meta_path, stats_path, stats_detailed_path, tool_color_path=NULL) {
  data <- list()
  data$tool_order <- NULL
  
  data$stats <- read_table(stats_path)
  data$stats_detailed <- read_table(stats_detailed_path)
  data$meta <- read_table(meta_path)
  
  if (!is.null(tool_color_path)) {
    data$tool_colors <- load_tool_colors(tool_color_path)
    data$tool_order <- load_tool_ordered(tool_color_path)
    
    missing_colors <- setdiff(unique(data$meta$Tool), names(data$tool_colors))
    for (tool in missing_colors) {
      data$tool_colors[[tool]] <- "#000000"
    }
    missing_tools <- setdiff(unique(data$meta$Tool), data$tool_order)
    data$tool_order <- c(data$tool_order, missing_tools)
  }
  
  data$meta$AvailableSpecies[data$meta$AvailableSpecies == "NA"] <- NA
  data$meta$GoldStdTree[data$meta$GoldStdTree == "NA"] <- NA
  
  
  path_df <- data$meta %>% select(Tool, AvailableSpecies) %>% filter(!is.na(AvailableSpecies)) %>% unique()
  path_list <- path_df$AvailableSpecies
  names(path_list) <- path_df$Tool
  data$spdict <- load_species_dict2(path_list)
  
  
  has_multiple_datasets <- length(unique(data$stats$dataset)) > 1
  data$has_multiple_datasets <- has_multiple_datasets
  if (has_multiple_datasets & !(ALL_DATASET %in% data$stats$dataset)) {
    datasets <- unique(data$stats$dataset)
    tmp <- data$stats
    tmp$dataset <- ALL_DATASET
    data$stats <- rbind(tmp, data$stats)
    data$stats$tool <- factor(data$stats$tool)
    data$stats$dataset <- factor(data$stats$dataset, levels=c(ALL_DATASET, datasets), ordered = TRUE)
    remove(tmp)
  } else {
    data$stats$dataset <- factor(data$stats$dataset, levels=c(unique(data$stats$dataset)), ordered = TRUE)
  }
  
  data$stats$sample <- factor(data$stats$sample)
  data$stats$rank <- ordered(data$stats$rank, levels=RANKS)
  data$stats$metric <- ordered(data$stats$metric, levels=METRICS)
  
  if (is.null(data$tool_order)) {
    data$stats$tool <- factor(data$stats$tool)
    data$stats_detailed$Tool <- factor(data$stats_detailed$Tool)
    data$meta$Tool <- factor(data$meta$Tool)
  } else {
    data$stats$tool <- ordered(data$stats$tool, levels=data$tool_order)
    data$stats_detailed$Tool <- ordered(data$stats_detailed$Tool, levels=data$tool_order)
    data$meta$Tool <- ordered(data$meta$Tool, levels=data$tool_order)
  }
  
  return(data)
}

generate_data <- function(data) {
  data$stats_detailed$Detectable <- TRUE
  for (tool in data$stats_detailed %>% filter(Tool %in% names(data$spdict)) %>% pull(Tool) %>% unique()) {
    spd <- data$spdict[[tool]]
    data$stats_detailed$Detectable[data$stats_detailed$Tool == tool & data$stats_detailed$Rank == "Species"] <- data$stats_detailed$Taxon[data$stats_detailed$Tool == tool & data$stats_detailed$Rank == "Species"] %in% spd
  }
  valid_taxa <- grepl("^((d|p|o|c|f|g)__[a-zA-Z0-9_-]+)|(s__[a-zA-Z0-9_-]+ [a-zA-Z0-9_-]+)$|^[0-9]+$", data$stats_detailed$Taxon)
  data$stats_detailed$Type[!valid_taxa & data$stats_detailed$Type == "FP"] <- "Unknown"
  
  new_details <- adjust_gtdb_performance(data$meta, data$stats_detailed, data$stats, 0.04, data$spdict)
  
  if (!is.null(new_details)) {
    new_stats <- details_to_binary_stats(new_details, data$stats %>% filter(dataset != ALL_DATASET)) %>% mutate(adjusted=TRUE, Rank="Species Adj")
    old_stats <- details_to_binary_stats(data$stats_detailed, data$stats %>% filter(dataset != ALL_DATASET)) %>% mutate(adjusted=FALSE)
    abundance_stats <- data$stats %>% filter(!metric %in% (old_stats$metric %>% unique()) & !is.na(metric) & dataset != ALL_DATASET) %>% mutate(adjusted=FALSE)
    colnames(new_stats)[colnames(new_stats) == "Rank"] <- "rank"
    colnames(old_stats)[colnames(old_stats) == "Rank"] <- "rank"
    
    all_stats <-rbind(
      new_stats,
      old_stats,
      abundance_stats
    )
    merged_details <- rbind(
      data$stats_detailed %>% filter(!Tool %in% unique(new_details$Tool)) %>% 
        mutate(Dist=NA, Detectable=NA),
      new_details %>% select(c(colnames(data$stats_detailed), "Dist", "Detectable")) %>%
        mutate(Rank="Species Adj")
    )
  } else {
    all_stats <- data$stats %>% mutate(adjusted=FALSE)
    merged_details <- data$stats_detailed %>% mutate(Dist=NA, Detectable=NA)
  }
  all_stats$rank <- ordered(all_stats$rank, levels=RANKS)
  merged_details$Rank <- ordered(merged_details$Rank, levels=RANKS)
  new_details$Rank <- ordered(new_details$Rank, levels=RANKS)
  all_stats$tool <- ordered(all_stats$tool, levels=data$tool_order)
  all_stats$metric <- ordered(all_stats$metric, levels=METRICS)
  
  data$stats_detailed_merged <- merged_details
  data$stats_detailed_new <- new_details
  data$stats_all <- all_stats
  
  
  
  if (is.null(data$tool_order)) {
    data$stats_all$tool <- factor(data$stats_all$tool)
    data$stats_detailed_new$tool <- factor(data$stats_detailed_new$tool)
    data$stats_detailed_merged$Tool <- factor(data$stats_detailed_merged$Tool)
  } else {
    data$stats_all$tool <- ordered(data$stats_all$tool, levels=data$tool_order)
    data$stats_detailed_new$Tool <- ordered(data$stats_detailed_new$Tool, levels=data$tool_order)
    data$stats_detailed_merged$Tool <- ordered(data$stats_detailed_merged$Tool, levels=data$tool_order)
  }
  
  print("explore fp neighborhood")
  data$fp_neighbors <- NULL
  fp_neighbors <- explore_fp_neighborhood(data$meta, data$stats_detailed_new)
  if (nrow(fp_neighbors) > 0) {
    data$fp_neighbors <- fp_neighbors %>%
      left_join(data$meta %>% select(ID, Dataset), by=join_by(Sample == ID))
  }
  
  if (length(unique(data$fp_neighbors$Dataset)) > 1) {
    data$fp_neighbors <- rbind(
      data$fp_neighbors,
      data$fp_neighbors %>% mutate(Dataset = ALL_DATASET)
    )
    data$stats_all <- rbind(
      data$stats_all,
      data$stats_all %>% mutate(dataset = ALL_DATASET)
    )
  }
  
  return(data)
}


load_tool_colors <- function(path) {
  ext <- file_ext(path)
  color.df <- NULL
  if (ext == "xlsx") {
    color.df <- read_xlsx(path)
  } else {
    color.df <-read.csv(path, header=TRUE, sep='\t')
  }
  color_dict <- color.df$Color
  names(color_dict) <- as.factor(color.df$Tool)
  color_dict
}

load_tool_ordered <- function(path) {
  ext <- file_ext(path)
  color.df <- NULL
  if (ext == "xlsx") {
    color.df <- read_xlsx(path)
  } else {
    color.df <-read.csv(path, header=TRUE, sep='\t')
  }
  color.df %>% arrange(Order) %>% pull(Tool)
}

prediction_colors_list <- list(
  TP="#003200",
  FP="#FF3030",
  FN="#FF9F00")
prediction_colors_vec <- c(
  "TP"="#003200",
  "FP"="#FF3030",
  "FN"="#FF9F00"
)

knit_plot_text <- function(p, width, height) {
  return(c(paste("```{r, echo=FALSE, results=\"asis\", fig.width=", width, ", fig.height=", height, ", warning=FALSE}", sep=''),
    "print(p)",
    "```"))
}

read_detailed_stats <- function(path) {
  ext <- file_ext(path)
  if (ext == 'xlsx') {
    return(read_xlsx(path))
  } else if (ext == 'tsv' | ext == 'csv') {
    return(read.table(path, sep = '\t', header=TRUE))
  }
}

plot.tree.abundance <- function(subtree, detail.df, available_species=NULL) {
  detail.df <- detail.df %>% filter(Type %in% c("TP", "FP", "FN"))
  detail.df$InDatabase <- TRUE
  if (!is.null(available_species)) {
    detail.df$InDatabase <-  ordered(detail.df$Taxon %in% available_species, levels=c(TRUE,FALSE))
  }
  
  detail.df$Type <- ordered(detail.df$Type, levels=c("TP", "FP", "FN"))
  tree.plot <- ggtree(subtree, layout="circular")  %<+% 
    (detail.df %>% mutate(ID=Taxon, Identity=Type, Abundance=round(abundance, 4), InDatabase=InDatabase) %>% select(ID, Identity, Abundance, InDatabase)) + 
    geom_tippoint(aes(color=Identity, shape=InDatabase), size=2) +
    scale_color_manual(values=prediction_colors_list)#sapply(unique(detail.df$Type), function(x) prediction_colors[[x]])
  
    
  detailplot2 <- tree.plot +
    geom_fruit(
      geom=geom_col,
      mapping=aes(y=ID, x=Abundance),  #The 'Abundance' of 'dat1' will be mapped to x
      pwidth=0.4,
      offset = 0,
      axis.params=list(
        axis="x", # add axis text of the layer.
        text.angle=-45, # the text size of axis.
        hjust=0  # adjust the horizontal position of text of axis.
      ),
      grid.params=list() # add the grid line of the external bar plot.
    ) + 
    geom_tiplab(aes(label=as.character(Abundance)), linesize=0, size=3, offset=0.1, align=TRUE) + 
    geom_tiplab(aes(color=Identity), linesize=0, size=3, offset=0.6, align=TRUE) + 
    theme(#legend.position=c(0.96, 0.5), # the position of legend.
      legend.background=element_rect(fill=NA), # the background of legend.
      legend.title=element_text(size=7), # the title size of legend.
      legend.text=element_text(size=6), # the text size of legend.
      legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
    )
  return(detailplot2)
}

load_tree <- function(nwk.path) {
  nwk.tree <- ape::read.tree(nwk.path)
  nwk.tree$tip.label <- gsub("'", "", nwk.tree$tip.label)
  nwk.tree$node.label <- gsub("'", "", nwk.tree$node.label)
  return(nwk.tree)
}


plot_all <- function(ds.all, title, tool_colors=NULL) {
  if (is.null(ds.all)) return(NULL)
  if (nrow(ds.all) == 0) return(NULL)
  # ds.plot <- ds.all %>% 
  #   ggplot(aes(x=reorder(tool, -value), y=value, fill=tool))
  ds.plot <- ds.all %>% 
    ggplot(aes(x=tool, y=value, fill=tool))
  ds.plot <- ds.plot +
    ggtitle(title) +
    geom_violin()+
    geom_boxplot(width=.2)+
    geom_jitter(color="black", size=0.4, alpha=0.2, height = 0) +
    facet_grid(cols = vars(metric), rows = vars(dataset), scales="free_y") +
    xlab("") + ylab("Value") + theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank())
  if (!is.null(tool_colors)) {
    ds.plot <- ds.plot + 
      scale_fill_manual(values = tool_colors)
  }
  
  return(ds.plot)
}

plot_abundance <- function(dssub) {
  ab <- dssub %>% mutate(Abundance=if_else(GOLD == 0, PRED, GOLD)) %>%
    ggplot(aes(x=Tool, y=Abundance, color=Type)) +
    geom_violin() +
    geom_jitter(alpha=0.8) + 
    coord_cartesian(ylim=c(0, 0.0005)) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.000001))
  return(ab)
}


plot_tool <- function(tool.df, title) {
  ds.plot <- tool.df %>% 
    dplyr::filter(metric %in% c("F1", "Sensitivity", "Precision")) %>% 
    ggplot(aes(x=rank, y=value, fill=rank)) +
    ggtitle(title) +
    geom_violin()+
    geom_boxplot(width=.2)+
    geom_jitter(color="black", size=0.4, alpha=0.2, height = 0) +
    facet_grid(cols = vars(metric), rows = vars(dataset), scales="free_y") +
    xlab("") + ylab("Value") + theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank())
  return(ds.plot)
}

evaluate_fp_fn <- function(sample.df, tree, dist_threshold, detectable_species=NULL) {
  test.fp <- sample.df %>% filter(Type == "FP")
  test.fn <- sample.df %>% filter(Type == "FN")
  
  sample.df$Dist <- 0
  
  if (nrow(test.fp) == 0 | nrow(test.fn) == 0) return(sample.df)
  #browser()
  taxa.fp <- intersect(test.fp$Taxon, tree$tip.label)
  taxa.fn <- intersect(test.fn$Taxon, tree$tip.label)
  subtree.cophen <- ape::cophenetic.phylo(tree)[taxa.fp, taxa.fn, drop=F]
  
  if (length(taxa.fn) == 0) return(sample.df)
  
  subtree_pairs <- data.frame(fn=colnames(subtree.cophen)[col(subtree.cophen)], fp=rownames(subtree.cophen)[row(subtree.cophen)], dist=c(subtree.cophen)) %>%
    group_by(fn) %>% 
    summarize(best_fp = fp[which.min(dist)], min_fp = min(dist)) %>% arrange(min_fp)
  used_fp <- c()
  
  
  for (fn_tax in subtree_pairs$fn) {
    sample.df$Detectable[sample.df$Taxon == fn_tax] <- fn_tax %in% detectable_species
    by_fp_list <- subtree_pairs %>% filter(fn == fn_tax & !(best_fp %in% used_fp))
    by_fp <- by_fp_list %>% pull(best_fp)
    
    if (nrow(by_fp_list) == 0) next
    
    at_dist <- subtree_pairs %>% filter(fn == fn_tax) %>% pull(min_fp)
    
    
    if (at_dist > dist_threshold) next
    
    # next
    if (!is.null(detectable_species)) {
      # print(sample.df$Tool[1])
      # print(paste(by_fp, by_fp %in% detectable_species))
      if (fn_tax %in% detectable_species) next
    }
    
    sample.df$Type[sample.df$Taxon == fn_tax] <- "TP"
    sample.df$PRED[sample.df$Taxon == fn_tax] <- sample.df$PRED[sample.df$Taxon == by_fp]
    sample.df$Type[sample.df$Taxon == by_fp] <- "FFP"
    sample.df$Dist[sample.df$Taxon == fn_tax] <- at_dist
    sample.df$Dist[sample.df$Taxon == by_fp] <- at_dist
    
    used_fp <- c(used_fp, by_fp)
  }
  return(sample.df)
}

sensitivity <- function(tp, fn) {
  return(tp/(tp+fn))
}
precision <- function(tp, fp) {
  return(tp/(tp+fp))
}
F1 <- function(tp, fp, fn) {
  return(tp/(tp + (fp+fn)/2))
}


binary_eval <- function(detail.df) {
  samples <- unique(detail.df$Sample)
  new.df <- detail.df %>% 
    aggregate(FP=count(Type == "FP"), FN=count(Type == "FN"), TP=count(Type == "TP"))
  new.df
}

adjust_gtdb_performance <- function(meta, ds.detail, ds.all, dist_threshold, available_species_dict=list()) {
  tree_samples <- meta$ID[!is.na(meta$GoldStdTree)]
  
  if (length(tree_samples) == 0) {
    return(NULL)
  }
  # if (nrow(ds.detail) == 0) {
  #   quit()
  # }
  
  ds.detail.tree <- ds.detail %>% filter(Sample %in% tree_samples & Rank == "Species")
  if (nrow(ds.detail.tree) == 0) {
    quit()
  }
  ds.detail.tree <- ds.detail.tree %>% left_join(meta %>% select(ID, GoldStdTree), by=c('Sample' = 'ID'))
  
  samples <- unique(ds.detail.tree$Sample)
  
  # if (length(samples) == 0) {
  #   quit()
  # }

  ignore_species <- c("s__", "")
  
  new.details <- lapply(samples, FUN = function(sample) {
    test.df <- ds.detail.tree %>% filter(Sample == sample)
    
    test.fp <- test.df %>% filter(Type == "FP")
    test.fn <- test.df %>% filter(Type == "FN")
    
    tool <- test.df$Tool[1]
    available_species <- available_species_dict[[as.character(tool)]] 
    # Important to use as.character. Indexing with factors will use integer and retrieve the wrong value
    # test.df$Detectable <- test.df$Taxon %in% available_species
    
    test.df$Type[test.df$Type == "FP" & !(startsWith(test.df$Taxon, "s__") | test.df$Taxon == "s__")] <- "TN"
    
    if (nrow(test.fp) == 0 | nrow(test.fn) == 0) {
      test.df$Dist <- 0
      return(test.df)
    }
    
    
    tree <- ape::read.tree(test.df$GoldStdTree[1])
    tree$tip.label <- gsub("'", "", tree$tip.label)
    
    shared_taxa <- intersect(test.df$Taxon, tree$tip.label)
    subtree <- ape::keep.tip(tree, shared_taxa)
    
    return(evaluate_fp_fn(test.df, subtree, dist_threshold, available_species))
  }) %>% bind_rows()
  
  
  return(new.details)
}

explore_fp_neighborhood <- function(meta, ds.detail) {
  tree_samples <- meta$ID[!is.na(meta$GoldStdTree)]
  
  if (length(tree_samples) == 0) {
    return(NULL)
  }
  # if (nrow(ds.detail) == 0) {
  #   quit()
  # }
  
  ds.detail.tree <- ds.detail %>% filter(Sample %in% tree_samples & Rank == "Species")
  if (nrow(ds.detail.tree) == 0) {
    quit()
  }
  if (!"GoldStdTree" %in% colnames(ds.detail)) {
    ds.detail <- ds.detail %>% left_join(meta %>% select(ID, GoldStdTree), by=c('Sample' = 'ID'))
  }
  
  samples <- unique(ds.detail$Sample)
  
  if (length(samples) == 0) {
    quit()
  }
  
  ignore_species <- c("s__", "")
  
  new.details <- lapply(samples, FUN = function(sample) {
    sample.df <- ds.detail %>% filter(Sample == sample)
    sample.fp <- sample.df %>% filter(Type == "FP")
    sample.fn <- sample.df %>% filter(Type == "FN")
    sample.tp <- sample.df %>% filter(Type == "TP")
    
    tool <- sample.df$Tool[1]
    
    if (nrow(sample.fp) == 0 | (nrow(sample.tp) == 0 & nrow(sample.fn) == 0)) {
      return(NULL)
    }
    
    tree <- ape::read.tree(sample.df$GoldStdTree[1])
    tree$tip.label <- gsub("'", "", tree$tip.label)
    
    # print(sample)
    shared_taxa <- intersect(sample.df$Taxon, tree$tip.label)
    subtree <- ape::keep.tip(tree, shared_taxa)
    
    num_fp <- sample.fp %>% pull(Taxon) %>% intersect(shared_taxa) %>% length()
    if (length(sample.fp %>% pull(Taxon) %>% intersect(shared_taxa)) == 0) {
      return(NULL)
    }
    shared_fp <- sample.fp %>% pull(Taxon) %>% intersect(shared_taxa)
    shared_fn <- sample.fn %>% pull(Taxon) %>% intersect(shared_taxa)
    shared_tp <- sample.tp %>% pull(Taxon) %>% intersect(shared_taxa)
    # Analyze FP and FN
    fn_df <- NULL
    tp_df <- NULL
    
    if (length(shared_fn) > 0) {
      subtree.cophen.fn <- ape::cophenetic.phylo(subtree)[shared_fp, shared_fn, drop=F]
      fn_df <- data.frame(neighbor=colnames(subtree.cophen.fn)[col(subtree.cophen.fn)], 
                          fp=rownames(subtree.cophen.fn)[row(subtree.cophen.fn)],
                          type="FN",
                          dist=c(subtree.cophen.fn))
    }
    if (length(shared_tp) > 0) {
      subtree.cophen.tp <- ape::cophenetic.phylo(subtree)[shared_fp, shared_tp, drop=F]
      tp_df <- data.frame(neighbor=colnames(subtree.cophen.tp)[col(subtree.cophen.tp)], 
                          fp=rownames(subtree.cophen.tp)[row(subtree.cophen.tp)],
                          type="TP",
                          dist=c(subtree.cophen.tp))
    }
  
    pairs <- rbind(fn_df, tp_df) %>%
      group_by(fp) %>% 
      summarize(closest_neighbor = neighbor[which.min(dist)], distance = min(dist), type=type[which.min(dist)]) %>% 
      arrange(distance) %>%
      left_join(sample.df %>% select(Taxon, GOLD, Sample, Tool), by=join_by(closest_neighbor == Taxon)) %>%
      left_join(sample.df %>% select(Taxon, PRED), by=join_by(fp == Taxon)) 
    # 
    # pairs.fn <- data.frame(fn=colnames(subtree.cophen)[col(subtree.cophen)], 
    #                     fp=rownames(subtree.cophen)[row(subtree.cophen)], 
    #                     dist=c(subtree.cophen)) %>%
    #   group_by(fp) %>% 
    #   summarize(best_fn = fn[which.min(dist)], min_fn = min(dist)) %>% arrange(min_fn)
    # 
    # # Analyse FP and TP
    # subtree.cophen <- ape::cophenetic.phylo(subtree)[sample.fp %>% pull(Taxon) %>% intersect(shared_taxa), sample.tp %>% pull(Taxon) %>% intersect(shared_taxa), drop=F]
    # pairs.tp <- data.frame(tp=colnames(subtree.cophen)[col(subtree.cophen)], 
    #                     fp=rownames(subtree.cophen)[row(subtree.cophen)], 
    #                     dist=c(subtree.cophen)) %>%
    #   group_by(fp) %>% 
    #   summarize(best_tp = tp[which.min(dist)], min_tp = min(dist)) %>% arrange(min_tp) %>%
    #   left_join(sample.tp %>% select(Taxon, GOLD, Sample, Tool), by=join_by(best_tp == Taxon)) 
    # 
    # pairs.tp <- pairs.tp %>% 
    #   left_join()
    # browser()
    
    return(pairs)
  }) %>% bind_rows()
  
  return(new.details)
}

details_to_binary_stats <- function(new.details, ds.all) {
  if (nrow(new.details) == 0) return(NULL)
  
  binary_eval <- new.details %>% 
    filter(Type %in% c("FP", "TP", "FN", "FFP")) %>% 
    count(Sample, Type, Rank) %>% 
    pivot_wider(names_from=Type, values_from = n, values_fill = 0)
  
  colnames(binary_eval)[colnames(binary_eval) == "Sample"] <- "sample"
  
  if (!"FFP" %in% colnames(binary_eval)) {
    binary_eval$FFP <- 0
  }
  if (!"FP" %in% colnames(binary_eval)) {
    binary_eval$FP <- 0
  }
  if (!"TP" %in% colnames(binary_eval)) {
    binary_eval$TP <- 0
  }
  if (!"FN" %in% colnames(binary_eval)) {
    binary_eval$TN <- 0
  }
  new_stats <- binary_eval %>% 
    mutate(Sensitivity=sensitivity(TP, FN), Precision=precision(TP, FP), F1=F1(TP, FP, FN)) %>%
    pivot_longer(c("Sensitivity", "Precision", "F1", "FP", "FN", "TP", "FFP"), names_to="metric") %>%
    left_join(ds.all %>% filter(dataset != ALL_DATASET) %>% select(sample, dataset, tool) %>% unique(), by="sample", relationship="many-to-many")
  
  return(new_stats)
}

details_to_binary_stats2 <- function(ds.details, ds.all) {
  ds.details %>% 
    filter(Type %in% c("FP", "TP", "FN") & Rank == "Species") %>% 
    count(Sample, Type, Threshold) %>% 
    pivot_wider(names_from=Type, values_from = n, values_fill = 0) %>%
    mutate(Sensitivity=sensitivity(TP, FN), Precision=precision(TP, FP), F1=F1(TP, FP, FN)) %>%
    pivot_longer(c("Sensitivity", "Precision", "F1", "FP", "FN", "TP"), names_to="metric") %>%
    left_join(ds.all %>% filter(rank == "Species") %>% select(sample, dataset, tool, rank) %>% unique(), by=join_by(Sample == sample), relationship="many-to-many")
}

apply_abundance_threshold <- function(ds.detail, threshold) {
  ds.detail$Type[ds.detail$Type == "FP" & ds.detail$PRED < threshold] <- "TN"
  ds.detail$Type[ds.detail$Type == "TP" & ds.detail$PRED < threshold] <- "FN"
  ds.detail$Threshold <- threshold
  return(ds.detail)
}

plot_abundance_threshold <- function(data) {
  data %>% 
    filter(metric %in% c("Sensitivity", "Precision", "F1")) %>%
    group_by(Threshold, metric, tool) %>%
    summarise(mean=mean(value), median=median(value), tool=head(tool,1)) %>%
    ggplot(aes(x=Threshold, y=mean, color=metric, group=metric)) + 
    geom_line() +
    ylab("Mean across datasets") +
    facet_wrap(~ tool) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

distances <- function(detailed, cophen) {
  dists <- data.frame(sample=character(0), species=character(0), fn_species=character(0), tp_species=character(0), distance_to_fn=numeric(0), distance_to_tp=numeric(0), GOLD=numeric(0), PRED=numeric(0))
  allowed <- rownames(cophen)

  result  <- list()
  for (sample in unique(detailed$Sample)) {
    # print(sample)
    dsub <- detailed %>% filter(Sample == sample)
    
    if (length(dsub %>% filter(Type == "FP")) == 0 | length(dsub %>% filter(Type == "FN")) == 0) next
    
    fps <- dsub %>% filter(Type == "FP") %>% pull(Taxon)
    fns <- dsub %>% filter(Type == "FN") %>% pull(Taxon)
    tps <- dsub %>% filter(Type == "TP") %>% pull(Taxon)
    
    
    if (length(fns) > 0 & length(fns %in% allowed) == 0) {
      if (is.null(result$samples[[sample]])) {
        result$samples[[sample]] <- list()
      }
      result$samples[[sample]]$fns <- fns
      result$samples[[sample]]$allowed <- allowed
      next
    }
    if (length(tps) > 0 & length(tps %in% allowed) == 0) {
      if (is.null(result$samples[[sample]])) {
        result$samples[[sample]] <- list()
      }
      result$samples[[sample]]$tps <- tps
      result$samples[[sample]]$allowed <- allowed
      next
    }
    
    for (fp in fps) {
      if (!(fp %in% colnames(cophen))) {
        next
      }
      if (!any(fns %in% allowed)) {
        next
      }
      mini.fn <- which.min(cophen[fp, fns[fns %in% allowed]])
      mini.tp <- which.min(cophen[fp, tps[tps %in% allowed]])
      min.fn <- min(cophen[fp, fns[fns %in% allowed]])
      min.tp <- min(cophen[fp, tps[tps %in% allowed]])
      
      if (is.null(names(mini.fn))) {
        next
      }
      
      dists[nrow(dists)+1,] <- c(sample, fp, names(mini.fn), names(mini.tp), as.numeric(min.fn), as.numeric(min.tp), dsub[dsub$Taxon == fp, "GOLD"], dsub[dsub$Taxon == fp, "PRED"])
    }
  }
  dists$distance_to_fn <- as.numeric(dists$distance_to_fn)
  dists$distance_to_tp <- as.numeric(dists$distance_to_tp)
  result$dist <- dists
  return(result)
}


load_species_dict <- function(filepaths) {
  dicts <- list()
  for (filepath in na.omit(unique(filepaths))) {
    dicts[[filepath]] <-read.csv(filepath, header=FALSE) %>% pull(V1)
  }
  return(dicts)
}

load_species_dict2 <- function(path_list) {
  dicts <- list()
  for (tool in names(path_list)) {
    dicts[[tool]] <-read.csv(path_list[[tool]], header=FALSE) %>% pull(V1)
  }
  return(dicts)
}

#, metric_col="metric"
summarize.df <- function(df, target_metrics, merge_datasets=FALSE) {
  # print(df %>% filter(metric %in% target_metrics) %>% pull(metric))
  return (agg.df <- df %>% filter(metric %in% target_metrics) %>%
    pivot_wider(names_from=metric, values_from=value) %>% 
    group_by(tool, dataset) %>%
    summarise(across(target_metrics, lst(mean, sd), .names = "{.col}_{.fn}"))  %>% 
    pivot_longer(cols = !c("tool", "dataset"), 
                names_to = c("metric", "type"),
                names_pattern = "(.*)_(.*)",
                values_to = "score") %>%
    pivot_wider(names_from = type, values_from = score) %>%
    mutate(summary_str=paste(round(mean, 3), "±", round(sd, 3), sep="")))
}

summarize.df.all <- function(df, target_metrics) {
  # print(df %>% filter(metric %in% target_metrics) %>% pull(metric))
  return (df %>% filter(metric %in% target_metrics) %>%
            select(-dataset) %>%
            group_by(tool, metric) %>%
            summarize(across(value, lst(mean, sd), .names = "{.fn}")) %>% as.data.frame() %>%
            mutate(summary_str=paste(round(mean, 3), "±", round(sd, 3), sep="")))
}

highlight_all <- function(all.plot, column_idx = 1) {
  g <- ggplot_gtable(ggplot_build(all.plot))
  strips <- which(grepl('strip-', g$layout$name))
  
  highlight_background <- "#555555"
  highlight_font <- "#FFFFFF"
  
  bold <- as.integer(2)
  
  i <- column_idx # Column number
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- highlight_background
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- highlight_font
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$font <- bold
  
  all.plot <- as.ggplot(g)
  return(all.plot)
}

















