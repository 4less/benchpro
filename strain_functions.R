library(vegan)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(phytools)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(knitr)
library(plotly)
library(reshape2)
library(phangorn)
library(cowplot)
library(tidyquant)

MAP_NAME_COL <<- "Name"
MAP_TREE_COL <<- "Tree"
MAP_GROUP_TREE_COL <<- "SpeciesTree"
MAP_MSA_COL <<- "MSA"
MAP_MSA_ANALYSIS_COL <<- "AnalysisMSA"
MAP_META_COL <<- "Meta"
META_GROUP_COL <<- "genome"



run_msage <- function(map) {
  python <- "/usr/users/QIB_fr017/fritsche/miniconda3/bin/python"
  msage_path <- "/usr/users/QIB_fr017/fritsche/Projects/labutils/scripts/tool_related/protal/msage"
  
  for (i in seq(nrow(map))) {
    msa.path <- map[i,] %>% pull(MAP_MSA_COL)
    meta.path <- map[i,] %>% pull(MAP_META_COL)
    name <- map[i,] %>% pull(MAP_NAME_COL)
    
    msa.analysis.path <- map[i,] %>% pull(MAP_MSA_ANALYSIS_COL)
    
    if (file.exists(msa.analysis.path) | !file.exists(msa.path)) {
      print(paste("Exists ", msa.analysis.path))
      next
    }
    
    cmd <- paste(
      python,
      msage_path,
      "analyse",
      "--meta", paste("\"", meta.path, "\"", sep=''),
      "--msa", paste("\"", msa.path, "\"", sep=''),
      "--output", msa.analysis.path, 
      "",
      sep=' '
    )
    cat(cmd)
    system(cmd)
  }
}

read_msage <- function(map) {
  result <- NULL
  cnames <- c("Genome", "GroupSize", "SNPs", "AlignmentLength", "AlignmentLength3", "ErrorRate")
  for (i in seq(nrow(map))) {
    msa.analysis.path <- map[i,] %>% pull(MAP_MSA_ANALYSIS_COL)
    name <- map[i,] %>% pull(MAP_NAME_COL)
    if (file.exists(msa.analysis.path)) {
      print(msa.analysis.path)
      df <- read.csv(msa.analysis.path, header=FALSE, sep='\t')
      ##colnames(df) <- c("Genome", "GroupSize", "SNPs", "AlignmentLength", "ErrorRate", "MinCov", "WholeGroup")
      colnames(df) <- cnames
      df$Name <- name
      #df$WholeGroup <- df$WholeGroup == "True"
      if (is.null(result)) {
        result <- df
      } else {
        result <- rbind(result, df)
      }
    }
  }
  return(result)
}

# Function definitions
load_map <- function(map_file) {
  csv <- read.csv(map_file, header=TRUE, sep='\t')
  rownames(csv) <- csv$Name
  csv$HasTree <- file.exists(csv$Tree)
  return(csv)
}

load_meta <- function(meta_file, id_col = "ID") {
  csv <- read.csv(meta_file, header=TRUE, sep='\t')
  rownames(csv) <- csv[,id_col]
  return(csv)
}

load_all_map <- function(map) {
  data <- list()
  data$map <- map
  meta.dict <- get_meta_dict(data$map)
  data$meta.dict <- meta.dict
  
  generate_gold_stds(data$map, data$meta.dict)
  data$tree.dict <- get_tree_dict(data$map)
  data$tree.species.dict <- get_tree_dict_gold(data$map)
  data$tree.dict.shared.tips <- get_tree_dict_shared_tips(data$map, data$tree.dict)
  data$tree.dict.pruned <- get_tree_dict_prune(data$map, data$meta.dict)
  data$tree.dict.shared.pruned <- get_tree_dict_shared_prune(data$map, data$meta.dict, data$tree.dict.shared.tips)
  data$shared_species <- species_for_tools(data$map)
  
  return(data)
}

get_tree_dict_prune <- function(map, meta.dict, midpoint.root=TRUE) {
  names <- map[, MAP_NAME_COL]
  paths <- map[, MAP_TREE_COL]
  paths_gold <- map[, MAP_GROUP_TREE_COL]
  metas <- map[, MAP_META_COL]
  tree_dict <- list()
  
  for (rid in 1:length(names)) {
    name <- names[rid]
    nwk.path <- paths[rid]
    nwk.gold.path <- paths_gold[rid]
    meta <- meta.dict[[metas[rid]]]
    
    if (!file.exists(nwk.path) & !file.exists(nwk.gold.path)) {
      tree_dict[[name]] <- NULL
      next
    }
    
    if (file.exists(nwk.gold.path)) {
      tree <- ape::read.tree(nwk.gold.path)
      tree <- midpoint.root(tree)
      tree_dict[[name]] <- tree
      next
    }
    
    tree <- ape::read.tree(nwk.path)
    
    if (midpoint.root) {
      tree <- midpoint.root(tree)
    }
    
    meta.sub <- meta[tree$tip.label,]
    tree.pruned <- prune_tree(tree, meta.sub)
    
    tree_dict[[name]] <- tree.pruned
  }
  return(tree_dict)
}

get_tree_dict_shared_prune <- function(map, meta.dict, tree.dict, midpoint.root=TRUE) {
  names <- map[, MAP_NAME_COL]
  paths <- map[, MAP_TREE_COL]
  paths_gold <- map[, MAP_GROUP_TREE_COL]
  metas <- map[, MAP_META_COL]
  tree_dict_prune <- list()
  
  for (rid in 1:length(names)) {
    name <- names[rid]
    nwk.path <- paths[rid]
    nwk.gold.path <- paths_gold[rid]
    meta <- meta.dict[[metas[rid]]]
    
    if (!file.exists(nwk.path) & !file.exists(nwk.gold.path)) {
      tree_dict_prune[[name]] <- NULL
      next
    }
    
    if (file.exists(nwk.gold.path)) {
      tree <- ape::read.tree(nwk.gold.path)
      tree <- midpoint.root(tree)
      tree_dict_prune[[name]] <- tree
      next
    }
    
    tree <- tree.dict[[name]]
    
    if (midpoint.root) {
      tree <- midpoint.root(tree)
    }
    
    meta.sub <- meta[tree$tip.label,]
    tree.pruned <- prune_tree(tree, meta.sub)
    
    tree_dict_prune[[name]] <- tree.pruned
  }
  return(tree_dict_prune)
}


get_tree_dict_gold <- function(map, midpoint.root=TRUE) {
  names <- map[, MAP_NAME_COL]
  paths <- map[, MAP_GROUP_TREE_COL]
  tree_dict <- list()
  
  for (rid in 1:length(names)) {
    name <- names[rid]
    nwk.path <- paths[rid]
    
    if (!file.exists(nwk.path)) {
      tree_dict[[name]] <- NULL
      next
    }
    tree <- ape::read.tree(nwk.path)
    if (midpoint.root) {
      tree <- midpoint.root(tree)
    }
    tree_dict[[name]] <- tree
  }
  return(tree_dict)
}



load_all <- function(map) {
  # map <- load_map(map.file)
  data <- load_all_map(map)
  return(data)
}

get_tree_dict <- function(map, midpoint.root=TRUE) {
  names <- map[, MAP_NAME_COL]
  paths <- map[, MAP_TREE_COL]
  tree_dict <- list()
  
  for (rid in 1:length(names)) {
    name <- names[rid]
    nwk.path <- paths[rid]
    
    # print(nwk.path)
    if (!file.exists(nwk.path)) {
      tree_dict[[name]] <- NULL
      next
    }
    tree <- ape::read.tree(nwk.path)
    if (midpoint.root) {
      tree <- midpoint.root(tree)
    }
    tree_dict[[name]] <- tree
  }
  return(tree_dict)
}

get_tree_dict_shared_tips <- function(map, tree.dict) {
  tree_dict_shared_tips <- list()
  for (species in unique(map$Species)) {
    names <- map %>% filter(Species == species) %>% pull(Name)
    tree.dict.sub <- Filter(function(x) !is.null(unlist(x)), tree.dict[names])
    shared_tips_dict <- shared_tips_tree_dict(tree.dict.sub)
    tree_dict_shared_tips <- append(tree_dict_shared_tips, shared_tips_dict)
  }
  return(tree_dict_shared_tips)
}

generate_gold_std_tree <- function(tree, mapping) {
  tag <- NA
  old_species <- NA
  for (rid in 1:nrow(mapping)) {
    new_label <- mapping[rid, 1]
    species <- mapping[rid, 2]
    if (is.na(old_species)) old_species <- species
    if (is.na(tag) || species != old_species) {
      tag <- species
      tree$tip.label[tree$tip.label == tag] <- new_label
    } else {
      index <- which(tree$tip.label == tag)
      #index <- grep(tag, tree$tip.label)
      if (length(index) > 1) {
        print("abort")
        browser()
        return()
      }
      tree <- bind.tip(tree, new_label, edge.length=0, where=index, pos=0)
    }
    tag <- new_label
    old_species <- species
  }
  print("return")
  return(tree)
}


get_meta_dict <- function(map) {
  paths <- map[, MAP_META_COL]
  meta_dict <- list()
  
  for (path in unique(paths)) {
    meta_dict[[path]] <- load_meta(path)
  }
  return(meta_dict)
}

quicktree <- function(nwk.file, meta=NA, group_var=NA, midpoint=TRUE) {
  tree <- nwk.file
  if (is.character(nwk.file)) {
    tree <- ape::read.tree(nwk.file)
  }
  if (midpoint) {
    tree <- midpoint.root(tree)
  }
  if (length(meta) != 1 || !is.na(meta)) {
    tree.plot <- ggtree(tree, layout='circular') %<+% meta  +
      geom_tippoint(aes(color=.data[[group_var]]), size=3) +
      geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  } else {
    tree.plot <- ggtree(tree, layout='circular') +
      geom_tippoint(size=3) +
      geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  }
  return(tree.plot)
}

quicktree <- function(tree, meta=NA, group_var=NA, midpoint=FALSE) {
  if (midpoint) {
    tree <- midpoint.root(tree)
  }
  if (length(meta) != 1 || !is.na(meta)) {  
    tree.plot <- ggtree(tree, layout='circular') %<+% meta  +
      geom_tippoint(aes(color=.data[[group_var]]), size=3) +
      geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  } else {
    tree.plot <- ggtree(tree, layout='circular') +
      geom_tippoint(size=3) +
      geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  }
  return(tree.plot)
}

shared_tips_tree_dict <- function(tree_dict) {
  new_tree_dict <- list()
  common.tips <- Reduce(intersect, lapply(tree_dict, FUN=function(x){return(x$tip.label)}))
  
  if (length(common.tips) == 0) {
    print("----------------")
    print(names(tree_dict))
  }
  
  for (name in names(tree_dict)) {
    tree <- tree_dict[[name]]
    new_tree_dict[[name]] <- ape::keep.tip(tree, common.tips)
  }
  return(new_tree_dict)
}

pairwise_mantel <- function(tree_dict, permutations=999, parallel=4) {
  tree.dict.sub <- Filter(function(x) !is.null(unlist(x)), tree_dict)
  names_order <- tree_dict[[1]]$tip.label
  
  cophenetic_list <- lapply(tree.dict.sub, FUN=function(x){return(cophenetic.phylo(x)[names_order, names_order])})
  
  df <- data.frame(df = data.frame(matrix(nrow = 0, ncol = 3)) )
  colnames(df) <- c("Name1", "Name2", "Cor")
  
  print(length(tree_dict))
  
  if (length(names(tree.dict.sub)) < 2) {
    return(df)
  }
    
    
  pairs <- t(combn(names(tree.dict.sub),2))
  
  for (rid in 1:nrow(pairs)) {
    name1 <- pairs[rid, 1]
    name2 <- pairs[rid, 2]
    mcor <- mantel(cophenetic_list[[name1]], cophenetic_list[[name2]], method="spearman", parallel=6, permutations=permutations)
    df[nrow(df) + 1,] <- c(name1, name2, mcor$statistic)
  }
  return(df)
}

prune_tree <- function(tree, meta) {
  new.tree <- tree
  
  for (group in unique(meta$genome)) {
    samples <- meta %>% 
      filter(genome == group) %>% 
      pull(ID) %>% 
      intersect(tree$tip.label)
    
    if (length(samples) == 0) next
    if (length(samples) > 1) {
      # for each group find LCA
      mrca <- ape::getMRCA(new.tree, samples)
      new.tree <- bind.tip(new.tree, group, edge.length=0, where=mrca)
    } else {
      new.tree$tip.label[new.tree$tip.label == samples[1]] <- group
    }
  }
  new.tree <- drop.tip(new.tree, tree$tip.label)
  return(new.tree)
}




pairwise_treedist <- function(tree_dict, rooted=FALSE) {
  tree.dict.sub <- Filter(function(x) !is.null(unlist(x)), tree_dict)
  names_order <- tree_dict[[1]]$tip.label
  
  df <- data.frame(df = data.frame(matrix(nrow = 0, ncol = 7)) )
  colnames(df) <- c("Name1", "Name2", "symmetric.difference", "branch.score.difference", "path.difference", "quadratic.path.difference", "weighted.rf")
  
  # print(length(names(tree.dict.sub)))
  # print(length(names(tree_dict)))
  if (length(names(tree.dict.sub)) < 2) {
    return(df)
  }
  
  pairs <- t(combn(names(tree.dict.sub),2))
  for (rid in 1:nrow(pairs)) {
    name1 <- pairs[rid, 1]
    name2 <- pairs[rid, 2]
    
    t1 <- tree.dict.sub[[name1]]
    t2 <- tree.dict.sub[[name2]]
    
    treedist <- phangorn::treedist(t1, t2)
    weighted.rf <- wRF.dist(t1, t2, rooted=rooted)
    df[nrow(df) + 1,] <- c(name1, name2, 
                           treedist[["symmetric.difference"]], 
                           treedist[["branch.score.difference"]], 
                           treedist[["path.difference"]], 
                           treedist[["quadratic.path.difference"]], 
                           weighted.rf)
  }
  return(df)
}


get_pairwise_tree_distances <- function(map, tree.dict) {
  tree_dist.df <- data.frame(
    Name1=character(0), 
    Name2=character(0), 
    symmetric.difference=double(0), 
    branch.score.difference=double(0), 
    path.difference=double(0), 
    quadratic.path.difference=double(0), 
    weighted.rf=double(0))
  
  for (species in unique(map$Species)) {
    names <- map %>% filter(Species == species) %>% pull(Name)
    tree.dict.sub <- Filter(function(x) !is.null(unlist(x)), tree.dict[names])
    treedists.species <- pairwise_treedist(tree.dict.sub)
    tree_dist.df <- rbind(tree_dist.df, treedists.species)
  }
  tree_dist.df$symmetric.difference <- as.numeric(tree_dist.df$symmetric.difference)
  tree_dist.df$branch.score.difference <- as.numeric(tree_dist.df$branch.score.difference)
  tree_dist.df$path.difference <- as.numeric(tree_dist.df$path.difference)
  tree_dist.df$quadratic.path.difference <- as.numeric(tree_dist.df$quadratic.path.difference)
  tree_dist.df$weighted.rf <- as.numeric(tree_dist.df$weighted.rf)
  
  return(gold_std_pairs_only(tree_dist.df))
}

get_pairwise_mantel <- function(map, tree.dict) {
  mantel.df <- data.frame(Name1=character(0), Name2=character(0), Cor=double(0))
  
  # Get pairwise mantel tests for each species.
  for (species in unique(map$Species)) {
    names <- map %>% filter(Species == species) %>% pull(Name)
    tree.dict.sub <- Filter(function(x) !is.null(unlist(x)), tree.dict[names])
    mantel.species <- pairwise_mantel(tree.dict.sub, parallel=8, permutations=200)
    mantel.df <- rbind(mantel.df, mantel.species)
  }
  mantel.df$Cor <- as.numeric(mantel.df$Cor)
  mantel.m <- gold_std_pairs_only(mantel.df)
  
  return(mantel.m)
}


intersect_trees <- function(tree1, tree2) {
  shared_tips <- intersect(tree1$tip.label, tree2$tip.label)
  result <- list()
  result$tree1 <- ape::keep.tip(tree1, shared_tips)
  result$tree2 <- ape::keep.tip(tree2, shared_tips)
  return(result)
}

dist_topo_shared <- function(tree1, tree2) {
  tree.intersection <- intersect_trees(tree1, tree2)
  return(ape::dist.topo(ape::unroot(tree.intersection$tree1), ape::unroot(tree.intersection$tree2)))
}

tokey <- function(key1, key2) {
  return(paste(min(key1, key2), max(key1, key2), collapse='_'))
}

pairwise_matrix_to_dict <- function(matrix, col1=1, col2=2, col_val=3) {
  return_dict <- list()
  
  for (i in 1:nrow(matrix)) {
    key1 <- matrix[i, col1]
    key2 <- matrix[i, col2]
    key <- tokey(key1, key2)
    value <- matrix[i, col_val]
    
    return_dict[[key]] <- value
  }
  
  return(return_dict)
}

gold_std_pairs_only <- function(mmatrix) {
  colnames_tail <- tail(colnames(mmatrix), -2)
  
  first <- mmatrix[grep("gold_std", mmatrix$Name1),  c("Name2", "Name1", colnames_tail)]
  second <- mmatrix[grep("gold_std", mmatrix$Name2),  c("Name1", "Name2", colnames_tail)]
  
  colnames(first) <- c("Name", "GoldStdName", colnames_tail)
  colnames(second) <- c("Name", "GoldStdName", colnames_tail)
  
  return(rbind(first, second))
}

get_tree_meta <- function(name, map, meta.dict) {
  md <- meta.dict[[map[map$Name == name, MAP_META_COL]]]
  return(md)
}

get_name_to_group <- function(name, map, meta.dict, group_var) {
  md <- get_tree_meta(name, map, meta.dict)
  group_dict <- setNames(as.list(md[, group_var]), md$ID)
  return(group_dict)
}

monophyly <- function(tree, tree.group) {
  monophyly.list <- c()
  groups <- unique(tree.group)
  working_groups <- c()
  for (group in groups) {
    # print(group)
    # browser()
    # Only take tips (samples) that are present in both tree and the given set of samples tree.group
    group.tips <- intersect(names(tree.group[tree.group %in% group]), tree$tip.label)
    if (length(group.tips) < 2) {
      next
    }
    
    
    lca <- getMRCA(tree, group.tips[group.tips %in% tree$tip.label])
    
    if (is.null(lca)) {
      next
    }
    
    working_groups <- c(working_groups, group)
    clade.tree <- extract.clade(tree, lca, collapse.singles = FALSE)
    monophyly.item <- length(group.tips)/length(clade.tree$tip.label)
    monophyly.list <- c(monophyly.list, monophyly.item)
  }
  
  return(setNames(as.list(monophyly.list), working_groups))
}

monophylies <- function(trees, map, meta.dict, group_name) {
  df <- data.frame(df = data.frame(matrix(nrow = 0, ncol = 6)) )
  colnames(df) <- c("Name", "Measure", "Group", "Species", "Value", "GroupSize")
  measure <- "Monophyly"
  
  for (name in names(trees)) {
    tree <- trees[[name]]
    if (is.null(tree)) next
    
    # print(name)
    tree.group <- get_name_to_group(name, map, meta.dict, group_name)
    monophyly.list <- monophyly(tree, tree.group)
    for (group in names(monophyly.list)) {
      group.tips <- intersect(names(tree.group[tree.group %in% group]), tree$tip.label)
      df[nrow(df)+1,] <- c(name, measure, group, map$Species[map$Name == name], monophyly.list[[group]], length(group.tips))
    }
  }
  df$Value <- as.numeric(df$Value)
  df$GroupSize <- as.numeric(df$GroupSize)
  return(df)
}

monophylies_per_coverage <- function(trees, map, meta.dict, group_name, coverages=c(50, 30, 20, 10, 5, 2)) {
  df <- data.frame(df = data.frame(matrix(nrow = 0, ncol = 6)) )
  colnames(df) <- c("Name", "Measure", "Group", "Species", "Value", "CoverageThreshold")
  measure <- "Monophyly"
  
  for (name in names(trees)) {
    tree <- trees[[name]]
    if (is.null(tree)) next
    
    
    # print(name)
    tree.group <<- get_name_to_group(name, map, meta.dict, group_name)
    tree.meta <- get_tree_meta(name, map, meta.dict)
    tree.meta$coverage_value <- as.numeric(sub("coverage", "", tree.meta$coverage))
    
    for (coverage_threshold in coverages) {
      allowed_labels <- tree.meta %>% filter(coverage_value >= coverage_threshold) %>% pull(ID) %>% intersect(tree$tip.label)
      if (length(allowed_labels) < 2) next
      subtree <- keep.tip(tree, allowed_labels)
      monophyly.list <- monophyly(subtree, tree.group[allowed_labels])
      
      for (group in names(monophyly.list)) {
        df[nrow(df)+1,] <- c(name, measure, group, map$Species[map$Name == name], monophyly.list[[group]], coverage_threshold)
      }
    }
  }
  df$Value <- as.numeric(df$Value)
  df$CoverageThreshold <- as.numeric(df$CoverageThreshold)
  df$CoverageThresholdF <- as.factor(df$CoverageThreshold)
  return(df)
}


get_list <- function(m, meta) {
  m[upper.tri(m)] <- NA
  m.melt <- melt(m) %>% filter(!is.na(value))
  colnames(m.melt) <- c("sample1", "sample2", "value")
  m.melt$identity <- apply(m.melt, 1, function(x) {
    return(meta[[x[1]]] == meta[[x[2]]])
  })
  return(m.melt)
}


# "sample1"  "sample2"  "value"    "identity" "Name"    
# also need genome
pairwise_distance <-  function(tree, tree.group) {
  tree.matrix <- cophenetic.phylo(tree)
  tree.melt <- get_list(tree.matrix, tree.group)
  
  tree.melt$genome1 <- unlist(tree.group[as.character(tree.melt$sample1)])
  tree.melt$genome2 <- unlist(tree.group[as.character(tree.melt$sample2)])
  
  tree.melt$value <- 1 - tree.melt$value
  return(tree.melt)
}

complete_pairwise <- function(pairwise) {
  result.tmp <- pairwise
  result.tmp$sample1 <- pairwise$sample2
  result.tmp$sample2 <- pairwise$sample1
  result.tmp$genome1 <- pairwise$genome2
  result.tmp$genome2 <- pairwise$genome1
  result <- rbind(pairwise, result.tmp)
  result$Group <- result$genome1
  return(result)
}


pairwise_distances <- function(trees, map, meta.dict, group_name) {
  result <- NA
  for (name in names(trees)) {
    tree <- trees[[name]]
    
    tree.group <- get_name_to_group(name, map, meta.dict, group_name)
    tree.melt <- pairwise_distance(tree, tree.group)
    
    tree.melt$Name <- name
    tree.melt$min <- min(cophenetic.phylo(tree))
    tree.melt$mean <- mean(cophenetic.phylo(tree))
    tree.melt$max <- max(cophenetic.phylo(tree))
    
    result <- if (is.null(nrow(result))) tree.melt else rbind(result, tree.melt)
  }
  return(result)
}


# More tree printing options
quicktreemeta <- function(map, tree.dict, meta.dict, tree.name) {
  tree <- tree.dict[[tree.name]]
  meta <- meta.dict[[map %>% filter(Name == tree.name) %>% pull(MAP_META_COL)]]
  
  tree.plot <- ggtree(tree, layout='circular') %<+% meta  +
    geom_tippoint(aes(color=.data[[META_GROUP_COL]]), size=3) +
    geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  print(tree.plot)
  
  return(tree.plot)
}

quicktree3 <- function(tree, meta) {
  tree.plot <- ggtree(tree, layout='circular') %<+% meta  +
    geom_tippoint(aes(color=.data[[META_GROUP_COL]]), size=3) +
    geom_tiplab(linesize=3, size=3, offset=0.001)+ theme(legend.position="None")
  print(tree.plot)
  
  return(tree.plot)
}

get_tree_meta_plot <- function(map, tree.dict, meta.dict, tree.name) {
  tree <- tree.dict[[tree.name]]
  meta <- meta.dict[[map %>% filter(Name == tree.name) %>% pull(MAP_META_COL)]]
  return(list(
    tree=tree,
    meta=meta,
    plot=quicktree3(tree, meta)
  ))
}

sample2group <- function(meta) {
  groups <- meta %>% pull(META_GROUP_COL)
  names(groups) <- meta %>% pull(ID)
  return(groups)
}

comsub<-function(x) {
  # sort the vector
  x<-sort(x)
  # split the first and last element by character
  d_x<-strsplit(x[c(1,length(x))],"")
  # search for the first not common element and so, get the last matching one
  der_com<-match(FALSE,do.call("==",d_x))-1
  # if there is no matching element, return an empty vector, else return the common part
  ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))
}


abbrev_species <- function(x) {
  spl <- unlist(str_split(gsub("s__", "", x), "_", 2))
  spl[1] <- substring(spl[1], 1, 1)
  spl[2] <- gsub("_", " ", spl[2])
  paste(spl, collapse=". ")
}

plot_tree_meta <- function(res, plot_group_size=TRUE, condense_labels=FALSE) {
  res$meta$tlab <- res$meta$genome
  if (condense_labels) {
    res$meta$shortlab <- gsub(comsub(res$meta$genome), "", res$meta$genome)
    res$meta$tlab <- res$meta$shortlab
  }
  if (plot_group_size) {
    s2g <- sample2group(res$meta)
    t2g <- s2g[res$tree$tip.label]
    g2counts <- table(t2g)
    res$meta$insampleg <- g2counts[res$meta$genome]
    res$meta$insampleg[is.na(res$meta$insampleg)] <- 0 
    res$meta$tlab <- paste(res$meta$tlab, "  (", res$meta$insampleg, ")", sep='')
  }
  
  tree.plot <- ggtree(res$tree, layout="circular")  %<+% res$meta + 
    geom_tippoint(aes(color=.data[[META_GROUP_COL]]), size=3)
  p <- tree.plot +
    geom_fruit(
      geom=geom_col,
      mapping=aes(y=ID, x=coverage),  #The 'Abundance' of 'dat1' will be mapped to x
      pwidth=0.4,
      offset = 0,
      axis.params=list(
        axis="x", # add axis text of the layer.
        text.angle=-45, # the text size of axis.
        hjust=0  # adjust the horizontal position of text of axis.
      ),
      grid.params=list() # add the grid line of the external bar plot.
    ) + 
    geom_tiplab(aes(label=as.character(round(coverage, digits=1))), linesize=0, size=2.5, offset=0.0002, align=TRUE) + 
    geom_tiplab(aes(label=tlab), linesize=0, size=2.5, offset=0.001, align=TRUE) + 
    theme(#legend.position=c(0.96, 0.5), # the position of legend.
      legend.background=element_rect(fill=NA), # the background of legend.
      legend.title=element_text(size=7), # the title size of legend.
      legend.text=element_text(size=6), # the text size of legend.
      legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
    )
  return(p)
}

generate_gold_stds <- function(map,  meta.dict) {
  for (treepath in map[map$Tool == "gold_std", MAP_GROUP_TREE_COL]) {
    if (!file.exists(treepath)) {
      next
    }
    gold.tree.path <- map %>% filter(SpeciesTree == treepath) %>% pull(Tree)
    
    if (file.exists(gold.tree.path)) {
      next
    }
    
    tree <- ape::read.tree(treepath)
    tip.map <- meta.dict[[map %>% filter(SpeciesTree == treepath) %>% pull(Meta)]] %>% select(ID, genome)
    tip.map <- tip.map[order(tip.map[,2],decreasing=FALSE),]
    gold_tree <- generate_gold_std_tree(tree, tip.map)
    
    ape::write.tree(gold_tree, gold.tree.path)
  }
}

rm_duplicate_suffix <- function(df, suffix="\\.x") {
  colnames(df) <- gsub(suffix, "", colnames(df))
  return(df)
}

species_for_tools <- function(df) {
  tools <- unique(df$Tool)
  tool2species <- list()
  for (tool in tools) {
    species <- unique(df %>% filter(Tool == tool & HasTree) %>% pull(Species))
    tool2species[[tool]] <- species
  }
  species.intersect <- c()
  for (tool in names(tool2species)) {
    species <- tool2species[[tool]]
    if (length(species.intersect) == 0) {
      species.intersect <- species
    } else {
      species.intersect <- intersect(species.intersect, species)
    }
  }
  return(species.intersect)
}

sliding_window_possible <- function(df, window_size=30) {
  has_one <- FALSE
  closest.df <- df %>% mutate(sim=rank(maxsim))
  for (tool in unique(df$Tool)) {
    tool.df <- closest.df %>% filter(Tool == tool) %>% arrange(sim)
    if (length(tool.df$Value) < window_size)
      next
    has_one <- TRUE
  }
  return(has_one)
}


sliding_window_monophyly <- function(df, window_size=30, theme=NULL) {
  browser()
  tool_list <- unique(df$Tool)
  closest.df <- df %>% mutate(sim=rank(maxsim))
  result.df <- NULL
  for (tool in tool_list) {
    tool.df <- closest.df %>% filter(Tool == tool) %>% arrange(sim)
    
    if (length(tool.df$Value) < window_size)
      next
    
    tool.df$sliding_window <- c(rep(NA, window_size/2), unlist(lapply(seq(nrow(tool.df)-window_size), FUN=function(i) {
      # Return the mean over a window
      start <- i
      end <- i + window_size
      start_maxsim <- tool.df[start, "maxsim"]
      end_maxsim <- tool.df[end, "maxsim"]
      window_mean <- mean(tool.df$Value[i:(i+window_size)])
      window_median <- mean(tool.df$Value[i:(i+window_size)])
      
      return(c(start, end, start_maxsim, end_maxsim, window_mean, ))
    })), rep(NA, window_size/2))
    
    if (all(is.na(result.df))) {
      result.df <- tool.df
    } else {
      result.df <- rbind(result.df, tool.df)
    }
  }
  if (is.null(result.df)) {
    return(NULL)
  }
  
  window.plot <- result.df %>%
    ggplot(aes(x=sim, y=sliding_window, color=Tool))
  if (!is.null(theme)) {
    window.plot <- window.plot + theme
  }
  window.plot <- window.plot +
    geom_line(key_glyph = "rect") +
    scale_x_continuous(labels=function(x){
      return(round(sort(result.df$maxsim)[x], digits=8))
    }, breaks=as.integer(seq(0,1, 0.1)*(nrow(result.df)-1))+1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    xlab("Similarity to closest neighboring genome in tree") +
    ylab("Monophyly per genome")
  return(window.plot)
}
