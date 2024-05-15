#' Human MGP
#' 
#' Plot genetic loadings from "runHumanMGP" output 
#' 
#' @param obj runHumanMGP object
#' @param comp Which PLS component to visualize (default = 1)
#' @param type Type of plot (options: "single", "sum1", "sum2", "weighted"; default = "single")
#' 
#' single = plot individual loadings per gene (or per genetic PC); sum1 = abs. sum of genetic PCs per gene (stacked); sum2 = abs. sum of genetic PCs per gene (individual); weighted = abs. sum of genetic PCs, weighted by number of PCs, per gene
#'
#' @author Hanne Hoskens
#' 
#' @export
plotMGPLoadings <- function(obj, comp, type, colormap){
  
  #  Required
  if (missing(obj)) {
    stop("Please provide runHumanMGP object")
  }
  
  # Optional
  if (missing(comp)) { comp = 1 }
  if (missing(type)) { type = "single" }
  
  max_comp = dim(obj$GeneLoadings)[2]
  if (length(comp)>1){
    stop("Please only select one PLS component at a time")
  }
  if (comp > max_comp) {
    cat("\033[33m", paste("Maximum number of components is ", max_comp, ". Plotting PLS1 instead.", sep = ""), "\033[0m", "\n")
    comp = 1
  }
  
  
  
  ngen = dim(obj$GeneLoadings)[1]
  tmp = unlist(strsplit(rownames(obj$GeneLoadings),"_",fixed=T))
  tmp = matrix(data=tmp, nrow=ngen, ncol=2, byrow=T)
  genes.name = tmp[,1]
  genes.pc = as.numeric(tmp[,2])
  

  if (type == "single"){
    # Remove "_1" (genes.pc) from gene names if npc==1
    if (sum(genes.pc==1)==ngen){ rownames(obj$GeneLoadings) = genes.name }
    
    pathway.loadings <- data.frame(gloadings = as.matrix(obj$GeneLoadings[,comp]), gnames = rownames(obj$GeneLoadings), pc = genes.pc)
      
    # Reorder
    bar_order <- pathway.loadings %>% 
      group_by(gnames) %>%
      summarise(idx = abs(gloadings)) %>%
      arrange(-idx) 
    
    pathway.loadings$gnames <- factor(x = pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)
    
  } else if (type == "sum1" | type == "sum2"){
    
    pathway.loadings <- data.frame(gloadings = as.matrix(obj$GeneLoadings[,comp]), gnames = genes.name, pc = genes.pc)
    
    # Sum of absolute values of gene loadings per gene
    genes = unique(genes.name)
    gloadings = matrix(0,length(genes),1)
    for (i in 1:length(genes)){
      ind = (genes.name %in% genes[i])
      gloadings[i] = sum(abs(obj$GeneLoadings[ind,comp]))
    }
    
    # Reorder
    ind = order(gloadings,decreasing = T)
    bar_order = tibble(gnames=genes[ind],test=gloadings[ind])
    
    pathway.loadings$gnames <- factor(x = pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)
    
  } else if (type == "weighted"){
    
    # Sum of absolute values of gene loading, weighted by number of PCs, per gene
    genes = unique(genes.name)
    genes.nPC = matrix(0,length(genes),1)
    gloadings = matrix(0,length(genes),1)
    for (i in 1:length(genes)){
      ind = (genes.name %in% genes[i])
      gloadings[i] = (sum(abs(obj$GeneLoadings[ind,comp]))) / (sum(ind))
      genes.nPC[i] = sum(ind)
    }
    
    pathway.loadings <- data.frame(gloadings = as.matrix(gloadings), gnames = genes, pc = genes.nPC)
    
    # Reorder
    bar_order <- pathway.loadings %>% 
      group_by(gnames) %>%
      summarise(ind = abs(gloadings)) %>%
      arrange(-ind) 
    
    pathway.loadings$gnames <- factor(x = pathway.loadings$gnames, levels = lapply(bar_order, as.character)$gnames)
  }
    
  
  # Plot gene loadings
  if (type == "single"){
    p = ggplot() +
      geom_bar(data = pathway.loadings,
             aes(x = gnames, y = abs(gloadings)),
             stat = "identity",
             width = .75,
             position=position_dodge())

  } else if (type == "sum1"){
    p = ggplot() + 
      geom_bar(data = pathway.loadings,
             aes(x = gnames, y = abs(gloadings), fill = as.factor(pc)),
             stat = "identity",
             width = .75,
             position = position_stack(reverse=T))
    
  } else if (type == "sum2"){
    p = ggplot() + 
      geom_bar(data = pathway.loadings,
             aes(x = gnames, y = abs(gloadings), fill = as.factor(pc)),
             stat = "identity",
             width = .75,
             position=position_dodge(preserve = "total"))
    
  } else if (type == "weighted"){
    p = ggplot() +
      geom_bar(data = pathway.loadings,
               aes(x = gnames, y = abs(gloadings), fill = as.factor(pc)),
               stat = "identity",
               width = .75,
               position=position_dodge())
  }
    
    
  # Add axes
  p = p +
    labs(x = "",
         y = "Gene loadings",
         #title = toupper(obj$GOterm),
         subtitle = paste(obj$Cohort, ", ", "PLS ", comp, sep = ""),
         fill = "number of composite PCs") +
    theme(text = element_text(size=6),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "plain"),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 8),
          axis.title.y = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "plain", hjust = .5),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10)) +
    scale_fill_manual(values=as.vector(ocean.delta(length(levels(as.factor(pathway.loadings$pc))))))

    if (length(obj$GOterm)>1){ p = p + labs(title = "Custom gene list") 
    } else { p = p + labs(title = toupper(obj$GOterm)) } 
  
    if (type == "single"){ p = p + theme(legend.position = "none") }

  p
      
  
}
