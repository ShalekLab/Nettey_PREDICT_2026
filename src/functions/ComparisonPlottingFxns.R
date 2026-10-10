require(Hmisc)
require(ggplot2)
require(dplyr)
require(tidyr)
require(forcats)
require(Seurat)
require(ggpubr)
require(ggrepel)

PlotBoxplotsProp <- function(df, ann.col, relative.to = NULL,
                             ncols = 4,
                             plotTitle = NULL, subTitle = NULL, 
                             filterCol = NULL, filterVal = NULL,
                             xLab = NULL, yLab = "Percentage of Cells",
                             sample.col = "biopsy",
                             showSubTitle = TRUE,
                             makeFacets = TRUE,
                             minDenomFreq = -Inf, denom.col = NULL, addCutoffToSubtitle = TRUE,
                             xVar = "diagSub", fillVar = "diagSub", fillVector = diagSubColVec,
                             doTests = TRUE, testMethod = "wilcox.test", testLabel = "p.signif", testComp = diagSub.comp, 
                             vjust = 1, testSize = 5,
                             dotSize = 2.5, jitterwidth = 0.75, seed.use = 220822,
                             scaleDotsBySize = TRUE, dotMin = 1, dotMax = 5, sizeTrans = "identity", alphaVal = NULL,
                             botFreqQ = 0, topFreqQ = 1.0,
                             showXVal = TRUE, showYVal = TRUE,
                             legendPos = "bottom", errorbarCol = "black",
                             gridSpace = "free_x", gridScales = "free_x", wrapScales = "free",
                             fillLegend = NULL, fillTitle = NULL, sizeLegend = TRUE, sizeTitle = NULL,
                             colorLegend = NULL,
                             showBoxplots = FALSE, stackedBar = FALSE,
                             titleSize = 18, subtitleSize = 14, 
                             axisXsize = 11, axisXangle = 0, axisXhjust = 0.5,
                             axisYsize = 12, titleYsize = 16, legendSize = 12, stripSize = 12,
                             errorWidth = 0.4, meanWidth = 0.6, meanLinewidth = 0.3,
                             exportDF = FALSE 
){
  
  df.orig <- df
  
  if(is.null(denom.col)){
    if(!is.null(relative.to)){
      denom.col <- relative.to
    } else {
      denom.col <- sample.col
    }
  }
  
  ann.denom.col <- paste(ann.col, denom.col, sep = "__")
  sample.denom.col <- paste(sample.col, denom.col, sep = "__")
  sample.ann.denom.col <- paste(sample.col, ann.col, denom.col, sep = "__")
  df <- df %>% 
    mutate(!!ann.denom.col:=paste(!!!rlang::syms(ann.col), !!!rlang::syms(denom.col),sep="__")) %>% 
    mutate(!!ann.denom.col:=as.factor(!!!rlang::syms(ann.denom.col))) %>%
    mutate(!!sample.denom.col:=paste(!!!rlang::syms(sample.col), !!!rlang::syms(denom.col),sep="__")) %>% 
    mutate(!!sample.denom.col:=as.factor(!!!rlang::syms(sample.denom.col))) 
  
  allxVars <- df %>% 
    distinct(!!!rlang::syms(sample.col), .keep_all = T) %>% 
    mutate(!!xVar := as.character(.[[xVar]])) %>%
    select({{xVar}}, {{sample.col}})
  xVarLevels <- df %>%
    pull({{xVar}}) %>%
    levels
  
  allfillVars <- df %>% 
    distinct(!!!rlang::syms(sample.col), .keep_all = T) %>% 
    mutate(!!fillVar := as.character(.[[fillVar]])) %>%
    select({{fillVar}}, {{sample.col}})
  fillVarLevels <- df %>%
    pull({{fillVar}}) %>%
    levels
  
  allAnnDenomVars <- df %>% 
    distinct(!!!rlang::syms(ann.denom.col)) %>% 
    pull({{ann.denom.col}})
  
  denomDF <- df %>% 
    add_count(!!!rlang::syms(sample.denom.col), name = "denominator.freq") %>% 
    complete(!!!rlang::syms(sample.col),  !!!rlang::syms(denom.col)) %>%
    distinct(!!!rlang::syms(sample.denom.col), .keep_all = T) %>% 
    mutate(!!sample.denom.col:=paste(!!!rlang::syms(sample.col), !!!rlang::syms(denom.col),sep="__")) %>% 
    mutate(!!sample.denom.col:=as.factor(!!!rlang::syms(sample.denom.col))) %>% 
    select({{sample.denom.col}}, denominator.freq)
  
  df <- df %>%
    add_count(!!!rlang::syms(sample.col), !!!rlang::syms(ann.col), name = "annotation.freq") %>% 
    group_by_at(vars(all_of(sample.col), all_of(relative.to)), .drop = FALSE) %>% 
    mutate(prop = annotation.freq/n()) %>%
    ungroup 
  
  if(!is.null(filterCol) & !is.null(filterVal)){
    df <- df %>% 
      filter(.[[filterCol]] %in% filterVal)
    if(!is.null(ann.col)){df <- df %>% mutate(!!ann.col := droplevels(.[[ann.col]]))}
  } else if ((!is.null(filterCol) & is.null(filterVal)) | (is.null(filterCol) & !is.null(filterVal)) ){
    error("Filtering not performed - must specify both the column to filter on and the value to use")
  }    
  
  df <- df %>% 
    complete(!!!rlang::syms(sample.col), !!!rlang::syms(ann.col), !!!rlang::syms(denom.col), fill = list(annotation.freq = 0, prop = 0)) %>%
    mutate(!!sample.ann.denom.col:=paste(!!!rlang::syms(sample.col), !!!rlang::syms(ann.col),!!!rlang::syms(denom.col),sep="__")) %>% 
    mutate(!!sample.ann.denom.col:=as.factor(!!!rlang::syms(sample.ann.denom.col))) %>% 
    mutate(!!sample.denom.col:=paste(!!!rlang::syms(sample.col), !!!rlang::syms(denom.col),sep="__")) %>% 
    mutate(!!sample.denom.col:=as.factor(!!!rlang::syms(sample.denom.col))) %>% 
    distinct(!!!rlang::syms(sample.ann.denom.col), .keep_all=T)
  
  df <- df %>%
    left_join(denomDF, by = sample.denom.col) %>%
    filter(denominator.freq >= minDenomFreq)
  
  df <- df %>%
    separate(col = {{sample.ann.denom.col}}, into = c("split.sample.col", "split.ann.col", "split.denom.col"), sep = "__") %>% 
    mutate(comb.ann.denom.col = paste(split.ann.col, split.denom.col,  sep = "__")) %>% 
    filter(comb.ann.denom.col %in% allAnnDenomVars)
  
  if(xVar == fillVar){
    df <- df %>%
      select(-!!xVar) %>% 
      left_join(allxVars, by = sample.col) 
    if(!is.null(xVarLevels)){df <- df %>% mutate(!!xVar := fct_relevel(.[[xVar]], xVarLevels))}
  } else {
    df <- df %>%
      select(-!!xVar, -!!fillVar) %>% 
      left_join(allxVars, by = sample.col) %>% 
      left_join(allfillVars, by = sample.col) 
    if(!is.null(xVarLevels)){df <- df %>% mutate(!!xVar := fct_relevel(.[[xVar]], xVarLevels))}
    if(!is.null(fillVarLevels)){df <- df %>% mutate(!!fillVar := fct_relevel(.[[fillVar]], fillVarLevels))}
    
  }
  
  qBot <- df %>% pull(annotation.freq) %>% quantile(botFreqQ)
  qTop <- df %>% pull(annotation.freq) %>% quantile(topFreqQ)
  
  if(exportDF){return(df %>% select(!!sample.col, !!ann.col, !!relative.to, !!filterCol, !!xVar, !!fillVar,
                                    annotation.freq, denominator.freq, prop))}
  
  if(stackedBar){
    {if(is.null(fillLegend)) fillLegend <- TRUE}
    ggplotobj <- ggplot(df, aes(x = get(sample.col), y = prop, fill = get(ann.col))) +
      geom_col(color = "black", position="fill", linewidth=0.25) +
      {if(makeFacets) facet_grid(~get(xVar), scales = gridScales, space = gridSpace)}
    
    if(is.null(subTitle) & !is.null(relative.to) & showSubTitle){
      if(!is.null(filterCol) & !is.null(filterVal)){
        subTitle <- sprintf("Relative to %s | Filtered to %s", relative.to, filterVal)
      } else {
        subTitle <- sprintf("Relative to %s", relative.to)
      }
      
    }
    
  } else {
    if(!showBoxplots & is.null(alphaVal)) {
      alphaVal <- 0.5
    } else {
      alphaVal <- alphaVal
    }
    {if(is.null(fillLegend)) fillLegend <- FALSE}
    {if(is.null(colorLegend)) colorLegend <- FALSE}
    ggplotobj <- ggplot(df, aes(x = get(xVar), y = prop, fill = get(fillVar))) +
      {if(showBoxplots) geom_boxplot(outlier.shape = NA)} +
      {if(!scaleDotsBySize) geom_point(size = dotSize, shape = 21, alpha = alphaVal,
                                       position = position_jitterdodge(jitter.width = jitterwidth, seed = seed.use))} +
      {if(scaleDotsBySize) geom_point(aes(size = annotation.freq), shape = 21, alpha = alphaVal,
                                      position = position_jitterdodge(jitter.width = jitterwidth, seed = seed.use))} +
      {if(!showBoxplots) 
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = errorWidth, color = errorbarCol) 
      } +
      {if(!showBoxplots) 
        stat_summary(fun = "mean", geom = "crossbar", width = meanWidth, linewidth = meanLinewidth, 
                     aes(color = get(fillVar)))
      } +
      facet_wrap(~get(ann.col), scales = wrapScales, ncol = ncols) +        
      scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax), limits = c(qBot, qTop)) +
      scale_fill_manual(values = fillVector) +
      scale_color_manual(values = fillVector) +
      {if(doTests) stat_compare_means(comparisons = testComp, 
                                      method = testMethod, label = testLabel, vjust = vjust, size = testSize)
      }
    
    if(is.null(subTitle) & !is.null(relative.to) & showSubTitle){
      subTitle <- sprintf("Relative to %s", relative.to)
    }
  }
  
  if(minDenomFreq != -Inf & addCutoffToSubtitle){
    cellStr <- ifelse(minDenomFreq == 1, "cell", "cells")
    if(is.null(subTitle)){
      subTitle <- sprintf("\u2265 %i %s", minDenomFreq, cellStr)
    } else {
      subTitle <- paste0(subTitle, sprintf("\n\u2265 %i %s", minDenomFreq, cellStr))
    }
  }    
  
  ggplotobj +
    scale_y_continuous(labels = percent) +
    {if(!showXVal) guides(x = "none")} +
    {if(!showYVal) guides(y = "none")} +
    {if(!fillLegend) guides(fill = "none")} +
    {if(!colorLegend) guides(color = "none")} +
    {if(!sizeLegend) guides(size = "none")} +
    labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, title = plotTitle, subtitle = subTitle) +
    theme_light() +
    theme(
      plot.title = element_text(size = titleSize, face = "bold"),
      plot.subtitle = element_text(size = subtitleSize, face = "bold"),
      axis.text.x = element_text(size = axisXsize, face = "bold", angle = axisXangle, hjust=axisXhjust),
      axis.text.y = element_text(size = axisYsize, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = titleYsize, face = "bold"),
      legend.position = legendPos,
      legend.text = element_text(size = legendSize, face = "bold"),
      legend.title = element_text(size = legendSize, face = "bold"),
      strip.text = element_text(size = stripSize, face = "bold", color = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    )
  
}


PlotBoxplotsFeatures <- function(df, features, relative.to = NULL,
                                 ncols = 4, slot = "data", assay = "RNA", 
                                 plotTitle = NULL, subTitle = NULL, 
                                 filterCol = NULL, filterVal = NULL,
                                 xLab = NULL, yLab = NULL,
                                 sample.col = "biopsy", observation.name = "barcode", seurat = NULL,
                                 minDenomFreq = -Inf,
                                 compare.to = NULL,
                                 stat.use = "mean", expression.threshold = 0,
                                 showCoef = TRUE, showReg = TRUE, groupReg = T, r2level = 0.999,
                                 xVar = "diagSub", fillVar = "diagSub", fillVector = diagSubColVec,
                                 order.by = TRUE,
                                 plotFacets = TRUE, facetType = "wrap", facet_scales = "free", facet_space = "fixed", swap_facets = F,
                                 doTests = TRUE, testMethod = "wilcox.test", testLabel = "p.signif", testComp = diagSub.comp, regSize = 4, facet_dir = "h",
                                 vjust = 1, testSize = 5,
                                 dotSize = 2.5, alphaVal = 0.25, jitterwidth = 0.75, seed.use = 220822,
                                 scaleDotsBySize = TRUE, dotMin = 1, dotMax = 5, sizeTrans = "identity", 
                                 botFreqQ = 0, topFreqQ = 1.0,
                                 fillLegend = FALSE, fillTitle = NULL, sizeLegend = TRUE, sizeTitle = NULL,
                                 showBoxplots = TRUE, 
                                 plotObservations = FALSE, removeZeros = FALSE,
                                 titleSize = 18, subtitleSize = 14, 
                                 axisXsize = 12, axisXangle = 0, axisXhjust = 0.5,
                                 axisYsize = 12, titleXsize = 16, titleYsize = 16, legendSize = 12, stripSize = 12,
                                 exportSampleDF = FALSE, exportObservationDF = FALSE
){
  
  if(!is.null(compare.to)){
    if(!(compare.to %in% features)){
      features <- c(features, compare.to) 
    }
  } 
  
  if(!(all(features %in% colnames(df)))){
    geneMat <- GetAssayData(get(seurat), slot = slot, assay = assay)       
    featuresInMat <- features[which(features %in% rownames(geneMat))]
    geneMatTrim <- geneMat[featuresInMat,]
    
    if(n_distinct(featuresInMat) == 1){
      geneDF <- tibble(!!observation.name := names(geneMatTrim))
      geneDF <- geneDF %>% mutate(!!featuresInMat := geneMatTrim)
    } else {
      geneDF <- geneMatTrim %>% as.matrix %>% t %>% as.data.frame %>% rownames_to_column(observation.name) 
    }
    
    if(!(observation.name %in% colnames(df))){df <- df %>% rownames_to_column(observation.name)}
    geneDF <- geneDF %>% filter(.[[observation.name]] %in% (df %>% pull({{observation.name}})))
    df <- left_join(geneDF, df, by = observation.name)
  }
  
  if(!is.null(filterCol) & !is.null(filterVal)){
    df <- df %>% 
      filter(.[[filterCol]] %in% filterVal)
    
    if(!is.null(relative.to)){
      df <- df %>% 
        mutate(!!relative.to := droplevels(.[[relative.to]]))
      
      if(relative.to == filterCol) {
        df <- df %>%
          mutate(!!relative.to := fct_relevel(.[[relative.to]], filterVal))
      }
      
    }
    
  } else if ((!is.null(filterCol) & is.null(filterVal)) | (is.null(filterCol) & !is.null(filterVal)) ){
    stop("Filtering not performed - must specify both the column to filter on and the value to use")
  }    
  
  df <- df %>% 
    add_count(!!!rlang::syms(sample.col), !!!rlang::syms(relative.to), name = "denominator.freq") %>%
    filter(denominator.freq > minDenomFreq)    
  
  qBot <-df %>% pull(denominator.freq) %>% quantile(botFreqQ)
  qTop <- df %>% pull(denominator.freq) %>% quantile(topFreqQ)
  
  df <- df %>%
    pivot_longer(cols = all_of(setdiff(features,compare.to)), names_to = "featureVar", values_to = "featureVal")
  
  {if(order.by) df <- df %>% mutate(featureVar = fct_relevel(featureVar, setdiff(features,compare.to)))}    
  
  if(exportObservationDF){return(df)}    
  
  if(plotObservations){
    if(!is.null(compare.to)){
      {if(removeZeros) df <- df %>% filter(.[[compare.to]] > 0 & featureVal > 0)}
      {if(is.null(xLab)) xLab <- sprintf("%s Expression", compare.to)}
      {if(is.null(yLab)) yLab <- "Expression"}
      ggplot(df, aes(x = get(compare.to), y = featureVal, fill = get(fillVar))) +
        geom_point(size = dotSize, alpha = alphaVal, shape = 21) +
        {if(showReg) stat_poly_line(method = "lm")} +
        {if(showCoef) stat_poly_eq(method = "lm", use_label(c("R2", "P")), size = regSize)} +
        {if(!is.null(relative.to)) facet_wrap(get(relative.to)~featureVar, scales = "free", ncol = ncols)} +
        {if(is.null(relative.to)) facet_wrap(~featureVar, scales = "free", ncol = ncols)} +
        {if(!is.null(fillVector)) scale_fill_manual(values = fillVector)} + 
        scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax), limits = c(qBot, qTop)) +
        {if(!fillLegend) guides(fill = "none")} +
        {if(!sizeLegend) guides(fill = "none")} +
        labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, title = plotTitle, subtitle = subTitle) +
        theme_light() +
        theme(
          plot.title = element_text(size = titleSize, face = "bold"),
          plot.subtitle = element_text(size = subtitleSize, face = "bold"),
          axis.text.x = element_text(size = axisXsize, face = "bold", angle = axisXangle, hjust = axisXhjust),
          axis.text.y = element_text(size = axisYsize, face = "bold"),
          axis.title.x = element_text(size = titleXsize, face = "bold"),
          axis.title.y = element_text(size = titleYsize, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = legendSize, face = "bold"),
          legend.text = element_text(size = legendSize, face = "bold"),
          strip.text = element_text(size = stripSize, face = "bold", color = "black"),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        )            
    } else {
      {if(removeZeros) df <- df %>% filter(featureVal > 0)}
      {if(is.null(yLab)) yLab <- "Expression"}
      ggplot(df, aes(x = get(xVar), y = featureVal, fill = get(fillVar))) +
        {if(showBoxplots) geom_boxplot(outlier.shape = NA)} +
        geom_point(size = dotSize, shape = 21, alpha = alphaVal, 
                   position = position_jitterdodge(jitter.width = jitterwidth, seed = seed.use)) +
        {if(!is.null(relative.to)) facet_wrap(get(relative.to)~featureVar, scales = "free", ncol = ncols)} +
        {if(is.null(relative.to)) facet_wrap(~featureVar, scales = "free", ncol = ncols)} +
        {if(!is.null(fillVector)) scale_fill_manual(values = fillVector)} + 
        scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax), limits = c(qBot, qTop)) +
        {if(!fillLegend) guides(fill = "none")} +
        {if(!sizeLegend) guides(fill = "none")} +
        labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, title = plotTitle, subtitle = subTitle) +
        {if(doTests) stat_compare_means(comparisons = testComp, 
                                        method = testMethod, label = testLabel, vjust = vjust, size = testSize)
        } +      
        theme_light() +
        theme(
          plot.title = element_text(size = titleSize, face = "bold"),
          plot.subtitle = element_text(size = subtitleSize, face = "bold"),
          axis.text.x = element_text(size = axisXsize, face = "bold", angle = axisXangle, hjust = axisXhjust),
          axis.text.y = element_text(size = axisYsize, face = "bold"),
          axis.title.x = element_text(size = titleXsize, face = "bold"),
          axis.title.y = element_text(size = titleYsize, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = legendSize, face = "bold"),
          legend.text = element_text(size = legendSize, face = "bold"),
          strip.text = element_text(size = stripSize, face = "bold", color = "black"),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        )
      
    }
    
  } else {
    
    
    if(!is.null(compare.to)){
      features <- features[-which(features %in% compare.to)]
      df <- df %>%
        mutate(percentTest = .[[compare.to]] > expression.threshold) %>%
        group_by_at(vars(all_of(sample.col), all_of(relative.to)), .drop = FALSE) %>%
        mutate(sumComp = sum(!!!rlang::syms(compare.to), na.rm = T)) %>%
        mutate(meanComp = mean(!!!rlang::syms(compare.to), na.rm = T)) %>%
        mutate(medianComp = median(!!!rlang::syms(compare.to), na.rm = T)) %>%
        mutate(percentComp = sum(percentTest)/n()) 
    }
    
    df <- df %>%
      group_by_at(vars(all_of(sample.col), all_of(relative.to), all_of("featureVar")), .drop = FALSE) %>%
      mutate(sumVal = sum(featureVal, na.rm = T)) %>%
      mutate(meanVal = mean(featureVal, na.rm = T)) %>%
      mutate(medianVal = median(featureVal, na.rm = T)) %>%
      mutate(percentVal = sum(featureVal > expression.threshold)/n()) %>%
      ungroup %>%
      distinct(!!!rlang::syms(sample.col), !!!rlang::syms(relative.to), featureVar, .keep_all=T)    
    
    if(exportSampleDF){return(df %>% relocate(!!sample.col, 
                                              !!relative.to, !!filterCol, !!xVar, !!fillVar,
                                              denominator.freq, meanVal, medianVal, .after = last_col()))}
    
    if(stat.use == "mean"){
      stat.val <- "meanVal"
      stat.comp <- "meanComp"
    } else if (stat.use == "sum"){
      stat.val <- "sumVal"
      stat.comp <- "sumComp"
    } else if (stat.use == "median"){
      stat.val <- "medianVal"
      stat.comp <- "medianComp"
    } else if (stat.use == "percent"){
      stat.val <- "percentVal"
      stat.comp <- "percentComp"
    } else {
      stop("Unknown Statistical Summary Method")
    }
    
    if(is.null(yLab)){yLab <- sprintf("%s Expression", stringr::str_to_title(stat.use))}
    
    if(!is.null(compare.to)){
      {if(is.null(xLab)) xLab <- sprintf("%s %s Expression", stringr::str_to_title(stat.use), compare.to)}
      
      if(!groupReg){
        p <- ggplot(df, aes(x = get(stat.comp), y = get(stat.val)))
      } else {
        p <- ggplot(df, aes(x = get(stat.comp), y = get(stat.val), group = get(fillVar)))
      }
      p +
        {if(!scaleDotsBySize) geom_point(size = dotSize, shape = 21, aes(fill = get(fillVar)))} +
        {if(scaleDotsBySize) geom_point(aes(size = denominator.freq, fill = get(fillVar)), shape = 21)} +
        {if(showReg) stat_poly_line(method = "lm")} +
        {if(showCoef) stat_poly_eq(method = "lm", use_label(c("R2", "P")), size = regSize, rsquared.conf.level = r2level)} +
        {if(!is.null(relative.to) & plotFacets & facetType == "wrap") facet_wrap(get(relative.to)~featureVar, scales = facet_scales, ncol = ncols, dir = facet_dir)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "wrap") facet_wrap(~featureVar, scales = facet_scales, ncol = ncols, dir = facet_dir)} +
        {if(!is.null(relative.to) & plotFacets & facetType == "grid" & !swap_facets) facet_grid(featureVar~get(relative.to), scales = facet_scales, space = facet_space)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "grid" & !swap_facets) facet_grid(~featureVar, scales = facet_scales, space = facet_space)} +   
        {if(!is.null(relative.to) & plotFacets & facetType == "grid" & swap_facets) facet_grid(get(relative.to) ~ featureVar, scales = facet_scales, space = facet_space)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "grid" & swap_facets) facet_grid(featureVar ~ ., scales = facet_scales, space = facet_space)} +   
        {if(!is.null(fillVector)) scale_fill_manual(values = fillVector)} + 
        scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax), limits = c(qBot, qTop)) +
        {if(stat.use == "percent") scale_x_continuous(labels = percent)} +
        {if(stat.use == "percent") scale_y_continuous(labels = percent)} +
        {if(!fillLegend) guides(fill = "none")} +
        {if(!sizeLegend) guides(size = "none")} +
        labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, title = plotTitle, subtitle = subTitle) +
        theme_light() +
        theme(
          plot.title = element_text(size = titleSize, face = "bold"),
          plot.subtitle = element_text(size = subtitleSize, face = "bold"),
          axis.text.x = element_text(size = axisXsize, angle = axisXangle, hjust = axisXhjust, face = "bold"),
          axis.text.y = element_text(size = axisYsize, face = "bold"),
          axis.title.x = element_text(size = titleXsize, face = "bold"),
          axis.title.y = element_text(size = titleYsize, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = legendSize, face = "bold"),
          legend.text = element_text(size = legendSize, face = "bold"),
          strip.text = element_text(size = stripSize, face = "bold", color = "black"),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        )            
    } else {
      
      ggplot(df, aes(x = get(xVar), y = get(stat.val), fill = get(fillVar))) +
        {if(showBoxplots) geom_boxplot(outlier.shape = NA)} +
        {if(!scaleDotsBySize) geom_point(size = dotSize, shape = 21, 
                                         position = position_jitterdodge(jitter.width = jitterwidth, seed = seed.use))} +
        {if(scaleDotsBySize) geom_point(aes(size = denominator.freq), shape = 21, 
                                        position = position_jitterdodge(jitter.width = jitterwidth, seed = seed.use))} +
        {if(!is.null(relative.to) & plotFacets & facetType == "wrap") facet_wrap(get(relative.to)~featureVar, scales = facet_scales, ncol = ncols, dir = facet_dir)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "wrap") facet_wrap(~featureVar, scales = facet_scales, ncol = ncols, dir = facet_dir)} +
        {if(!is.null(relative.to) & plotFacets & facetType == "grid" & !swap_facets) facet_grid(featureVar~get(relative.to), scales = facet_scales, space = facet_space)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "grid" & !swap_facets) facet_grid(~featureVar, scales = facet_scales, space = facet_space)} +   
        {if(!is.null(relative.to) & plotFacets & facetType == "grid" & swap_facets) facet_grid(get(relative.to) ~ featureVar, scales = facet_scales, space = facet_space)} +
        {if((is.null(relative.to) | !plotFacets) & facetType == "grid" & swap_facets) facet_grid(featureVar ~ ., scales = facet_scales, space = facet_space)} +   
        {if(!is.null(fillVector)) scale_fill_manual(values = fillVector)} + 
        scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax), limits = c(qBot, qTop)) +
        {if(stat.use == "percent") scale_y_continuous(labels = percent)} +
        {if(!fillLegend) guides(fill = "none")} +
        {if(!sizeLegend) guides(fill = "none")} +
        labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, title = plotTitle, subtitle = subTitle) +
        {if(doTests) stat_compare_means(comparisons = testComp, 
                                        method = testMethod, label = testLabel, vjust = vjust, size = testSize)
        } +      
        theme_light() +
        theme(
          plot.title = element_text(size = titleSize, face = "bold"),
          plot.subtitle = element_text(size = subtitleSize, face = "bold"),
          axis.text.x = element_text(size = axisXsize, face = "bold", angle = axisXangle, hjust = axisXhjust),
          axis.text.y = element_text(size = axisYsize, face = "bold"),
          axis.title.x = element_text(size = titleXsize, face = "bold"),
          axis.title.y = element_text(size = titleYsize, face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = legendSize, face = "bold"),
          legend.text = element_text(size = legendSize, face = "bold"),
          strip.text = element_text(size = stripSize, face = "bold", color = "black"),
          strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
        )
      
    }
    
  }
  
}


SenderReceiver <- function(df.orig, 
                           features.sender = NULL,
                           features.receiver = NULL,
                           senderCol = NULL, senderVal = NULL,
                           receiverCol = NULL, receiverVal = NULL,
                           ncols = 4, slot = "data", assay = "RNA", 
                           plotTitle = NULL, subTitle = NULL, 
                           filterCol = NULL, filterVal = NULL,
                           senderProp = FALSE, senderRelative = NULL, senderDenom = NULL, senderMinDenomFreq = -Inf, 
                           receiverProp = FALSE, receiverRelative = NULL, receiverDenom = NULL, receiverMinDenomFreq = -Inf,
                           xLab = NULL, yLab = NULL,
                           sample.col = "biopsy", observation.name = "barcode", seurat = "seuobj",
                           minSenderFreq = -Inf, minReceiverFreq = -Inf,
                           stat.use = "mean", expression.threshold = 0,
                           showCoef = TRUE, showReg = TRUE, regSE = TRUE, regColour = "#3366FF", matchRegColour = FALSE, 
                           regSize = 4, fillReg = TRUE, groupReg = TRUE, r2level = 0.999,
                           facetType = "grid", facet_scales = "free", facet_space = "fixed", facet_dir = "h", swap_facets = F,
                           swapSenderOrder = FALSE, swapReceiverOrder = FALSE,
                           xVar = "diagSub", fillVar = "diagSub", fillVector = diagSubColVec,
                           order.by = TRUE,
                           doTests = TRUE, testMethod = "wilcox.test", testLabel = "p.signif", testComp = diagSub.comp, 
                           vjust = 1, testSize = 5,
                           dotSize = 2.5, alphaVal = 0.25, jitterwidth = 0.75, seed.use = 220822,
                           scaleDotsBySize = TRUE, dotMin = 1, dotMax = 5, sizeTrans = "identity", 
                           botFreqQ = 0, topFreqQ = 1.0,
                           fillLegend = TRUE, fillTitle = NULL, sizeLegend = TRUE, sizeTitle = NULL, colorLegend = TRUE, colorTitle = NULL,
                           legendPosition = "right",
                           showBoxplots = TRUE, 
                           labelPoints = FALSE, labelSize = 3, labelFace = "bold",
                           plotObservations = FALSE, removeZeros = FALSE,
                           titleSize = 18, subtitleSize = 14, 
                           axisXsize = 12, axisYsize = 12, titleXsize = 16, titleYsize = 16, legendSize = 12, stripSize = 12,
                           exportSampleDF = FALSE, exportObservationDF = FALSE
){
  
  # Create new dataframe to work on
  df <- df.orig
  
  # Remove accidentally provided features, if proportion is desired
  if(senderProp){features.sender <- NULL}
  if(receiverProp){features.receiver <- NULL}
  
  # Ensure that the correct filter variables are set
  if (is.null(senderCol) | is.null(senderVal) | is.null(receiverCol) | is.null(receiverVal) ){
    stop("Sender/Receiver filter not performed - must specify both the column to filter on and the value to use")
  } else {
    
    # Add the cell ID to the dataframe if not already present, (assumes cell IDs are rownames, i.e. Seurat-style)
    if(!(observation.name %in% colnames(df))){df <- df %>% rownames_to_column(observation.name)}
    
    # Subset to only celltypes of interest.
    df <- df %>% 
      filter(.[[senderCol]] %in% senderVal | .[[receiverCol]] %in% receiverVal) %>%
      mutate(!!senderCol := droplevels(.[[senderCol]])) %>%
      mutate(!!senderCol := fct_relevel(.[[senderCol]], senderVal)) %>%
      mutate(!!receiverCol := droplevels(.[[receiverCol]])) %>%
      mutate(!!receiverCol := fct_relevel(.[[receiverCol]], receiverVal)) %>% 
      add_count(!!!rlang::syms(sample.col), !!!rlang::syms(senderCol), name = "sender.freq") %>%
      add_count(!!!rlang::syms(sample.col), !!!rlang::syms(receiverCol), name = "receiver.freq")
  }
  
  # Collect all features of interest, if present
  features.all <- unique(c(features.sender, features.receiver))
  
  if(!is.null(features.all)){
    
    # Some features will be in the meta data. Others will be in the assay data. This will look for those in assay data
    if(!(all(features.all %in% colnames(df)))){
      geneMat <- GetAssayData(get(seurat), slot = slot, assay = assay)       
      featuresInMat <- features.all[which(features.all %in% rownames(geneMat))]
      geneMatTrim <- geneMat[featuresInMat,]
      
      # Must handle things slightly different if only provided 1 feature (creates a vector) vs. multiple features (data frame)
      if(n_distinct(featuresInMat) == 1){
        geneDF <- tibble(!!observation.name := names(geneMatTrim))
        geneDF <- geneDF %>% mutate(!!featuresInMat := geneMatTrim)
      } else {
        geneDF <- geneMatTrim %>% as.matrix %>% t %>% as.data.frame %>% rownames_to_column(observation.name) 
      }
      
      # Add columns with feature values to main dataframe (right join to get feature values to beginning of dataframe, but subsetted cells only)
      df <- right_join(geneDF, df, by = observation.name)
    }
  }
  
  # If desired, export individual observation (cell) based dataframe with features included
  if(exportObservationDF){return(df)}
  
  # set which statistic will be used for the x and y summation, if applicable
  if(stat.use == "mean"){
    stat.val <- "meanVal"
  } else if (stat.use == "sum"){
    stat.val <- "sumVal"
  } else if (stat.use == "median"){
    stat.val <- "medianVal"
  } else if (stat.use == "percent"){
    stat.val <- "percentVal"
  } else {
    stop("Unknown Statistical Summary Method")
  }
  
  
  if(!is.null(features.sender)){
    
    df.sender <- df %>%
      filter(.[[senderCol]] %in% senderVal & sender.freq >= minSenderFreq) %>%
      pivot_longer(cols = all_of(features.sender), names_to = "featureVar", values_to = "featureVal") %>%
      group_by_at(vars(all_of(sample.col), all_of(senderCol), all_of("featureVar")), .drop = FALSE) %>%
      mutate(sumVal = sum(featureVal, na.rm = T)) %>%
      mutate(meanVal = mean(featureVal, na.rm = T)) %>%
      mutate(medianVal = median(featureVal, na.rm = T)) %>%
      mutate(percentVal = sum(featureVal > expression.threshold)/n()) %>%
      ungroup %>%
      distinct(!!!rlang::syms(sample.col), !!!rlang::syms(senderCol), featureVar, .keep_all=T) %>%
      pivot_wider(
        id_cols = c(all_of(sample.col), all_of(xVar), all_of(fillVar)),
        names_from = c(all_of(senderCol), featureVar),
        values_from = all_of(stat.val) 
      )
    
    senderColName <- as.vector(outer(senderVal, features.sender, paste, sep="_"))
    
  } else {
    
    df.sender <- df %>% 
      filter(.[[senderCol]] %in% senderVal & sender.freq >= minSenderFreq) %>%
      ungroup %>%
      distinct(!!!rlang::syms(sample.col), !!!rlang::syms(senderCol), .keep_all=T) %>% 
      select(c(all_of(sample.col), all_of(xVar), all_of(fillVar), all_of(senderCol))) %>% 
      rename_with(.cols = all_of(c(senderCol)), .fn = ~paste0(., ".sender"))
    
    senderColName <- paste0(senderCol, ".sender")
    
  }
  
  if(!is.null(features.receiver)){
    
    df.receiver <- df %>%
      filter(.[[receiverCol]] %in% receiverVal & receiver.freq >= minReceiverFreq) %>%
      pivot_longer(cols = all_of(features.receiver), names_to = "featureVar", values_to = "featureVal") %>%
      group_by_at(vars(all_of(sample.col), all_of(receiverCol), all_of("featureVar")), .drop = FALSE) %>%
      mutate(sumVal = sum(featureVal, na.rm = T)) %>%
      mutate(meanVal = mean(featureVal, na.rm = T)) %>%
      mutate(medianVal = median(featureVal, na.rm = T)) %>%
      mutate(percentVal = sum(featureVal > expression.threshold)/n()) %>%
      ungroup %>%
      distinct(!!!rlang::syms(sample.col), !!!rlang::syms(receiverCol), featureVar, .keep_all=T) %>%
      pivot_wider(
        id_cols = c(all_of(sample.col), all_of(xVar), all_of(fillVar)),
        names_from = c(all_of(receiverCol), featureVar),
        values_from = all_of(stat.val) 
      )
    
    receiverColName <- as.vector(outer(receiverVal, features.receiver, paste, sep="_"))
    
  } else {
    
    df.receiver <- df %>% 
      filter(.[[receiverCol]] %in% receiverVal & receiver.freq >= minReceiverFreq) %>%
      ungroup %>%
      distinct(!!!rlang::syms(sample.col), !!!rlang::syms(receiverCol), .keep_all=T) %>% 
      select(c(all_of(sample.col), all_of(xVar), all_of(fillVar), all_of(receiverCol))) %>% 
      rename_with(.cols = all_of(c(receiverCol)), .fn = ~paste0(., ".receiver"))
    
    receiverColName <- paste0(receiverCol, ".receiver")
    
  }
  
  # if(exportSampleDF){return(list(df.sender, df.receiver))}
  # if(senderProp & receiverProp){message("Expected many-to-many warning because of frequency vs. frequency comparison:")}
  df.joined <- inner_join(df.sender, df.receiver %>% select(-all_of(xVar), -all_of(fillVar)), by = sample.col)
  
  df.long <- df.joined %>%
    pivot_longer(cols = any_of(senderColName), names_to = "senderName", values_to = "senderExp") %>% 
    pivot_longer(cols = any_of(receiverColName), names_to = "receiverName", values_to = "receiverExp") 
  
  
  xPlotVar <- "senderName"
  xPlotVal <- "senderExp"
  yPlotVar <- "receiverName"
  yPlotVal <- "receiverExp"
  
  if(senderProp){
    
    if(is.null(senderRelative)){senderRelative <- sample.col}
    
    if(is.null(senderDenom)){
      if(!is.null(senderRelative)){
        senderDenom <- senderRelative
      } else {
        senderDenom <- sample.col
      }
    }
    
    df.long <- left_join(df.long,
                         PlotBoxplotsProp(
                           df.orig,
                           ann.col = senderCol,
                           filterCol = senderCol,
                           filterVal = senderVal,
                           relative.to = senderRelative,
                           sample.col = sample.col,
                           denom.col = senderDenom,
                           minDenomFreq = senderMinDenomFreq,
                           xVar = xVar,
                           fillVar = fillVar,
                           exportDF = T
                         ) %>% 
                           select(-any_of(c(xVar, fillVar, paste(sample.col, senderDenom, sep = "__")))) %>% 
                           rename_with(.cols = any_of(c(senderCol, "annotation.freq", "denominator.freq", "prop")),
                                       .fn = ~paste0(., ".sender"))
                         ,
                         by = c(sample.col, senderExp = paste0(senderCol, ".sender"))
    )
    
    xPlotVar <- "senderExp"
    xPlotVal <- "prop.sender"
  }
  
  if(receiverProp){
    
    if(is.null(receiverRelative)){receiverRelative <- sample.col}
    
    if(is.null(receiverDenom)){
      if(!is.null(receiverRelative)){
        receiverDenom <- receiverRelative
      } else {
        receiverDenom <- sample.col
      }
    }
    
    df.long <- left_join(df.long,
                         PlotBoxplotsProp(
                           df.orig,
                           ann.col = receiverCol,
                           filterCol = receiverCol,
                           filterVal = receiverVal,
                           relative.to = receiverRelative,
                           sample.col = sample.col,
                           denom.col = receiverDenom,
                           minDenomFreq = receiverMinDenomFreq,
                           xVar = xVar,
                           fillVar = fillVar,
                           exportDF = T
                         ) %>% 
                           select(-any_of(c(xVar, fillVar, paste(sample.col, receiverDenom, sep = "__")))) %>% 
                           rename_with(.cols = any_of(c(receiverCol, "annotation.freq", "denominator.freq", "prop")),
                                       .fn = ~paste0(., ".receiver"))
                         ,
                         by = c(sample.col, receiverExp = paste0(receiverCol, ".receiver"))
    )
    
    yPlotVar <- "receiverExp"
    yPlotVal <- "prop.receiver"
  }
  
  {if(removeZeros & !senderProp) df.long <- df.long %>% filter(senderExp > 0)}
  {if(removeZeros & !receiverProp) df.long <- df.long %>% filter(receiverExp > 0)}
  
  # Reset factor levels to input for features
  # https://stackoverflow.com/questions/16143700/pasting-two-vectors-with-combinations-of-all-vectors-elements
  
  if(!is.null(features.sender)){
    if(!swapSenderOrder){
      # this goes through each denominator in a feature, then switches features (better for comparing features among denominators)
      df.long <- df.long %>% mutate(senderName = fct_relevel(senderName, as.vector(outer(senderVal, features.sender, paste, sep="_"))))
    } else {
      # this goes through each feature in a denominator, then switches denominators (better for comparing features within a single denominator)
      df.long <- df.long %>% mutate(senderName = fct_relevel(senderName, paste(rep(senderVal, each = length(features.sender)), features.sender, sep = "_")))
    }
  }
  
  if(!is.null(features.receiver)){
    if(!swapReceiverOrder){
      # this goes through each denominator in a feature, then switches features (better for comparing features among denominators)
      df.long <- df.long %>% mutate(receiverName = fct_relevel(receiverName, as.vector(outer(receiverVal, features.receiver, paste, sep="_"))))
    } else {
      # this goes through each feature in a denominator, then switches denominators (better for comparing features within a single denominator)
      df.long <- df.long %>% mutate(receiverName = fct_relevel(receiverName, paste(rep(receiverVal, each = length(features.receiver)), features.receiver, sep = "_")))
    }
  }
  
  if(exportSampleDF){return(df.long)}
  
  if(is.null(xLab)){
    if(senderProp){
      xLab <- sprintf("Cell Type Frequency\n(relative to %s)", stringr::str_to_title(senderRelative))
    } else {
      xLab <- sprintf("%s Expression", stringr::str_to_title(stat.use))
    }
  }
  
  if(is.null(yLab)){
    if(receiverProp){
      yLab <- sprintf("Cell Type Frequency\n(relative to %s)", stringr::str_to_title(receiverRelative))
    } else {
      yLab <- sprintf("%s Expression", stringr::str_to_title(stat.use))
    }
  }
  
  if(!groupReg){
    p <- ggplot(df.long, aes(x = get(xPlotVal), y = get(yPlotVal)))
  } else {
    p <- ggplot(df.long, aes(x = get(xPlotVal), y = get(yPlotVal), group = get(xVar)))
  }
  p +
    geom_point(size = dotSize, shape = 21, aes(fill = get(fillVar))) +
    {if(showReg & !matchRegColour & !fillReg & doTests) stat_poly_line(method = "lm", se = regSE, color = regColour)} +
    {if(showReg & matchRegColour & !fillReg & doTests) stat_poly_line(method = "lm", se = regSE, aes(color = get(fillVar)))} +
    {if(showReg & !matchRegColour & fillReg & doTests) stat_poly_line(method = "lm", se = regSE, color = regColour, aes(fill = get(fillVar)))} +
    {if(showReg & matchRegColour & fillReg & doTests) stat_poly_line(method = "lm", se = regSE, aes(color = get(fillVar), fill = get(fillVar)))} +
    {if(showCoef & doTests) stat_poly_eq(method = "lm", use_label(c("R2", "P")), size = regSize, rsquared.conf.level = r2level)} +
    {if(facetType == "wrap") facet_wrap(get(xPlotVar)~get(yPlotVar), scales = facet_scales, ncol = ncols, dir = facet_dir)} +
    {if(facetType == "grid" & !swap_facets) facet_grid(get(yPlotVar)~get(xPlotVar), scales = facet_scales, space = facet_space)} +
    {if(facetType == "grid" & swap_facets) facet_grid(get(xPlotVar)~get(yPlotVar), scales = facet_scales, space = facet_space)} +
    scale_fill_manual(values = fillVector) + 
    scale_color_manual(values = fillVector) + 
    scale_size_continuous(trans = sizeTrans, range = c(dotMin, dotMax)) +
    {if(senderProp) scale_x_continuous(labels = scales::percent)} +
    {if(receiverProp) scale_y_continuous(labels = scales::percent)} +
    {if(stat.use == "percent") scale_x_continuous(labels = percent)} +
    {if(stat.use == "percent") scale_y_continuous(labels = percent)} +
    {if(!fillLegend) guides(fill = "none")} +
    {if(!sizeLegend) guides(fill = "none")} +
    {if(!colorLegend) guides(fill = "none")} +
    {if(labelPoints) geom_text_repel(aes(label = get(sample.col)), size = labelSize, fontface = labelFace)} +
    labs(x = xLab, y = yLab, fill = fillTitle, size = sizeTitle, color = colorTitle, title = plotTitle, subtitle = subTitle) +
    theme_light() +
    theme(
      plot.title = element_text(size = titleSize, face = "bold"),
      plot.subtitle = element_text(size = subtitleSize, face = "bold"),
      axis.text.x = element_text(size = axisXsize, face = "bold"),
      axis.text.y = element_text(size = axisYsize, face = "bold"),
      axis.title.x = element_text(size = titleXsize, face = "bold"),
      axis.title.y = element_text(size = titleYsize, face = "bold"),
      legend.position = legendPosition,
      legend.text = element_text(size = legendSize, face = "bold"),
      legend.title = element_text(size = legendSize, face = "bold"),
      strip.text = element_text(size = stripSize, face = "bold", color = "black"),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    )            
  
}