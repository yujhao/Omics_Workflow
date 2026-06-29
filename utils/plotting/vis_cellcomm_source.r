ParseCpdb <- function(
  cpdb,
  pvalue = 0.05
){
  deconvoluted <- read.table(glue::glue("{cpdb}/deconvoluted.txt"),
                            header=T, stringsAsFactors = F,
                            sep="\t", comment.char = '', check.names=F)
  raw_pval <- read.table(glue::glue("{cpdb}/pvalues.txt"), header=T,
                        stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  raw_pval$interacting_pair = apply(raw_pval,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2], ignore.case = T)
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })
  raw_means <- read.table(glue::glue("{cpdb}/means.txt"), header=T,
                         stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
  raw_means$interacting_pair = apply(raw_means,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2])
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })
  raw_sig_means <- read.table(glue::glue("{cpdb}/significant_means.txt"),
                         header=T, stringsAsFactors = F, sep="\t",
                         comment.char = '', check.names=F)
  raw_sig_means$interacting_pair = apply(raw_sig_means,1, function(x){
                      if(grepl("^complex:", x[3], perl = T)[1]){
                        pt_a = gsub("^complex:", "", x[3])
                        x[2]= gsub(glue::glue("{pt_a}_"), glue::glue("{pt_a}|"), x[2])
                      }else if ( grepl("^complex:", x[4], perl = T)[1] ){
                        pt_b = gsub("^complex:", "", x[4])
                        x[2] = gsub(glue::glue("_{pt_b}"), glue::glue("|{pt_b}"), x[2])
                      } else if ( grepl("^simple:",x[3], perl = T) && grepl("^simple:", x[4],perl =T)){
                        x[2] = gsub("^([a-zA-Z0-9-]+)_([a-zA-Z0-9-]+)$", "\\1|\\2", x[2])
                      }
    })

  desired_cols  <- colnames(raw_pval)
  sig_gene_df <- raw_pval %>% dplyr::select(desired_cols) %>%
                  dplyr::select(1:2,12:dim(.)[2]) %>%
                  tidyr::gather( "cell_pair","pval", 3:dim(.)[2]) %>%
                  distinct() %>% filter( pval < pvalue ) 
  sig_gene_pair = sig_gene_df %>% pull(interacting_pair)

  desired_means <- raw_means %>%
                dplyr::filter( interacting_pair %in% sig_gene_pair ) %>%
                dplyr::select( desired_cols)
  desired_pval <- raw_pval %>%
                dplyr::filter( interacting_pair %in% sig_gene_pair ) %>%
                dplyr::select( desired_cols )

  desired_sig_pvalx <- desired_pval %>%
                        dplyr::select( id_cp_interaction,interacting_pair, receptor_a, receptor_b, secreted, is_integrin, 12:dim(.)[2] ) %>%
                        tidyr::gather( "cell_pair","pval", 7:dim(.)[2]) %>% distinct()
  desired_sig_meansx <- desired_means %>%
                        dplyr::select( id_cp_interaction,interacting_pair, 12:dim(.)[2] ) %>%
                        tidyr::gather( "cell_pair","expr", 3:dim(.)[2]) %>% distinct()
  merged_df <- dplyr::full_join(desired_sig_meansx, desired_sig_pvalx,
            by =c("id_cp_interaction" = "id_cp_interaction","interacting_pair" = "interacting_pair", "cell_pair" = "cell_pair") ) %>%
            tidyr::separate( cell_pair, into = c("part_a_cell", "part_b_cell"), sep = "\\|") %>%
            tidyr::separate(interacting_pair, into = c("part_a_gene", "part_b_gene"), sep ="\\|")
   lr_df <- future.apply::future_lapply(1:nrow(merged_df), function(i){
        if ( merged_df[i,"receptor_a"]== "False" & merged_df[i,"receptor_b"] == "True" ){
             a <- merged_df[i,c(1,3,2,5,4,6,9,10,11)]
            names(a) <- c("id_cp_interaction","receptor", "ligand", "receptor_cell", "ligand_cell", "expr","secreted","is_integrin","pval")
        }else{
            a <- merged_df[i,c(1,2,3,4,5,6,9,10,11)]
            names(a) <- c("id_cp_interaction","receptor", "ligand", "receptor_cell", "ligand_cell", "expr","secreted","is_integrin","pval")
        }
        a
    })
    lr_df <- do.call(rbind, lr_df)

# 将新列添加到数据框
lr_df$receptor_expr=""
lr_df$ligand_expr=""

for(i in c(1:nrow(lr_df))){
    tmp_df = deconvoluted[which(deconvoluted$id_cp_interaction == lr_df[i, "id_cp_interaction"]),]
    lr_df[i, "receptor_expr"] = min(tmp_df[apply(tmp_df[, c("gene_name", "complex_name")], MARGIN = 1, function(y) any(grepl(paste0("\\b", lr_df[i, "receptor"], "\\b"), y))), lr_df[i, "receptor_cell"]])
    lr_df[i, "ligand_expr"] = min(tmp_df[apply(tmp_df[, c("gene_name", "complex_name")], MARGIN = 1, function(y) any(grepl(paste0("\\b", lr_df[i, "ligand"], "\\b"), y))), lr_df[i, "ligand_cell"]])
}
lr_df$receptor_expr = log2(as.numeric(lr_df$receptor_expr) + 0.0001)
lr_df$ligand_expr = log2(as.numeric(lr_df$ligand_expr) + 0.0001)

  cellphonedb <- list(raw.pvalues=raw_pval, raw.means = raw_means,
                      ligrec = lr_df,desired_means = desired_means, desired_pval = desired_pval,
                     raw.sigmean = raw_sig_means, deconvoluted = deconvoluted)
  return(cellphonedb)
}
LRDotplot <- function(
  data,
  is_onlySig = T,
  topn = 5,
  xangle = 90,
  xsize = 10,
  remove.isolate = TRUE,
  palette = c("black", "blue", "yellow", "red")
){
  if( class(data) == "data.frame" ){
      all_data = data[,c("receptor","ligand","receptor_cell","ligand_cell","pval","expr","receptor_expr","ligand_expr")]
  }else{
      stop("NO cell communication matrix is procided!")
  }

  all_data = all_data %>% unite( "pair", receptor, ligand, sep = "|") %>% unite( "clusters", receptor_cell, ligand_cell, sep = "|")

  if ( is_onlySig ){
    desired_pair = ""
    for (i in unique(all_data$clusters)) {
      desired_pair2 <- all_data %>% dplyr::filter( pval < 0.05 & receptor_expr != 0 & ligand_expr != 0 ) %>% arrange( desc(expr)) %>% filter( clusters == i ) %>% pull( pair )
      desired_pair <- c(desired_pair,unique(desired_pair2)[1:topn])
    }
    filter_data = all_data %>% dplyr::filter( pair %in% unique(desired_pair))
  }else{
    filter_data = all_data
  }

  filter_data$pval[filter_data$pval==0] = 0.0009
  while ( length(filter_data$expr[filter_data$expr==0]) >0 ) {
          filter_data$expr[filter_data$expr==0][1] = (2^filter_data[filter_data$expr==0,][1,]$receptor_expr+2^filter_data[filter_data$expr==0,][1,]$ligand_expr)/2 
  }
  filter_data$mean = as.numeric(log2(filter_data$expr))

  plot.data = filter_data[,c("pair","clusters","pval","mean")]
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  # 为横坐标排序
  clusters_pairs = as.character(unique(plot.data$clusters))
  # 提取横坐标括号外和括号内的字符串
  clusters_df <- data.frame(name = clusters_pairs, pair = stringr::str_extract(clusters_pairs, "^[^\\(]+"), 
                            group = stringr::str_extract(clusters_pairs, "(?<=\\()[^\\)]+")) 
  # 根据细胞类型首字符和data文件中的分组顺序排序
  clusters_df <- clusters_df[order(clusters_df$pair,match(clusters_df$group, unique(data$group))), ]  
  if (! is.na(unique(clusters_df$group)) & ! remove.isolate ){
    #如果保留所有受配体和分组列
    clusters_levels = NULL
    for (i in unique(clusters_df$pair)){ clusters_levels = c(clusters_levels,paste0(i , "(",unique(clusters_df$group), ")" ))}
    diff_clusters = setdiff(clusters_levels, unique(plot.data$clusters)) 
    if(length(diff_clusters) != 0 ){ #所有分组受配体都一致就不用添加na_df
      na_df = data.frame(pair = plot.data$pair[1:length(diff_clusters)], clusters = diff_clusters, pvalue = NA, mean = NA)
      plot.data = rbind(plot.data, na_df)
    }
  } else { 
    clusters_levels = as.character(clusters_df$name) 
  }
  # 给横坐标加上levels
  plot.data$clusters = factor(plot.data$clusters, levels = clusters_levels)
  plot.data$pvalue=as.numeric(plot.data$pvalue) 
  plot.data$mean=as.numeric(plot.data$mean)
  my_palette <- colorRampPalette(palette, alpha=TRUE)(n=399)
  dotplot <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() + #coord_fixed() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=xsize, face = "bold", colour = "black"),
          axis.text.x = element_text(size = xsize, face = "bold", angle = xangle, hjust = 1),
          axis.text.y = element_text(size = xsize, face = "bold", colour = "black"),
          axis.title=element_blank(),
          text = element_text(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
# 如果有分组，添加文字颜色和竖线
  if( ! is.na(unique(clusters_df$group))){
    clusters_df <- data.frame(name = clusters_levels, pair = stringr::str_extract(clusters_levels, "^[^\\(]+"), 
                            group = stringr::str_extract(clusters_levels, "(?<=\\()[^\\)]+")) 
    dataset.name.order <- clusters_df$group
    # dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
    # dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order) - 1)
    color <- SelectColors(palette = "tableau10medium",length(unique(dataset.name.order)))
    names(color) <- unique(dataset.name.order)
    xtick.color <- color[dataset.name.order]
    # 计算各受配体细胞类型关系对的长度（lengths表示各元素差异序列）
    diffs <- rle(as.vector(clusters_df$pair))$lengths
    # cumsum计算每个子序列的结束位置
    xintercept = cumsum(diffs) + 0.5
    dotplot <- dotplot + theme(axis.text.x = element_text(colour = xtick.color)) + 
                        geom_vline(xintercept = xintercept[-length(xintercept)], linetype = "dashed", color = "grey60", size = 0.4)
  }
  out=list(dotplot,plot.data)
  return(out)
}

LRNetwork<-function(
  data,
  col,
  label=TRUE,
  edge.curved=0.2,
  shape='circle',
  layout=nicely(),
  vertex.size=20,
  margin=0.4,
  vertex.label.cex=1.5,
  vertex.label.color='black',
  arrow.width=1.5,
  edge.label.color='black',
  edge.label.cex=1,
  edge.max.width=10,
  edge.arrow.size=1
){
  set.seed(1234)
  net<-data %>% group_by(ligand_cell, receptor_cell) %>% dplyr::summarize(n=n())
  net<-as.data.frame(net,stringsAsFactors=FALSE)
  g<-graph.data.frame(net,directed=TRUE)
  edge.start <- ends(g, es=E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  loop.angle<-ifelse(coords_scale[V(g),1]>0,-
                       atan(coords_scale[V(g),2]/coords_scale[V(g),1]),
                     pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  vertex.size <- data %>% group_by(ligand_cell) %>% dplyr::summarize(n=n())
  vertex.size <- vertex.size$n*0.07
  V(g)$size<-vertex.size
  V(g)$color<-col[V(g)]
  #V(g)$label.color<-vertex.label.color
  #V(g)$label.cex<-vertex.label.cex
  V(g)$frame.color <- NA
  #V(g)$label.font <- 2
  V(g)$label <- NA
  if(label){
    E(g)$label<-E(g)$n
  }
  if(max(E(g)$n)==min(E(g)$n)){
    E(g)$width<-2
  }else{
    E(g)$width<-1+edge.max.width/(max(E(g)$n)-min(E(g)$n))*(E(g)$n-min(E(g)$n))
  }
  E(g)$arrow.width<-arrow.width
  E(g)$label.color<-edge.label.color
  E(g)$label.cex<-edge.label.cex
#  E(g)$color<-V(g)$color[edge.start[,1]]
  E(g)$color<-"#BDBDBD"
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=layout.circle,margin=margin,edge.arrow.size=edge.arrow.size)

  # Create a combined legend (ligand + receptor cells) with corresponding colors
  if(is.factor(data$ligand_cell)){
    legend_labels <- unique(c(levels(data$ligand_cell), levels(data$receptor_cell)))
  }else{
    ligand_legend <- unique(data$ligand_cell)
    receptor_legend <- unique(data$receptor_cell)
    legend_labels <- unique(c(ligand_legend, receptor_legend))
  }
  legend_colors <- c(col[legend_labels])
  # Add a legend with colors matching the nodes
  legend("right", 
         legend = legend_labels, 
         fill = legend_colors, # Use the colors of ligands and receptors from the `col` vector
         title = expression(bold("Ligand and Receptor Cells")),
         cex = 1.2, 
         bty = "n")
  return(g)
}
LRCircos <- function(
  data,
  screenvar="expr",
  topn=5,
  pval=0.05,
  gap.degree=3,
  color_expr=c("green","yellow","red"),
#  expr_range=c(-15,-6,3),
  cell_col=NULL,
  link.lty=NULL,
  labels.cex=0.58,
  arr.length=0.2,
  link.col=NULL,
  link.lwd=2,
  link.arr.width=NULL,
  link.arr.type=NULL,
  facing='clockwise',
  track.height_1=uh(8,'mm'),
  track.height_2=uh(12,'mm'),
  text.vjust=0.5,
  ...
  ){
  if(class(data) == "data.frame"){
    data = data
  }else if(file.exists(data)){
    data <- read.csv(data, header = T, sep = '\t')
  }else{
    stop("No cell communication result is provided!")
  }

  if(is.null(screenvar)){
    screenvar = "expr"
  }

  if(is.null(topn)){
    topn = 5
  }

  data <- data %>% group_by(receptor_cell,ligand_cell) %>% filter( pval < 0.05 ) %>% top_n(topn, !!ensym(screenvar))
  data <- as.data.frame(data)
  data <- data %>% dplyr::mutate(lr = 'ligand', lr2 = 'receptor')
  ldf <- data %>% dplyr::select(ligand_cell, lr, ligand, ligand_expr,
                receptor_cell, lr2, receptor, receptor_expr, everything()) %>%
          dplyr::rename("cell" = "ligand_cell", "gene" = "ligand", "expr1" = "ligand_expr",
                "cell2" = "receptor_cell", "gene2" = "receptor", "expr2" = "receptor_expr")
  rdf <- data %>% dplyr::select(receptor_cell, lr2, receptor, receptor_expr,
                ligand_cell, lr, ligand, ligand_expr, everything()) %>%
          dplyr::rename("cell" = "receptor_cell", "lr2" = "lr", "gene" = "receptor", "expr1" = "receptor_expr",
                "cell2" = "ligand_cell", "lr" = "lr2", "gene2" = "ligand", "expr2" = "ligand_expr")
  df <- rbind(ldf,rdf)
  df <- df[order(df$cell),]
  df <- df %>% dplyr::mutate(gene_id = NA, gene_id2 = NA)
  df$gene_id <- paste("geneid", 1:nrow(df))
  d1 <- data.frame(index = paste(df$cell,df$lr,df$gene,df$cell2,df$lr2,df$gene2,sep = '_'), gene_id = df$gene_id)
  d2 <- data.frame(index = paste(df$cell2,df$lr2,df$gene2,df$cell,df$lr,df$gene,sep = '_'), gene_id = df$gene_id2)
  for(i in 1:nrow(d1)){
    index = which(d2$index %in% d1[i,1])
    for(ind in index){
      d2[ind,2] = as.character(d1[i,2])
    }
  }
  df$gene_id2 <- d2$gene_id
  df <- df %>% dplyr::rename("ligand_cell" = "cell", "ligand" = "gene", "ligand_expr" = "expr1",
               "receptor_cell" = "cell2", "receptor" = "gene2", "receptor_expr" = "expr2") %>%
        select(ligand_cell, lr, ligand, gene_id, ligand_expr, receptor_cell, lr2,
               receptor, gene_id2, receptor_expr, everything())
  lrid <- df[which(df$lr == 'ligand'),]
  
  max_range = max(ceiling(max(lrid$ligand_expr)),ceiling(max(lrid$receptor_expr)))
  min_range = min(floor(min(lrid$ligand_expr)), floor(min(lrid$receptor_expr)))

  expr_range = c(min_range, (max_range+min_range)/2, max_range)

  if(is.null(color_expr)){
    stop("No color set of gene expression in the plot is provide!")
  }else{
    color_expr = color_expr
  }

  if(is.null(pval)){
    pval = 0.05
  }else{
    pval = pval
  }
  if(is.null(link.lty)){
    link.lty = structure(ifelse(lrid$pval < pval, 'solid', 'dashed'),
                             names = paste(lrid$ligand_cell,lrid$receptor))
  }
  if(is.null(lrid$comm_type)){
    lrid <-lrid %>% dplyr::mutate(link_col = 'black')
  }else{
    if(is.null(link.col)){
      comm_col <- structure(brewer.pal(length(unique(lrid$comm_type)),'Paired'),
                            names = as.character(unique(lrid$comm_type)))
      lrid <-lrid %>% dplyr::mutate(link_col = as.character(comm_col[as.character(lrid$comm_type)]))
    }else{
      lrid$link_col = link.col
    }
  }
  if(is.null(cell_col)){
    cell_col<-structure(c("#7FC97F","#BEAED4","#FDC086","#FBB4AE","#1B9E77","#FFFF99","#386CB0","#F0027F",
                          "#666666","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666",
                          "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                          "#CAB2D6","#6A3D9A","#FFFF99","#B15928","#B3CDE3","#CCEBC5","#DECBE4","#BF5B17",
                          "#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8",
                          "#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC","#E41A1C","#377EB8","#4DAF4A",
                          "#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
                          "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3",
                          "#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD",
                          "#CCEBC5", "#FFED6F"),names=as.character(unique(df$receptor_cell)))
#    cell_col <- structure(CustomCol(1:length(unique(df$receptor_cell))),
#                          names = as.character(unique(df$receptor_cell)))
  }
  if(is.null(link.lwd)){
    link.lwd = 2
  }else{
    link.lwd = link.lwd
  }
  if(is.null(labels.cex)){
    labels.cex = 0.58
  }else{
    labels.cex = labels.cex
  }
  if(is.null(arr.length)){
    arr.length = 0.2
  }else{
    arr.length = arr.length
  }

  circle_size = unit(1, "snpc")
  fa = df$gene_id
  fa = factor(fa,levels = fa)
  circos.par(canvas.xlim = c(-1.1,1.1), canvas.ylim = c(-1.1,1.1),
             cell.padding = c(0, 0, 0, 0), gap.degree = gap.degree)
  circos.initialize(factors = fa, xlim = c(0,1))


  col_fun <- colorRamp2(expr_range, color_expr)
  circos.trackPlotRegion(
    ylim = c(0,1), track.height = track.height_1, bg.border = 'black',
    bg.col = col_fun(df$ligand_expr),
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for (i in 1:nrow(df)) {
    circos.axis(sector.index = df[i,4],direction = "outside",
                labels = FALSE,col = 'black', minor.ticks = 0, major.at = seq(0.5, length(df$ligand)/2))
    circos.text(
      x = 0.5,                     # x轴位置
      y = 1.2,                     # 控制距离圆弧的距离（越大越外）
      labels = df[i, 3],           # 显示的文字
      sector.index = df[i, 4],
      facing = "clockwise",        # 顺时针方向
      niceFacing = TRUE,           # 让文字自动调整朝向
      cex = labels.cex * 1.5,      # 字体大小（略微放大）
      font = 2,                    # 加粗（相当于 fontface="bold"）
      adj = c(0, 0.5)             # 对齐方式
    )
  }


  circos.trackPlotRegion(
    ylim = c(0,1), track.height = track.height_2, bg.border = NA,
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for(i in unique(df$ligand_cell)){
    num1 = min(which(df$ligand_cell == i))
    num2 = max(which(df$ligand_cell == i))
    highlight.sector(as.character(df$gene_id[num1:num2]),
                     track.index = 2, text = i,
                     niceFacing = T, font = 2,cex=1,
                     col = cell_col[which(unique(df$ligand_cell) == i)])
  }


  circos.trackPlotRegion(
    ylim = c(0,1),track.height = track.height_1, bg.border = NA,
    panel.fun = function(x, y){
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  for(i in unique(df$ligand_cell)){
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "ligand"))) > 0 ){
      num1 = min(intersect(which(df$ligand_cell == i),which(df$lr == "ligand")))
    }else{
      num1 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "ligand"))) > 0 ){
      num2 = max(intersect(which(df$ligand_cell == i),which(df$lr == "ligand")))
    }else{
      num2 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "receptor"))) > 0 ){
      num3 = min(intersect(which(df$ligand_cell == i),which(df$lr == "receptor")))
    }else{
      num3 = NULL
    }
    if ( length(intersect(which(df$ligand_cell == i),which(df$lr == "receptor"))) > 0 ){
      num4 = max(intersect(which(df$ligand_cell == i),which(df$lr == "receptor")))
    }else{
      num4 = NULL
    }
    if ( !(is.null(num1) | is.null(is.integer(num2))) ) {
      highlight.sector(as.character(df$gene_id[num1:num2]),
                       track.index = 3, text = 'L', cex = 0.85,
                       text.vjust = text.vjust,
                       text.col = 'white', niceFacing = T,
                       col = cell_col[which(unique(df$ligand_cell) == i)])
    }
    if ( !(is.null(num3) | is.null(is.integer(num4))) ) {
    highlight.sector(as.character(df$gene_id[num3:num4]),
                     track.index = 3, text = 'R',
                     text.vjust = text.vjust,
                     text.col = 'white', niceFacing = T, cex = 0.85,
                     col = cell_col[which(unique(df$ligand_cell) == i)])
    }
  }


  for(i in 1:nrow(lrid)){
    circos.link(sector.index1 = lrid[i,4], point1 = 0.5,
                sector.index2 = lrid[i,9], point2 = 0.5,
                directional = 1, h = 0.85, lwd = link.lwd,
                col = lrid[i,'link_col'],
                lty = ifelse(length(link.lty) == 1, link.lty, link.lty[i]),
                arr.length = 0.2 , arr.col = col)
  }


  lgd_expr <- Legend(title = 'gene expression', at = expr_range,
                     col_fun = col_fun, title_position = 'topleft',
                    title_gp = gpar(fontsize = 21, fontface = "bold"),  # 增大图例标题字体
                    title_gap = unit(3, "mm"),
                    labels_gp = gpar(fontsize = 19) )
#  lgd_pval = Legend(title = "ligand_receptor pval", at = 1:2,
#                    labels = c(paste('lr_pval < ', pval), paste('lr_pval >= ', pval)),
#                    type = "lines", legend_gp = gpar(lty = 1:2))
  lgd_pval = Legend(title = "ligand_receptor pval", at = 1:2,
                    labels = c(paste('lr_pval < ', pval)),
                    type = "lines", legend_gp = gpar(lty = 1:1),
                    title_gp = gpar(fontsize = 21, fontface = "bold"),  # 增大图例标题字体
                    title_gap = unit(3, "mm"),
                    labels_gp = gpar(fontsize = 19) )
  lgd_cell = Legend(title = "cell type", at = 1:length(unique(df$ligand_cell)),
                    labels = as.character(unique(df$ligand_cell)),
                    grid_height = unit(ifelse(length(unique(df$ligand_cell))/2<6,6,length(unique(df$ligand_cell))/2), "mm"),
                    legend_gp = gpar(fill = cell_col[1:length(unique(df$ligand_cell))]),
                    title_position = "topleft",
                    title_gp = gpar(fontsize = 21, fontface = "bold"),  # 增大图例标题字体
                    title_gap = unit(2, "mm"),
                    labels_gp = gpar(fontsize = 19) )
  if(is.null(lrid$comm_type)){
    lgd_list_vertical = packLegend(lgd_expr, lgd_cell, lgd_pval,row_gap = unit(4, "mm"))
  }else{
    lgd_comm <- Legend(title = "communication type", at = 1:length(unique(lrid$comm_type)),
                       labels = as.character(unique(lrid$comm_type)), type = 'lines',
                       legend_gp = gpar( col = comm_col[1:length(unique(lrid$comm_type))]),
                       title_position = "topleft",
                       title_gp = gpar(fontsize = 21, fontface = "bold"),  # 增大图例标题字体
                       labels_gp = gpar(fontsize = 19) )
    lgd_list_vertical = packLegend(lgd_expr, lgd_cell, lgd_pval, lgd_comm,row_gap = unit(4, "mm"))
  }
  pushViewport(viewport(x = 0.8, y = 0.5))
  draw(lgd_list_vertical, x = circle_size, just = "left")
  upViewport()

  circos.clear()
}

LRChorddiagram <- function(
  data,
  grid_orbit_h=0.02,
  label_h=0.04,
  diffHeight=2,
  output_dir="./",
  colx = colx,
  extra_legend = TRUE,
  ...
){
  ## significant summary   
  data <- data %>% group_by(ligand_cell,receptor_cell) %>%
      filter( pval < 0.05 ) %>% 
      dplyr::summarize(n=n()) %>% 
      dplyr::rename(significant_pairs =  n )
  #write.table(data,file.path(output, "cell_comm_chorddiagram_summary.xls"),quote=F,sep="\t",row.names=F)

  ## col set
  #如果原始数据有levels保持原来的顺序，如果没有按照字母排序
  if(is.factor(data$receptor_cell)){
    celltypes = unique(c(levels(data$receptor_cell),levels(data$ligand_cell)))
  }else{
    celltypes = sort( unique(c(data$receptor_cell, data$ligand_cell)) )
  }
  #colx = SelectColors(palette = opt$colorschema, length(celltypes))
  #names(colx) <- celltypes
  if (extra_legend != TRUE){
    par(mar = c(20, 20, 20, 20),  # 内边距：上、右、下、左（数值越大留白越多）
        xpd = TRUE)            # 允许文字绘制到边距区域（避免标签被裁剪）
    }
  ## plot 
  chordDiagram(data, 
              grid.col=colx,
              directional = 1,
              diffHeight = -uh(diffHeight, "mm"), 
              direction.type = c("diffHeight", "arrows"),
              link.arr.type = "big.arrow",
              annotationTrack = c("grid"), 
              annotationTrackHeight = c(label_h, grid_orbit_h),
              preAllocateTracks = 1
              )
  if (extra_legend != TRUE){
    # 在空轨迹中放置文本标签
    circos.track(
      track.index = 1, panel.fun = function(x, y) {
        circos.text(
          CELL_META$xcenter, CELL_META$ylim[1], 
          CELL_META$sector.index,  facing = "clockwise", 
          niceFacing = TRUE, adj = c(0,0.5),cex=1.5,rotation = 45
        )
      }, bg.border = NA
    )
    }
# # 单独添加图例
  if (extra_legend == TRUE){
  legend("right", 
        legend = celltypes,  # 图例内容为细胞类型
        fill = colx[celltypes],  # 使用对应的颜色
        title = expression(bold("Ligand & Receptor Cells")), 
        cex = 1.3,               # 图例文字大小
        bty = "n"               # 不绘制图例边框
        )
    }

}

LRHeatmap <- function(
  data,
  ...
){
  ## significant summary
  sum_data = data %>% filter(pval < 0.05) %>% group_by(ligand_cell, receptor_cell) %>%
        dplyr::summarize(n=n()) %>% tidyr::spread(receptor_cell,n)
#  write.table(sum_data, file.path(output_dir, "cell_comm_heatmap_summary.xls"), quote=F,sep="\t", row.names=F)

  plot_data = sum_data %>% tibble::column_to_rownames(var="ligand_cell")

  ## plot
  Heatmap(as.matrix(plot_data),
          col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(50),
          rect_gp = gpar(col = "white"),
          row_title = "Cell types for ligands",
          row_title_side = c("left"),
          row_title_gp = gpar(fontsize = 22, fontface = "bold"),
          column_title = "Cell types for receptors",
          column_title_side = c("top"),
          column_title_gp = gpar(fontsize = 22, fontface = "bold"),
          cluster_rows = F,
          cluster_columns = F,
          row_labels = rownames(plot_data),
          row_names_side = c("right"),
          show_row_names = TRUE,
          column_labels = colnames(plot_data),
          column_names_side = c("bottom"),
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 20, fontface = "bold"),
          column_names_gp = gpar(fontsize = 20, fontface = "bold"),
          name = "Number",
          row_names_max_width = unit(8, "cm"),
          column_names_max_height = unit(8, "cm"),
          heatmap_legend_param = list(
          title_gp = gpar(fontsize = 14, fontface = "bold"),   # 图例标题字体
          labels_gp = gpar(fontsize = 14, fontface = "bold"),  # 图例刻度字体
          legend_height = unit(4, "cm")
          )
  )
#  return(p)
}

LRBarplot <- function(
  data,
  bar_width = 0.6,
  colx = colx,
  ...
) {

  ## significant summary
  plot_data = data %>% filter(pval < 0.05) %>% group_by(ligand_cell, receptor_cell) %>%
        dplyr::summarize(n=n()) %>% group_by(ligand_cell) %>% arrange(ligand_cell) %>% mutate(cumsum = sum(n))

  ## sort the ligand_cell
  #如果之前没有指定顺序需要重新指定，否则保持原来的顺序即可
  if(!is.factor(plot_data$ligand_cell)){
    celltype = unique(plot_data[order(as.matrix(plot_data)[,4],decreasing=T),]$ligand_cell)
    plot_data$ligand_cell = factor(plot_data$ligand_cell,levels=celltype)
  }

  ## plot
  p = ggplot(plot_data,
              aes(fill = receptor_cell, y = n, x = ligand_cell)) +
      geom_bar(position="stack", stat="identity", width = bar_width) +
      labs(x="Cell types of ligands", y = "Number of interactions")+
#      scale_fill_discrete(name="Cell types of receptors") +
      scale_fill_manual("Cell types of receptors", values = colx) +
      theme(panel.background = element_rect(fill = "transparent", colour = NA), 
            panel.grid = element_blank(), 
            axis.title = element_text(size = 18, face = "bold"),
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.line = element_line(size=0.5, colour = "black"),
            legend.title = element_text(size = 15, face = "bold", color = "black"),   # 图例标题字体
            legend.text = element_text(size = 14, color = "black"),
            plot.margin = unit(c(3,3,3,3), "cm"))
  return(p)
}
