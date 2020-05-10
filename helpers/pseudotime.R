euclidean.mat <- function(A,B){
  x <- do.call(cbind,rep(list(colSums(A^2)),ncol(B)))
  y <- do.call(rbind,rep(list(colSums(B^2)),ncol(A)))
  suppressWarnings(res <- sqrt(x + y - 2*crossprod(A,B)))
  res[is.na(res) | is.nan(res)] <- 0
  return(res)
}


project.point.to.emb <- function(points,X,emb){
  R <- euclidean.mat(t(points), t(X))
  colnames(R)=rownames(X)
  R <- exp(-R/.1)
  R[is.na(R) | is.nan(R)] <- 0;R <- t(R/rowSums(R))
  return(t(t(emb[rownames(R),]) %*% R)/colSums(R))
}


plot.graph.fates.to.umap <- function(X,NodePositions,Edges,umap,spectrum=NULL,labels=F){
  umap1=umap[rownames(X),]
  
  Nu=project.point.to.emb(NodePositions,X,umap)
  
  segments=cbind(Nu[Edges[,1],],Nu[Edges[,2],])
  
  pl=ggplot()+geom_point(data=data.frame(umap),aes(x=UMAP1,y=UMAP2),color="grey",size=1)+
    geom_point(data=data.frame(umap1),aes(x=UMAP1,y=UMAP2),color="black",size=1.5)+theme(aspect.ratio = 1)
  
  
  pl = pl+ geom_point(data=data.frame(umap1),aes(x=UMAP1,y=UMAP2), color="white",size=1)
  
  if (!is.null(spectrum)){
    umap2=umap[names(spectrum),]
    pl = pl+ geom_point(data=data.frame(umap2),aes(x=UMAP1,y=UMAP2,col=spectrum),size=1)+
      theme(legend.position =c(.9,.2))+
      scale_color_distiller(palette="Spectral", name="fates",trans = "reverse",breaks = c(-1,1),labels=c("endpoint","root cell"))
  } 
  
  
  pl=pl+geom_point(aes(x=Nu[,1],y=Nu[,2]),color="black",size=2)+theme_void()+
    geom_segment(aes(x=segments[,1],y=segments[,2],xend=segments[,3],yend=segments[,4]))
  
  if (labels){
    pl=pl+geom_point(aes(x=Nu[,1],y=Nu[,2]),color="grey",size=6)+
      geom_text(aes(x=Nu[,1],y=Nu[,2]),label=1:nrow(Nu),size=5)
  }
  pl
  
}

remove.nodes <- function(NodePositions,Edges,Nodes){
  keepnodes=(1:nrow(NodePositions))[-Nodes]
  
  keepnodes2=logical(dim(NodePositions)[1])
  keepnodes2[keepnodes]=T
  
  keepedges=!apply(Edges,1,function(check) any(check%in%Nodes))
  
  NodePositions=NodePositions[keepnodes,]
  Edges=Edges[keepedges,]
  
  for (i in 1:nrow(Edges)){
    for (j in 1:2){
      Edges[i,j]=match(Edges[i,j],keepnodes)
    }
  }
  
  return(list(NodePositions,Edges))
}


graph.to.ppt <- function(NodePositions,Edges=NULL,X,emb=NULL,plot=T){
  
  if (is.null(emb)){emb=X[,1:2]}
  
  Xe=t(X)
  F.mat <- t(NodePositions)
  rownames(F.mat) <- NULL
  colnames(F.mat) <- NULL
  
  
  #R <- euclidean.mat(t(X),F.mat)
  R <- euclidean.mat(t(X),F.mat)
  
  R <- exp(-R/.1)
  d <- euclidean.mat(F.mat, F.mat)
  R[is.na(R) | is.nan(R)] <- 0;R <- R/rowSums(R)
  
  if (!is.null(Edges)){
    Net <- igraph::graph.empty(n = max(Edges), 
                               directed = FALSE)
    igraph::V(Net)$name <- paste(1:max(Edges))
    Net <- igraph::add.edges(graph = Net, as.character(t(Edges)))
    B <- as.matrix(get.adjacency(Net))
  }else{
    bt <- minimum.spanning.tree(graph.adjacency(as.matrix(d), 
                                                weighted = T, mode = "undirected"))
    B <- as.matrix(get.adjacency(bt))
  }
  
  D <- diag(nrow(B)) * rowSums(B)
  L <- D - B
  M <- L + diag(ncol(R)) * colSums(R)
  
  g = graph.adjacency(B, mode = "undirected")
  tips = V(g)[degree(g) == 1]
  forks = V(g)[degree(g) > 2]
  
  score = 42
  
  
  
  z=list(score = score, F = F.mat, B = B, R = R, L = L, 
         DT = d, lambda = 1, sigma = .1, n.steps = 1, 
         metrics = "euclidean", M = 101, tips = tips, forks = forks)
  
  if (plot){plotppt(z,emb,tips=F,forks=F,cex.tree = 0.2,lwd.tree = 2)}
  return(z)
}

ppt.to.dyno <- function(z,emb,X=NULL,simplify=T) {
  require(dyno)
  if (is.null(X)){X=emb}
  g = graph.adjacency(z$B, mode = "undirected")
  ml = igraph::as_data_frame(g,"edges")
  ml$from=as.numeric(ml$from)
  ml$to=as.numeric(ml$to)
  
  PartStruct <- PartitionData(
    X,
    NodePositions =t(z$F)
  )
  
  ProjStruct <- project_point_onto_graph(
    X,
    NodePositions = t(z$F),
    Edges = as.matrix(ml),
    Partition = PartStruct$Partition
  )
  
  milestone_network <- ProjStruct$Edges %>%
    dplyr::as_data_frame()%>%
    mutate(
      from = paste0("M", from),
      to = paste0("M", to),
      length = ProjStruct$EdgeLen,
      directed = FALSE
    )
  
  progressions <- tibble(cell_id = rownames(X), edge_id = ProjStruct$EdgeID) %>%
    left_join(milestone_network %>% dplyr::select(from, to) %>% mutate(edge_id = row_number()), "edge_id") %>%
    dplyr::select(-edge_id) %>%
    mutate(percentage = pmin(1, pmax(0, ProjStruct$ProjectionValues)))
  
  traj=wrap_data(
    cell_id = rownames(X),
  ) %>% 
    add_trajectory(
      milestone_network = milestone_network,
      progressions = progressions # either milestone_percentages or progressions have to be provided
    ) %>% add_dimred(emb)
  
  if (!is.null(z$root)){traj=add_root(traj,root_milestone_id = paste0("M",z$root))}
  
  if (simplify){return(simplify_trajectory(traj))}else{return(traj)}
}

scde.process.dataset <- function(dat,name,env=go.env,batch=NULL,k=min(20,ncol(cd)/n.groups),max.model.plots=50,cd=NULL,varinfo=NULL,knn=NULL,max.adj.var=5,skip.pca=FALSE,min.size.entries=1e3,control.for.depth.variation=TRUE,max.quantile=1,seed=0) {
  cat("processing ",name,": ");
  if(is.null(cd)) {
    cd <- dat;
    
    CairoPNG(file=paste(name,"reads.per.cell.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(colSums(dat)/1e6,col="wheat",xlab="reads per cell (M)",main="Read counts across cells")
    abline(v=1e5/1e6,lty=2,col=2)
    dev.off()
    table(colSums(dat)>=min.cell.reads)
    
    CairoPNG(file=paste(name,"reads.per.gene.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(log10(rowSums(dat)+1),col="wheat",xlab="reads per gene (log10)",main="Read counts across genes")
    abline(v=log10(min.gene.reads+1),lty=2,col=2)
    dev.off()
    
    CairoPNG(file=paste(name,"genes.per.cell.png",sep="."),width=350,height=350)
    par(mfrow=c(1,1), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 0.9);
    hist(log10(colSums(dat>0)+1),col="wheat",xlab="genes per cell (log10)",main="Gene counts across cells")
    abline(v=log10(min.cell.genes+1),lty=2,col=2)
    dev.off()
    
    # filter out low-gene cells
    vi <- colSums(cd)>min.cell.reads; table(vi)
    cd <- cd[,vi];
    
    # remove genes that don't have many reads
    vi <- rowSums(cd)>min.gene.reads; table(vi)
    cd <- cd[vi,];
    
    # remove genes that are not seen in a sufficient number of cells
    vi <- rowSums(cd>0)>min.gene.cells; table(vi)
    cd <- cd[vi,];
  }  
  cat("proceeding with ",nrow(cd)," genes across ",ncol(cd)," cells ");
  
  set.seed(seed);
  
  if(is.null(knn) || is.null(varinfo)) {
    knn <- knn.error.models(cd,groups=as.factor(rep(name,ncol(cd))),k=k,n.cores=n.cores,min.count.threshold=1,min.nonfailed=min.nonfailed,verbose=0,max.model.plots=max.model.plots,min.size.entries=min.size.entries)
    cat("models ")
    prior <- scde.expression.prior(models=knn,counts=cd,length.out=400,show.plot=F,max.quantile=max.quantile)
    pdf(file=paste(name,"varnorm.pdf",sep="."),height=4,width=8)
    varinfo <- pagoda.varnorm(knn,counts=cd,trim=trim/ncol(cd),plot=T,verbose=1,prior=prior,max.adj.var=max.adj.var,weight.df.power=1,batch=batch)
    dev.off();
    cat("varinfo ")
    if(control.for.depth.variation) {
      varinfo <- pagoda.subtract.aspect(varinfo,colSums(cd[,rownames(knn)]>0))
    }
  }
  
  if(!skip.pca) {
    pwpca <- pagoda.pathway.wPCA(varinfo,env,n.components=1,n.cores=n.cores,n.internal.shuffles=0,verbose=1,n.randomizations=5)
    cat("pathways ")
    pdf(file=paste(name,"clvar.pdf",sep="."),height=4,width=8)
    clpca <- pagoda.gene.clusters(varinfo,trim=(trim+5)/ncol(varinfo$mat),n.clusters=150,n.cores=n.cores,verbose=1,plot=T)
    dev.off();
    cat("clusters\n")
  } else {
    pwpca <- clpca <- NULL;
  }
  return(list(cd=cd,knn=knn,prior=prior,varinfo=varinfo,pwpca=pwpca,clpca=clpca,batch=batch))
}


plot.tree <- function(ppt,dyn,emb,subtree=NULL,proj_sd=.8){
  emb2=emb[!rownames(emb)%in%rownames(ppt$cell.summary),]
  pl=plot_dimred(dyn,trajectory_projection_sd = proj_sd,size_trajectory = 2)+geom_point(data=data.frame(emb2),aes(x=emb2[,1],y=emb2[,2]),color="grey")+theme(aspect.ratio = 1)
  pl=gginnards::move_layers(pl,idx=7,position = "bottom")
  pl=gginnards::move_layers(pl,idx=8,position = "bottom")
  pl$layers[[5]]$aes_params$size=2
  pl$layers[[6]]$geom_params$lineend="round"
  pl$data$color=ppt$cell.summary[pl$data$cell_id,]$color
  if (!is.null(subtree)){
    pl$data[!pl$data$cell_id%in%rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% subtree$segs],]$color=
      "grey"
  }
  
  return(pl)
}
