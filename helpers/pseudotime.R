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
