pc.select <- function(p2,plt=F,elbow=T){
  
  x <- cbind(1:length(p2$misc$PCA$d), p2$misc$PCA$d)
  line <- x[c(1, nrow(x)),]
  proj <- princurve::project_to_curve(x, line)
  return(which.max(proj$dist_ind))
  
}

doUMAP <- function(PCA,n_neighbors,min_dist,max_dim=2,seed.use=42){
  require(reticulate)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(n_neighbors = as.integer(x = n_neighbors), 
                           n_components = as.integer(x = max_dim), metric = "correlation", 
                           min_dist = min_dist)
  
  umap_output <- umap$fit_transform(as.matrix(x = PCA))
  rownames(umap_output)=rownames(PCA)
  colnames(umap_output)=paste0("UMAP",1:max_dim)
  
  return(umap_output)
}

doPalantir <- function(PCA,n_neighbors,min_dist,n_eig=NULL,seed.use=42){
  library(reticulate)
  
  palantir=import("palantir")
  pd=import("pandas")
  umap=import("umap")
  
  pca_py=pd$DataFrame(r_to_py(PCA))
  cat("Running diffusion maps... ")
  dm_res=palantir$utils$run_diffusion_maps(pca_py)
  cat("done\n")
  if (!is.null(n_eig)){
    ms_data = palantir$utils$determine_multiscale_space(dm_res,n_eigs=as.integer(n_eig))
  } else {
    ms_data = palantir$utils$determine_multiscale_space(dm_res)
  }
  
  ms_data=as.matrix(ms_data);
  rownames(ms_data)=rownames(PCA);colnames(ms_data)=paste0("Dim",1:ncol(ms_data))
  
  
  set.seed(seed = seed.use)
  py_set_seed(seed = seed.use)
  
  cat("Running UMAP... ")
  
  fit=umap$UMAP(n_neighbors=as.integer(n_neighbors),min_dist=min_dist)
  u=fit$fit_transform(ms_data)
  
  cat("done\n")
  
  rownames(u)=rownames(ms_data)
  return(list(ms_data=ms_data,umap=u))
}


p2.wrapper <- function(counts,n_neighbors=30,min_dist=.3,...) {
  rownames(counts) <- make.unique(rownames(counts))
  p2 <- Pagoda2$new(counts,n.cores=parallel::detectCores()/2,...)
  p2$adjustVariance(plot=T,gam.k=10)
  p2$calculatePcaReduction(nPcs=100,n.odgenes=NULL,maxit=1000)
  
  opt=pc.select(p2);cat(paste0(opt," PCs retained\n"))
  
  p2$reductions$PCA=p2$reductions$PCA[,1:opt]
  cat("Computing UMAP... ")
  p2$embeddings$PCA$UMAP=doUMAP(p2$reductions$PCA,n_neighbors,min_dist)
  
  cat("done\n")
  p2$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
  p2$getKnnClusters(method=conos::leiden.community,type='PCA',name = "leiden")
  invisible(p2)
}


locplot <- function(emb,p2w,reads.filtered,pos=c(.1,.9)){
  emb=data.frame(emb)
  ggplot(emb)+geom_point(aes(x=UMAP1,y=UMAP2),size=1.5,color="black")+
    geom_point(aes(x=UMAP1,y=UMAP2,col=factor(reads.filtered$timepoints[rownames(emb)])),size=1)+
    scale_color_manual(values=p2w$cellmetadata$tp$palette[p2w$cellmetadata$tp$levels%in%levels(factor(reads.filtered$timepoints[rownames(emb)]))],
                       labels=p2w$cellmetadata$tp$levels[p2w$cellmetadata$tp$levels%in%levels(factor(reads.filtered$timepoints[rownames(emb)]))],
                       name="Locations")+theme_void()+theme(aspect.ratio = 1,legend.position=pos)
}


timplot <- function(emb,p2w,reads.filtered,pos=c(.1,.9),reverse=F){
  emb=data.frame(emb)
  emb$devtime=as.character(p2w$cellmetadata$time$data)
  colnames(emb)[1:2]=c("Dim1","Dim2")
  emb=emb[order(emb$devtime,decreasing = reverse),]
  ggplot(emb)+geom_point(aes(x=Dim1,y=Dim2),size=1.5,color="black")+
    geom_point(aes(x=Dim1,y=Dim2,col=devtime),size=1)+
    scale_color_manual(values=p2w$cellmetadata$time$palette,
                       labels=p2w$cellmetadata$time$levels,
                       name="Developmental time")+theme_void()+theme(aspect.ratio = 1,legend.position=pos)
}


cluplot <- function(emb,p2=NULL,p2w=NULL,pos=c(.1,.9)){
  emb=data.frame(emb)
  colnames(emb)=c("Dim1","Dim2")
  if (is.null(p2)){p2=p2w$originalP2object}
  pl=ggplot(emb)+geom_point(aes(x=Dim1,y=Dim2),size=1.5,color="black")+
    geom_point(aes(x=Dim1,y=Dim2,col=p2$clusters$PCA$leiden),size=1)+theme_void()+theme(aspect.ratio = 1,legend.position=pos)
  
  if(!is.null(p2w)){
    pl=pl+ scale_color_manual(values=p2w$cellmetadata$leiden$palette,
                              labels=p2w$cellmetadata$leiden$levels,
                              name="Leiden clusters")
  }
  return(pl)
  
}
