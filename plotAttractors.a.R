

plotAttractors.a<-function(data,rowOrder=FALSE,rowIndexes,
                           colOrder=FALSE,colIndexes){
  data.u<-data
  colnames(data.u)<-make.unique(colnames(data))
  rownames(data.u)<-make.unique(rownames(data))
  if ( colOrder == TRUE ) {
    data.u<-data.u[,rev(colIndexes)]
  }
  if ( rowOrder == TRUE ) {
    data.u<-data.u[rowIndexes,]
  }
  data.m<-melt(data.u)
  data.m$X2<-factor(data.m$X2,
                                levels = unique(data.m$X2),
                                ordered = TRUE)
  data.m$X1<-factor(data.m$X1,
                                levels = unique(data.m$X1),
                                ordered = TRUE)
  gg<-ggplot(data = data.m, aes(x=X1, y=X2, colours=value,fill=value )) 
  gg<-gg + geom_tile(color="white",size=0.1)  
  gg<-gg + scale_fill_gradient(low = "black", high = "gray", limit = c(0,1))
  gg<-gg + theme(text=element_text(size = 16),axis.text.x=element_text(angle = 90,vjust = 0.4,colour = "black",face = "bold",size = 15))
  gg<-gg + labs(x="",y="")
  gg<-gg + coord_equal()
  plot(gg)
}

#plotAttractors.a(binarizeMatrix(mutants,0) )