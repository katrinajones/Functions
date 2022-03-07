#' Concatenate data in a geomorph data frame
#'
#' @param gm.data.frame geomorph data frame
#' @param by ID variable with which to concatenate e.g., ID variable
#'
#' @return concatenated data frame
#' @export
#'
#'
concatgeom<-function(gm.data.frame,by){
  for(i in 1:length(gm.data.frame)){
    item<-gm.data.frame[[i]]

 ##For landmarks/coords
       if((length(dim(item))==3)){
             new<-NULL
      for(a in 1:length(levels(by))){
      IDvar<-levels(by)[a]
      dat<-item[,,which(by==IDvar)]#grab data
      con<-NULL
            for(b in 1:dim(dat)[3]){ #bind together
        con<-rbind(con,dat[,,b])
      }
      new<-abind::abind(new,con,along=3)#Make new array
    }
dimnames(new)[[3]]<-levels(by)
assign(attributes(gm.data.frame)$names[i],new)##Create unique identifier
    }

###For PC scores
    if((length(dim(item))==2)){
    new<-NULL
    for(a in 1:length(levels(by))){
      IDvar<-levels(by)[a]
      dat<-item[which(by==IDvar),]#grab data
      con<-utils::stack(as.data.frame(t(dat)))[,1]
      new<-rbind(new,con)
    }
    dimnames(new)[[1]]<-levels(by)
    assign(attributes(gm.data.frame)$names[i],new)##Create unique identifier
    }

####For factors
    if(is.factor(item)==T){
      new<-NULL
       for(a in 1:length(levels(by))){
         IDvar<-levels(by)[a]
         dat<-item[which(by==IDvar)]#grab data
         con<-dat[1]
        new<-c(new,paste(con))
       }
         new<-as.factor(new)
         new<-droplevels(new)
         assign(attributes(gm.data.frame)$names[i],new)##Create unique identifier
      }

  ####For vectors
  if(is.vector(item)==T){
    new<-NULL
    for(a in 1:length(levels(by))){
      IDvar<-levels(by)[a]
      dat<-item[which(by==IDvar)]#grab data
      con<-mean(dat)
      new<-c(new,con)
    }
    assign(attributes(gm.data.frame)$names[i],new)##Create unique identifier
  }
    }

      ##Create new geomorph data frame
  new.df<-NULL
  name<-attributes(gm.data.frame)$names[1]
  obj<-get(name)
  new.df<-geomorph::geomorph.data.frame(obj)
    for(i in 2:length(gm.data.frame)){
    name<-attributes(gm.data.frame)$names[i]
    obj<-get(name)
    new.df<-geomorph::geomorph.data.frame(new.df,obj)
  }
    attributes(new.df)$names<-attributes(gm.data.frame)$names
    return(new.df)
      }

#' Calculate mean shape from concatenated shape
#'
#' @param coords Concatenated shape
#' @param nobj number of objects
#'
#' @return mean shape
#' @export
#'
#'
meanconc<-function(coords,nobj){

  nlands<-nrow(coords)/nobj

  newcoords<-array(data=NA, dim=c(nlands, 3, nobj))
  for(i in 1:nobj){
    start<-((nlands*i)-nlands)+1
    end<-(start + nlands) -1
    newcoords[,,i]<-coords[start:end,]
  }

  mconc<-geomorph::mshape(newcoords)
  rownames(mconc)<-rownames(coords)[1:nlands]

  return(list(coords=newcoords, mean=mconc))
}

#Blake functions

#' Project a new set of coordinates into morphospace
#'
#'   From Blake Dickson
#'
#'
#' @param coords 2d or 3d coordinates of new shape
#' @param pca pca object from prcomp
#' @param gpa gpa object or geomorph data frame
#'
#' @return principal components scores
#' @export
#'
Coords2PC<-function(coords,pca,gpa){

  k <- dim(gpa$coords)[2]
  p <- dim(gpa$coords)[1]
  n <- dim(gpa$coords)[3]

  Mshape <- mshape(gpa$coords)
  Rotation <- pca$rotation

  x <- array(data= NA,dim=c(p,k,1))

  if (dim(coords)[2]==k){
    if (length(dim(coords))==2) {

      dim(coords)[3]<-1
      x <- two.d.array(coords)

    }

    if (length(dim(coords))==3) {

      x <- two.d.array(coords)

    }
  }
  if (dim(coords)[2]>k) x <- coords

  w <- x %*% Rotation
  return(w)

}

#' Calculate shapes from morphospace
#'
#' from Blake Dickson
#'
#' @param pc pc scores
#' @param pca pca object from prcomp
#' @param gpa gpa object
#'
#' @return 2d or 3d coordinates
#' @export
#'
PC2Coords<-function(pc,pca,gpa){

  Mshape <- mshape(gpa$coords)
  Rotation <- pca$rotation
  k <- dim(gpa$coords)[2]
  p <- dim(gpa$coords)[1]
  n <- dim(gpa$coords)[3]

  w <- as.matrix( (pc %*%  t(Rotation)) + pca$center)
  w <- arrayspecs(w,p,k)
  w[1:3,1:3,1]
  #coords[1:3,1:3,1]

  return(w)
}
