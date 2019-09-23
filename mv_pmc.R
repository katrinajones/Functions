#' Multivariate phylogenetic montecarlo test
#'
#' Monte Carlo simulation for hypothesis testing of evolutionary models on multivariate data. Function based on \code{\link{pmc}} from pmc package but utilizing \code{\link{mvMORPH}} for multivariate data.
#'
#' Tests null hypothesis A against hypothesis B using Monte carlo simulations. The likelihood ratio of the original data given hypothesis A and B is calculated. Next data are simulated nsim times under both hypotheses. Null and test distributions are generated as the likelihood ratio of hypothesis B given simulation A and hypothesis A given simulation A; and hypothesis B given simulation B and hypothesis A given simulation B. This allows one to test the power and significance of the hypotheses.
#'
#' This function differs from the original \code{\link{pmc}} in several ways. Firstly, it tests multivariate data using the \code{\link{mvMORPH}} package. The 'func' and 'param' argument are passed directly to mvMORPH. Secondly, four separate trees are accepted. The tree for hypothesis A and hypothesis B will be the same, but may have different regime mapping if testing between multiregime models. The 'sim' and 'test' trees may differ if the user wishes to incorporate the effect of sampling into the hypothesis test. In this case, the data may be simulated over a larger tree (for example for a fossil analysis where sampling may be challenging), which will then be subset to the just the taxa in the test tree (and original data) for testing. This simulates the effect of subsampling on evolutionary patterns.
#'
#'
#' \code{\link{mv_pmc_multiPhylo}} enables incorporation of error due to phylogenetic uncertainty (e.g., branch length estimation). Instead of a single tree, and object of class "multiPhylo" is provided for both test and simulation trees. The likelihood ratio of the original data is calculated across the set of trees and the mean returned for the test. For each simulation, a tree is drawn at random from the population of simulation trees, which is tested against the equivalent test tree. It is assumed that the test trees are a taxonomic subset of the simulation tree population and are in the same order.
#'
#'
#' @param tree.simA Simulation tree for model A
#' @param tree.testA Test tree for model A
#' @param tree.simB Simulation tree for model B
#' @param tree.testB Test tree for model B
#' @param data data
#' @param funcA mvMORPH function
#' @param modelA mvMORPH model
#' @param funcB mvMORPH function
#' @param modelB mvMORPH model
#' @param nsim number of simulations
#' @param paramA list passed to param arguement
#' @param paramB list passed to param arguement
#' @param no_cores number of cores to be used with doParallel
#' @param jack conduct jacknife analysis?
#'
#' @return a list containing:
#' \item{lr}{Likelihood ratio based on original data}
#' \item{null}{distribution of likelihood ratios under null hypothesis}
#' \item{test}{distribution of likelihood ratios under test hypothesis}
#' \item{A}{original model A}
#' \item{B}{original model B}
#' \item{stats}{list of pvalue and power}
#'
#' @export
#'
#'@references Boettiger, Carl, Graham Coop, and Peter Ralph. "Is your phylogeny informative? Measuring the power of comparative methods." Evolution: International Journal of Organic Evolution 66.7 (2012): 2240-2251.
#'@references Clavel, Julien, Gilles Escarguel, and Gildas Merceron. "mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data." Methods in Ecology and Evolution 6.11 (2015): 1311-1319.
#'
#' @seealso \code{\link{mvBM}}, \code{\link{mvOU}}, \code{\link{mvSHIFT}}
#'
#' @import mvMORPH
#' @import ape
#' @import doParallel
#' @import parallel
#' @import foreach
#'
#' @examples
#' ## Simulated dataset
#' set.seed(14)
#' # Generating a random tree with 50 species
#' tree<-pbtree(n=50)
#'
#' #subsample to 25 taxa
#' tree.sub<-drop.tip(tree, tree$tip.label[sample(1:50,25)])
#'
#' #paint regimes
#' tree<-paintSubTree(tree, node=ape::getMRCA(tree, c("t26","t28")), state="2")
#' tree.sub<-paintSubTree(tree.sub, node=ape::getMRCA(tree.sub, c("t26","t28")),state="2")
# '
#' ## Simulate trait evolution according to a bivariate "BMM" model
#' # Number of traits
#' ntraits<-2
#' # Rates matrices for the "Forest" and the "Savannah" regimes
#' sigma<-list(Forest=matrix(c(1,0.5,0.5,1),2), Savannah=matrix(c(10,3,3,10),2))
#' # ancestral states for each traits
#' theta<-c(0,0)
#'
#' # Simulate
#' data<-mvSIM(tree.sub,nsim=1, model="BMM",param=list(ntraits=ntraits,sigma=sigma,theta=theta))
#'
#' #Hypothesis test
#' out.BMM<-myfunctions::mv_pmc(tree.simA=tree, tree.testA = tree.sub, tree.simB=tree, tree.testB = tree.sub, data=data,funcA="mvBM", modelA="BM1",funcB="mvBM", modelB="BMM", nsim = 10)
#'
#' #plot
#' plot(out.BMM)
#'
#' # Simulate null
#' data<-mvSIM(tree.sub,nsim=1, model="BM1",param=out.BMM$A)
#'
#' #Hypothesis test
#' out.BM<-myfunctions::mv_pmc(tree.simA=tree, tree.testA = tree.sub, tree.simB=tree, tree.testB = tree.sub, data=data,funcA="mvBM", modelA="BM1",funcB="mvBM", modelB="BMM", nsim = 10)
#'
#' #plot
#' plot(out.BM)
#'
#' #Calculate over multiple trees to account for uncertainty
#' mtree<-as.multiPhylo(c(tree,tree))
#' mtree.sub<-as.multiPhylo(c(tree.sub,tree.sub))
#'
#' #Hypothesis test
#' out.BM_multi<-mv_pmc_multiPhylo(tree.simA=mtree, tree.testA = mtree.sub, tree.simB=mtree, tree.testB = mtree.sub, data=data,funcA="mvBM", modelA="BM1",funcB="mvBM", modelB="BMM", nsim = 10)
#'
#' #get summary statistics for parameters over multiple trees
#' summary(out.BM_multi)#No variation because trees identical
#' summary(out.BM_multi, jack=T)#calculate over jacknifed data

mv_pmc<-function (tree.simA, tree.testA,tree.simB, tree.testB, data, funcA="mvBM", modelA="BM1",
                  funcB="mvOU", modelB="OU1", nsim = 2, paramA = list(),
                  paramB = list(),  no_cores = parallel::detectCores(), jack=T)

{
  if(length(setdiff(rownames(data), tree.testA$tip.label))>0|length(setdiff(rownames(data), tree.testB$tip.label))>0){stop("rownames of data must match tip labels")}

  #sort data - order of tips must match
  dataA<-data[tree.testA$tip.label,]#sort
  dataB<-data[tree.testB$tip.label,]#sort

  #set up the parallelization
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl, cores=no_cores)

  #fit original data to get observed likelihood ratio
  fit_A<-do.call(funcA, list(tree.testA, dataA, model=modelA, param = paramA, diagnostic=F, echo=F))
  fit_B<-do.call(funcB, list(tree.testB, dataB, model=modelB, param = paramB, diagnostic=F, echo=F))
  #calculate likelihood
  lr_orig <- -2 * (logLik(fit_A) - logLik(fit_B))

  #jack-knife
  if(isTRUE(jack)){
    #leave one out analysis
    datajack<-as.list(c(1:nrow(data)))
    row<-0
    datajack<-lapply(datajack, function(x) {row<<-row+1;x<-data[-row,]})

    jackoutA<-foreach(x=1:nrow(data),.packages = c("mvMORPH")) %dopar%
    {
      tree.jackA<-ape::drop.tip(tree.testA,tree.testA$tip.label[!(tree.testA$tip.label%in%rownames(datajack[[x]]))])#drop tips
      datajack[[x]]<-datajack[[x]][tree.jackA$tip.label,]
      do.call(funcA, list(tree.jackA, datajack[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
    }
    jackoutB<-foreach(x=1:nrow(data),.packages = c("mvMORPH")) %dopar%
    {
      tree.jackB<-ape::drop.tip(tree.testB,tree.testB$tip.label[!(tree.testB$tip.label%in%rownames(datajack[[x]]))])
      datajack[[x]]<-datajack[[x]][tree.jackB$tip.label,]
      do.call(funcB, list(tree.jackB, datajack[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))
    }

    jack_dist = (sapply(jackoutB, function(x) x$AICc))-(sapply(jackoutA, function(x) x$AICc))
    jack_lik_dist = -2*((sapply(jackoutA, function(x) x$LogLik))-(sapply(jackoutB, function(x) x$LogLik)))


    AICcdiff=median(jack_dist)
    Likjack=median(jack_lik_dist)
    CI.lik<-quantile(jack_lik_dist, c(0.025,0.975))
    paic<-length(which(jack_dist>=0))/length(jack_dist)
  }else{AICcdiff=NULL; paic=NULL; Likjack=NULL; CI.lik=NULL;jackoutA=NULL; jackoutB=NULL}

  #Now simulate to produce test distribution
  A_sims<-mvSIM(tree.simA, nsim=nsim,  param=fit_A)
  B_sims<-mvSIM(tree.simB, nsim=nsim,  param=fit_B)

  #subset the sims to the test tree
  A_sims<-lapply(A_sims, function(x){x<-x[rownames(x)%in%tree.testA$tip.label,]})
  B_sims<-lapply(B_sims, function(x){x<-x[rownames(x)%in%tree.testB$tip.label,]})


  #begin simulations
  AA <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
    do.call(funcA, list(tree.testA, A_sims[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
  AB <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
    do.call(funcB, list(tree.testB, A_sims[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))
  BA <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
    do.call(funcA, list(tree.testA, B_sims[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
  BB <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
    do.call(funcB, list(tree.testB, B_sims[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))

  #end cluster
  stopCluster(cl)

  #calculate test values
  null_dist = -2 * (sapply(AA, logLik) - sapply(AB, logLik))
  test_dist = -2 * (sapply(BA, logLik) - sapply(BB, logLik))

  #stats
  pval<-length(which(null_dist>=lr_orig))/length(null_dist)
  if(!is.null(Likjack)){
  pval_jack<-length(which(null_dist>=Likjack))/length(null_dist)
  }else{pval_jack=NULL}
  cutoff<-quantile(null_dist,probs=c(0.95))
  power<-length(which(test_dist>=cutoff))/length(test_dist)
  stats<-list(lnlik=lr_orig,pval=pval, cutoff=cutoff, power=power, AICcdiff=AICcdiff, paic=paic,Likjack=Likjack, pval_jack=pval_jack, CI.lik=CI.lik)

  #output
  out <- list(lr = lr_orig, null = null_dist, test = test_dist,
              A = fit_A, B = fit_B, stats=stats, jackA=jackoutA, jackB=jackoutB)
  class(out) <- "pmc"
  return(out)
}


#' @describeIn mv_pmc Apply to multiple trees
#'
#' @export

mv_pmc_multiPhylo<-function (tree.simA, tree.testA,tree.simB, tree.testB, data, funcA="mvBM", modelA="BM1",
                  funcB="mvOU", modelB="OU1", nsim = 2, paramA = list(),
                  paramB = list(),  no_cores = parallel::detectCores(), jack=T)

{

  if(!class(tree.testA)=="multiPhylo"){stop("tree of class multiPhylo required for tree.testA")}
  if(!class(tree.testB)=="multiPhylo"){stop("tree of class multiPhylo required for tree.testB")}
  if(!class(tree.simA)=="multiPhylo"){stop("tree of class multiPhylo required for tree.simA")}
  if(!class(tree.simB)=="multiPhylo"){stop("tree of class multiPhylo required for tree.simB")}

  if(length(setdiff(rownames(data), tree.testA[[1]]$tip.label))>0|length(setdiff(rownames(data), tree.testB[[1]]$tip.label))>0){stop("rownames of data must match tip labels")}

  #sort data - order of tips must match
  dataA<-data[tree.testA[[1]]$tip.label,]#sort
  dataB<-data[tree.testB[[1]]$tip.label,]#sort

  #set up the parallelization
  cl <- makeCluster(no_cores)
  doParallel::registerDoParallel(cl, cores=no_cores)

  #index trees to sims
  ind<-sample(1:length(tree.simA), nsim, replace=T)

  #original data
      fit_A <- foreach(x=1:length(tree.testA),.packages = c("mvMORPH")) %dopar%
      do.call(funcA, list(tree.testA[[x]], data, model=modelA, param = paramA, diagnostic=F, echo=F))
    fit_B <- foreach(x=1:length(tree.testB),.packages = c("mvMORPH")) %dopar%
      do.call(funcB, list(tree.testB[[x]], data, model=modelB, param = paramB, diagnostic=F, echo=F))

    lr_orig <- median(-2 * ((sapply(fit_A, logLik))-(sapply(fit_B, logLik))))


    #jack-knife
    if(isTRUE(jack)){
      #leave one out analysis
      datajack<-as.list(c(1:nrow(data)))
      row<-0
      datajack<-lapply(datajack, function(x) {row<<-row+1;x<-data[-row,]})
      #index trees to jacks
      indj<-sample(1:length(tree.simA), nrow(data), replace=T)

      jackoutA<-foreach(x=1:nrow(data),.packages = c("mvMORPH")) %dopar%
      {
        tree.jackA<-ape::drop.tip(tree.testA[[indj[x]]],tree.testA[[indj[x]]]$tip.label[!(tree.testA[[indj[x]]]$tip.label%in%rownames(datajack[[x]]))])#drop tips
        datajack[[x]]<-datajack[[x]][tree.jackA$tip.label,]
        do.call(funcA, list(tree.jackA, datajack[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
      }
      jackoutB<-foreach(x=1:nrow(data),.packages = c("mvMORPH")) %dopar%
      {
        tree.jackB<-ape::drop.tip(tree.testB[[indj[x]]],tree.testB[[indj[x]]]$tip.label[!(tree.testB[[indj[x]]]$tip.label%in%rownames(datajack[[x]]))])
        datajack[[x]]<-datajack[[x]][tree.jackB$tip.label,]
        do.call(funcB, list(tree.jackB, datajack[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))
      }

      jack_dist = (sapply(jackoutB, function(x) x$AICc))-(sapply(jackoutA, function(x) x$AICc) )
      jack_lik_dist = -2* ((sapply(jackoutA, function(x) x$LogLik))-(sapply(jackoutB, function(x) x$LogLik) ))
      AICcdiff<-median(jack_dist)
      Likjack<-median(jack_lik_dist)
      CI.lik<-quantile(jack_lik_dist, c(0.025,0.975))
      paic<-length(which(jack_dist>=0))/length(jack_dist)
    }else{AICcdiff=NULL; paic=NULL; Likjack=NULL; CI.lik=NULL; jackoutA=NULL; jackoutB=NULL}

 #simulations
    A_sims <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
      do.call("mvSIM", list(tree.simA[[ind[x]]], nsim=1,  param=fit_A[[ind[x]]]))
    B_sims <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
      do.call("mvSIM", list(tree.simB[[ind[x]]], nsim=1,  param=fit_B[[ind[x]]]))

    #subset the sims to the test tree
    A_sims<-lapply(A_sims, function(x){x<-x[rownames(x)%in%tree.testA[[1]]$tip.label,]})
    B_sims<-lapply(B_sims, function(x){x<-x[rownames(x)%in%tree.testB[[1]]$tip.label,]})


      #begin simulations
    AA <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
     do.call(funcA, list(tree.testA[[ind[x]]], A_sims[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
    AB <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
      do.call(funcB, list(tree.testB[[ind[x]]], A_sims[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))
    BA <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
      do.call(funcA, list(tree.testA[[ind[x]]], B_sims[[x]], model=modelA, param = paramA, diagnostic=F, echo=F))
    BB <- foreach(x=1:nsim,.packages = c("mvMORPH")) %dopar%
      do.call(funcB, list(tree.testB[[ind[x]]], B_sims[[x]], model=modelB, param = paramB, diagnostic=F, echo=F))


  #end cluster
  stopCluster(cl)

  #calculate test values
  null_dist = -2 * (sapply(AA, logLik) - sapply(AB, logLik))
  test_dist = -2 * (sapply(BA, logLik) - sapply(BB, logLik))

  #stats
  pval<-length(which(null_dist>=lr_orig))/length(null_dist)
  if(!is.null(Likjack)){
    pval_jack<-length(which(null_dist>=Likjack))/length(null_dist)
  }else{pval_jack=NULL}
  cutoff<-quantile(null_dist,probs=c(0.95))
  power<-length(which(test_dist>=cutoff))/length(test_dist)
  stats<-list(lnlik=lr_orig,pval=pval, cutoff=cutoff, power=power, AICcdiff=AICcdiff, paic=paic,Likjack=Likjack, pval_jack=pval_jack, CI.lik=CI.lik)

  #output
  out <- list(lr = lr_orig, null = null_dist, test = test_dist,
              A = fit_A, B = fit_B, stats=stats, jackA=jackoutA, jackB=jackoutB)
  class(out) <- "pmc"
  return(out)
}

#' Plot for pmc object
#'
#'
#' @param pmcobj Object generated by pmc or mv_pmc
#'
#' @return plot
#' @export
#'
#'
plot.pmc<-function(pmcobj, ...){

  if(!is.null(pmcobj$stats$Likjack)){
  lr<-pmcobj$stats$Likjack
  pval<-pmcobj$stats$pval_jack}else{
    lr<-pmcobj$lr
    pval<-pmcobj$pval
  }

  d1<-density(pmcobj$null)
  d2<-density(pmcobj$test)
  rangex<-range(c(d1$x,d2$x))
  rangey<-range(c(d1$y,d2$y))
  plot(d1,xlim=rangex,ylim=rangey, col="red",xlab=paste0("lnLik. diff.:",round(lr, digits=1), "   P-val:", round(pval, digits=2), "   Power:", round(pmcobj$stats$power, digits=2)), ...)
  polygon(d1, col=rgb(1,0,0,alpha=0.3))
  lines(d2, col="blue")
  polygon(d2, col=rgb(0,0,1,alpha=0.3))

  #add true value
  if(!is.null(pmcobj$stats$Likjack)){
    abline(v=pmcobj$stats$Likjack)
    abline(v=pmcobj$stats$CI.lik[1],lty=3)
    abline(v=pmcobj$stats$CI.lik[2],lty=3)
  }else{
    abline(v=pmcobj$lr)}
  rect(rangex[1],0,pmcobj$stats$cutoff,rangey[2], col=rgb(0,0,0,alpha=0.1), border=NA)

}

#' Summarize mv_pmc object
#'
#' Calculates median and ranges of parameters across trees
#'
#' @param pmcobj mv_pmc object
#' @param const list of constraints indicating models to be excluded
#' @param jack calculate over jacknifed data?
#'
#' @return median and ranges
#' @export
#'
#'
summary.pmc<-function(pmcobj, const=list(), jack=F){

  sumout<-list()

  for(hyp in 1:2){

    hypname<-c("A","B")[hyp]
    if(isTRUE(jack)){
      hypname<-c("jackA","jackB")[hyp]
    }

    x<-pmcobj[hypname][[1]]
    if(is.null(names(x[[1]]))){next()}#skip if only one

    elems<-names(x[[1]])
    elems<-elems[!elems%in%c("param", "llik", "convergence", "hess.values")]
    output<-list()

    for(b in 1:length(elems)){

      i<-which(names(x[[1]])==elems[b])

      if(is.matrix(x[[1]][[i]])==F){
        val <- sapply(1:length(x), function(a) x[[a]][[i]])
        #apply constraints
        if(elems[b]%in%names(const)){
          ind<-const[which(names(const)==elems[b])][[1]]
          keep<-which(ind[1]>val|val>ind[2])
          val<-val[-keep]
        }
        val<-c(median=median(val), range=range(val))
        output[[i]]<-val
      }
      if(is.matrix(x[[1]][[i]])==T){
        val<-lapply(1:length(x), function(a) x[[a]][[i]])
        #apply constraints
        if(elems[b]%in%names(const)){
          ind<-const[which(names(const)==elems[b])][[1]]
          keep<-lapply(val, function(x) !any(ind[1]>x|x>ind[2]))
          val<-val[which(keep==T)]
        }
        #calculate median
        val<-array(unlist(val), dim = c(nrow(val[[1]]), ncol(val[[1]]), length(val)), dimnames = list(rownames(val[[1]]), colnames(val[[1]]), c(1:length(val))))
        val<-list(median=apply(val,c(1,2), median), range1=apply(val,c(1,2), range)[1,,],range2=apply(val,c(1,2), range)[2,,], CI1=apply(val,c(1,2), quantile, probs=c(0.025)),CI2=apply(val,c(1,2), quantile, probs=c(0.975)))
        output[[i]]<-val
      }
    }
    names(output)<-elems

    sumout[hypname]<-list(output)
  }

  return(sumout)

}
