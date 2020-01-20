msefun <- function(yhat,y) {
    (y - yhat)^2
}

binfun <- function(yhat, y) {
    prob_min = 1e-05
    prob_max = 1 - prob_min
    yhat = pmin(pmax(yhat, prob_min), prob_max)
    y <- matrix(c(1 - y, y), ncol = 2)
    lp = y[, 1] * log(1 - yhat) + y[, 2] * log(yhat)
    ly = log(y)
    ly[y == 0] = 0
    ly = drop((y * ly) %*% c(1, 1))
    2 * (ly - lp)
}

poifun <- function(eta, y) {
    deveta = y * eta - exp(eta)
    devy = y * log(y) - y
    devy[y == 0] = 0
    2 * (devy - deveta)
}

# helper function for plotting CV curve: draws the error bars
error.bars <- function(x, upper, lower, width = 0.02, ...) {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

# compute gradient for cox model
coxgrad=function(f,time,d,w,eps=0.00001){
    ### f is fitted function from glmnet at a particular lambda
    ### time is death or censoring time
    ### d is death indicator; d=0 means censored, d=1 means death
    ### w is a weight vector of non-negative weights, which will be normalized to sum to 1
    if(missing(w))w=rep(1,length(f))
    w=w/sum(w)
    f=scale(f,TRUE,FALSE)#center f so exponents are not too large
    time=time-d*eps#break ties between death times and non death times, leaving tied death times tied
    o=order(time)
    ef=exp(f)[o]
    time=time[o]
    d=d[o]
    w=w[o]
    rskden=rev(cumsum(rev(ef*w))) ##reverse order inside;last guy is in all the risk sets
    ### See if there are dups in death times
    dups=fid(time[d==1],seq(length(d))[d==1])
    dd=d
    ww=w
    ### next code replaces each sequence of tied death indicators by a new
    ### sequence where only the first is a 1 and the rest are zero. This
    ### makes the accounting in the following step work properly we also
    ### sums the weights in each of the tied death sets, and assign that
    ### weight to the first
    if(!is.null(ties<-dups$index_ties)){
        dd[unlist(ties)]=0
        dd[dups$index_first]=1
        wsum=sapply(ties,function(i,w)sum(w[i]),ww)
        tie1=sapply(ties,function(i)i[1])
        ww[tie1]=wsum
    }
    ### Get counts over risk sets at each death time
    rskcount=cumsum(dd)#this says how many of the risk sets each observation is in; 0 is none
    ### We now form partial sums of the 1/den just at the risk sets
    rskdeninv=cumsum((ww/rskden)[dd==1])
    ### pad with a zero, so we can index it
    rskdeninv=c(0,rskdeninv)
    ### compute gradient for each obs
    #   grad=(d-rskdeninv[rskcount+1]*ef)*w   # rob changed this
    grad=(d-rskdeninv[rskcount+1]*ef)*w*length(f)
    grad[o]=grad
    grad
}

# helper function for coxgrad
fid <- function(x, index) {
    ### Input:
    ### x is a sorted vector of death times
    ### index is vector of indices of this set
    ### Output:
    ### index of first member of every death set as they appear in sorted list
    ### list of ties for each element of index, in the case of two or more ties;
    ## if no ties, this list is NULL
    idup <- duplicated(x)
    if(!any(idup)) list(index_first=index,index_ties=NULL)  # no ties
    else {
        ndup=!idup  # the first for each death time
        xu=x[ndup] # first death times
        index_first=index[ndup]
        ities=match(x,xu)
        index_ties=split(index,ities)
        nties=sapply(index_ties,length)
        list(index_first=index_first,index_ties=index_ties[nties>1])
    }
}
