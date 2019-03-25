# The code below was written by S. Guindon

###################################################################################################
###################################################################################################

## Calculate expected posterior error for all ambiguity levels. Each level has a
## corresponding value of alpha
which.min.error <- function(alpha,p)
{
    n = length(p);

    sorted.p = sort(p,decreasing=TRUE);
    cdf.sorted.p = cumsum(sorted.p);
    pme = 1:(n-1);
    
    for(i in 1:(n-1)) ## for the different levels of ambiguity
    {
        a = alpha[i];
        b = (n-1)/(n-i) - a * i/(n-i);
        pme[i] = a + (b-a)*(1-cdf.sorted.p[i]);                
    }
    
    idx.min = which.min(pme); ## ambiguity level that minimizes posterior error
    
    ## binary encoding of state that minimises posterior error
    ordered.p = order(p,decreasing=TRUE);
    best.state.int = 0;
    for(i in 1:idx.min)
    {
        best.state.int = best.state.int + 2^(ordered.p[i]-1);
    }
    return(best.state.int);
}

###################################################################################################
###################################################################################################
## Returns binary vector corresponding to integer i
int.to.bit <- function(i)
{
    return(as.integer(intToBits(i)));
}

###################################################################################################
###################################################################################################
## Returns set of character states where bit > 0 (i.e. bit==1)
bit.to.state <- function(bit,n)
{
    if(n == 4)
    {
        states = c("A","C","G","T");
    }
    else
    {
        states = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V");        
    }
    return(states[bit[1:n]>0]);
}

###################################################################################################
###################################################################################################

state.to.bit <- function(state,n)
{
    # state is a string
    if(n == 4)
    {
        return(as.numeric(lapply(c("A","C","G","T"),grepl,state)));
    }
    else
    {
        return(as.numeric(lapply(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),grepl,state)));
    }
}

###################################################################################################
###################################################################################################


mpee <- function(p,mesh.size)
{
    n = length(p);
    alpha = 1:n;
    best.state = 1:mesh.size;
    
    for(k in 0:mesh.size)
    {
        for(i in 1:n) { alpha[i] = (i-1)/i * (k/mesh.size); }
        best.state[k+1] = which.min.error(alpha,p);
    }
 
    z=as.integer(names(sort(table(best.state),decreasing=TRUE)[1]));
    return(paste(sort(bit.to.state(int.to.bit(z),n),decreasing=FALSE),collapse=""));
}

###################################################################################################
###################################################################################################

brier.score <- function(est.true.p,n)
{
                                        # est and true are strings
    est = as.character(est.true.p[1]);
    true = as.character(est.true.p[2]);
    p = as.numeric(as.vector(est.true.p[3:(3+n-1)]));
    
    p.true = state.to.bit(true,n);

    p.est = state.to.bit(est,n);
    p.est = p.est * p;
    p.est = p.est / sum(p.est);

    
    diff = (p.est - p.true)^2;

    ## return(sqrt(sum(diff)));
    return(sum(diff));
    
}

###################################################################################################
###################################################################################################

brier <- function(p)
{
    n = length(p);
    brier.score = 1:n;
    sorted.p = sort(p,decreasing=TRUE);

    
    for(k in 1:n) # calculate Brier score for each ambiguity level
    {
        brier.score[k] = 0.0;
        for(i in 1:k)
        {
            brier.score[k] = brier.score[k] + (sorted.p[i] - 1./k)^2;
        }
        if(k < n)
        {
            for(i in (k+1):n)
            {
                brier.score[k] = brier.score[k] + sorted.p[i]^2;
            }
        }
    }

    idx.min = which.min(brier.score); ## ambiguity level that minimizes Brier score
    
    ordered.p = order(p,decreasing=TRUE);
    ## binary encoding of state that minimises Brier score
    best.state.int = 0;
    for(i in 1:idx.min) { best.state.int = best.state.int + 2^(ordered.p[i]-1); }

    return(paste(sort(bit.to.state(int.to.bit(best.state.int),n),decreasing=FALSE),collapse=""));
}

###################################################################################################
###################################################################################################

map <- function(p)
{
    n = length(p);
    return(idx.to.state(which.max(p),n));
}

###################################################################################################
###################################################################################################
# Diff criterion
map.plus <- function(p,eps)
{
    n = length(p);

    sort.p = sort(p,decreasing=TRUE);

    ordered.p = order(p,decreasing=TRUE);
    
    ## binary encoding of state that minimises map.plus score
    best.state.int = 2^(ordered.p[1]-1);
    i = 2;
    while(((sort.p[i-1] - sort.p[i]) < eps) && i <= n)
    {
        best.state.int = best.state.int + 2^(ordered.p[i]-1);
        i = i+1;
    }
    
    return(paste(sort(bit.to.state(int.to.bit(best.state.int),n),decreasing=FALSE),collapse=""));
}

###################################################################################################
###################################################################################################
# CumProb criterion
cum.prob <- function(p,eps)
{
    n = length(p);
    

    sort.p = sort(p,decreasing=TRUE);
    cum.sort.p = cumsum(sort.p);
    
    ordered.p = order(p,decreasing=TRUE);
    
    ## binary encoding of state that minimises map.plus score
    best.state.int = 0;
    i = 1;
    repeat
    {
        best.state.int = best.state.int + 2^(ordered.p[i]-1);
        if(cum.sort.p[i] > eps || i >= n) { break; }
        i = i+1;
    }
    
    return(paste(sort(bit.to.state(int.to.bit(best.state.int),n),decreasing=FALSE),collapse=""));

}

###################################################################################################
###################################################################################################
# Thresh criterion
prob.thresh <- function(p,thresh)
{
    n = length(p);
    
    
    sort.p    = sort(p,decreasing=TRUE);    
    ordered.p = order(p,decreasing=TRUE);
    
    best.state.int = 0;
    i = 1;
    while(sort.p[i] > thresh && i <= n)
    {
        best.state.int = best.state.int + 2^(ordered.p[i]-1);
        i = i+1;
    }
    
    return(paste(sort(bit.to.state(int.to.bit(best.state.int),n),decreasing=FALSE),collapse=""));

}

###################################################################################################
###################################################################################################
## true is the vector of all true ancestral states and est is the corresponding vector of ancestral
## states estimated with a given method.
accuracy <- function(true,est)
{
    z=apply(cbind(true,est),1,match.length);
    return(z);
}
match.length <- function(true.est)
{
    length(grep(true.est[1],true.est[2]));
}
precision <- function(est)
{
    z=apply(as.array(est),1,nchar);
    return(z);
}
brier.score.all <- function(true,est,p,n)
{
    z = apply(matrix(c(est,true,p),ncol=2+n),1,brier.score,n);
    return(z);
}
