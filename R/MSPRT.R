
if(getRversion() >= "2.15.1"){
  utils::globalVariables("k")
}

      
compmerge.list = function(l1, l2){
  l = vector("list",length(l1))
  for (i in 1:length(l)) {
    l[[i]] = c(l1[[i]], l2[[i]])
  }
  return(l)
}


### objective function for Bernoulli family

##  arguments
# p     : value in a region to be optimized over
# delta : threshold of bayes factor to find UMPBT
# n.obs : sample size
# p0    : null value

## returns value of the objective function at p

obj.func.ber = function( p, delta, n.obs, p0){
  num = log(delta) - (n.obs*(log(1-p) - log(1- p0)))
  den = log(p/(1-p)) - log(p0 /(1- p0))
  t = num/den
  
  return(t)
}


### evaluates (LHS - c_0) as a function of \gamma as in appendix of sequential.pdf

##  arguments
# delta         : threshold for bayes factor to find UMPBT
# side          : direction of H1(right/left)
# n.obs         : sample size
# p0            : null value
# opt.interval  : interval containing p under H1
# root          : c_0 to be specific as in supplemental file

## returns (LHS - c_0) as a function of \gamma as in supplemental file

find.thresold.ber = function( delta, side = "right", n.obs, p0 = 0.5,
                              opt.interval, root){
  
  if(side=="right"){
    
    if(missing(opt.interval)==T){
      opt.interval = c(p0, 1)
    }
    
    out = optimize( f=obj.func.ber, interval = opt.interval,
                    delta = delta, n.obs = n.obs, p0 = p0)
  }else if(side=="left"){
    
    if(missing(opt.interval)==T){
      opt.interval = c(0, p0)
    }
    
    out = optimize( f=obj.func.ber, interval = opt.interval, maximum = T,
                    delta = delta, n.obs = n.obs, p0 = p0)
  }
  
  return(out$objective -root)
}



### solving \gamma by matching UMP rejection region to find umpbt

##  arguments
# side  : direction of H1(right/left)
# type1 : specified type 1 error
# n.obs : sample size
# p0    : null value

## returns \gamma for which ump rejection region is matched


ump.match.ber = function( side = "right", type1 = 0.005, n.obs, p0 = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n.obs, p0) +1
    
    solve.out = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = "right", p0=p0,
                         opt.interval = c(p0, 1), root = (c0 -1) )
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n.obs, prob = p0)
    quant.prob = pbinom( q = c0, size = n.obs, prob = p0)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    
    solve.out = nleqslv( x=3, fn = find.thresold.ber, n.obs = n.obs, side = side, p0 = p0,
                         opt.interval = c(0, p0), root = (c0 +1) )
  }
  
  return(solve.out$x)
}



### finding mixture UMPBT alternative prior for exact proportion test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n.obs   : sample size
# null    : null value

## returns UMPBT alternative prior by matching UMP rejection region

point.umpbt.ber = function( side = "right", type1 = 0.005, n.obs, null = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n.obs, null) +1
    solve.out = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                         opt.interval = c(null, 1), root = (c0 -1) )
    out = optimize( f = obj.func.ber, interval = c(null, 1),
                    delta = solve.out$x, n.obs = n.obs, p0 = null)
    umpbt = round(out$minimum, 4)
    
    return(umpbt)
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n.obs, prob = null)
    quant.prob = pbinom( q = c0, size = n.obs, prob = null)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    
    solve.out = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                         opt.interval = c(0, null), root = (c0 +1) )
    out = optimize( f = obj.func.ber, interval = c(0, null), maximum = T,
                    delta = solve.out$x, n.obs = n.obs, p0 = null)
    umpbt = round(out$maximum,4)
    
    return(umpbt)
  }
  
}

find.umpbt.ber = function( side = "right", type1 = 0.005, n.obs, null = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n.obs, null) +1
    quant.prob = 1- pbinom( q = (c0 -1), size = n.obs, prob = null)
    
    solve.out.r = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                           opt.interval = c(null, 1), root = (c0 -1) )
    out.r = optimize( f = obj.func.ber, interval = c(null, 1),
                      delta = solve.out.r$x, n.obs = n.obs, p0 = null)
    umpbt.r = round(out.r$minimum, 4)
    
    solve.out.l = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                           opt.interval = c(null, 1), root = (c0 -2) )
    out.l = optimize( f = obj.func.ber, interval = c(null, 1),
                      delta = solve.out.l$x, n.obs = n.obs, p0 = null)
    umpbt.l = round(out.l$minimum, 4)
    
    psi = (type1 - quant.prob)/dbinom(x=(c0 -1), size = n.obs, prob = null)
    psi = round(psi, 5)
    
    return(list("u"=c(umpbt.l, umpbt.r), "psi" = psi))
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n.obs, prob = null)
    quant.prob = pbinom( q = c0, size = n.obs, prob = null)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    quant.prob = pbinom( q = c0, size = n.obs, prob = null)
    
    
    solve.out.l = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                           opt.interval = c(0, null), root = (c0 +1) )
    out.l = optimize( f = obj.func.ber, interval = c(0, null), maximum = T,
                      delta = solve.out.l$x, n.obs = n.obs, p0 = null)
    umpbt.l = round(out.l$maximum, 4)
    
    solve.out.r = nleqslv( x=3, fn = find.thresold.ber, n.obs=n.obs, side = side, p0=null,
                           opt.interval = c(0, null), root = (c0 +2) )
    out.r = optimize( f = obj.func.ber, interval = c(0, null), maximum = T,
                      delta = solve.out.r$x, n.obs = n.obs, p0 = null)
    umpbt.r = round(out.r$maximum, 4)
    
    psi = (type1 - quant.prob)/dbinom(x=(c0 +1), size = n.obs, prob = null)
    psi = round(psi, 5)
    
    return(list("u"=c(umpbt.r, umpbt.l), "psi" = psi))
  }
  
}



### finding UMPBT for z-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n.obs   : sample size
# null    : null value
# sigma0  : known s.d.

## returns UMPBT alternative by matching UMP rejection region

find.umpbt.norm = function( side = "right", type1 = 0.005, n.obs, null = 0, sigma0 = 1){
  
  z.alpha = qnorm( p = type1, lower.tail = F)
  
  if(side=="right"){
    umpbt = null + ((sigma0 * z.alpha)/(sqrt(n.obs)))
  }else if(side=="left"){
    umpbt = null - ((sigma0 * z.alpha)/(sqrt(n.obs)))
  }
  
  umpbt = round(umpbt, 4)
  
  return(umpbt)
}


### finding UMPBT for t-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n.obs   : sample size
# null    : null value
# s       : sample s.d.

## returns UMPBT alternative by matching UMP rejection region

find.umpbt.t = function( side = "right", type1 = 0.005, n.obs, null = 0, obs, s){
  
  if(missing(s)==T){
    s = sd(obs)
  }
  
  t.alpha = qt( type1, df = (n.obs -1), lower.tail = F)
  
  if(side=="right"){
    umpbt = null + ((s * t.alpha)/(sqrt(n.obs)))
  }else if(side=="left"){
    umpbt = null - ((s * t.alpha)/(sqrt(n.obs)))
  }
  
  umpbt = round(umpbt, 4)
  
  return(umpbt)
}



### likelihood ratio(LR) for proportion test

##  arguments
# m         : sample size
# suff.stat : value of sufficient statistic
# null      : null value
# alt       : alt value

## retuns LR

LR.ber = function( m, suff.stat, null = 0.5, alt){
  out = (((1- alt)/(1- null))^m) * (((alt*(1- null))/(null*(1- alt)))^suff.stat)
  return(out)
}


### likelihood ratio(LR) for z-test

##  arguments
# m         : sample size
# suff.stat : value of sufficient statistic
# null      : null value
# alt       : alt value
# sigma0    : known s.d.

## retuns LR

LR.norm = function( m, suff.stat, null = 0, alt, sigma0 = 1){
  t1 = (suff.stat*(alt - null)) - ((m*((alt^2) - (null^2)))/2)
  out = exp(t1/(sigma0^2))
  return(out)
}


### likelihood ratio(LR) for t-test

##  arguments
# m         : sample size
# suff.stat : value of sufficient statistic
# null      : null value
# alt       : alt value
# s         : sample s.d.

## retuns LR

LR.t = function( m, suff.stat, null = 0, alt, s){
  t0 = (((suff.stat/m) -null)/s)
  t1 = (((suff.stat/m) -alt)/s)
  t = ((1+ ((m/(m-1))*(t0 ^2)) )/(1+ ((m/(m-1))*(t1 ^2)) ))^(m/2)
  return(t)
}



### computes (type-2 error -root) for proportion test

##  arguments
# alt       : alt value
# side      : direction of H1(right, left or both)
# null.val  : null value
# N         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.ber = function( alt, side = "right", null.val = 0.5, N, type1 = 0.005, root = 0){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), N, null.val) +1
    t = pbinom( q= (c0 -1), size = N, prob = alt)
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = N, prob = null.val)
    quant.prob = pbinom( q = c0, size = N, prob = null.val)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    t = 1- pbinom( q= c0, size = N, prob = alt)
    
  }
  
  return((t- root))
}



### computes (type-2 error -root) for z-test

##  arguments
# side      : direction of H1(right, left or both)
# alt       : alt value
# null.val  : null value
# sigma0    : known s.d.
# N         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.norm = function( alt, side = "right", null.val = 0, sigma0 = 1, N, type1 = 0.005, root = 0){
  
  if(side=="right"){
    z.alpha = qnorm( p=type1, lower.tail = F)
    c0 = null.val + ((z.alpha*sigma0)/sqrt(N))
    t = pnorm( c0, mean = alt, sd = (sigma0/sqrt(N)))
  }else if(side=="left"){
    z.alpha = qnorm( p=type1, lower.tail = F)
    c0 = null.val - ((z.alpha*sigma0)/sqrt(N))
    t = pnorm( c0, mean = alt, sd = (sigma0/sqrt(N)), lower.tail = F)
  }
  
  return(t -root)
}


### computes (type-2 error -root) for t-test

##  arguments
# side      : direction of H1(right, left or both)
# alt       : alt value
# null.val  : null value
# N         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.t = function( alt, side = "right", null.val = 0, N, type1 = 0.005, root = 0){
  
  if(side=="right"){
    t.alpha = qt( type1, df = (N -1), lower.tail = F)
    t = pt( t.alpha, df = (N-1), ncp = (sqrt(N)*(alt-null.val)) )
  }else if(side=="left"){
    t.alpha = qt( type1, df = (N -1), lower.tail = F)
    t = pt( -t.alpha, df = (N-1), ncp = (sqrt(N)*(alt-null.val)), lower.tail = F )
  }
  
  return(t -root)
}


### finds alternative given type-2 error

##  arguments
# test.type   : test type; prop.test, z.test or t.test
# side        : direction of H1(right, left or both)
# null        : null value
# size        : direction of H1
# type1       : specified Type-1 error
# type2       : specified Type-2 error
# sigma0      : known s.d.

## returns {\theta}_{\beta} corresponding to the specified \beta

find.alt = function( test.type, side = "right", null, size, type1 = 0.005, type2 = 0.2, sigma0 = 1){
  
  if((test.type!="prop.test") & (test.type!="z.test") & (test.type!="t.test")){
    return(print("unidentified test type! Has to be one of 'prop.test', 'z.test' or 't.test' "))
  }
  
  if(test.type=="prop.test"){
    
    if(missing(null)==T){
      null = 0.5
    }
    
    if(side=="right"){
      
      solve.out = uniroot( f = type2.error.ber, interval = c( null,1), side = side,
                           null.val=null, N=size, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
      
    }else if(side=="left"){
      
      solve.out = uniroot( f = type2.error.ber, interval = c( 0,null), side = side,
                           null.val=null, N=size, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
    }
    
  }else if(test.type=="z.test"){
    
    if(missing(null)==T){
      null = 0
    }
    
    z.alpha = qnorm( p = type1, lower.tail = F)
    z.beta = qnorm( p = type2, lower.tail = F)
    
    if(side=="right"){
      
      alt = null + ((sigma0*(z.alpha+z.beta))/sqrt(size))
      
    }else if(side=="left"){
      
      alt = null - ((sigma0*(z.alpha+z.beta))/sqrt(size))
    }
    
  }else if(test.type=="t.test"){
    
    if(missing(null)==T){
      null = 0
    }
    
    if(side=="right"){
      
      solve.out = uniroot( f = type2.error.t, interval = c( null, .Machine$integer.max),
                           null.val=null, side=side, N=size, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
      
    }else if(side=="left"){
      
      solve.out = uniroot( f = type2.error.t, interval = c( -.Machine$integer.max, null),
                           null.val=null, side=side, N=size, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
    }
    
  }
  return(alt)
}



### a general checking function for comparing a sequence of statistic with a given sequence of threshold
### after inconclusion at N, either accept at H_0 or checks bayes factor >= threshold to reject H_0

##  arguments
# statistic : a sequence containing values of statistic
# upper     : a sequence of upper thresholds
# lower     : a sequence of lower thresholds
# batch.seq : a sequence of sample sizes where we observe data and make comparison
# threshold : \gamma in ModSPRT

## returns decision accept or reject or continue, and sample size needed for coming to a decision


check = function( test.type, statistic, upper, lower, batch.seq, threshold){
  
  # checking the lengths of upper, lower and batch.seq
  if((length(upper)!=length(batch.seq)) || (length(lower)!=length(batch.seq))){
    
    return("length of 'batch.seq', 'upper' and 'lower' doesn't match!")
  }
  
  if(length(which(statistic>=upper[1:length(statistic)]))==0){
    
    if(length(which(statistic<=lower[1:length(statistic)]))==0){
      
      if(length(statistic)==length(batch.seq)){
        
        n = batch.seq[length(batch.seq)]
        
        if(statistic[length(batch.seq)] >= threshold){
          decision = "reject"
        }else{
          decision = "accept"
        }
        
      }else{
        n = batch.seq[length(statistic)]
        decision = "continue"
      }
      
      
    }else{
      n = batch.seq[min( which(statistic<=lower[1:length(statistic)]) )]
      decision = "accept"
    }
  }else{
    
    if(length(which(statistic<=lower[1:length(statistic)]))==0){
      n = batch.seq[min( which(statistic>=upper[1:length(statistic)]) )]
      decision = "reject"
    }else{
      n = batch.seq[min(union(which(statistic>=upper[1:length(statistic)]),
                              which(statistic<=lower[1:length(statistic)])))]
      
      if(min(which(statistic>=upper[1:length(statistic)])) <  min(which(statistic<=lower[1:length(statistic)]))){
        decision = "reject"
      }else{decision = "accept"}
    }
  }
  
  return(list("decision"=decision, "n"=n))
}




#### computes estimated type1 & type2 error, avg n0 and n1 given boundaries for MSPRT and SUMPBT

### for a single replication

##  arguments
# test.type   : test type; prop.test, z.test or t.test
# side        : direction of H1
# batch.seq   : a sequence of sample sizes where we observe data and make comparison
# type1       : specified Type-1 error
# type2       : specified Type-2 error
# null        : null value
# gen.par     : \theta at which observations are generated from
# other.par   : parameter value other than hypothesized parameter
# alt.LR      : alternative at which LR is computed, needed only for prop.test and z.test
# upper.cut   : sequence of upper threshold
# lower.cut   : sequence of lower threshold
# N           : max affordable sample size
# delta       : \gamma as in ModSPRT

## given a \delta, returns a list containing counts and an inconclusive vector. count is used to estimate
## P( \tau <= N) in SPRT. the vector is the vector of LR values at N which corresponds to \tau > N in SPRT


ovr.repl.ber = function( error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
                         up, low, N, seed){
  
  n = incr.count = 0
  inconclusive = numeric()
  
  # data generation under gen.par
  set.seed(seed)
  x.seq = rbinom(N, 1, gen.par)
  
  indx = 1
  while((indx<=length(batch.seq)) && (n==0)){
    
    cumsum.obs = sum(x.seq[1:batch.seq[indx]])
    LR1 = (((1- alt.LR[1])/(1- null))^batch.seq[indx]) * (((alt.LR[1]*(1- null))/(null*(1- alt.LR[1])))^cumsum.obs)
    LR2 = (((1- alt.LR[2])/(1- null))^batch.seq[indx]) * (((alt.LR[2]*(1- null))/(null*(1- alt.LR[2])))^cumsum.obs)
    LR = (alt.psi*LR1) + ((1- alt.psi)*LR2)
    
    if((error.type=="type1") && (LR>=up)){
      n = batch.seq[indx]
      incr.count = 1
    }else if((error.type=="type2") && (LR>=up)){
      n = batch.seq[indx]
    }else if((error.type=="type1") && (LR<=low)){
      n = batch.seq[indx]
    }else if((error.type=="type2") && (LR<=low)){
      n = batch.seq[indx]
      incr.count = 1
    }
    indx = indx +1
  }
  
  if((indx>length(batch.seq)) && (n==0)){
    inconclusive = LR
    n = batch.seq[length(batch.seq)]
  }
  
  return(list(incr.count,inconclusive, n))
}


ovr.repl.norm = function( error.type, batch.seq, null, gen.par, alt.LR,
                          up, low, N, seed){
  
  n = incr.count = 0
  inconclusive = numeric()
  
  # data generation under gen.par
  set.seed(seed)
  x.seq = rnorm(N, mean = gen.par[1], sd = gen.par[2])
  
  indx = 1
  while((indx<=length(batch.seq)) && (n==0)){
    
    cumsum.obs = sum(x.seq[1:batch.seq[indx]])
    LR = (cumsum.obs*(alt.LR - null)) - ((batch.seq[indx]*((alt.LR^2) - (null^2)))/2)
    LR = exp(LR/(gen.par[2]^2))
    
    if((error.type=="type1") && (LR>=up)){
      n = batch.seq[indx]
      incr.count = 1
    }else if((error.type=="type2") && (LR>=up)){
      n = batch.seq[indx]
    }else if((error.type=="type1") && (LR<=low)){
      n = batch.seq[indx]
    }else if((error.type=="type2") && (LR<=low)){
      n = batch.seq[indx]
      incr.count = 1
    }
    indx = indx +1
  }
  
  if((indx>length(batch.seq)) && (n==0)){
    inconclusive = LR
    n = batch.seq[length(batch.seq)]
  }
  
  return(list(incr.count,inconclusive, n))
}


ovr.repl.t = function( side, error.type, batch.seq, type1, null, gen.par,
                       up, low, N, seed){
  
  n = incr.count = 0
  inconclusive = numeric()
  
  # data generation under gen.par
  set.seed(seed)
  x.seq = rnorm(N, mean = gen.par, sd = 1)
  
  indx = 1
  while((indx<=length(batch.seq)) && (n==0)){
    
    cumsum.obs = sum(x.seq[1:batch.seq[indx]])
    sx = sd(x.seq[1:batch.seq[indx]])
    alt.LR = find.umpbt.t(side = side, type1 = type1, n.obs = N, null = null, s = sx)
    t0 = (((cumsum.obs/batch.seq[indx]) -null)/sx)
    t1 = (((cumsum.obs/batch.seq[indx]) -alt.LR)/sx)
    LR = ((1+ ((batch.seq[indx]/(batch.seq[indx]-1))*(t0 ^2)) )/(1+ ((batch.seq[indx]/(batch.seq[indx]-1))*(t1 ^2)) ))^(batch.seq[indx]/2)
    
    if((error.type=="type1") && (LR>=up)){
      n = batch.seq[indx]
      incr.count = 1
    }else if((error.type=="type2") && (LR>=up)){
      n = batch.seq[indx]
    }else if((error.type=="type1") && (LR<=low)){
      n = batch.seq[indx]
    }else if((error.type=="type2") && (LR<=low)){
      n = batch.seq[indx]
      incr.count = 1
    }
    indx = indx +1
  }
  
  if((indx>length(batch.seq)) && (n==0)){
    inconclusive = LR
    n = batch.seq[length(batch.seq)]
  }
  
  return(list(incr.count,inconclusive, n))
}


overshoot.ber = function( error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
                          up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.ber, seed = ((k-1)*1e+5 + seq(1e+5)), 
                      MoreArgs = list(error.type=error.type, batch.seq=batch.seq,
                                      null=null, gen.par=gen.par, alt.LR=alt.LR, alt.psi = alt.psi,
                                      up=up, low=low, N=N))
    
    count.t = count.t + sum(as.numeric(out.temp[1,]))
    inc = as.numeric(out.temp[2,])
    vec.t = c(vec.t,inc[!is.na(inc)])
    n = as.numeric(out.temp[3,])
    n.vec.t = c(n.vec.t,n)
    
    t = list( count.t, vec.t, n.vec.t)
  }
  
  inc = as.numeric(out[[2]])
  vec = c(vec,inc[!is.na(inc)])
  
  if(return.n==T){
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]])))
  }else{
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = as.numeric(out[[3]])))
  }
}



overshoot.norm = function( error.type, batch.seq, null, gen.par, alt.LR,
                           up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.norm, seed = ((k-1)*1e+5 + seq(1e+5)), 
                      MoreArgs = list(error.type=error.type, batch.seq=batch.seq,
                                      null=null, gen.par=gen.par, alt.LR=alt.LR,
                                      up=up, low=low, N=N))
    
    count.t = count.t + sum(as.numeric(out.temp[1,]))
    inc = as.numeric(out.temp[2,])
    vec.t = c(vec.t,inc[!is.na(inc)])
    n = as.numeric(out.temp[3,])
    n.vec.t = c(n.vec.t,n)
    
    t = list( count.t, vec.t, n.vec.t)
  }
  
  inc = as.numeric(out[[2]])
  vec = c(vec,inc[!is.na(inc)])
  
  if(return.n==T){
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]])))
  }else{
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = as.numeric(out[[3]])))
  }
}



overshoot.t = function( side, error.type, batch.seq, type1, null, gen.par,
                        up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.t, seed = ((k-1)*1e+5 + seq(1e+5)), 
                      MoreArgs = list(side = side, error.type=error.type, batch.seq=batch.seq,
                                      type1 = type1, null=null, gen.par=gen.par, 
                                      up=up, low=low, N=N))
    
    count.t = count.t + sum(as.numeric(out.temp[1,]))
    inc = as.numeric(out.temp[2,])
    vec.t = c(vec.t,inc[!is.na(inc)])
    n = as.numeric(out.temp[3,])
    n.vec.t = c(n.vec.t,n)
    
    t = list( count.t, vec.t, n.vec.t)
  }
  
  inc = as.numeric(out[[2]])
  vec = c(vec,inc[!is.na(inc)])
  
  if(return.n==T){
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]])))
  }else{
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = as.numeric(out[[3]])))
  }
}




### combined for all replications

##  arguments

# type1       : specified Type-1 error
# delta       : \gamma as in MSPRT or SUMPBT
# R           : no. of replications

## given a \delta, based on R replications this returns (estimated type1 error - specified type1)


error.summary = function( error.type, delta, root, count, inconclusive.vec, R, type1, type2){
  
  if(error.type=="type1"){
    est = (count + sum(inconclusive.vec>=delta))/R
    
    k.dec = abs(ceiling(log10(sqrt((type1*(1- type1))/R))))
    est = floor(est*10^k.dec)/(10^k.dec)
    
  }else if(error.type=="type2"){
    est = (count + sum(inconclusive.vec<delta))/R
    
    k.dec = abs(ceiling(log10(sqrt((type2*(1- type2))/R))))
    est = floor(est*10^k.dec)/(10^k.dec)
  }
  
  return((est - root))
}



### finding termination threshold in MSPRT, and it's performance

##  arguments
# test.type   : test type; prop.test, z.test or t.test
# side        : direction of H1
# batch.seq   : a sequence of sample sizes where we observe data and make comparison
# type1       : specified Type-1 error
# type2       : specified Type-2 error
# null        : null value
# sigma0      : known s.d. for z-test
# N.max       : max affordable sample size
# repl        : no. of replications
# seed        : seed for randomization

## based on R replications, this returns optimal delta(\gamma) and a bunch of other estimated quantities
## reflecting it's performance

design.MSPRT = function( test.type, side, batch.seq, type1 =.005, type2 = .2,
                         null, sigma0 = 1, N.max, alt.comp,
                         repl, verbose=T, core.no){
  
  # checking if N.max is provided
  if(missing(N.max)==T){
    return("Maximum budget on sample size is not provided")
  }
  
  
  ## msg
  if(verbose==T){
    print("-------------------------------------------------------------------------")
    if(test.type=="prop.test"){
      print("Designing the MSPRT in case of a test for a binomial proportion")
    }else if(test.type=="z.test"){
      print("Designing the MSPRT in case of a Z-test")
    }else if(test.type=="t.test"){
      print("Designing the MSPRT in case of a T-test")
    }
    print("-------------------------------------------------------------------------")
  }
  
  
  ## setting default arguments
  
  # direction of H1
  if(missing(side)==T){
    side = "right"
    
    ## msg
    if(verbose==T){
      print("Working with side = right (default)")
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with side = ",side," (user specified)", sep = ""))
    }
  }
  
  
  # default batch sequence
  if(missing(batch.seq)==T){
    if(test.type!="t.test"){
      batch.seq = 1:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",1,":N.max (default)", sep = ""))
      }
    }else{
      batch.seq = 2:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",2,":N.max (default)", sep = ""))
      }
    }
    
  }else{
    
    if((batch.seq[1]<2) && (test.type=="t.test")){
      return(print("batch.seq[1]<2 . Need at least 2 samples to compute sample standard deviation"))
    }else{
      
      ##msg
      if(verbose==T){
        print("Working with batch.seq (user specified)")
      }
    }
  }
  
  
  # type1
  if(missing(type1)==T){
    
    ## msg
    if(verbose==T){
      print(paste("Working with type1 = ",type1," (default)", sep = ""))
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with user specified type1 = ",type1," (user specified)", sep = ""))
    }
  }
  
  
  # type2
  if(missing(type2)==T){
    
    ## msg
    if(verbose==T){
      print(paste("Working with type2 = ",type2," (default)", sep = ""))
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with type2 = ",type2," (user specified)", sep = ""))
    }
  }
  
  
  ##
  if(missing(repl)==T){
    if(test.type=="prop.test"){
      repl = 2e+6
    }else{repl = 1e+6}
  }else{repl = repl}
  k.dec = floor(log(sqrt(repl)))
  
  # default number of cores
  if(missing(core.no)==T){
    c = detectCores()
    core.no = max(1,(c-1))
  }else{core.no=core.no}
  
  
  ## proptest
  if(test.type=="prop.test"){
    
    # setting default null
    if(missing(null)==T){
      null = 0.5
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    umpbt.out = find.umpbt.ber( side = side, type1 = type1, n.obs = N.max, null = null)
    alt.psi = umpbt.out$psi
    alt.LR = umpbt.out$u
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      if(length(alt.LR)==1){
        print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      }else if(length(alt.LR)>1){
        print(paste("The UMPBT alternative has been obtained: ", round(alt.LR[1],3),
                    " & ",round(alt.LR[2],3)," with probabilities ", round(alt.psi,3),
                    " & ",(1- round(alt.psi,3)),", respectively", sep = "" ))
      }
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    # calculating overshooting summary of \tau <= N.max in case of SPRT
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Searching for Termination Threshold")
    }
    
    # calculating overshooting probability in SPRT under null
    registerDoParallel(cores = core.no)
    
    count = 0
    vec = n.vec = numeric()
    k.rep = repl/(1e+5)
    
    out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
      
      count.t = 0
      vec.t = n.vec.t = numeric()
      out.temp = mapply(FUN = ovr.repl.ber, seed = ((k-1)*1e+5 + seq(1e+5)), 
                        MoreArgs = list(error.type="type1", batch.seq=batch.seq,
                                        null=null, gen.par=null, alt.LR=alt.LR, alt.psi = alt.psi,
                                        up=wald.up, low=wald.low, N=N.max))
      
      count.t = count.t + sum(as.numeric(out.temp[1,]))
      inc = as.numeric(out.temp[2,])
      vec.t = c(vec.t,inc[!is.na(inc)])
      n = as.numeric(out.temp[3,])
      n.vec.t = c(n.vec.t,n)
      
      t = list( count.t, vec.t, n.vec.t)
    }
    
    inc = as.numeric(out[[2]])
    vec = inc[!is.na(inc)]
    
    overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
    
    
    # increasingly ordered inconclusive vector
    vec = sort(overshoot.summ$inconclusive.vec)
    
    vec.unique = unique(vec)
    
    # setting default starting delta
    delta.opt = ump.match.ber( side, type1, N.max, null)
    
    if(sum(vec.unique>delta.opt)==0){
      
      i.opt = length(vec.unique)
    }else if(sum(vec.unique<delta.opt)==0){
      
      i.opt = 1
    }else{
      i.opt = min(which(vec.unique>=delta.opt))
    }
    
    fnew = error.summary( error.type = "type1", delta= vec.unique[i.opt], root = type1, count= overshoot.summ$count,
                          inconclusive.vec= vec, R= repl, type1 = type1)
    
    if(fnew<=0){
      
      if(i.opt==1){
        
        delta.opt.t = (floor(wald.low*1000))/1000
        delta.opt = (floor(delta.opt.t*100))/100
        while((delta.opt - delta.opt.t)<=(10^-k.dec)){
          delta.opt = delta.opt +.01
        }
        
      }else{
        
        i.opt = i.opt - 1
        while( (fnew<=0) && (i.opt>0)){
          
          fnew = error.summary( error.type = "type1", delta= vec.unique[i.opt], root = type1, count= overshoot.summ$count,
                                inconclusive.vec= vec, R= repl, type1 = type1)
          i.opt = i.opt - 1
        }
        
        if( (fnew>0) && (i.opt>0) ){
          
          i.opt = i.opt + 1
          
          i.opt.u = i.opt
          delta.opt.t.u = (floor(vec.unique[i.opt.u]*1000))/1000
          delta.opt.u = (floor(delta.opt.t.u*100))/100
          while((delta.opt.u - delta.opt.t.u)<=(10^-k.dec)){
            delta.opt.u = delta.opt.u +.01
          }
          
          i.opt.l = i.opt -1
          delta.opt.t.l = (floor(vec.unique[i.opt.l]*1000))/1000
          delta.opt.l = (floor(delta.opt.t.l*100))/100
          while((delta.opt.l - delta.opt.t.l)<=(10^-k.dec)){
            delta.opt.l = delta.opt.l +.01
          }
          
          quant.prob = error.summary( error.type = "type1", delta= vec.unique[i.opt+1], root = 0, 
                                      count= overshoot.summ$count, inconclusive.vec= vec, R= repl, type1 = type1)
          psi = ((type1 - quant.prob)*repl)/sum(vec==vec.unique[i.opt])
          delta.opt = c(delta.opt.l, delta.opt.u)
          
        }else if( (fnew<=0) && (i.opt==0) ){
          
          delta.opt.t = (floor(wald.low*1000))/1000
          delta.opt = (floor(delta.opt.t*100))/100
          while((delta.opt - delta.opt.t)<=(10^-k.dec)){
            delta.opt = delta.opt +.01
          }
          
        }else if( (fnew>0) && (i.opt==0) ){
          
          i.opt = i.opt + 1
          
          i.opt.u = i.opt
          delta.opt.t.u = (floor(vec.unique[i.opt.u]*1000))/1000
          delta.opt.u = (floor(delta.opt.t.u*100))/100
          while((delta.opt.u - delta.opt.t.u)<=(10^-k.dec)){
            delta.opt.u = delta.opt.u +.01
          }
          
          delta.opt.t.l = (floor(wald.low*1000))/1000
          delta.opt.l = (floor(delta.opt.t.l*100))/100
          while((delta.opt.l - delta.opt.t.l)<=(10^-k.dec)){
            delta.opt.l = delta.opt.l +.01
          }
          
          quant.prob = error.summary( error.type = "type1", delta= vec.unique[i.opt+1], root = 0, 
                                      count= overshoot.summ$count, inconclusive.vec= vec, R= repl, type1 = type1)
          psi = ((type1 - quant.prob)*repl)/sum(vec==vec.unique[i.opt])
          delta.opt = c(delta.opt.l, delta.opt.u)
          
        }
      }
    }else if(fnew>0){
      
      if(i.opt==length(vec.unique)){
        
        delta.opt.t = (floor(vec.unique[length(vec.unique)]*1000))/1000
        delta.opt = (floor(delta.opt.t*100))/100
        while((delta.opt - delta.opt.t)<=(10^-k.dec)){
          delta.opt = delta.opt +.01
        }
        
      }else{
        
        i.opt = i.opt + 1
        while( (fnew>0) && (i.opt<(length(vec.unique) +1)) ){
          
          fnew = error.summary( error.type = "type1", delta= vec.unique[i.opt], root = type1, count= overshoot.summ$count,
                                inconclusive.vec= vec, R= repl, type1 = type1)
          i.opt = i.opt + 1
        }
        
        if( (fnew<=0) && (i.opt<(length(vec.unique) +1)) ){
          
          i.opt = i.opt -2
          
          i.opt.u = i.opt
          delta.opt.t.u = (floor(vec.unique[i.opt.u]*1000))/1000
          delta.opt.u = (floor(delta.opt.t.u*100))/100
          while((delta.opt.u - delta.opt.t.u)<=(10^-k.dec)){
            delta.opt.u = delta.opt.u +.01
          }
          
          i.opt.l = i.opt -1
          delta.opt.t.l = (floor(vec.unique[i.opt.l]*1000))/1000
          delta.opt.l = (floor(delta.opt.t.l*100))/100
          while((delta.opt.l - delta.opt.t.l)<=(10^-k.dec)){
            delta.opt.l = delta.opt.l +.01
          }
          
          quant.prob = error.summary( error.type = "type1", delta= vec.unique[i.opt+1], root = 0, 
                                      count= overshoot.summ$count, inconclusive.vec= vec, R= repl, type1 = type1)
          psi = ((type1 - quant.prob)*repl)/sum(vec==vec.unique[i.opt])
          delta.opt = c(delta.opt.l, delta.opt.u)
          
          
        }else if( (fnew>0) && (i.opt==(length(vec.unique) +1)) ){
          
          delta.opt.t = (floor(vec.unique[length(vec.unique)]*1000))/1000
          delta.opt = (floor(delta.opt.t*100))/100
          while((delta.opt - delta.opt.t)<=(10^-k.dec)){
            delta.opt = delta.opt +.01
          }
          
        }else if( (fnew<=0) && ( i.opt==(length(vec.unique) +1) ) ){
          
          i.opt = i.opt -2
          
          i.opt.u = i.opt
          delta.opt.t.u = (floor(vec.unique[i.opt.u]*1000))/1000
          delta.opt.u = (floor(delta.opt.t.u*100))/100
          while((delta.opt.u - delta.opt.t.u)<=(10^-k.dec)){
            delta.opt.u = delta.opt.u +.01
          }
          
          i.opt.l = i.opt -1
          delta.opt.t.l = (floor(vec.unique[i.opt.l]*1000))/1000
          delta.opt.l = (floor(delta.opt.t.l*100))/100
          while((delta.opt.l - delta.opt.t.l)<=(10^-k.dec)){
            delta.opt.l = delta.opt.l +.01
          }
          
          quant.prob = error.summary( error.type = "type1", delta= vec.unique[i.opt+1], root = 0, 
                                      count= overshoot.summ$count, inconclusive.vec= vec, R= repl)
          psi = ((type1 - quant.prob)*repl)/sum(vec==vec.unique[i.opt], type1 = type1)
          delta.opt = c(delta.opt.l, delta.opt.u)
          
        }
      }
      
    }
    
    
    ## msg
    if(verbose==T){
      
      print(paste("Termination threshold is ", delta.opt.u, sep = ""))
      print("Done")
    }
    
    
    # performance under null
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("OC of the obtained MSPRT under the null is being calculated:")
    }
    
    type1.est = error.summary( error.type = "type1", delta= delta.opt.u, root = 0, count= overshoot.summ$count,
                               inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type1 = type1)
    
    var.n0 = var(overshoot.summ$n.vec)
    k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
    avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
    
    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)
    
    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ", type1.est, sep = ""))
      print("Done")
    }
    
    
    # finding alternative corresponding to fixed type1, type2, N
    
    if(missing(alt.comp)==T){
      
      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }
      
      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0, 
                   "umpbt.alt"=alt.LR, "psi.umpbt"=alt.psi,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"= delta.opt.u ))
      
      
    }else{
      
      if(alt.comp==T){
        
        alt.comp = find.alt(test.type=test.type, side = side, null = null, size = N.max, type1 = type1, type2 = type2)
        
        ## msg
        if(verbose==T){
          print("-------------------------------------------------------------------------")
          print(paste("Working with the fixed design alternative ",round(alt.comp,3), sep = ""))
        }
        
      }else if(class(alt.comp)=="numeric"){
        
        alt.comp = alt.comp
        
        ## msg
        if(verbose==T){
          print("-------------------------------------------------------------------------")
          print(paste("Working with the user specified alternative ",round(alt.comp,3), sep = ""))
        }
      }
      
      
      # actual type 2 error probability at umpbt and fixed design
      alt.type2 = type2.error.ber( side = side, alt = alt.comp, null.val = null, N=N.max, type1 = type1, root = 0)
      
      
      # performance under the corresponding alternative
      
      ## msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("OC of the obtained MSPRT under the user desired alternative is being calculated:")
      }
      
      
      # calculating overshooting probability in SPRT uneder alternative
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.ber, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type2", batch.seq=batch.seq,
                                          null=null, gen.par=alt.comp, alt.LR=alt.LR, alt.psi = alt.psi,
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = inc[!is.na(inc)]
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= delta.opt.u, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ",type2.est,sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      
      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1,
                   "umpbt.alt"=alt.LR, "psi.umpbt"= alt.psi,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"= delta.opt.u ))
      
    }
  }
  
  
  ## z test
  if(test.type=="z.test"){
    
    # setting default null
    if(missing(null)==T){
      null = 0
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      null = null
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    if(missing(sigma0)==T){
      sigma0 = 1
      
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{
      
      sigma0 = sigma0
      
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    alt.LR = find.umpbt.norm( side = side, type1 = type1, n.obs = N.max, null = null,
                              sigma0 = sigma0)
    
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      print(paste("UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    # calculating overshooting summary of \tau <= N.max in case of SPRT
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Searching for Termination Threshold")
    }
    
    # calculating overshooting probability in SPRT under null
    registerDoParallel(cores = core.no)
    
    count = 0
    vec = n.vec = numeric()
    k.rep = repl/(1e+5)
    
    out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
      
      count.t = 0
      vec.t = n.vec.t = numeric()
      out.temp = mapply(FUN = ovr.repl.norm, seed = ((k-1)*1e+5 + seq(1e+5)), 
                        MoreArgs = list(error.type="type1", batch.seq=batch.seq,
                                        null=null, gen.par=c(null,sigma0), alt.LR=alt.LR, 
                                        up=wald.up, low=wald.low, N=N.max))
      
      count.t = count.t + sum(as.numeric(out.temp[1,]))
      inc = as.numeric(out.temp[2,])
      vec.t = c(vec.t,inc[!is.na(inc)])
      n = as.numeric(out.temp[3,])
      n.vec.t = c(n.vec.t,n)
      
      t = list( count.t, vec.t, n.vec.t)
    }
    
    inc = as.numeric(out[[2]])
    vec = c(vec,inc[!is.na(inc)])
    
    overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
    
    
    # increasingly ordered inconclusive vector
    vec = sort(overshoot.summ$inconclusive.vec)
    
    
    root.out = uniroot(f=error.summary, interval = c(wald.low,wald.up), error.type = "type1",
                       root = type1, count=overshoot.summ$count, inconclusive.vec=vec, R=repl,
                       type1 = type1, extendInt = "downX")
    delta.opt = floor((root.out$root*10^3))/(10^3)
    
    
    ## msg
    if(verbose==T){
      print(paste("Termination threshold is ", delta.opt, sep = ""))
      print("Done")
    }
    
    
    # performance under null
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("OC of the obtained MSPRT under the null is being calculated:")
    }
    
    k.dec.type1 = abs(ceiling(log10(sqrt((type1*(1- type1))/repl))))
    type1.est = round( ((overshoot.summ$count + sum(overshoot.summ$inconclusive.vec>=delta.opt))/repl), (k.dec.type1 -1))
    
    var.n0 = var(overshoot.summ$n.vec)
    k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
    avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
    
    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)
    
    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ",type1.est,sep = ""))
      print("Done")
    }
    
    
    # finding alternative corresponding to fixed type1, type2, N
    
    if(missing(alt.comp)==T){
      
      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }
      
      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0, 
                   "umpbt.alt"=alt.LR,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))
      
    }else{
      
      if(alt.comp==T){
        
        alt.comp = find.alt(test.type=test.type, side = side, null = null, size = N.max, type1 = type1, type2 = type2,
                            sigma0 = sigma0)
        
        ## msg
        if(verbose==T){
          print("-------------------------------------------------------------------------")
          print(paste("Working with fixed design alternative ",round(alt.comp,3), sep = ""))
        }
        
      }else if(class(alt.comp)=="numeric"){
        
        alt.comp = alt.comp
        
        ## msg
        if(verbose==T){
          print("-------------------------------------------------------------------------")
          print(paste("Working with user specified alternative ",round(alt.comp,3), sep = ""))
        }
      }
      
      
      # actual type 2 error probability at umpbt and fixed design
      alt.type2 = type2.error.norm( side = side, alt = alt.comp, null.val = null, sigma0 = sigma0, N=N.max,
                                    type1 = type1, root = 0)
      
      
      # performance under the corresponding alternative
      
      ## msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("OC of the obtained MSPRT under the user desired alternative is being calculated:")
      }
      
      
      # calculating overshooting probability in SPRT under alternative
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.norm, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type2", batch.seq=batch.seq,
                                          null=null, gen.par=c(alt.comp,sigma0), alt.LR=alt.LR, 
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1, "umpbt.alt"=alt.LR,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))
      
    }
  }
  
  
  ## t test
  if(test.type=="t.test"){
    
    # setting default null
    if(missing(null)==T){
      null = 0
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      null = null
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    # calculating overshooting summary of \tau <= N.max in case of SPRT
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Searching for Termination Threshold")
    }
    
    
    # calculating overshooting probability in SPRT under null
    registerDoParallel(cores = core.no)
    
    count = 0
    vec = n.vec = numeric()
    k.rep = repl/(1e+5)
    
    out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
      
      count.t = 0
      vec.t = n.vec.t = numeric()
      out.temp = mapply(FUN = ovr.repl.t, seed = ((k-1)*1e+5 + seq(1e+5)), 
                        MoreArgs = list(side = side, error.type="type1", batch.seq=batch.seq,
                                        type1 = type1, null=0, gen.par=0, 
                                        up=wald.up, low=wald.low, N=N.max))
      
      count.t = count.t + sum(as.numeric(out.temp[1,]))
      inc = as.numeric(out.temp[2,])
      vec.t = c(vec.t,inc[!is.na(inc)])
      n = as.numeric(out.temp[3,])
      n.vec.t = c(n.vec.t,n)
      
      t = list( count.t, vec.t, n.vec.t)
    }
    
    inc = as.numeric(out[[2]])
    vec = c(vec,inc[!is.na(inc)])
    
    overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
    
    ## msg
    if(verbose==T){
      print("UMPBT alternatives obtained")
    }
    
    
    root.out = uniroot(f=error.summary, interval = c(wald.low,wald.up), error.type = "type1",
                       root = type1, count=overshoot.summ$count, inconclusive.vec=overshoot.summ$inconclusive.vec,
                       R=repl, type1 = type1, extendInt = "downX")
    delta.opt = floor(root.out$root*10^3)/(10^3)
    
    
    ## msg
    if(verbose==T){
      print(paste("Termination threshold is ", delta.opt, sep = ""))
      print("Done")
    }
    
    
    # performance under null
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("OC of the obtained MSPRT under the null is being calculated:")
    }
    
    type1.est = error.summary( error.type = "type1", delta= delta.opt, root = 0, count= overshoot.summ$count,
                               inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type1 = type1)
    
    var.n0 = var(overshoot.summ$n.vec)
    k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
    avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
    
    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)
    
    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ", type1.est, sep = ""))
      print("Done")
    }
    
    
    # finding alternative corresponding to fixed type1, type2, N
    
    if(missing(alt.comp)==T){
      
      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }
      
      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))
      
      
    }else{
      
      if(class(alt.comp)=="numeric"){
        alt.comp = alt.comp
        
        ## msg
        if(verbose==T){
          print("-------------------------------------------------------------------------")
          print(paste("Working with user specified point alternative ",round(alt.comp,3), sep = ""))
          
        }else if(alt.comp==T){
          
          alt.comp = find.alt(test.type=test.type, side = side, null = null, size = N.max, type1 = type1, type2 = type2)
          
          ## msg
          if(verbose==T){
            print("-------------------------------------------------------------------------")
            print(paste("Working with fixed design alternative ",round(alt.comp,3), sep = ""))
          }
          
        }
      }
      
      
      # actual type 2 error probability at umpbt and fixed design
      alt.type2 = type2.error.t( side = side, alt = alt.comp, null.val = null, N=N.max, type1 = type1, root = 0)
      
      
      # performance under the corresponding alternative
      alt.comp.shifted = alt.comp - null
      null.shifted = 0
      
      ## msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("OC of the obtained MSPRT under the user desired alternative is being calculated:")
      }
      
      
      # calculating overshooting probability in SPRT under alternative
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.t, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(side = side, error.type="type2", batch.seq=batch.seq,
                                          type1 = type1, null=null.shifted, gen.par=alt.comp.shifted,
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
      
      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))
      
    }
  }
}



### finding effective N in proportion test

effective.N = function(N.max, side = "right", type1 = 0.005, null.val = 0.5, plot.it = T){
  
  N.seq = seq(1, N.max, by = 1)
  
  # finding umpbt seq
  umpbt.seq = mapply(FUN = point.umpbt.ber, n.obs = N.seq,
                     MoreArgs = list(side = side, type1 = type1, null = null.val))
  
  if(side=="right"){
    
    u = umpbt.seq[1]
    N.eff = N.seq[1]
    i = 2
    while (i<(length(N.seq)+1)) {
      if(umpbt.seq[i]<u[length(u)]){
        u = c(u,umpbt.seq[i])
        N.eff = c(N.eff, N.seq[i])
      }else{
        u = c(u,u[length(u)])
        N.eff = c(N.eff, N.eff[length(N.eff)])
      }
      i = i+1
    }
    
    if(plot.it==T){
      
      layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
             widths=c(1,1), heights=c(3,1))
      plot(umpbt.seq~N.seq, type="b", pch=16, cex=.7, cex.main=1.5,
           xlab = "sample size", ylab = paste("The point UMPBT alternative at",type1),
           main="Finding the 'max sample size' for the MSPRT in a proportion test")
      points(u~N.seq, type="b", pch=16, cex=.7, col=2)
      points(u[which(N.seq==N.eff)]~N.seq[which(N.seq==N.eff)], cex=1.5, col="forestgreen", lwd=1.5)
      abline(v=N.eff[length(N.seq)], lty=2)
      text(N.eff[length(N.seq)], mean(range(umpbt.seq)), labels = paste("N =",N.eff[length(N.seq)]), pos = 2)
      plot.new()
      legend( "center", xpd = TRUE, xjust = .5, yjust = 1, lwd = 2, bty = "n", cex= .8,
              legend=c( "Original point UMPBT alternatives", "The decreasing point UMPBT alternatives",
                        "Possible max sample size values", "The chosen max sample size"),
              lty = c(NA,NA,NA,2), col = c(1,2,"forestgreen",1), pch = c(16,16,1,NA), merge = TRUE, ncol = 2)
      layout(1)
    }
    
  }else if(side=="left"){
    
    u = umpbt.seq[1]
    N.eff = N.seq[1]
    i = 2
    while (i<(length(N.seq)+1)) {
      if(umpbt.seq[i]>u[length(u)]){
        u = c(u,umpbt.seq[i])
        N.eff = c(N.eff, N.seq[i])
      }else{
        u = c(u,u[length(u)])
        N.eff = c(N.eff, N.eff[length(N.eff)])
      }
      i = i+1
    }
    
    
    if(plot.it==T){
      
      layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
             widths=c(1,1), heights=c(3,1))
      plot(umpbt.seq~N.seq, type="b", pch=16, cex=.7, cex.main=1.5,
           xlab = "sample size", ylab = paste("The point UMPBT alternative at",type1),
           main="Finding the 'max sample size' for the MSPRT in a proportion test")
      points(u~N.seq, type="b", pch=16, cex=.7, col=2)
      points(u[which(N.seq==N.eff)]~N.seq[which(N.seq==N.eff)], cex=1.5, col="forestgreen", lwd=1.5)
      abline(v=N.eff[length(N.seq)], lty=2)
      text(N.eff[length(N.seq)], mean(range(umpbt.seq)), labels = paste("N =",N.eff[length(N.seq)]), pos = 2)
      plot.new()
      legend( "center", xpd = TRUE, xjust = .5, yjust = 1, lwd = 2, bty = "n", cex= 1,
              legend=c( "Original point UMPBT alternatives", "The decreasing point UMPBT alternatives",
                        "Possible max sample size values", "The chosen max sample size"),
              lty = c(NA,NA,NA,2), col = c(1,2,"forestgreen",1), pch = c(16,16,1,NA), merge = TRUE, ncol = 2)
      layout(1)
    }
  }
  
  return(N.eff[length(N.seq)])
}



### finding sample size for higher significance

find.samplesize = function( test.type, N, lower.signif = 0.05, higher.signif = 0.005, null.val, side = "right",
                            pow = 0.8, alt, sigma0 = 1, n.seq, plot.it=T){
  
  # defaults
  if(missing(n.seq)==T){
    n.seq = seq(N, (4*N), by=1)
  }
  
  type2 = 1- pow
  if(test.type=="prop.test"){
    
    if(missing(null.val)==T){
      null.val = .5
    }
    
    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null.val, size = N,
                     type1 = lower.signif, type2 = type2)
    }
    
    t.seq = mapply(FUN = type2.error.ber, N = n.seq,
                   MoreArgs = list(alt = alt, side = side, null.val = null.val, type1 = higher.signif))
    N.star = max(n.seq[which(t.seq>type2)]) +1
    
  }else if(test.type=="z.test"){
    
    if(missing(null.val)==T){
      null.val = 0
    }
    
    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null.val, size = N,
                     sigma0 = sigma0, type1 = lower.signif, type2 = type2)
    }
    
    t.seq = mapply(FUN = type2.error.norm, N = n.seq,
                   MoreArgs = list(side = side, alt = alt, null.val = null.val, sigma0 = sigma0,
                                   type1 = higher.signif))
    N.star = max(n.seq[which(t.seq>type2)]) +1
    
  }else if(test.type=="t.test"){
    
    if(missing(null.val)==T){
      null.val = 0
    }
    
    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null.val, size = N,
                     type1 = lower.signif, type2 = type2)
    }
    
    t.seq = mapply(FUN = type2.error.t, N = n.seq,
                   MoreArgs = list(side = side, alt = alt, null.val = null.val, type1 = higher.signif))
    N.star = max(n.seq[which(t.seq>type2)]) +1
  }
  
  
  if(plot.it==T){
    
    layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), 
           widths=c(1,1), heights=c(3,1))
    plot((1- t.seq)~n.seq, type="l", cex.main = 1.5,
         main = paste("Sample size required for the higher significance at ", round(alt,4)),
         xlab = "Sample size", ylab = paste("Power at ", round(alt,4)) )
    points((1- t.seq)~n.seq, pch = 16, cex=.8)
    abline(h=pow)
    points( (1- t.seq[max(which(t.seq>type2)) +1])~N.star, col=2, pch=16)
    text( floor(.7*n.seq[length(n.seq)]), labels = paste((pow*100), "% power", sep=""), pos = 3)
    abline(v=N.star, lty = 2)
    text(N.star, mean(range((1- t.seq))), labels = paste("N =",N.star), pos = 4)
    plot.new()
    legend("center", xpd = TRUE, pch = c(16,16,NA,NA), col=c(1,2,1,1), lwd=2, lty = c(NA,NA,2,1), seg.len = 1,
           legend=c( "Power for increasing sample size", "Power for the chosen sample size",
                     "The chosen sample size", paste("Desired power at", round(alt, 4)) ),
           xjust = .5, yjust = 1, bty = "n", merge = T, x.intersp = .4, y.intersp = 1, ncol=2, cex= .8)
  }
  return(N.star)
}




### carrys out squential comparison given a set of observations and two boundaries
### also returns a plot

##  arguments
# obs         : observed sequential data
# test.type   : test type; prop.test, z.test or t.test
# side        : direction of H1
# batch.seq   : a sequence of sample sizes where we observe data and make comparison
# type1       : specified Type-1 error
# type2       : specified Type-2 error
# null        : null value
# sigma0      : known s.d. for z-test
# delta.opt   : termination threshold
# N.max       : max affordable sample size
# plot.it     : logical; asks whether a plot should be returned
# standardized.compare  : logical; asks whether comparison should be made based on standardized data

implement.MSPRT = function( obs, test.type, side, batch.seq, type1= .005, type2= .2,
                            null, sigma0, term.thresh, N.max,
                            plot.it=T, verbose=T){
  
  # checking if N.max is provided
  if(missing(N.max)==T){
    return("Maximum available sample size is not provided")
  }
  
  if(test.type=="prop.test"){
    
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Implementing the MSPRT in case of a test for a binomial proportion")
      print("-------------------------------------------------------------------------")
    }
    
    # direction of H1
    if(missing(side)==T){
      side = "right"
    }
    
    #default batch sequence
    if(missing(batch.seq)==T){
      batch.seq = 1:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",1,":N.max (default)", sep = ""))
      }
    }else{
      
      ##msg
      if(verbose==T){
        print("Working with batch.seq (user specified)")
      }
    }
    
    # type1
    if(missing(type1)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type1 = ",type1," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with user specified type1 = ",type1," (user specified)", sep = ""))
      }
    }
    
    
    # type2
    if(missing(type2)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (user specified)", sep = ""))
      }
    }
    
    
    # setting default null
    if(missing(null)==T){
      null = 0.5
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    umpbt.out = find.umpbt.ber( side = side, type1 = type1, n.obs = N.max, null = null)
    alt.psi = umpbt.out$psi
    alt.LR = umpbt.out$u
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      if(length(alt.LR)==1){
        print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      }else if(length(alt.LR)>1){
        print(paste("The UMPBT alternative has been obtained: ", round(alt.LR[1],3),
                    " & ",round(alt.LR[2],3)," with probabilities ", round(alt.psi,3),
                    " & ",(1- round(alt.psi,3)), ", respectively", sep = "" ))
      }
      
    }
    
    # computing wald's thresholds for LR
    
    # upper threshold
    wald.up = (1- type2)/type1
    upper.threshold = rep( wald.up, length(batch.seq))
    
    
    # lower threshold
    wald.low = type2/(1- type1)
    lower.threshold = rep( wald.low, length(batch.seq))
    
    
    
    # finding until which batch to compute
    stage = max(which(batch.seq<=length(obs)))
    
    
    # computing sequence of sufficient statistic
    t = cumsum(obs)
    cumsum.obs = t[batch.seq[1:stage]]
    
    
    # computing sequence of bayes factor or LR
    if(length(alt.LR)==1){
      LR.seq = mapply(FUN = LR.ber, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                      MoreArgs = list( null= null, alt= alt.LR))
    }else if(length(alt.LR)>1){
      
      LR.seq1 = mapply(FUN = LR.ber, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                       MoreArgs = list( null= null, alt= alt.LR[1]))
      LR.seq2 = mapply(FUN = LR.ber, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                       MoreArgs = list( null= null, alt= alt.LR[2]))
      LR.seq = (alt.psi*LR.seq1) + ((1- alt.psi)*LR.seq2)
    }
    
    
    # comparing sequence of bayes factor with thresholds
    if(length(term.thresh)==1){
      out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq, 
                   upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)
    }else{
      out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq, 
                   upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)
    }
    
    
    
    ## plots names
    testname = "Proportion test"
    ylabname = "Wtd. likelihood ratio in favor of the UMPBT alternative"
    
  }else if(test.type=="z.test"){
    print("-------------------------------------------------------------------------")
    print("Implementing the MSPRT in case of a Z-test")
    print("-------------------------------------------------------------------------")
    
    
    # direction of H1
    if(missing(side)==T){
      side = "right"
    }
    
    #default batch sequence
    if(missing(batch.seq)==T){
      batch.seq = 1:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",1,":N.max (default)", sep = ""))
      }
    }else{
      
      ##msg
      if(verbose==T){
        print("Working with batch.seq (user specified)")
      }
    }
    
    # type1
    if(missing(type1)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type1 = ",type1," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with user specified type1 = ",type1," (user specified)", sep = ""))
      }
    }
    
    
    # type2
    if(missing(type2)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (user specified)", sep = ""))
      }
    }
    
    # setting default null
    if(missing(null)==T){
      null = 0
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    if(missing(sigma0)==T){
      sigma0 = 1
      
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    
    # upper threshold
    wald.up = (1- type2)/type1
    upper.threshold = rep( wald.up, length(batch.seq))
    
    
    # lower threshold
    wald.low = type2/(1- type1)
    lower.threshold = rep( wald.low, length(batch.seq))
    
    
    alt.LR = find.umpbt.norm( side = side, type1 = type1, n.obs = N.max, null = null,
                              sigma0 = sigma0)
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      
    }
    
    
    # finding until which batch to compute
    stage = max(which(batch.seq<=length(obs)))
    
    
    # computing sequence of sufficient statistic
    t = cumsum(obs)
    cumsum.obs = t[batch.seq[1:stage]]
    
    
    # computing sequence of LR
    LR.seq = mapply(FUN = LR.norm, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                    MoreArgs = list( null= null, alt= alt.LR, sigma0=sigma0))
    
    
    # comparing sequence of bayes factor with thresholds
    
    out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq, 
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)
    
    
    
    ## plots names
    testname = "Z-test"
    ylabname = "Likelihood ratio in favor of the UMPBT alternative"
    
  }else if(test.type=="t.test"){
    print("-------------------------------------------------------------------------")
    print("Implementing the MSPRT in case of a T-test")
    print("-------------------------------------------------------------------------")
    
    
    # direction of H1
    if(missing(side)==T){
      side = "right"
    }
    
    #default batch sequence
    if(missing(batch.seq)==T){
      batch.seq = 2:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",2,":N.max (default)", sep = ""))
      }
    }else{
      
      ##msg
      if(verbose==T){
        print("Working with batch.seq (user specified)")
      }
    }
    
    # type1
    if(missing(type1)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type1 = ",type1," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with user specified type1 = ",type1," (user specified)", sep = ""))
      }
    }
    
    
    # type2
    if(missing(type2)==T){
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with type2 = ",type2," (user specified)", sep = ""))
      }
    }
    
    # setting default null
    if(missing(null)==T){
      null = 0
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified))", sep = ""))
      }
    }
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    # shifting the observation
    obs = obs - null
    null.shifted = 0
    
    
    # upper threshold
    wald.up = (1- type2)/type1
    upper.threshold = rep( wald.up, length(batch.seq))
    
    
    # lower threshold
    wald.low = type2/(1- type1)
    lower.threshold = rep( wald.low, length(batch.seq))
    
    
    
    # finding until which batch to compute
    stage = max(which(batch.seq<=length(obs)))
    
    
    # computing sequence of sufficient statistic
    t = cumsum(obs)
    cumsum.obs = t[batch.seq[1:stage]]
    
    
    # computing sequence of bayes factor
    LR.seq = alt.LR = numeric(stage)
    for (i in seq(stage)) {
      sx = sd(obs[1:batch.seq[i]])
      alt.LR[i] = find.umpbt.t(side = side, type1 = type1, n.obs = N.max, null = null.shifted, s = sx)
      LR.seq[i] = LR.t(m= batch.seq[i], suff.stat= cumsum.obs[i], null= null.shifted, alt= alt.LR[i], s=sx)
    }
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      print(paste("The UMPBT alternatives have been obtained"))
      
    }
    
    # comparing sequence of bayes factor with thresholds
    out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq, 
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)
    
    alt.LR = alt.LR + null
    ## plots names
    testname = "T-test"
    ylabname = "Bayes factor in favor of the UMPBT alternative"
  }
  
  if(verbose==T){
    print("-------------------------------------------------------------------------")
    print("Sequential comparison done!")
    print("-------------------------------------------------------------------------")
  }
  
  # plotting sequence of bayes factor together with thresholds
  
  if(plot.it==T){
    
    ylow = 0
    if(out$decision=="accept"){
      plot.title= paste(testname,": Accept Null (n =",out$n,")")
      yup = max(LR.seq)
    }else if( (out$decision=="reject") && (out$n<batch.seq[length(batch.seq)]) ){
      plot.title= paste(testname,": Reject Null (n =",out$n,")")
      yup = max(LR.seq)
    }else if( (out$decision=="reject") && (out$n==batch.seq[length(batch.seq)]) ){
      plot.title= paste(testname,": Reject Null (n =",out$n,")")
      yup = wald.up
    }else if(out$decision=="continue"){
      plot.title= paste(testname,": Continue sampling (n =",out$n,")")
      yup = wald.up
    }
    
    
    legend.label = c("Rejection Threshold","Acceptance Threshold","Likelihood ratio","Termination Threshold")
    
    
    plot.df = data.frame(xcol=batch.seq,ycol=upper.threshold,group=rep("a",length(batch.seq)))
    plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch.seq,ycol=lower.threshold,group=rep("b",length(batch.seq))))
    plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch.seq[1:stage],ycol=LR.seq,group=rep("c",stage)))
    plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch.seq[length(batch.seq)],ycol=term.thresh,group="s"))
    
    p = ggplot( plot.df, aes_string(x="xcol", y="ycol", color="group")) +
      scale_color_manual(labels=legend.label, values=c("firebrick1","chartreuse4","dodgerblue1","black")) +
      geom_segment(aes(x = batch.seq[length(batch.seq)], y = wald.low,
                       xend = batch.seq[length(batch.seq)], yend = term.thresh), color="chartreuse4", size=1.1) +
      geom_segment(aes(x = batch.seq[length(batch.seq)], y = wald.up,
                       xend = batch.seq[length(batch.seq)], yend = term.thresh), color="firebrick1", size=1.1) +
      geom_point( size = 3) +
      geom_line(size=1.1) +
      theme(plot.title = element_text(size=14, face="bold"),
            axis.title.x = element_text(size=14, face="bold"),
            axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),
            panel.background = element_rect(fill = "ivory2", colour = "ivory2"),
            legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
      guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,0), shape=16))) +
      ggtitle(plot.title) +
      xlab("Sample size") + ylab(ylabname) +
      ylim(ylow,yup)
    
    suppressWarnings(print(p))
  }
  
  if(test.type=="prop.test"){
    return(list("decision"= out$decision, "n" =out$n, "lhood.ratio"=LR.seq, "rej.threshold"=wald.up, "acc.threshold"=wald.low,
                "umpbt.alt" = alt.LR, "psi.umpbt" = alt.psi ))
  }else{
    return(list("decision"= out$decision, "n" =out$n, "lhood.ratio"=LR.seq, "rej.threshold"=wald.up, "acc.threshold"=wald.low,
                "umpbt.alt" = alt.LR ))
  }
  
}



##  arguments
# test.type         : test type; prop.test, z.test or t.test
# side              : direction of H1
# batch.seq         : a sequence of sample sizes where we observe data and make comparison
# gen.par           : vector of \theta at which observations are generated from
# null              : null value
# upper.threshold   : sequence of upper threshold
# lower.threshold   : sequence of lower threshold
# delta             : termination threshold
# type1             : specified Type-1 error
# type2             : specified Type-2 error
# N.max             : max affordable sample size
# repl              : no. of replications
# seed              : seed for randomization

## given a \gamma, based on R replications this returns estimated Type1/Type2 error and required avg sample size
## with a bunch of other quantities

OC.MSPRT = function( test.type, side, batch.seq, null, term.thresh, theta, sigma0,
                     type1= .005, type2= .2, N.max, verbose = T, repl, core.no){
  
  
  # checking if N.max is provided
  if(missing(N.max)==T){
    return("Maximum budget on sample size is not provided")
  }
  
  
  ## msg
  if(verbose==T){
    print("-------------------------------------------------------------------------")
    if(test.type=="prop.test"){
      print("MSPRT in case of a test for a binomial proportion")
    }else if(test.type=="z.test"){
      print("MSPRT in case of a Z-test")
    }else if(test.type=="t.test"){
      print("MSPRT in case of a T-test")
    }
    print("-------------------------------------------------------------------------")
  }
  
  
  ## setting default arguments
  
  # direction of H1
  if(missing(side)==T){
    side = "right"
    
    ## msg
    if(verbose==T){
      print("Working with side = right (default)")
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with side = ",side," (user specified)", sep = ""))
    }
  }
  
  
  # default batch sequence
  if(missing(batch.seq)==T){
    if(test.type!="t.test"){
      batch.seq = 1:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",1,":N.max (default)", sep = ""))
      }
    }else{
      batch.seq = 2:N.max
      
      ## msg
      if(verbose==T){
        print(paste("Working with batch.seq = ",2,":N.max (default)", sep = ""))
      }
    }
    
  }else{
    
    if((batch.seq[1]<2) && (test.type=="t.test")){
      return(print("batch.seq[1]<2 . Need at least 2 samples to compute sample standard deviation"))
    }else{
      
      ##msg
      if(verbose==T){
        print("Working with batch.seq (user specified)")
      }
    }
  }
  
  
  # type1
  if(missing(type1)==T){
    
    ## msg
    if(verbose==T){
      print(paste("Working with type1 = ",type1," (default)", sep = ""))
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with user specified type1 = ",type1," (user specified)", sep = ""))
    }
  }
  
  
  # type2
  if(missing(type2)==T){
    
    ## msg
    if(verbose==T){
      print(paste("Working with type2 = ",type2," (default)", sep = ""))
    }
  }else{
    
    ## msg
    if(verbose==T){
      print(paste("Working with type2 = ",type2," (user specified)", sep = ""))
    }
  }
  
  
  ##
  if(missing(repl)==T){
    if(test.type=="prop.test"){
      repl = 2e+6
    }else{repl = 1e+6}
  }else{repl = repl}
  # decimal accuracy
  k.dec = floor(log10(sqrt(repl)))
  
  
  # default number of cores
  if(missing(core.no)==T){
    c = detectCores()
    core.no = max(1,(c-1))
  }else{core.no=core.no}
  
  
  ## proptest
  if(test.type=="prop.test"){
    
    # setting default null
    if(missing(null)==T){
      null = 0.5
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    if(missing(theta)==T){
      theta = null
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (default)", sep = ""))
      }
    }else{
      theta = theta
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (user specified)", sep = ""))
      }
    }
    
    gen.par = theta
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    umpbt.out = find.umpbt.ber( side = side, type1 = type1, n.obs = N.max, null = null)
    alt.psi = umpbt.out$psi
    alt.LR = umpbt.out$u
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      if(length(alt.LR)==1){
        print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      }else if(length(alt.LR)>1){
        print(paste("Then UMPBT alternative has been obtained: ", round(alt.LR[1],3),
                    " & ",round(alt.LR[2],3)," with probabilities ", round(alt.psi,3),
                    " & ",(1- round(alt.psi,3)), " respectively", sep = ""))
      }
      
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    # performance at generating parameter
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Operating characteristics of the MSPRT at the desired parameter value is being calculated:")
    }
    
    
    if( ((side=="right") & (gen.par[1]<=null)) || ((side=="left") & (gen.par[1]>=null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.ber, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type1", batch.seq=batch.seq,
                                          null=null, gen.par=gen.par, alt.LR=alt.LR, 
                                          alt.psi = alt.psi, up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = inc[!is.na(inc)]
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type1 = type1)
      
      var.n0 = var(overshoot.summ$n.vec)
      k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
      avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
      
      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }else if( ((side=="right") & (gen.par[1]>null)) || ((side=="left") & (gen.par[1]<null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.ber, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type2", batch.seq=batch.seq,
                                          null=null, gen.par=gen.par, alt.LR=alt.LR, 
                                          alt.psi = alt.psi, up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = inc[!is.na(inc)]
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }
    
  }
  
  
  ## z test
  if(test.type=="z.test"){
    
    # setting default null
    if(missing(null)==T){
      null = 0
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    # setting default gen.par
    if(missing(theta)==T){
      
      theta = null
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (user specified)", sep = ""))
      }
    }
    
    if(missing(sigma0)==T){
      
      sigma0 = 1
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("Working with sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }
    
    gen.par = c(theta, sigma0)
    
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    
    alt.LR = find.umpbt.norm( side = side, type1 = type1, n.obs = N.max, null = null,
                              sigma0 = gen.par[2])
    
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      
      print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
      
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    # calculating overshooting summary of \tau <= N.max in case of SPRT
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Operating characteristics of the MSPRT at the desired parameter value is being calculated:")
    }
    
    
    if( ((side=="right") & (gen.par[1]<=null)) || ((side=="left") & (gen.par[1]>=null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.norm, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type1", batch.seq=batch.seq,
                                          null=null, gen.par=gen.par, alt.LR=alt.LR, 
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type1 = type1)
      
      var.n0 = var(overshoot.summ$n.vec)
      k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
      avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
      
      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }else if( ((side=="right") & (gen.par[1]>null)) || ((side=="left") & (gen.par[1]<null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.norm, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(error.type="type2", batch.seq=batch.seq,
                                          null=null, gen.par=gen.par, alt.LR=alt.LR, 
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }
    
  }
  
  
  
  ## t test
  if(test.type=="t.test"){
    
    # setting default null
    if(missing(null)==T){
      
      null = 0
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (default)", sep = ""))
      }
    }else{
      null = null
      ## msg
      if(verbose==T){
        print(paste("Working with null = ",null," (user specified)", sep = ""))
      }
    }
    
    
    if(missing(theta)==T){
      theta = null
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (default)", sep = ""))
      }
    }else{
      theta = theta
      ## msg
      if(verbose==T){
        print(paste("Working with theta = ",theta," (user specified)", sep = ""))
      }
    }
    
    # shifting the null
    gen.par = theta - null
    null = 0
    
    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }
    
    # computing wald's thresholds for LR
    
    wald.up = (1- type2)/type1
    wald.low = type2/(1- type1)
    
    
    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Operating characteristics of the MSPRT at the desired parameter value is being calculated:")
    }
    
    
    if( ((side=="right") & (gen.par[1]<=null)) || ((side=="left") & (gen.par[1]>=null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.t, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(side = side, error.type="type1", batch.seq=batch.seq,
                                          type1 = type1, null=null, gen.par=gen.par, 
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      ## msg
      if(verbose==T){
        print("The UMPBT alternatives has been obtained")
      }
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type1 = type1)
      
      var.n0 = var(overshoot.summ$n.vec)
      k.dec.n0 = abs(ceiling(log10(sqrt(var.n0/repl))))
      avg.n0 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n0)/(10^k.dec.n0)
      
      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }else if( ((side=="right") & (gen.par[1]>null)) || ((side=="left") & (gen.par[1]<null)) ){
      
      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)
      
      count = 0
      vec = n.vec = numeric()
      k.rep = repl/(1e+5)
      
      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
        
        count.t = 0
        vec.t = n.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.t, seed = ((k-1)*1e+5 + seq(1e+5)), 
                          MoreArgs = list(side = side, error.type="type2", batch.seq=batch.seq,
                                          type1 = type1, null=null, gen.par=gen.par, 
                                          up=wald.up, low=wald.low, N=N.max))
        
        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n = as.numeric(out.temp[3,])
        n.vec.t = c(n.vec.t,n)
        
        t = list( count.t, vec.t, n.vec.t)
      }
      
      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])
      
      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, "n.vec"=as.numeric(out[[3]]))
      
      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, R= repl, type2 = type2)
      
      var.n1 = var(overshoot.summ$n.vec)
      k.dec.n1 = abs(ceiling(log10(sqrt(var.n1/repl))))
      avg.n1 = floor(mean(overshoot.summ$n.vec)*10^k.dec.n1)/(10^k.dec.n1)
      
      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)
      
      ## msg
      if(verbose==T){
        print("Done")
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("-------------------------------------------------------------------------")
      }
      
    }
    
  }
  
  return(out)
}






