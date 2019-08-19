

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
# n     : sample size
# p0    : null value

## returns value of the objective function at p

objfunc.oneProp = function( p, delta, n, p0){
  
  num = log(delta) - (n*(log(1-p) - log(1- p0)))
  den = log(p/(1-p)) - log(p0 /(1- p0))
  t = num/den
  
  return(t)
}


### evaluates (LHS - c_0) as a function of \gamma as in appendix of sequential.pdf

##  arguments
# delta         : threshold for bayes factor to find UMPBT
# side          : direction of H1(right/left)
# n             : sample size
# p0            : null value
# opt.interval  : interval containing p under H1
# root          : c_0 to be specific as in supplemental file

## returns (LHS - c_0) as a function of \gamma as in supplemental file

find.threshold.oneProp = function( delta, side = "right", n, p0 = 0.5,
                                  opt.interval, root){
  
  if(side=="right"){
    
    if(missing(opt.interval)==T){
      opt.interval = c(p0, 1)
    }
    
    out = optimize( f=objfunc.oneProp, interval = opt.interval,
                    delta = delta, n = n, p0 = p0)
  }else if(side=="left"){
    
    if(missing(opt.interval)==T){
      opt.interval = c(0, p0)
    }
    
    out = optimize( f=objfunc.oneProp, interval = opt.interval, maximum = T,
                    delta = delta, n = n, p0 = p0)
  }
  
  return(out$objective -root)
}



### solving \gamma by matching UMP rejection region to find umpbt

##  arguments
# side  : direction of H1(right/left)
# type1 : specified type 1 error
# n     : sample size
# p0    : null value

## returns \gamma for which ump rejection region is matched


ump.match.oneProp = function( side = "right", type1 = 0.005, n, p0 = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n, p0) +1
    
    solve.out = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = "right", p0=p0,
                         opt.interval = c(p0, 1), root = (c0 -1) )
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n, prob = p0)
    quant.prob = pbinom( q = c0, size = n, prob = p0)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    
    solve.out = nleqslv( x=3, fn = find.threshold.oneProp, n = n, side = side, p0 = p0,
                         opt.interval = c(0, p0), root = (c0 +1) )
  }
  
  return(solve.out$x)
}



### finding mixture UMPBT alternative prior for exact proportion test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n       : sample size
# null    : null value

## returns UMPBT alternative prior by matching UMP rejection region

point.umpbt.oneProp = function( side = "right", type1 = 0.005, n, null = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n, null) +1
    solve.out = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                         opt.interval = c(null, 1), root = (c0 -1) )
    out = optimize( f = objfunc.oneProp, interval = c(null, 1),
                    delta = solve.out$x, n = n, p0 = null)
    umpbt = round(out$minimum, 4)
    
    return(umpbt)
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n, prob = null)
    quant.prob = pbinom( q = c0, size = n, prob = null)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    
    solve.out = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                         opt.interval = c(0, null), root = (c0 +1) )
    out = optimize( f = objfunc.oneProp, interval = c(0, null), maximum = T,
                    delta = solve.out$x, n = n, p0 = null)
    umpbt = round(out$maximum,4)
    
    return(umpbt)
  }
  
}

umpbt.oneProp = function( side = "right", type1 = 0.005, n, null = 0.5){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n, null) +1
    quant.prob = 1- pbinom( q = (c0 -1), size = n, prob = null)
    
    solve.out.r = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                           opt.interval = c(null, 1), root = (c0 -1) )
    out.r = optimize( f = objfunc.oneProp, interval = c(null, 1),
                      delta = solve.out.r$x, n = n, p0 = null)
    umpbt.r = round(out.r$minimum, 4)
    
    solve.out.l = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                           opt.interval = c(null, 1), root = (c0 -2) )
    out.l = optimize( f = objfunc.oneProp, interval = c(null, 1),
                      delta = solve.out.l$x, n = n, p0 = null)
    umpbt.l = round(out.l$minimum, 4)
    
    psi = (type1 - quant.prob)/dbinom(x=(c0 -1), size = n, prob = null)
    psi = round(psi, 5)
    
    return(list("u"=c(umpbt.l, umpbt.r), "psi" = psi))
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n, prob = null)
    quant.prob = pbinom( q = c0, size = n, prob = null)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    quant.prob = pbinom( q = c0, size = n, prob = null)
    
    
    solve.out.l = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                           opt.interval = c(0, null), root = (c0 +1) )
    out.l = optimize( f = objfunc.oneProp, interval = c(0, null), maximum = T,
                      delta = solve.out.l$x, n = n, p0 = null)
    umpbt.l = round(out.l$maximum, 4)
    
    solve.out.r = nleqslv( x=3, fn = find.threshold.oneProp, n=n, side = side, p0=null,
                           opt.interval = c(0, null), root = (c0 +2) )
    out.r = optimize( f = objfunc.oneProp, interval = c(0, null), maximum = T,
                      delta = solve.out.r$x, n = n, p0 = null)
    umpbt.r = round(out.r$maximum, 4)
    
    psi = (type1 - quant.prob)/dbinom(x=(c0 +1), size = n, prob = null)
    psi = round(psi, 5)
    
    return(list("u"=c(umpbt.r, umpbt.l), "psi" = psi))
  }
  
}



### finding UMPBT for one-sample z-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n       : sample size
# null    : null value
# sigma0  : known pop s.d.

## returns UMPBT alternative by matching UMP rejection region

umpbt.oneZ = function( side = "right", type1 = 0.005, n, null = 0, sigma0 = 1){
  
  z.alpha = qnorm( p = type1, lower.tail = F)
  
  if(side=="right"){
    umpbt = null + ((sigma0 * z.alpha)/(sqrt(n)))
  }else if(side=="left"){
    umpbt = null - ((sigma0 * z.alpha)/(sqrt(n)))
  }
  
  umpbt = round(umpbt, 4)
  
  return(umpbt)
}



### finding UMPBT for t-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n       : sample size
# null    : null value
# obs     : vector of observations
# s       : sample s.d.

## returns UMPBT alternative by matching UMP rejection region

umpbt.oneT = function( side = "right", type1 = 0.005, n, null = 0, obs, s){
  
  if(missing(s)==T){
    
    if(missing(obs)==T){
      return("Need to provide either 's' or 'obs'.")
    }else{
      s = sd(obs)
    }
    
  }else{
    if(missing(obs)==F) print("Ignoring 'obs'. Specifying 's' is enough.")
  }
  
  t.alpha = qt( type1, df = (n -1), lower.tail = F)
  
  if(side=="right"){
    umpbt = null + ((s * t.alpha)/(sqrt(n)))
  }else if(side=="left"){
    umpbt = null - ((s * t.alpha)/(sqrt(n)))
  }
  
  umpbt = round(umpbt, 4)
  
  return(umpbt)
}


### finding UMPBT for two-sample z-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n1      : sample size for 1st group
# n2      : sample size for 2nd group
# sigma0  : known common population s.d. for both groups

## returns UMPBT alternative by matching UMP rejection region

umpbt.twoZ = function( side = "right", type1 = 0.005, n1, n2, sigma0 = 1){
  
  z.alpha = qnorm( p = type1, lower.tail = F)
  log.gamma = (z.alpha^2)/2
  
  alt = sigma0*sqrt((2*(n1+n2)*log.gamma)/(n1*n2))
  
  if(side=="right"){
    umpbt = alt
  }else if(side=="left"){
    umpbt = -alt
  }
  
  umpbt = round(umpbt, 4)
  
  return(umpbt)
}


### finding UMPBT for two-sample t-test

##  arguments
# side    : direction of H1(right/left)
# type1   : specified type 1 error
# n1      : sample size for 1st group
# n2      : sample size for 2nd group
# obs1    : vector of observations for 1st group
# obs2    : vector of observations for 2nd group
# s       : pooled sample s.d. for both groups

## returns UMPBT alternative by matching UMP rejection region

umpbt.twoT = function( side = "right", type1 = 0.005, n1, n2,
                       obs1, obs2, s){
  
  if(missing(s)==T){
    
    n1.obs = length(obs1)
    n2.obs = length(obs2)
    
    s = sqrt((((n1.obs-1)*var(obs1))+((n2.obs-1)*var(obs2)))/(n1.obs+ n2.obs -2))
  }
  
  nu = n1 + n2 -2
  t.alpha = qt( type1, df = nu, lower.tail = F)
  # gamma.thresh = (((t.alpha^2)/(n1 + n2 -2)) +1)^((n1+n2)/2)
  # gamma.star = gamma.thresh^(2/(n1 + n2 -1)) -1
  
  # alt = s*sqrt(((n1+n2)*gamma.star*(n1+n2-2))/(n1*n2))
  
  gamma.star = ((((t.alpha^2)/nu) +1)^((n1 +n2)/(n1 +n2 -1))) -1
  alt = s*sqrt(((n1+ n2)*gamma.star*nu)/(n1*n2))
  
  if(side=="right"){
    umpbt = alt
  }else if(side=="left"){
    umpbt = -alt
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

LR.oneProp = function( m, suff.stat, null = 0.5, alt){
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

LR.oneZ = function( m, suff.stat, null = 0, alt, sigma0 = 1){
  
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

LR.oneT = function( m, suff.stat, null = 0, alt, s){
  t0 = (((suff.stat/m) -null)/s)
  t1 = (((suff.stat/m) -alt)/s)
  t = ((1+ ((m/(m-1))*(t0 ^2)) )/(1+ ((m/(m-1))*(t1 ^2)) ))^(m/2)
  return(t)
}


### likelihood ratio(LR) for two sample z-test

##  arguments
# m1          : sample size for group 1
# m2          : sample size for group 2
# suff.stat1  : value of mean for group 1
# suff.stat2  : value of mean for group 2
# alt         : alternative value
# sigma0      : known population s.d.

## retuns LR

LR.twoZ = function( m1, m2, suff.stat1, suff.stat2, alt, sigma0 = 1){
  
  z.stat = sqrt((m1*m2)/(m1+m2))*((suff.stat2 - suff.stat1)/sigma0)
  t1 = z.stat*(alt/sigma0)*sqrt((m1*m2)/(m1+m2)) - ((alt/sigma0)^2)*((m1*m2)/(2*(m1+m2)))
  out = exp(t1)
  
  return(out)
}


### likelihood ratio(LR) for two sample t-test

##  arguments
# m1          : sample size for group 1
# m2          : sample size for group 2
# suff.stat1  : value of mean for group 1
# suff.stat2  : value of mean for group 2
# alt         : alternative value
# s           : pooled sample s.d.

## retuns LR

LR.twoT = function( m1, m2, suff.stat1, suff.stat2, alt, s){
  
  t.stat = sqrt((m1*m2)/(m1+m2))*((suff.stat2 - suff.stat1)/s)
  t = (alt/s)*sqrt((m1*m2)/(m1+m2))
  out = (((m1+ m2 -2) + (t.stat^2))/((m1+ m2 -2) + ((t.stat- t)^2)))^((m1+m2)/2)
  
  return(out)
}


### computes (type-2 error -root) for proportion test

##  arguments
# alt       : alt value
# side      : direction of H1(right, left or both)
# null      : null value
# n         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.oneProp = function( alt, side = "right", null = 0.5, n, 
                                type1 = 0.005, root = 0){
  
  if(side=="right"){
    
    c0 = qbinom( (1- type1), n, null) +1
    t = pbinom( q= (c0 -1), size = n, prob = alt)
    
  }else if(side=="left"){
    
    c0 = qbinom( p = type1, size = n, prob = null)
    quant.prob = pbinom( q = c0, size = n, prob = null)
    if(quant.prob>type1){
      c0 = c0 -1
    }else if(quant.prob==type1){
      c0 = c0
    }
    t = 1- pbinom( q= c0, size = n, prob = alt)
    
  }
  
  return((t- root))
}



### computes (type-2 error -root) for z-test

##  arguments
# alt       : alt value
# side      : direction of H1(right, left or both)
# null      : null value
# sigma0    : known pop s.d.
# n         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.oneZ = function( alt, side = "right", null = 0, sigma0 = 1, n, 
                             type1 = 0.005, root = 0){
  
  if(side=="right"){
    z.alpha = qnorm( p=type1, lower.tail = F)
    c0 = null + ((z.alpha*sigma0)/sqrt(n))
    t = pnorm( c0, mean = alt, sd = (sigma0/sqrt(n)))
  }else if(side=="left"){
    z.alpha = qnorm( p=type1, lower.tail = F)
    c0 = null - ((z.alpha*sigma0)/sqrt(n))
    t = pnorm( c0, mean = alt, sd = (sigma0/sqrt(n)), lower.tail = F)
  }
  
  return(t -root)
}


### computes (type-2 error -root) for t-test

##  arguments
# alt       : alt value
# side      : direction of H1(right, left or both)
# null      : null value
# n         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.oneT = function( alt, side = "right", null = 0, n, 
                             type1 = 0.005, root = 0){
  
  if(side=="right"){
    t.alpha = qt( type1, df = (n -1), lower.tail = F)
    t = pt( t.alpha, df = (n-1), ncp = (sqrt(n)*(alt-null)) )
  }else if(side=="left"){
    t.alpha = qt( type1, df = (n -1), lower.tail = F)
    t = pt( -t.alpha, df = (n-1), ncp = (sqrt(n)*(alt-null)), lower.tail = F )
  }
  
  return(t -root)
}


### computes (type-2 error -root) for two-sample z-test

##  arguments
# alt       : alt value
# side      : direction of H1(right, left or both)
# sigma0    : known s.d.
# n1        : sample size for Group-1
# n2        : sample size for Group-2
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.twoZ = function( alt, side = "right", sigma0 = 1, n1, n2, 
                             type1 = 0.005, root = 0){
  
  z.alpha = qnorm( p=type1, lower.tail = F)
  
  if(side=="right"){
    
    t = pnorm((z.alpha - sqrt((n1*n2)/(n1+n2))*(alt/sigma0)))
  }else if(side=="left"){
    
    t = pnorm((z.alpha + sqrt((n1*n2)/(n1+n2))*(alt/sigma0)))
  }
  
  return(t -root)
}


### computes (type-2 error -root) for two sample t-test

##  arguments
# side      : direction of H1(right, left or both)
# alt       : alt value
# null      : null value
# N         : sample size
# type1     : specified Type-1 error
# root      : any value \in (0,1)

## retuns (type-2 error -root)

type2.error.twoT = function( alt, side = "right", n1, n2, type1 = 0.005, root = 0){
  
  t.alpha = qt( type1, df = (n1+n2-2), lower.tail = F)
  
  if(side=="right"){
    
    t = pt( t.alpha, df = (n1+n2-2), ncp = (sqrt((n1*n2)/(n1+n2))*alt) )
    
  }else if(side=="left"){
    
    t = pt( -t.alpha, df = (n1+n2-2), ncp = (sqrt((n1*n2)/(n1+n2))*alt), lower.tail = F)
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

find.alt = function( test.type, side = "right", null, n, n1, n2, 
                     type1 = 0.005, type2 = 0.2, sigma0 = 1){
  
  if((test.type!="oneProp") & (test.type!="oneZ") & (test.type!="oneT") &
     (test.type!="twoZ") & (test.type!="twoT")){
    return(print("'test type' can't be identified. Has to be one of 'oneProp', 'oneZ', 'oneT',
                 'twoZ' or 'twoT'."))
  }
  
  if(test.type=="oneProp"){
    
    if(missing(null)==T){
      null = 0.5
    }
    
    # ignoring n1
    if(missing(n1)==F) print("'n1' is ignored. Not required in a one-sample test.")
    
    # ignoring n2
    if(missing(n2)==F) print("'n2' is ignored. Not required in a one-sample test.")
    
    
    if(side=="right"){
      
      solve.out = uniroot( f = type2.error.oneProp, interval = c( null,1), side = side,
                           null=null, n=n, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
      
    }else if(side=="left"){
      
      solve.out = uniroot( f = type2.error.oneProp, interval = c( 0,null), side = side,
                           null=null, n=n, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
    }
    
  }else if(test.type=="oneZ"){
    
    if(missing(null)==T){
      null = 0
    }
    
    # ignoring n1
    if(missing(n1)==F) print("'n1' is ignored. Not required in a one-sample test.")
    
    # ignoring n2
    if(missing(n2)==F) print("'n2' is ignored. Not required in a one-sample test.")
    
    
    z.alpha = qnorm( p = type1, lower.tail = F)
    z.beta = qnorm( p = type2, lower.tail = F)
    
    if(side=="right"){
      
      alt = null + ((sigma0*(z.alpha+z.beta))/sqrt(n))
      
    }else if(side=="left"){
      
      alt = null - ((sigma0*(z.alpha+z.beta))/sqrt(n))
    }
    
  }else if(test.type=="oneT"){
    
    if(missing(null)==T){
      null = 0
    }
    
    # ignoring n1
    if(missing(n1)==F) print("'n1' is ignored. Not required in a one-sample test.")
    
    # ignoring n2
    if(missing(n2)==F) print("'n2' is ignored. Not required in a one-sample test.")
    
    
    if(side=="right"){
      
      solve.out = uniroot( f = type2.error.oneT, interval = c( null, .Machine$integer.max),
                           null=null, side=side, n=n, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
      
    }else if(side=="left"){
      
      solve.out = uniroot( f = type2.error.oneT, interval = c( -.Machine$integer.max, null),
                           null=null, side=side, n=n, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
    }
    
  }else if(test.type=="twoZ"){
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in two-sample tests.")
    
    # ignoring n
    if(missing(n)==F) print("'n' is ignored. Not required in a two-sample test.")
    
    
    z.alpha = qnorm( p = type1, lower.tail = F)
    z.beta = qnorm( p = type2, lower.tail = F)
    
    if(side=="right"){
      
      alt = sigma0*(z.alpha+z.beta)*sqrt((n1+n2)/(n1*n2))
      
    }else if(side=="left"){
      
      alt = -(sigma0*(z.alpha+z.beta)*sqrt((n1+n2)/(n1*n2)))
    }
    
  }else if(test.type=="twoT"){
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in two-sample tests.")
    
    # ignoring n
    if(missing(n)==F) print("'n' is ignored. Not required in a two-sample test.")
    
    
    if(side=="right"){
      
      solve.out = uniroot( f = type2.error.twoT, interval = c( 0, .Machine$integer.max),
                           side=side, n1=n1, n2=n2, type1 = type1, root = type2)
      alt = round(solve.out$root, 4)
      
    }else if(side=="left"){
      
      solve.out = uniroot( f = type2.error.twoT, interval = c( -.Machine$integer.max, 0),
                           side=side, n1=n1, n2=n2, type1 = type1, root = type2)
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


check = function( test.type, statistic, upper, lower, 
                  batch.seq, batch1.seq, batch2.seq, threshold){
  
  exit.stage = 0
  
  if( (test.type=="oneProp") || (test.type=="oneZ") || (test.type=="oneT") ){
    
    # checking the lengths of upper, lower and batch.seq
    if((length(upper)!=length(batch.seq)) || (length(lower)!=length(batch.seq))){
      
      return("length of 'batch.seq', 'upper' and 'lower' don't match!")
    }
    
    if(length(which(statistic>=upper[1:length(statistic)]))==0){
      
      if(length(which(statistic<=lower[1:length(statistic)]))==0){
        
        if(length(statistic)==length(batch.seq)){
          
          exit.stage = length(batch.seq)
          n = batch.seq[exit.stage]
          
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
        exit.stage = min( which(statistic<=lower[1:length(statistic)]) )
        n = batch.seq[exit.stage]
        
        decision = "accept"
      }
    }else{
      
      if(length(which(statistic<=lower[1:length(statistic)]))==0){
        
        exit.stage = min( which(statistic>=upper[1:length(statistic)]) )
        n = batch.seq[exit.stage]
        decision = "reject"
        
      }else{
        exit.stage = min(union(which(statistic>=upper[1:length(statistic)]),
                               which(statistic<=lower[1:length(statistic)])))
        n = batch.seq[exit.stage]
        
        if(min(which(statistic>=upper[1:length(statistic)])) <  min(which(statistic<=lower[1:length(statistic)]))){
          decision = "reject"
        }else{decision = "accept"}
      }
    }
    
    return(list("decision"=decision, "n"=n, "exit.stage" = exit.stage))
    
  }else if( (test.type=="twoZ") || (test.type=="twoT") ){
    
    # checking the lengths of batch1.seq & batch2.seq
    if(length(batch1.seq)!=length(batch2.seq)){
      
      return("length of 'batch1.seq' and 'batch2.seq' doesn't match!")
    }
    
    # checking the lengths of upper & lower with batch1.seq & batch2.seq
    if(length(upper)!=length(batch1.seq)){
      
      return("length of 'batch1.seq' and 'upper' doesn't match!")
    }
    if(length(lower)!=length(batch1.seq)){
      
      return("length of 'batch1.seq' and 'lower' doesn't match!")
    }
    if(length(upper)!=length(batch2.seq)){
      
      return("length of 'batch2.seq' and 'upper' doesn't match!")
    }
    if(length(lower)!=length(batch2.seq)){
      
      return("length of 'batch2.seq' and 'lower' doesn't match!")
    }
    
    if(length(which(statistic>=upper[1:length(statistic)]))==0){
      
      if(length(which(statistic<=lower[1:length(statistic)]))==0){
        
        if(length(statistic)==length(batch1.seq)){
          
          exit.stage = length(batch1.seq)
          n1 = batch1.seq[exit.stage]
          n2 = batch2.seq[exit.stage]
          
          if(statistic[length(batch1.seq)] >= threshold){
            decision = "reject"
          }else{
            decision = "accept"
          }
          
        }else{
          exit.stage = length(statistic)
          n1 = batch1.seq[exit.stage]
          n2 = batch2.seq[exit.stage]
          decision = "continue"
        }
        
        
      }else{
        exit.stage = min( which(statistic<=lower[1:length(statistic)]) )
        n1 = batch1.seq[exit.stage]
        n2 = batch2.seq[exit.stage]
        decision = "accept"
      }
    }else{
      
      if(length(which(statistic<=lower[1:length(statistic)]))==0){
        
        exit.stage = min( which(statistic>=upper[1:length(statistic)]) )
        n1 = batch1.seq[exit.stage]
        n2 = batch2.seq[exit.stage]
        decision = "reject"
      }else{
        
        exit.stage = min(union(which(statistic>=upper[1:length(statistic)]),
                               which(statistic<=lower[1:length(statistic)])))
        n1 = batch1.seq[exit.stage]
        n2 = batch2.seq[exit.stage]
        
        if(min(which(statistic>=upper[1:length(statistic)])) <  min(which(statistic<=lower[1:length(statistic)]))){
          decision = "reject"
        }else{decision = "accept"}
      }
    }
    
    return(list("decision"=decision, "n1"=n1, "n2"=n2, "exit.stage" = exit.stage))
  }
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


ovr.repl.oneProp = function( error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
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


ovr.repl.oneZ = function( error.type, batch.seq, null, gen.par, alt.LR,
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


ovr.repl.oneT = function( side, error.type, batch.seq, type1, null, gen.par,
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
    alt.LR = umpbt.oneT(side = side, type1 = type1, n = N, null = null, s = sx)
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



ovr.repl.twoZ = function( error.type, batch1.seq, batch2.seq, alt.LR,
                          gen.par, up, low, N1, N2, seed){
  
  n1 = n2 = incr.count = 0
  inconclusive = numeric()
  
  # data generation under gen.par
  set.seed(seed)
  x1.seq = rnorm(N1, mean = -gen.par[1], sd = gen.par[2])
  x2.seq = rnorm(N2, mean = gen.par[1], sd = gen.par[2])
  
  indx = 1
  while((indx<=length(batch1.seq)) && (n1==0) && (n2==0)){
    
    mean.obs1 = mean(x1.seq[1:batch1.seq[indx]])
    mean.obs2 = mean(x2.seq[1:batch2.seq[indx]])
    m1 = batch1.seq[indx]
    m2 = batch2.seq[indx]
    z.stat = sqrt((m1*m2)/(m1+m2))*((mean.obs2 - mean.obs1)/gen.par[2])
    t1 = z.stat*(alt.LR/gen.par[2])*sqrt((m1*m2)/(m1+m2)) - ((alt.LR/gen.par[2])^2)*((m1*m2)/(2*(m1+m2)))
    LR = exp(t1)
    
    if((error.type=="type1") && (LR>=up)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
      incr.count = 1
    }else if((error.type=="type2") && (LR>=up)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
    }else if((error.type=="type1") && (LR<=low)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
    }else if((error.type=="type2") && (LR<=low)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
      incr.count = 1
    }
    indx = indx +1
  }
  
  if((indx>length(batch1.seq)) && (n1==0) && (n2==0)){
    inconclusive = LR
    n1 = batch1.seq[length(batch1.seq)]
    n2 = batch2.seq[length(batch2.seq)]
  }
  
  return(list(incr.count,inconclusive, n1, n2))
}


ovr.repl.twoT = function( side, error.type, batch1.seq, batch2.seq,
                          type1, gen.par, up, low, N1, N2, seed){
  
  n1 = n2 = incr.count = 0
  inconclusive = numeric()
  
  # data generation under gen.par
  set.seed(seed)
  x1.seq = rnorm(N1, mean = -gen.par, sd = 1)
  x2.seq = rnorm(N2, mean = gen.par, sd = 1)
  
  indx = 1
  while((indx<=length(batch1.seq)) && (n1==0) && (n2==0)){
    
    mean.obs1 = mean(x1.seq[1:batch1.seq[indx]])
    mean.obs2 = mean(x2.seq[1:batch2.seq[indx]])
    m1 = batch1.seq[indx]
    m2 = batch2.seq[indx]
    s.pooled = sqrt((((m1-1)*var(x1.seq[1:batch1.seq[indx]]))+((m2-1)*var(x2.seq[1:batch2.seq[indx]])))/(m1+ m2 -2))
    t.stat = sqrt((m1*m2)/(m1+m2))*((mean.obs2 - mean.obs1)/s.pooled)
    alt.LR = umpbt.twoT( side = side, type1 = type1, n1 = N1, n2 = N2, s = s.pooled)
    t = (alt.LR/s.pooled)*sqrt((m1*m2)/(m1+m2))
    LR = (((m1+ m2 -2) + (t.stat^2))/((m1+ m2 -2) + ((t.stat- t)^2)))^((m1+m2)/2)
    
    if((error.type=="type1") && (LR>=up)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
      incr.count = 1
    }else if((error.type=="type2") && (LR>=up)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
    }else if((error.type=="type1") && (LR<=low)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
    }else if((error.type=="type2") && (LR<=low)){
      n1 = batch1.seq[indx]
      n2 = batch2.seq[indx]
      incr.count = 1
    }
    indx = indx +1
  }
  
  if((indx>length(batch1.seq)) && (n1==0) && (n2==0)){
    inconclusive = LR
    n1 = batch1.seq[length(batch1.seq)]
    n2 = batch2.seq[length(batch2.seq)]
  }
  
  return(list(incr.count,inconclusive, n1, n2))
}



# combined for all replications

overshoot.oneProp = function( error.type, batch.seq, null, gen.par, alt.LR, alt.psi,
                              up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.oneProp, seed = ((k-1)*1e+5 + seq(1e+5)), 
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



overshoot.oneZ = function( error.type, batch.seq, null, gen.par, alt.LR,
                           up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.oneZ, seed = ((k-1)*1e+5 + seq(1e+5)), 
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


overshoot.oneT = function( side, error.type, batch.seq, type1, null, gen.par,
                           up, low, N, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.oneT, seed = ((k-1)*1e+5 + seq(1e+5)), 
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



overshoot.twoZ = function( error.type, batch1.seq, batch2.seq, gen.par, alt.LR,
                           up, low, N1, N2, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n1.vec = n2.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n1.vec.t = n2.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.twoZ, seed = ((k-1)*1e+5 + seq(1e+5)), 
                      MoreArgs = list(error.type=error.type, batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                      alt.LR=alt.LR,  gen.par=gen.par, up=up, low=low,
                                      N1=N1, N2=N2))
    
    count.t = count.t + sum(as.numeric(out.temp[1,]))
    inc = as.numeric(out.temp[2,])
    vec.t = c(vec.t,inc[!is.na(inc)])
    n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
    n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))
    
    t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
  }
  
  inc = as.numeric(out[[2]])
  vec = c(vec,inc[!is.na(inc)])
  
  if(return.n==T){
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) ))
  }else{
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = as.numeric(out[[3]])))
  }
}



overshoot.twoT = function( side, error.type, batch1.seq, batch2.seq, type1, gen.par,
                           up, low, N1, N2, R, core.no, return.n = T){
  
  registerDoParallel(cores = core.no)
  
  count = 0
  vec = n1.vec = n2.vec = numeric()
  k.rep = R/(1e+5)
  
  out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {
    
    count.t = 0
    vec.t = n1.vec.t = n2.vec.t = numeric()
    out.temp = mapply(FUN = ovr.repl.twoT, seed = ((k-1)*1e+5 + seq(1e+5)), 
                      MoreArgs = list(side = side, error.type=error.type, batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                      type1 = type1, gen.par=gen.par, up=up, low=low, N1=N1, N2=N2))
    
    count.t = count.t + sum(as.numeric(out.temp[1,]))
    inc = as.numeric(out.temp[2,])
    vec.t = c(vec.t,inc[!is.na(inc)])
    n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
    n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))
    
    t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
  }
  
  inc = as.numeric(out[[2]])
  vec = c(vec,inc[!is.na(inc)])
  
  if(return.n==T){
    return(list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec, 
                "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) ))
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
    
    k.dec=1
    while(type1!=round(type1,k.dec)){
      k.dec = k.dec+1
    }
    k.dec = k.dec+1
    est = round(est, k.dec)
    
  }else if(error.type=="type2"){
    
    est = (count + sum(inconclusive.vec<delta))/R
    
    k.dec = 4
    est = round(est, k.dec)
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

design.MSPRT = function( test.type, side, batch.seq, batch1.seq, batch2.seq,
                         type1 =.005, type2 = .2,
                         null, sigma0 = 1, N.max, N1.max, N2.max, alt.comp,
                         repl, verbose=T, core.no){


  ## proptest
  if(test.type=="oneProp"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Designing the group sequential MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Designing the MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    

    # setting default null
    if(missing(null)==T){
      null = 0.5

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }


    ##
    if(missing(repl)==T){
      repl = 2e+6
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))
    # k.dec=1
    # while(type1!=round(type1,k.dec)){
    #   k.dec = k.dec+1
    # }
    # k.dec = k.dec+1
    

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    umpbt.out = umpbt.oneProp( side = side, type1 = type1, n = N.max, null = null)
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
      out.temp = mapply(FUN = ovr.repl.oneProp, seed = ((k-1)*1e+5 + seq(1e+5)),
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
    err.seq = mapply(FUN = error.summary, delta=vec.unique, 
                     MoreArgs = list(error.type = "type1", count= overshoot.summ$count, 
                                     inconclusive.vec= vec, root=0, R= repl, 
                                     type1 = type1))
    if(sum(err.seq<=type1)==0){
      
      delta.opt = floor(vec.unique[length(vec.unique)]*100)/100 + 0.01
      
    }else if(sum(err.seq>type1)==0){
      
      delta.opt = floor(wald.low*100)/100 + 0.01
      
    }else{
      
      i.opt = min(which(err.seq<=type1))
      delta.opt = delta.opt = floor(vec.unique[i.opt -1]*100)/100 + 0.01
    }
    

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

    type1.est = error.summary( error.type = "type1", delta = delta.opt, root = 0,
                               count= overshoot.summ$count, 
                               inconclusive.vec = overshoot.summ$inconclusive.vec,
                               R = repl, type1 = type1)
    avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)

    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ", type1.est, sep = ""))
      print("Done")
    }


    # finding alternative corresponding to fixed type1, type2, N

    if( (missing(alt.comp)==T) || (alt.comp==F) ){

      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0,
                   "umpbt.alt"=alt.LR, "psi.umpbt"=alt.psi,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"= delta.opt ))


    }else{

      if(alt.comp==T){

        alt.comp = find.alt(test.type=test.type, side = side, null = null, n = N.max,
                            type1 = type1, type2 = type2)

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
      alt.type2 = round(type2.error.oneProp( side = side, alt = alt.comp, null = null, 
                                             n=N.max, type1 = type1, root = 0), 4)


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
        out.temp = mapply(FUN = ovr.repl.oneProp, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta = delta.opt, root = 0,
                                 count= overshoot.summ$count, inconclusive.vec = overshoot.summ$inconclusive.vec,
                                 R = repl, type2 = type2)
      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ",type2.est,sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }


      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1,
                   "umpbt.alt"=alt.LR, "psi.umpbt"= alt.psi,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"= delta.opt ))

    }
  }


  ## z test
  if(test.type=="oneZ"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Designing the group sequential MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Designing the MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    

    # setting default null
    if(missing(null)==T){
      null = 0

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      null = null
      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }

    if(missing(sigma0)==T){
      sigma0 = 1

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{

      sigma0 = sigma0

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }


    ##
    if(missing(repl)==T){
      repl = 1e+6
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))
    # k.dec=1
    # while(type1!=round(type1,k.dec)){
    #   k.dec = k.dec+1
    # }
    # k.dec = k.dec+1

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    alt.LR = umpbt.oneZ( side = side, type1 = type1, n = N.max, null = null,
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
      out.temp = mapply(FUN = ovr.repl.oneZ, seed = ((k-1)*1e+5 + seq(1e+5)),
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

    type1.est = ((overshoot.summ$count + sum(overshoot.summ$inconclusive.vec>=delta.opt))/repl)
    avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)

    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ",type1.est,sep = ""))
      print("Done")
    }


    # finding alternative corresponding to fixed type1, type2, N

    if( (missing(alt.comp)==T) || (alt.comp==F) ){

      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0,
                   "umpbt.alt"=alt.LR,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }else{
      
      if(alt.comp==T){

        alt.comp = find.alt(test.type=test.type, side = side, null = null, n = N.max,
                            type1 = type1, type2 = type2, sigma0 = sigma0)

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
      alt.type2 = round(type2.error.oneZ( side = side, alt = alt.comp, null = null, sigma0 = sigma0,
                                          n=N.max, type1 = type1, root = 0), 4)


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
        out.temp = mapply(FUN = ovr.repl.oneZ, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)
      
      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1, "umpbt.alt"=alt.LR,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }
  }


  ## t test
  if(test.type=="oneT"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored and is set at 0. This is required in one-sample T-test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq[-1]!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Designing the group sequential MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Designing the MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 2:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 2, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch.seq[1]<2) return("Error! First batch size should be at least 2")
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    
    
    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    

    # setting default null
    null = 0
    
    ## msg
    if(verbose==T){
      print(paste("             null = ",null, sep = ""))
    }


    ##
    if(missing(repl)==T){
      repl = 1e+6
    }else{repl = repl}
    # k.dec=1
    # while(type1!=round(type1,k.dec)){
    #   k.dec = k.dec+1
    # }
    # k.dec = k.dec+1

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


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
      out.temp = mapply(FUN = ovr.repl.oneT, seed = ((k-1)*1e+5 + seq(1e+5)),
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

    type1.est = error.summary( error.type = "type1", delta= delta.opt, root = 0, 
                               count= overshoot.summ$count,
                               inconclusive.vec= overshoot.summ$inconclusive.vec, 
                               R= repl, type1 = type1)

    avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

    out.null.opt = list("type1.est"= type1.est, "avg.n0"= avg.n0)

    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ", type1.est, sep = ""))
      print("Done")
    }


    # finding alternative corresponding to fixed type1, type2, N

    if( (missing(alt.comp)==T) || (alt.comp==F) ){

      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "avg.n0"= out.null.opt$avg.n0,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }else{
      
      if(alt.comp==T){
        
        alt.comp = find.alt(test.type=test.type, side = side, null = null, n = N.max,
                            type1 = type1, type2 = type2)
        
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
      alt.type2 = round(type2.error.oneT( side = side, alt = alt.comp, null = null, n=N.max,
                                          type1 = type1, root = 0), 4)


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
        out.temp = mapply(FUN = ovr.repl.oneT, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)

      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out.alt.opt = list("type2.est"=type2.est, "avg.n1"=avg.n1)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }


      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n0"= out.null.opt$avg.n0, "avg.n1"=out.alt.opt$avg.n1,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }
  }



  ## two sample z test
  if(test.type=="twoZ"){
    
    
    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("Designing the group sequential MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Designing the MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 1:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 1, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 1:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 1, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # null
    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }

    if(missing(sigma0)==T){
      sigma0 = 1

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{

      sigma0 = sigma0

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }


    ##
    if(missing(repl)==T){
      repl = 1e+6
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))
    # k.dec=1
    # while(type1!=round(type1,k.dec)){
    #   k.dec = k.dec+1
    # }
    # k.dec = k.dec+1

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    alt.LR = umpbt.twoZ( side = side, type1 = type1, n1 = N1.max, n2 = N2.max, sigma0 = sigma0)


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
    vec = n1.vec = n2.vec = numeric()
    k.rep = repl/(1e+5)

    out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

      count.t = 0
      vec.t = n1.vec.t = n2.vec.t = numeric()
      out.temp = mapply(FUN = ovr.repl.twoZ, seed = ((k-1)*1e+5 + seq(1e+5)),
                        MoreArgs = list(error.type="type1", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                        alt.LR=alt.LR,  gen.par=c(0,sigma0), up=wald.up, low=wald.low,
                                        N1=N1.max, N2=N2.max))

      count.t = count.t + sum(as.numeric(out.temp[1,]))
      inc = as.numeric(out.temp[2,])
      vec.t = c(vec.t,inc[!is.na(inc)])
      n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
      n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

      t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
    }

    inc = as.numeric(out[[2]])
    vec = c(vec,inc[!is.na(inc)])

    overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                          "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )


    # increasingly ordered inconclusive vector
    vec = sort(overshoot.summ$inconclusive.vec)


    root.out = uniroot(f=error.summary, interval = c(wald.low, wald.up), error.type = "type1",
                       root = type1, count = overshoot.summ$count, inconclusive.vec = vec, R = repl,
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

    type1.est = error.summary( error.type = "type1", delta= delta.opt, root = 0, 
                               count= overshoot.summ$count,
                               inconclusive.vec= overshoot.summ$inconclusive.vec, 
                               R= repl, type1 = type1)

    avg.n1_0 = floor(mean(overshoot.summ$n1.vec)*100)/100
    avg.n2_0 = floor(mean(overshoot.summ$n2.vec)*100)/100

    out.null.opt = list("type1.est"= type1.est, "avg.n1_0"= avg.n1_0, "avg.n2_0"= avg.n2_0)

    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ",type1.est,sep = ""))
      print("Done")
    }


    # finding alternative corresponding to fixed type1, type2, N

    if( (missing(alt.comp)==T) || (alt.comp==F) ){

      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est,
                   "avg.n1_0"= out.null.opt$avg.n1_0, "avg.n2_0"= out.null.opt$avg.n2_0,
                   "umpbt.alt"=alt.LR,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }else{

      if(alt.comp==T){

        alt.comp = find.alt(test.type=test.type, side = side, n1 = N1.max, n2 = N2.max,
                            type1 = type1, type2 = type2, sigma0 = sigma0)

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
      alt.type2 = round(type2.error.twoZ( alt = alt.comp, side = side, sigma0 = sigma0,
                                          n1=N1.max, n2=N2.max, type1 = type1, root = 0), 4)


      # performance under the corresponding alternative

      ## msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("OC of the obtained MSPRT under the user desired alternative is being calculated:")
      }


      # calculating overshooting probability in SPRT under alternative
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoZ, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(error.type="type2", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          alt.LR=alt.LR,  gen.par=c((alt.comp/2),sigma0), up=wald.up, low=wald.low,
                                          N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )


      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)

      avg.n1_1 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_1 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out.alt.opt = list("type2.est"= type2.est, "avg.n1_1"= avg.n1_1, "avg.n2_1"= avg.n2_1)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n1_0"= out.null.opt$avg.n1_0, "avg.n2_0"= out.null.opt$avg.n2_0,
                   "avg.n1_1"= out.alt.opt$avg.n1_1, "avg.n2_1"= out.alt.opt$avg.n2_1,
                   "umpbt.alt"=alt.LR,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }
  }


  ## two sample T-test
  if(test.type=="twoT"){
    

    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq[-1]!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq[-1]!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("Designing the group sequential MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Designing the MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 2:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 2, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch1.seq[1]<2) return("Error! First batch size for Group-1 should be at least 2")
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 2:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 2, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch2.seq[1]<2) return("Error! First batch size for Group-2 should be at least 2")
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # null
    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }


    ##
    if(missing(repl)==T){
      repl = 1e+6
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))
    # k.dec=1
    # while(type1!=round(type1,k.dec)){
    #   k.dec = k.dec+1
    # }
    # k.dec = k.dec+1
    

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


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
    vec = n1.vec = n2.vec = numeric()
    k.rep = repl/(1e+5)

    out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

      count.t = 0
      vec.t = n1.vec.t = n2.vec.t = numeric()
      out.temp = mapply(FUN = ovr.repl.twoT, seed = ((k-1)*1e+5 + seq(1e+5)),
                        MoreArgs = list(side = side, error.type="type1", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                        type1 = type1, gen.par=0, up=wald.up, low=wald.low,
                                        N1=N1.max, N2=N2.max))

      count.t = count.t + sum(as.numeric(out.temp[1,]))
      inc = as.numeric(out.temp[2,])
      vec.t = c(vec.t,inc[!is.na(inc)])
      n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
      n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

      t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
    }

    inc = as.numeric(out[[2]])
    vec = c(vec,inc[!is.na(inc)])

    overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                          "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )


    # increasingly ordered inconclusive vector
    vec = sort(overshoot.summ$inconclusive.vec)


    root.out = uniroot(f=error.summary, interval = c(wald.low, wald.up), error.type = "type1",
                       root = type1, count = overshoot.summ$count, inconclusive.vec = vec, R = repl,
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

    type1.est = error.summary( error.type = "type1", delta= delta.opt, root = 0, 
                               count= overshoot.summ$count,
                               inconclusive.vec= overshoot.summ$inconclusive.vec, 
                               R= repl, type1 = type1)

    avg.n1_0 = floor(mean(overshoot.summ$n1.vec)*100)/100
    avg.n2_0 = floor(mean(overshoot.summ$n2.vec)*100)/100

    out.null.opt = list("type1.est"= type1.est, "avg.n1_0"= avg.n1_0, "avg.n2_0"= avg.n2_0)

    ## msg
    if(verbose==T){
      print(paste("Type 1 error probability is ",type1.est,sep = ""))
      print("Done")
    }


    # finding alternative corresponding to fixed type1, type2, N

    if( (missing(alt.comp)==T) || (alt.comp==F) ){

      if(verbose==T){
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est,
                   "avg.n1_0"= out.null.opt$avg.n1_0, "avg.n2_0"= out.null.opt$avg.n2_0,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }else{

      if(alt.comp==T){

        alt.comp = find.alt(test.type=test.type, side = side, n1 = N1.max, n2 = N2.max,
                            type1 = type1, type2 = type2)

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
      alt.type2 = round(type2.error.twoT( alt = alt.comp, side = side,
                                          n1=N1.max, n2=N2.max, type1 = type1, root = 0), 4)


      # performance under the corresponding alternative

      ## msg
      if(verbose==T){
        print("-------------------------------------------------------------------------")
        print("OC of the obtained MSPRT under the user desired alternative is being calculated:")
      }


      # calculating overshooting probability in SPRT under alternative
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoT, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(side = side, error.type="type2", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          type1 = type1, gen.par=(alt.comp/2), up=wald.up, low=wald.low,
                                          N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )


      type2.est = error.summary( error.type = "type2", delta= delta.opt, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)

      avg.n1_1 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_1 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out.alt.opt = list("type2.est"= type2.est, "avg.n1_1"= avg.n1_1, "avg.n2_1"= avg.n2_1)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }

      return(list( "type1.est"= out.null.opt$type1.est, "type2.est"=out.alt.opt$type2.est,
                   "avg.n1_0"= out.null.opt$avg.n1_0, "avg.n2_0"= out.null.opt$avg.n2_0,
                   "avg.n1_1"= out.alt.opt$avg.n1_1, "avg.n2_1"= out.alt.opt$avg.n2_1,
                   "alt"=alt.comp, "alt.type2"=alt.type2,
                   "rej.threshold" = wald.up, "acc.threshold" = wald.low, "term.thresh"=delta.opt))

    }
  }
}




### finding effective N in proportion test

effective.N = function(N, side = "right", type1 = 0.005, null = 0.5, plot.it = T){

  N.seq = seq(1, N, by = 1)

  # finding umpbt seq
  umpbt.seq = mapply(FUN = point.umpbt.oneProp, n = N.seq,
                     MoreArgs = list(side = side, type1 = type1, null = null))

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
           main="'Effective max sample size' for the MSPRT \n in one-sample proportion tests")
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
           main="'Effective max sample size' for the MSPRT \n in one-sample proportion tests")
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

find.samplesize = function( test.type, N, lower.signif = 0.05, higher.signif = 0.005,
                            null, side = "right", pow = 0.8, alt, sigma0 = 1,
                            n.seq, verbose=T, plot.it=T){

  # defaults
  if(missing(n.seq)==T){
    n.seq = seq(N, (4*N), by=1)
  }

  type2 = 1- pow
  if(test.type=="oneProp"){

    if(missing(null)==T){
      null = .5
    }

    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null, n = N,
                     type1 = lower.signif, type2 = type2)
    }

    t.seq = mapply(FUN = type2.error.oneProp, n = n.seq,
                   MoreArgs = list(alt = alt, side = side, null = null, 
                                   type1 = higher.signif))
    N.star = max(n.seq[which(t.seq>type2)]) +1

  }else if(test.type=="oneZ"){

    if(missing(null)==T){
      null = 0
    }

    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null, n = N,
                     sigma0 = sigma0, type1 = lower.signif, type2 = type2)
    }

    t.seq = mapply(FUN = type2.error.oneZ, n = n.seq,
                   MoreArgs = list(side = side, alt = alt, null = null,
                                   sigma0 = sigma0, type1 = higher.signif))
    N.star = max(n.seq[which(t.seq>type2)]) +1

  }else if(test.type=="oneT"){

    if(missing(null)==T){
      null = 0
    }

    if(missing(alt)==T){
      alt = find.alt(test.type = test.type, side = side, null = null, n = N,
                     type1 = lower.signif, type2 = type2)
    }

    t.seq = mapply(FUN = type2.error.oneT, n = n.seq,
                   MoreArgs = list(side = side, alt = alt, null = null,
                                   type1 = higher.signif))
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
  
  # msg
  if(verbose==T){
    
    print(paste("Using ", N.star, " samples at the significance level ",
                higher.signif,", the achieved power at ", round(alt,4), " is ",
                round(1- t.seq[max(which(t.seq>type2)) +1], 4) ))
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

implement.MSPRT = function( test.type, obs, obs1, obs2, side, batch.seq, batch1.seq, batch2.seq,
                            type1= .005, type2= .2, null, sigma0, term.thresh,
                            N.max, N1.max, N2.max, plot.it=T, verbose=T){

  if(test.type=="oneProp"){
    
    
    # ignoring obs1 & obs2
    if(missing(obs1)==F) print("'obs1' is ignored. Not required in a one-sample test.")
    if(missing(obs2)==F) print("'obs2' is ignored. Not required in a one-sample test.")
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Implementing the group sequential MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Implementing the MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
      }
    }
    

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

    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # setting default null
    if(missing(null)==T){
      null = 0.5

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    umpbt.out = umpbt.oneProp( side = side, type1 = type1, n = N.max, null = null)
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
      LR.seq = mapply(FUN = LR.oneProp, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                      MoreArgs = list( null= null, alt= alt.LR))
    }else if(length(alt.LR)>1){

      LR.seq1 = mapply(FUN = LR.oneProp, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                       MoreArgs = list( null= null, alt= alt.LR[1]))
      LR.seq2 = mapply(FUN = LR.oneProp, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                       MoreArgs = list( null= null, alt= alt.LR[2]))
      LR.seq = (alt.psi*LR.seq1) + ((1- alt.psi)*LR.seq2)
    }


    # comparing sequence of bayes factor with thresholds
    out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq,
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)


    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Sequential comparison done!")
      print("-------------------------------------------------------------------------")
    }


    ## plots names
    testname = "One-sample Proportion test"
    ylabname = "Wtd. likelihood ratio in favor of the UMPBT alternative"

    # plotting sequence of bayes factor together with thresholds

    if(plot.it==T){

      # ylow = 0
      # if(out$decision=="accept"){
        # plot.title= paste(testname,": Accept Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n<batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n==batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = wald.up
      # }else if(out$decision=="continue"){
      #   plot.title= paste(testname,": Continue sampling (n =",out$n,")")
      #   yup = wald.up
      # }
      
      if( (out$decision=="accept") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = min(LR.seq)
        yup = max(LR.seq)
      }else if( (out$decision=="accept") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if( (out$decision=="reject") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = max(LR.seq)
      }else if( (out$decision=="reject") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if(out$decision=="continue"){
        plot.title= paste(testname,": Continue sampling (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }


      legend.label = c("Rejection Threshold","Acceptance Threshold",
                       "Wtd. likelihood ratio","Termination Threshold")


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

    return(list("decision"= out$decision, "n" =out$n, "lhood.ratio"=LR.seq,
                "rej.threshold"=wald.up, "acc.threshold"=wald.low,
                "umpbt.alt" = alt.LR, "psi.umpbt" = alt.psi ))

  }else if(test.type=="oneZ"){
    
    
    # ignoring obs1 & obs2
    if(missing(obs1)==F) print("'obs1' is ignored. Not required in a one-sample test.")
    if(missing(obs2)==F) print("'obs2' is ignored. Not required in a one-sample test.")
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Implementing the group sequential MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Implementing the MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
    }


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

    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }


    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }

    # setting default null
    if(missing(null)==T){
      null = 0

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }

    if(missing(sigma0)==T){
      sigma0 = 1

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
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


    alt.LR = umpbt.oneZ( side = side, type1 = type1, n = N.max, null = null,
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
    LR.seq = mapply(FUN = LR.oneZ, m = batch.seq[1:stage], suff.stat = cumsum.obs,
                    MoreArgs = list( null= null, alt= alt.LR, sigma0=sigma0))


    # comparing sequence of bayes factor with thresholds

    out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq,
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)


    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Sequential comparison done!")
      print("-------------------------------------------------------------------------")
    }


    ## plots names
    testname = "One-sample Z-test"
    ylabname = "Likelihood ratio in favor of the UMPBT alternative"

    # plotting sequence of bayes factor together with thresholds

    if(plot.it==T){

      # ylow = 0
      # if(out$decision=="accept"){
      #   plot.title= paste(testname,": Accept Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n<batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n==batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = wald.up
      # }else if(out$decision=="continue"){
      #   plot.title= paste(testname,": Continue sampling (n =",out$n,")")
      #   yup = wald.up
      # }
      
      if( (out$decision=="accept") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = min(LR.seq)
        yup = max(LR.seq)
      }else if( (out$decision=="accept") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if( (out$decision=="reject") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = max(LR.seq)
      }else if( (out$decision=="reject") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if(out$decision=="continue"){
        plot.title= paste(testname,": Continue sampling (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }


      legend.label = c("Rejection Threshold","Acceptance Threshold",
                       "Likelihood ratio","Termination Threshold")


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

    return(list("decision"= out$decision, "n" =out$n, "lhood.ratio"=LR.seq, "rej.threshold"=wald.up, "acc.threshold"=wald.low,
                "umpbt.alt" = alt.LR ))

  }else if(test.type=="oneT"){
    
    
    # ignoring obs1 & obs2
    if(missing(obs1)==F) print("'obs1' is ignored. Not required in a one-sample test.")
    if(missing(obs2)==F) print("'obs2' is ignored. Not required in a one-sample test.")
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq[-1]!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("Implementing the group sequential MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Implementing the MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
    }


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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 2:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 2, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch.seq[1]<2) return("Error! First batch size should be at least 2")
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    
    
    # setting default null
    if(missing(null)==T){
      null = 0

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified))", sep = ""))
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
      alt.LR[i] = umpbt.oneT(side = side, type1 = type1, n = N.max, null = null.shifted, s = sx)
      LR.seq[i] = LR.oneT(m= batch.seq[i], suff.stat= cumsum.obs[i], null= null.shifted, alt= alt.LR[i], s=sx)
    }

    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")

      print(paste("The UMPBT alternatives have been obtained"))

    }

    # comparing sequence of Likelihood ratio with thresholds
    out = check( test.type = test.type, statistic = LR.seq, batch.seq = batch.seq,
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)

    alt.LR = alt.LR + null


    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Sequential comparison done!")
      print("-------------------------------------------------------------------------")
    }


    ## plots names
    testname = "One-sample T-test"
    ylabname = "Likelihood ratio in favor of the UMPBT alternative"

    # plotting sequence of bayes factor together with thresholds

    if(plot.it==T){

      # ylow = 0
      # if(out$decision=="accept"){
      #   plot.title= paste(testname,": Accept Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n<batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$n==batch.seq[length(batch.seq)]) ){
      #   plot.title= paste(testname,": Reject Null (n =",out$n,")")
      #   yup = wald.up
      # }else if(out$decision=="continue"){
      #   plot.title= paste(testname,": Continue sampling (n =",out$n,")")
      #   yup = wald.up
      # }
      
      if( (out$decision=="accept") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = min(LR.seq)
        yup = max(LR.seq)
      }else if( (out$decision=="accept") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Accept Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if( (out$decision=="reject") && (out$exit.stage<length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = max(LR.seq)
      }else if( (out$decision=="reject") && (out$exit.stage==length(batch.seq)) ){
        plot.title= paste(testname,": Reject Null (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }else if(out$decision=="continue"){
        plot.title= paste(testname,": Continue sampling (n =",out$n,")")
        ylow = 0
        yup = wald.up
      }


      legend.label = c("Rejection Threshold","Acceptance Threshold",
                       "Likelihood ratio","Termination Threshold")


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

    return(list("decision"= out$decision, "n" =out$n, "lhood.ratio"=LR.seq, "rej.threshold"=wald.up, "acc.threshold"=wald.low,
                "umpbt.alt" = alt.LR ))

  }else if(test.type=="twoZ"){
    
    
    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring obs
    if(missing(obs)==F) print("'obs' is ignored. Not required in a two-sample test.")
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("Implementing the group sequential MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Implementing the MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 1:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 1, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 1:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 1, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # null
    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }

    if(missing(sigma0)==T){
      sigma0 = 1

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{

      sigma0 = sigma0

      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }



    # upper threshold
    wald.up = (1- type2)/type1
    upper.threshold = rep( wald.up, length(batch1.seq))


    # lower threshold
    wald.low = type2/(1- type1)
    lower.threshold = rep( wald.low, length(batch1.seq))


    alt.LR = umpbt.twoZ( side = side, type1 = type1, n1 = N1.max, n2 = N2.max,
                               sigma0 = sigma0)

    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")

      print(paste("The UMPBT alternative ",round(alt.LR,3)," has been obtained"))
    }


    # finding until which batch to compute
    stage = max(which(batch1.seq<=length(obs1)))


    # computing sequence of sufficient statistic
    t1 = cumsum(obs1)
    t2 = cumsum(obs2)

    cumsum.obs1 = t1[batch1.seq[1:stage]]
    cumsum.obs2 = t2[batch2.seq[1:stage]]

    mean.obs1 = cumsum.obs1/batch1.seq[1:stage]
    mean.obs2 = cumsum.obs2/batch2.seq[1:stage]


    # computing sequence of LR
    LR.seq = mapply(FUN = LR.twoZ, m1 = batch1.seq[1:stage], m2 = batch2.seq[1:stage],
                    suff.stat1 = mean.obs1, suff.stat2 = mean.obs2,
                    MoreArgs = list( alt = alt.LR, sigma0 = sigma0))


    # comparing sequence of bayes factor with thresholds

    out = check( test.type = test.type, statistic = LR.seq,
                 batch1.seq = batch1.seq, batch2.seq = batch2.seq,
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)


    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Sequential comparison done!")
      print("-------------------------------------------------------------------------")
    }


    ## plots names
    testname = "Two-sample Z-test"
    ylabname = "Likelihood ratio in favor of the UMPBT alternative"

    # plotting sequence of bayes factor together with thresholds

    if(plot.it==T){

      # ylow = 0
      # if(out$decision=="accept"){
      #   plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$exit.stage<length(batch1.seq)) ){
      #   plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$exit.stage==length(batch1.seq)) ){
      #   plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = wald.up
      # }else if(out$decision=="continue"){
      #   plot.title= paste(testname,": Continue sampling (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = wald.up
      # }
      
      
      if( (out$decision=="accept") && (out$exit.stage<length(batch1.seq)) ){
        plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = min(LR.seq)
        yup = max(LR.seq)
      }else if( (out$decision=="accept") && (out$exit.stage==length(batch1.seq)) ){
        plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }else if( (out$decision=="reject") && (out$exit.stage<length(batch1.seq)) ){
        plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = max(LR.seq)
      }else if( (out$decision=="reject") && (out$exit.stage==length(batch1.seq)) ){
        plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }else if(out$decision=="continue"){
        plot.title= paste(testname,": Continue sampling (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }


      legend.label = c("Rejection Threshold","Acceptance Threshold",
                       "Likelihood ratio","Termination Threshold")


      plot.df = data.frame(xcol = seq(length(batch1.seq)), ycol = upper.threshold, group = rep("a",length(batch1.seq)))
      plot.df = rbind.data.frame(plot.df, data.frame(xcol = seq(length(batch1.seq)), ycol = lower.threshold,
                                                     group = rep("b",length(batch1.seq))))
      plot.df = rbind.data.frame(plot.df, data.frame(xcol = seq(stage), ycol = LR.seq, group = rep("c",stage)))
      plot.df = rbind.data.frame(plot.df, data.frame(xcol = length(batch1.seq), ycol = term.thresh, group="s"))

      p = ggplot( plot.df, aes_string(x="xcol", y="ycol", color="group")) +
        scale_color_manual(labels=legend.label, values=c("firebrick1","chartreuse4","dodgerblue1","black")) +
        geom_segment(aes(x = length(batch1.seq), y = wald.low, xend = length(batch1.seq), yend = term.thresh),
                     color="chartreuse4", size=1.1) +
        geom_segment(aes(x = length(batch1.seq), y = wald.up, xend = length(batch1.seq), yend = term.thresh),
                     color="firebrick1", size=1.1) +
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

    return(list("decision"= out$decision, "n1" = out$n1, "n2" = out$n2, "lhood.ratio" = LR.seq,
                "rej.threshold" = wald.up, "acc.threshold" = wald.low, "umpbt.alt" = alt.LR ))

  }else if(test.type=="twoT"){

    # checking if N1.max is provided
    if(missing(N1.max)==T){
      return("Maximum budget on sample size for Group-1 is not provided")
    }

    # checking if N2.max is provided
    if(missing(N2.max)==T){
      return("Maximum budget on sample size for Group-2 is not provided")
    }

    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring obs
    if(missing(obs)==F) print("'obs' is ignored. Not required in a two-sample test.")
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq[-1]!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq[-1]!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("Implementing the group sequential MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("Implementing the MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 2:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 2, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch1.seq[1]<2) return("Error! First batch size for Group-1 should be at least 2")
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 2:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 2, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch2.seq[1]<2) return("Error! First batch size for Group-2 should be at least 2")
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # null
    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }

    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    # upper threshold
    wald.up = (1- type2)/type1
    upper.threshold = rep( wald.up, length(batch1.seq))


    # lower threshold
    wald.low = type2/(1- type1)
    lower.threshold = rep( wald.low, length(batch1.seq))


    # finding until which batch to compute
    stage = max(which(batch1.seq<=length(obs1)))


    # computing sequence of sufficient statistic
    t1 = cumsum(obs1)
    t2 = cumsum(obs2)

    cumsum.obs1 = t1[batch1.seq[1:stage]]
    cumsum.obs2 = t2[batch2.seq[1:stage]]

    mean.obs1 = cumsum.obs1/batch1.seq[1:stage]
    mean.obs2 = cumsum.obs2/batch2.seq[1:stage]


    # computing sequence of bayes factor
    LR.seq = alt.LR = numeric(stage)
    for (i in seq(stage)) {

      m1 = batch1.seq[i]
      m2 = batch2.seq[i]
      s.pooled = sqrt((((m1-1)*var(obs1[1:batch1.seq[i]]))+((m2-1)*var(obs2[1:batch2.seq[i]])))/(m1+ m2 -2))
      t.stat = sqrt((m1*m2)/(m1+m2))*((mean.obs2[i] - mean.obs1[i])/s.pooled)
      alt.LR[i] = umpbt.twoT( side = side, type1 = type1, n1 = N1.max, n2 = N2.max, s = s.pooled)
      t = (alt.LR[i]/s.pooled)*sqrt((m1*m2)/(m1+m2))
      LR.seq[i] = (((m1+ m2 -2) + (t.stat^2))/((m1+ m2 -2) + ((t.stat- t)^2)))^((m1+m2)/2)
    }


    ## msg
    if(verbose==T){
      print("-------------------------------------------------------------------------")

      print(paste("UMPBT alternatives and likelihood ratios have been obtained"))

    }

    # comparing sequence of bayes factor with thresholds
    out = check( test.type = test.type, statistic = LR.seq,
                 batch1.seq = batch1.seq, batch2.seq = batch2.seq,
                 upper = upper.threshold, lower = lower.threshold, threshold = term.thresh)


    if(verbose==T){
      print("-------------------------------------------------------------------------")
      print("Sequential comparison done!")
      print("-------------------------------------------------------------------------")
    }


    ## plots names
    testname = "Two-sample T-test"
    ylabname = "Likelihood ratio in favor of the UMPBT alternative"

    # plotting sequence of bayes factor together with thresholds

    if(plot.it==T){

      # ylow = 0
      # if( (out$decision=="accept") && (out$exit.stage<length(batch1.seq)) ){
      #   plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$exit.stage<length(batch1.seq)) ){
      #   plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = max(LR.seq)
      # }else if( (out$decision=="reject") && (out$exit.stage==length(batch1.seq)) ){
      #   plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = wald.up
      # }else if(out$decision=="continue"){
      #   plot.title= paste(testname,": Continue sampling (n1 =",out$n1,", n2 =",out$n2,")")
      #   yup = wald.up
      # }


      legend.label = c("Rejection Threshold","Acceptance Threshold",
                       "Likelihood ratio","Termination Threshold")


      # plot.df = data.frame(xcol = seq(length(batch1.seq)), ycol = upper.threshold, group = rep("a",length(batch1.seq)))
      # plot.df = rbind.data.frame(plot.df, data.frame(xcol = seq(length(batch1.seq)), ycol = lower.threshold,
      #                                                group = rep("b",length(batch1.seq))))
      # # plot.df = rbind.data.frame(plot.df, data.frame(xcol = seq(stage), ycol = LR.seq, group = rep("c",stage)))
      # plot.df = rbind.data.frame(plot.df, data.frame(xcol = seq(length(LR.seq)), ycol = LR.seq, group = rep("c",stage)))
      # plot.df = rbind.data.frame(plot.df, data.frame(xcol = length(batch1.seq), ycol = term.thresh, group="s"))
      # 
      # p = ggplot( plot.df, aes_string(x="xcol", y="ycol", color="group")) +
      #   scale_color_manual(labels=legend.label, values=c("firebrick1","chartreuse4","dodgerblue1","black")) +
      #   geom_segment(aes(x = length(batch1.seq), y = wald.low, xend = length(batch1.seq), yend = term.thresh),
      #                color="chartreuse4", size=1.1) +
      #   geom_segment(aes(x = length(batch1.seq), y = wald.up, xend = length(batch1.seq), yend = term.thresh),
      #                color="firebrick1", size=1.1) +
      #   geom_point( size = 3) +
      #   geom_line(size=1.1) +
      #   theme(plot.title = element_text(size=14, face="bold"),
      #         axis.title.x = element_text(size=14, face="bold"),
      #         axis.title.y = element_text(size=14, face="bold"),legend.title=element_blank(),
      #         panel.background = element_rect(fill = "ivory2", colour = "ivory2"),
      #         legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue")) +
      #   guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,0), shape=16))) +
      #   ggtitle(plot.title) +
      #   xlab("Comparison stages") + ylab(ylabname) +
      #   ylim(ylow,yup)
      
      
      if( (out$decision=="accept") && (out$exit.stage<length(batch1.seq)) ){
        plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = min(LR.seq)
        yup = max(LR.seq)
      }else if( (out$decision=="accept") && (out$exit.stage==length(batch1.seq)) ){
        plot.title= paste(testname,": Accept Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }else if( (out$decision=="reject") && (out$exit.stage<length(batch1.seq)) ){
        plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = max(LR.seq)
      }else if( (out$decision=="reject") && (out$exit.stage==length(batch1.seq)) ){
        plot.title= paste(testname,": Reject Null (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }else if(out$decision=="continue"){
        plot.title= paste(testname,": Continue sampling (n1 =",out$n1,", n2 =",out$n2,")")
        ylow = 0
        yup = wald.up
      }
      
      plot.df = data.frame(xcol=batch1.seq,ycol=upper.threshold,group=rep("a",length(batch1.seq)))
      plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch1.seq,ycol=lower.threshold,group=rep("b",length(batch1.seq))))
      plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch1.seq[1:stage],ycol=LR.seq,group=rep("c",stage)))
      plot.df = rbind.data.frame(plot.df,data.frame(xcol=batch1.seq[length(batch1.seq)],ycol=term.thresh,group="s"))
      
      p = ggplot( plot.df, aes_string(x="xcol", y="ycol", color="group")) +
        scale_color_manual(labels=legend.label, values=c("firebrick1","chartreuse4","dodgerblue1","black")) +
        geom_segment(aes(x = batch1.seq[length(batch1.seq)], y = wald.low,
                         xend = batch1.seq[length(batch1.seq)], yend = term.thresh), color="chartreuse4", size=1.1) +
        geom_segment(aes(x = batch1.seq[length(batch1.seq)], y = wald.up,
                         xend = batch1.seq[length(batch1.seq)], yend = term.thresh), color="firebrick1", size=1.1) +
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

    return(list("decision"= out$decision, "n1" = out$n1, "n2" = out$n2, "lhood.ratio" = LR.seq,
                "rej.threshold" = wald.up, "acc.threshold" = wald.low, "umpbt.alt" = alt.LR ))

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

OC.MSPRT = function( test.type, side, batch.seq, batch1.seq, batch2.seq, null, term.thresh,
                     theta, sigma0, type1= .005, type2= .2, N.max, N1.max, N2.max,
                     verbose = T, repl, core.no){


  ## proptest
  if(test.type=="oneProp"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("The group sequential MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("The MSPRT in case of a one-sample proportion test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    

    # setting default null
    if(missing(null)==T){
      null = 0.5

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }


    if(missing(theta)==T){
      theta = null
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (default)", sep = ""))
      }
    }else{
      theta = theta
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (user specified)", sep = ""))
      }
    }

    gen.par = theta


    # replication number
    if(missing(repl)==T){
      if(test.type=="prop.test"){
        repl = 2e+6
      }else{repl = 1e+6}
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    umpbt.out = umpbt.oneProp( side = side, type1 = type1, n = N.max, null = null)
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
        out.temp = mapply(FUN = ovr.repl.oneProp, seed = ((k-1)*1e+5 + seq(1e+5)),
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
      
      # k.dec=1
      # while(type1!=round(type1,k.dec)){
      #   k.dec = k.dec+1
      # }
      # k.dec = k.dec+1
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0,
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type1 = type1)

      avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("Done")
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
        out.temp = mapply(FUN = ovr.repl.oneProp, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)

      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }
    }
  }


  ## z test
  if(test.type=="oneZ"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("The group sequential MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("The MSPRT in case of a one-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 1:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 1, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    

    # setting default null
    if(missing(null)==T){
      null = 0

      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{

      null = null
      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified)", sep = ""))
      }
    }

    # setting default gen.par
    if(missing(theta)==T){

      theta = null
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (user specified)", sep = ""))
      }
    }

    if(missing(sigma0)==T){

      sigma0 = 1
      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }

    gen.par = c(theta, sigma0)


    # replication number
    if(missing(repl)==T){
      if(test.type=="prop.test"){
        repl = 2e+6
      }else{repl = 1e+6}
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    alt.LR = umpbt.oneZ( side = side, type1 = type1, n = N.max, null = null,
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
        out.temp = mapply(FUN = ovr.repl.oneZ, seed = ((k-1)*1e+5 + seq(1e+5)),
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
      
      # k.dec=1
      # while(type1!=round(type1,k.dec)){
      #   k.dec = k.dec+1
      # }
      # k.dec = k.dec+1
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0,
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type1 = type1)
      
      avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("Done")
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
        out.temp = mapply(FUN = ovr.repl.oneZ, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)
      
      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }
    }
  }



  ## t test
  if(test.type=="oneT"){
    
    
    # ignoring batch1.seq & batch2.seq
    if(missing(batch1.seq)==F) print("'batch1.seq' is ignored. Not required in a one-sample test.")
    if(missing(batch2.seq)==F) print("'batch2.seq' is ignored. Not required in a one-sample test.")
    
    # ignoring N1.max & N2.max
    if(missing(N1.max)==F) print("'N1.max' is ignored. Not required in a one-sample test.")
    if(missing(N2.max)==F) print("'N2.max' is ignored. Not required in a one-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if((missing(batch.seq)==F) && (sum(batch.seq[-1]!=1)>0)){
        
        print("-------------------------------------------------------------------------")
        print("The group sequential MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("The MSPRT in case of a one-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch sequence and N.max
    if(missing(batch.seq)==T){
      
      if(missing(N.max)==T){
        return("Error! Need to specify at least batch.seq (batch sizes) or N.max (Maximum available sample size)")
      }else{
        
        batch.seq = 2:N.max
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print(paste("             batch.seq = ", 2, ":N.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch.seq[1]<2) return("Error! First batch size should be at least 2")
      
      if(missing(N.max)==T){
        
        N.max = sum(batch.seq)
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch.seq)!=N.max) return("Error! Sum of batch sizes should add up to N.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N.max = ", N.max, sep = ""))
          print("             batch.seq (user specified)")
        }
      }
      
      batch.seq = cumsum(batch.seq)
    }
    

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # setting default null
    if(missing(null)==T){
      null = 0
      
      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (default)", sep = ""))
      }
    }else{
      
      ## msg
      if(verbose==T){
        print(paste("             null = ",null," (user specified))", sep = ""))
      }
    }


    if(missing(theta)==T){
      theta = null
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (default)", sep = ""))
      }
    }else{
      theta = theta
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (user specified)", sep = ""))
      }
    }

    # shifting the null
    gen.par = theta - null
    null = 0


    # replication number
    if(missing(repl)==T){
      if(test.type=="prop.test"){
        repl = 2e+6
      }else{repl = 1e+6}
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


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
        out.temp = mapply(FUN = ovr.repl.oneT, seed = ((k-1)*1e+5 + seq(1e+5)),
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
      # k.dec=1
      # while(type1!=round(type1,k.dec)){
      #   k.dec = k.dec+1
      # }
      # k.dec = k.dec+1
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0,
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type1 = type1)
      
      avg.n0 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type1.est"= type1.est, "avg.n0"= avg.n0, "n0.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("Done")
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
        out.temp = mapply(FUN = ovr.repl.oneT, seed = ((k-1)*1e+5 + seq(1e+5)),
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

      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)
      
      avg.n1 = floor(mean(overshoot.summ$n.vec)*100)/100

      out = list("type2.est"=type2.est, "avg.n1"=avg.n1, "n1.vec" = overshoot.summ$n.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }
    }
  }


  ## two sample z test
  if(test.type=="twoZ"){
    

    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq[-1]!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq[-1]!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("The group sequential MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("The MSPRT in case of a two-sample Z-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 1:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 1, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 1:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 1, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }
    


    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }


    # null
    null = 0

    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }

    # setting default gen.par
    if(missing(theta)==T){

      theta = null
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (user specified)", sep = ""))
      }
    }

    if(missing(sigma0)==T){

      sigma0 = 1
      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("             sigma0 = ",sigma0," (user specified)", sep = ""))
      }
    }


    # replication number
    if(missing(repl)==T){
      if(test.type=="prop.test"){
        repl = 2e+6
      }else{repl = 1e+6}
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}



    ## msg
    if(verbose==T){
      print("Required parameters are taken care of!")
    }


    alt.LR = umpbt.twoZ( side = side, type1 = type1, n1 = N1.max, n2 = N2.max, sigma0 = sigma0)


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
      print("Operating characteristics of the MSPRT at the desired parameter value is being calculated:")
    }


    if( ((side=="right") & (theta<=null)) || ((side=="left") & (theta>=null)) ){

      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoZ, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(error.type="type1", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          alt.LR=alt.LR,  gen.par=c((theta/2), sigma0),
                                          up=wald.up, low=wald.low, N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )
      # k.dec=1
      # while(type1!=round(type1,k.dec)){
      #   k.dec = k.dec+1
      # }
      # k.dec = k.dec+1
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0,
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type1 = type1)

      avg.n1_0 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_0 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out = list("type1.est"= type1.est, "avg.n1_0"= avg.n1_0, "avg.n2_0"= avg.n2_0,
                 "n1_0.vec" = overshoot.summ$n1.vec, "n2_0.vec" = overshoot.summ$n2.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }

    }else if( ((side=="right") & (theta>null)) || ((side=="left") & (theta<null)) ){

      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoZ, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(error.type="type2", batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          alt.LR=alt.LR,  gen.par=c((theta/2), sigma0),
                                          up=wald.up, low=wald.low, N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )

      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)
      
      avg.n1_1 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_1 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out = list("type2.est"= type2.est, "avg.n1_1"= avg.n1_1, "avg.n2_1"= avg.n2_1,
                 "n1_1.vec" = overshoot.summ$n1.vec, "n2_1.vec" = overshoot.summ$n2.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }
    }
  }



  ## two sample t test
  if(test.type=="twoT"){
    
    
    # checking if length(batch1.seq) and length(batch2.seq) are equal
    if( (missing(batch1.seq)==F) && (missing(batch2.seq)==F) && (length(batch1.seq)!=length(batch2.seq)) ){
      return("Lenghts of batch1.seq and batch2.seq are not matching. They should be same.")
    }
    
    # ignoring batch.seq
    if(missing(batch.seq)==F) print("'batch.seq' is ignored. Not required in a two-sample test.")
    
    # ignoring null
    if(missing(null)==F) print("'null' is ignored. This package only works with null=0 in a two-sample test.")
    
    # ignoring N.max
    if(missing(N.max)==F) print("'N.max' is ignored. Not required in a two-sample test.")
    
    
    
    ## msg
    if(verbose==T){
      
      if(((missing(batch1.seq)==F) && (sum(batch1.seq[-1]!=1)>0))||((missing(batch2.seq)==F) && (sum(batch2.seq[-1]!=1)>0))){
        
        print("-------------------------------------------------------------------------")
        print("The group sequential MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
        
      }else{
        
        print("-------------------------------------------------------------------------")
        print("The MSPRT in case of a two-sample T-test:")
        print("-------------------------------------------------------------------------")
      }
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
    
    
    #default batch1 sequence and N1.max
    if(missing(batch1.seq)==T){
      
      if(missing(N1.max)==T){
        return("Error! Need to specify at least batch1.seq (batch sizes for Group-1) or N1.max (Maximum available sample size for Group-1)")
      }else{
        
        batch1.seq = 2:N1.max
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print(paste("             batch1.seq = ", 2, ":N1.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch1.seq[1]<2) return("Error! First batch size for Group-1 should be at least 2")
      
      if(missing(N1.max)==T){
        
        N1.max = sum(batch1.seq)
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch1.seq)!=N1.max) return("Error! Sum of batch sizes for Group-1 should add up to N1.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N1.max = ", N1.max, sep = ""))
          print("             batch1.seq (user specified)")
        }
      }
      
      batch1.seq = cumsum(batch1.seq)
    }
    
    
    #default batch2 sequence and N2.max
    if(missing(batch2.seq)==T){
      
      if(missing(N2.max)==T){
        return("Error! Need to specify at least batch2.seq (batch sizes for Group-2) or N2.max (Maximum available sample size for Group-2)")
      }else{
        
        batch2.seq = 2:N2.max
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print(paste("             batch2.seq = ", 2, ":N2.max (default)", sep = ""))
        }
      }
      
    }else{
      
      if(batch2.seq[1]<2) return("Error! First batch size for Group-2 should be at least 2")
      
      if(missing(N2.max)==T){
        
        N2.max = sum(batch2.seq)
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
        
      }else{
        
        if(sum(batch2.seq)!=N2.max) return("Error! Sum of batch sizes for Group-2 should add up to N2.max")
        
        ## msg
        if(verbose==T){
          print(paste("             N2.max = ", N2.max, sep = ""))
          print("             batch2.seq (user specified)")
        }
      }
      
      batch2.seq = cumsum(batch2.seq)
    }

    # type1
    if(missing(type1)==T){

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type1 = ",type1," (user specified)", sep = ""))
      }
    }


    # type2
    if(missing(type2)==T){

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (default)", sep = ""))
      }
    }else{

      ## msg
      if(verbose==T){
        print(paste("             type2 = ",type2," (user specified)", sep = ""))
      }
    }
    


    # null
    null = 0

    ## msg
    if(verbose==T){
      print(paste("             null = ",0, sep = ""))
    }


    # setting default gen.par
    if(missing(theta)==T){

      theta = null

      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (default)", sep = ""))
      }
    }else{
      ## msg
      if(verbose==T){
        print(paste("             theta = ",theta," (user specified)", sep = ""))
      }
    }


    # replication number
    if(missing(repl)==T){
      if(test.type=="prop.test"){
        repl = 2e+6
      }else{repl = 1e+6}
    }else{repl = repl}
    # k.dec = floor(log(sqrt(repl)))

    # default number of cores
    if(missing(core.no)==T){
      c = detectCores()
      core.no = max(1,(c-1))
    }else{core.no=core.no}


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
      print("Operating characteristics of the MSPRT at the desired parameter value is being calculated:")
    }


    if( ((side=="right") & (theta<=null)) || ((side=="left") & (theta>=null)) ){

      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoT, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(side = side, error.type="type1",
                                          batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          type1 = type1, gen.par=theta/2, up=wald.up, low=wald.low,
                                          N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )

      # k.dec=1
      # while(type1!=round(type1,k.dec)){
      #   k.dec = k.dec+1
      # }
      # k.dec = k.dec+1
      
      type1.est = error.summary( error.type = "type1", delta= term.thresh, root = 0,
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type1 = type1)
      
      avg.n1_0 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_0 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out = list("type1.est"= type1.est, "avg.n1_0"= avg.n1_0, "avg.n2_0"= avg.n2_0,
                 "n1_0.vec" = overshoot.summ$n1.vec, "n2_0.vec" = overshoot.summ$n2.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 1 error probability is ", type1.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }

    }else if( ((side=="right") & (theta>null)) || ((side=="left") & (theta<null)) ){

      # calculating overshooting probability in SPRT
      registerDoParallel(cores = core.no)

      count = 0
      vec = n1.vec = n2.vec = numeric()
      k.rep = repl/(1e+5)

      out = foreach( k = seq(k.rep), .combine = compmerge.list) %dopar% {

        count.t = 0
        vec.t = n1.vec.t = n2.vec.t = numeric()
        out.temp = mapply(FUN = ovr.repl.twoT, seed = ((k-1)*1e+5 + seq(1e+5)),
                          MoreArgs = list(side = side, error.type="type2",
                                          batch1.seq=batch1.seq, batch2.seq=batch2.seq,
                                          type1 = type1, gen.par=theta/2, up=wald.up, low=wald.low,
                                          N1=N1.max, N2=N2.max))

        count.t = count.t + sum(as.numeric(out.temp[1,]))
        inc = as.numeric(out.temp[2,])
        vec.t = c(vec.t,inc[!is.na(inc)])
        n1.vec.t = c(n1.vec.t, as.numeric(out.temp[3,]))
        n2.vec.t = c(n2.vec.t, as.numeric(out.temp[4,]))

        t = list( count.t, vec.t, n1.vec.t, n2.vec.t)
      }

      inc = as.numeric(out[[2]])
      vec = c(vec,inc[!is.na(inc)])

      overshoot.summ = list("count" = sum(as.numeric(out[[1]])), "inconclusive.vec" = vec,
                            "n1.vec"=as.numeric(out[[3]]), "n2.vec"=as.numeric(out[[4]]) )

      type2.est = error.summary( error.type = "type2", delta= term.thresh, root = 0, 
                                 count= overshoot.summ$count,
                                 inconclusive.vec= overshoot.summ$inconclusive.vec, 
                                 R= repl, type2 = type2)
      
      avg.n1_1 = floor(mean(overshoot.summ$n1.vec)*100)/100
      avg.n2_1 = floor(mean(overshoot.summ$n2.vec)*100)/100

      out = list("type2.est"= type2.est, "avg.n1_1"= avg.n1_1, "avg.n2_1"= avg.n2_1,
                 "n1_1.vec" = overshoot.summ$n1.vec, "n2_1.vec" = overshoot.summ$n2.vec)

      ## msg
      if(verbose==T){
        print(paste("Type 2 error probability is ", type2.est, sep = ""))
        print("Done")
        print("-------------------------------------------------------------------------")
      }
    }
  }

  return(out)
}






