library('tidyverse')
library('stringr')
library('nloptr')
library('extRemes')

##### Functions necessary to initialize modified profile likelihood #####

ascovar = function(l, p, g, s, t){
  covar = matrix(nrow = 4, ncol = 4)
  covar[1,] = c(l/t, 0, 0, 0)
  covar[2,] = c(0, p*(1-p)/(l*t), 0, 0)
  covar[3,] = c(0, 0, ((1+g)^2)/(l*t*p), (s*(1+g))/(l*t*p))
  covar[4,] = c(0, 0, (s*(1+g))/(l*t*p), 2*((1+g)*s^2)/(l*t*p))
  return(covar)
}

pois.est = function(data, thresh, measure, x){
  thresh = -thresh
  z = -x-thresh
  gpddata = data[data[measure] != -1,]
  gpddata = gpddata[gpddata[measure] < -thresh,]
  
  b00 = as.POSIXct(paste(str_sub(data[['Date']][1], 1, 10), "00:00:01"),
                   tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  bend = as.POSIXct(paste(str_sub(data[['Date']][length(data[,1])],1,10), "23:59:59"),
                    tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  
  t = as.numeric(difftime(bend, b00, units = "hours"))
  
  l = length(data[,1])/t
  p = sum(-gpddata[[measure]] > thresh)/length(data[,1])
  gpdobj = fevd(-gpddata[[measure]], threshold = thresh, type = "GP", method = "MLE")
  g = gpdobj$results$par[['shape']]
  s = gpdobj$results$par[['scale']]
  lowend = -thresh+s/g
  estimates = matrix(c(l, p, g, s, lowend), nrow = 1, ncol = 5)
  colnames(estimates) = c("lambda", "Pi", "Gamma", "Sigma", "Lower endpoint")
  
  
  #standard errors and confidence intervals for parameters
  c = ascovar(l, p, g, s, t)
  std.err = matrix(c(sqrt(c[1,1]), sqrt(c[2,2]), sqrt(c[3,3]), sqrt(c[4,4])),
                   nrow = 1, ncol = 4)
  colnames(std.err) = c("lambda", "Pi", "Gamma", "Sigma")
  
  
  confint = matrix(nrow = 4, ncol = 2)
  confint[1,] = c(l-1.96*std.err["lambda"], l+1.96*std.err["lambda"])
  confint[2,] = c(p-1.96*std.err["Pi"], p+1.96*std.err["Pi"])
  confint[3,] = c(g-1.96*std.err["Gamma"], g+1.96*std.err["Gamma"])
  confint[4,] = c(s-1.96*std.err["Sigma"], s+1.96*std.err["Sigma"])
  
  # Estimate of lambda for crashes and confidence interval using delta method
  lcrash = l*t*p*(1+g*(z)/s)^(-1/g)
  
  h = matrix(nrow = 1, ncol = 4)
  h[1] = t*p*(1+g*(z/s))^(-1/g)
  h[2] = l*t*(1+g*(z/s))^(-1/g)
  h[3] = l*t*p*((1/(g^2))*log(1+g*z/s)-(z/(s*g+z*g^2)))*(1+g*z/s)^(-1/g)
  h[4] = l*t*p*(z/(s^2))*(1+g*z/s)^(-1/g - 1)
  
  vlcrash = h%*%c%*%t(h)
  
  l.lcrash = lcrash - 1.96*sqrt(vlcrash)
  u.lcrash = lcrash + 1.96*sqrt(vlcrash)
  lcrashconf = matrix(c(l.lcrash, u.lcrash), nrow = 1, ncol = 2)
  
  l.pcrash = qpois(0.025, lambda = lcrash)
  u.pcrash = qpois(0.975, lambda = lcrash)
  pcrash = matrix(c(l.pcrash, u.pcrash), nrow = 1, ncol = 2)
  
  out = list("par" = estimates, "std. error" = std.err, "Crash intensity" = lcrash, 
             "CI Crash Intensity" = lcrashconf, "PI Crash" = pcrash, "time" = t)
  return(out)
  
}

### confidence intervals for pi_ne - pi_dk based on modified profile likelihood ####


mpl_diffpi = function(dataNE, dataDK, uNE, uDK, q, measure1){
  thresh1 = uDK
  thresh2 = uNE
  x1 = q
  
  N1 = length(dataDK[[measure1]])
  data1 = dataDK[dataDK[[measure1]] != -1, ]
  data.s = data1[data1[[measure1]]<thresh1, ]
  # z1 = -data.s[[measure1]]+thresh1
  q1 = -x1+thresh1
  # y1 = length(z1)
  
  poisson_time_dk = as.POSIXct(dataDK[["Date"]], tz = "", 
                               tryFormats = "%Y-%m-%d %H:%M:%OS")
  b00.dk = as.POSIXct(paste(str_sub(dataDK[['Date']][1], 1, 10), "00:00:01"),
                      tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  poisson.time.dk = as.numeric(difftime(poisson_time_dk, b00.dk, units = "hours"))
  t.final.dk = as.numeric(difftime(
    as.POSIXct(paste(str_sub(dataDK[['Date']][length(dataDK[,1])],1,10), "23:59:59"),
               tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS"), 
    b00.dk, units = "hours"))
  
  dups.dk = duplicated(poisson.time.dk)
  for(i in 1:length(poisson.time.dk)){
    if(dups.dk[i]==TRUE){
      poisson.time.dk[i] = poisson.time.dk[i-1]+0.0000001
    } else {
      poisson.time.dk[i] = poisson.time.dk[i]
    }
  }
  
  n.dk = matrix(nrow = N1, ncol = 1)
  z1 = matrix(nrow = N1, ncol = 1)
  for(i in 1:N1){
    if(dataDK[[measure1]][i] == -1|dataDK[[measure1]][i]>thresh1){
      n.dk[i]=0
      z1[i] = NaN
    } else{
      n.dk[i]=1
      z1[i] = -dataDK[[measure1]][i]+thresh1
    }
  }
  
  N2 = length(dataNE[[measure1]])
  data2 = dataNE[dataNE[[measure1]] != -1, ]
  data.s = data2[data2[[measure1]]<thresh2, ]
  # z2 = -data.s[[measure1]]+thresh2
  q2 = -x1+thresh2
  # y2 = length(z2)
  
  poisson_time_ne = as.POSIXct(dataNE[["Date"]], tz = "", 
                               tryFormats = "%Y-%m-%d %H:%M:%OS")
  b00.ne = as.POSIXct(paste(str_sub(dataNE[['Date']][1], 1, 10), "00:00:01"),
                      tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  poisson.time.ne = as.numeric(difftime(poisson_time_ne, b00.ne, units = "hours"))
  
  dups.ne = duplicated(poisson.time.ne)
  for(i in 1:length(poisson.time.ne)){
    if(dups.ne[i]==TRUE){
      poisson.time.ne[i] = poisson.time.ne[i-1]+0.0000001
    } else {
      poisson.time.ne[i] = poisson.time.ne[i]
    }
  }
  
  t.final.ne = as.numeric(difftime(
    as.POSIXct(paste(str_sub(dataNE[['Date']][length(dataNE[,1])],1,10), "23:59:59"),
               tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS"), 
    b00.ne, units = "hours"))
  
  n.ne = matrix(nrow = N2, ncol = 1)
  z2 = matrix(nrow = N2, ncol = 1)
  for(i in 1:N2){
    if(dataNE[[measure1]][i] == -1|dataNE[[measure1]][i]>thresh2){
      n.ne[i]=0
      z2[i] = NaN
    } else{
      n.ne[i]=1
      z2[i] = -dataNE[[measure1]][i]+thresh2
    }
  }
  
  
  ### ML estimates for DK and NE ###
  
  ml.est.dk = pois.est(data = dataDK, measure = measure1, x = q,
                       thresh = uDK)$par
  ml.est.ne = pois.est(data = dataNE, measure = measure1, x = q,
                       thresh = uNE)$par
  
  ml.pi.dk = ml.est.dk[2]*(1+ml.est.dk[3]*q1/ml.est.dk[4])^(-1/ml.est.dk[3])
  ml.pi.ne = ml.est.ne[2]*(1+ml.est.ne[3]*q2/ml.est.ne[4])^(-1/ml.est.ne[3])
  
  ### Setup matrix computations ###
  
  hess.1 = function(l,p,g,s,q,zu, t1, t0){
    z = zu
    cg = matrix(c(-0.1e1 / l ^ 2,0,0,0,2 / g ^ 3 * log1p(g * q / s) - 2 / g ^ 2 * q / s / (1 + g * q / s) - 1 / g * q ^ 2 / s ^ 2 / (1 + g * q / s) ^ 2 - 2 / g ^ 3 * log1p(g * z / s) + 2 / g ^ 2 * z / s / (1 + g * z / s) - (-1 / g - 1) * z ^ 2 / s ^ 2 / (1 + g * z / s) ^ 2,q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 - 1 / g * z / s ^ 2 / (1 + g * z / s) - (-1 / g - 1) * z / s ^ 2 / (1 + g * z / s) + (-1 / g - 1) * z ^ 2 / s ^ 3 / (1 + g * z / s) ^ 2 * g,0,q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 - 1 / g * z / s ^ 2 / (1 + g * z / s) - (-1 / g - 1) * z / s ^ 2 / (1 + g * z / s) + (-1 / g - 1) * z ^ 2 / s ^ 3 / (1 + g * z / s) ^ 2 * g,2 * q / s ^ 3 / (1 + g * q / s) - q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 * g + (2 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) + 4 / s ^ 4 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s) + 1 / s ^ 5 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) ^ 2 * g ^ 2 * z ^ 2 / (1 + g * z / s) ^ 2 - 1 / s ^ 5 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g ^ 2 * z ^ 2 / (1 + g * z / s) ^ 2) * s / (1 + g * z / s) ^ (-1 / g - 1) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) / (1 + g * z / s) ^ (-1 / g - 1) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) / s / (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)),nrow=3,ncol=3)
    return(cg)
  }
  
  hess.2 = function(l,p,g,s,q,zu, t1, t0){
    z = zu
    cg1 = matrix(c(-0.1e1 / l ^ 2,0,0,0,-p * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) - p * (1 + g * q / s) ^ (1 / g) * (2 / g ^ 3 * log1p(g * q / s) - 2 / g ^ 2 * q / s / (1 + g * q / s) - 1 / g * q ^ 2 / s ^ 2 / (1 + g * q / s) ^ 2) / (1 - p * (1 + g * q / s) ^ (1 / g)) - p ^ 2 * ((1 + g * q / s) ^ (1 / g)) ^ 2 * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) ^ 2,p * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) / (1 - p * (1 + g * q / s) ^ (1 / g)) - p * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) + p ^ 2 * ((1 + g * q / s) ^ (1 / g)) ^ 2 * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) / (1 - p * (1 + g * q / s) ^ (1 / g)) ^ 2 * q / s ^ 2 / (1 + g * q / s),0,p * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) / (1 - p * (1 + g * q / s) ^ (1 / g)) - p * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) + p ^ 2 * ((1 + g * q / s) ^ (1 / g)) ^ 2 * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) / (1 - p * (1 + g * q / s) ^ (1 / g)) ^ 2 * q / s ^ 2 / (1 + g * q / s),-p * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) - 2 * p * (1 + g * q / s) ^ (1 / g) * q / s ^ 3 / (1 + g * q / s) / (1 - p * (1 + g * q / s) ^ (1 / g)) + p * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) * g - p ^ 2 * ((1 + g * q / s) ^ (1 / g)) ^ 2 * q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 / (1 - p * (1 + g * q / s) ^ (1 / g)) ^ 2),nrow=3,ncol=3)
    return(cg1)
  }
  
  likdiff.1 = function(l,p,g,s,q,zu, t1, t0){
    z = zu
    diff.l = -t1 + t0 + 1/l
    diff.g = -1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s) + 1 / g ^ 2 * log1p(g * z / s) + (-1 / g - 1) * z / s / (1 + g * z / s)
    diff.s = -q / s ^ 2 / (1 + g * q / s) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) * s / (1 + g * z / s) ^ (-1 / g - 1)
    return(matrix(c(diff.l, diff.g, diff.s), nrow = 3, ncol = 1))
  }
  
  likdiff.2 = function(l,p,g,s,q,zu,t1,t0){
    z = zu
    diff.l = -t1 + t0 + 0.1e1 / l
    diff.g = p * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) / (1 - p * (1 + g * q / s) ^ (1 / g))
    diff.s = p * (1 + g * q / s) ^ (0.1e1 / g) * q / s ^ 2 / (1 + g * q / s) / (1 - p * (1 + g * q / s) ^ (0.1e1 / g))
    return(matrix(c(diff.l, diff.g, diff.s), nrow = 3, ncol = 1))
  }
  
  
  obs.mat = function(theta.dk, theta.ne){
    l.dk = theta.dk[1]; p.dk = theta.dk[2]; g.dk = theta.dk[3]; s.dk = theta.dk[4]
    l.ne = theta.ne[1]; p.ne = theta.ne[2]; g.ne = theta.ne[3]; s.ne = theta.ne[4]
    mat = list()
    for(i in 1:max(N1,N2)){
      mat.in = matrix(rep(0,36), nrow = 6, ncol = 6)
      if(i == 1){
        t0.ne = 0
        t0.dk = 0
      } else {
        t0.ne = poisson.time.ne[i-1]
        t0.dk = poisson.time.dk[i-1]
      }
      
      if(i <= min(N1, N2)){
        if(n.ne[i] == 1){
          ne = hess.1(l = l.ne, p = p.ne, g=g.ne, s = s.ne, q = q2, zu = z2[i], 
                      t1 = poisson.time.ne[i], t0 = t0.ne)
        } else{
          ne = hess.2(l = l.ne, p = p.ne, g=g.ne, s = s.ne, q = q2, zu = z2[i], 
                      t1 = poisson.time.ne[i], t0 = t0.ne)
        }
        
        if(n.dk[i]==1){
          dk = hess.1(l = l.dk, p = p.dk, g = g.dk, s=s.dk, q = q1, zu = z1[i],
                      t1 = poisson.time.dk[i], t0 = t0.dk)
        } else {
          dk = hess.2(l = l.dk, p = p.dk, g = g.dk, s=s.dk, q = q1, zu = z1[i],
                      t1 = poisson.time.dk[i], t0 = t0.dk)
        }
        mat.in[1:3, 1:3] = ne
        mat.in[4:6, 4:6] = dk
        mat[[i]] = mat.in
      } else if(N1>N2) {
        if(n.dk[i]==1){
          dk = hess.1(l = l.dk, p = p.dk, g = g.dk, s=s.dk, q = q1, zu = z1[i],
                      t1 = poisson.time.dk[i], t0 = t0.dk)
          ne = matrix(rep(0,9), nrow = 3, ncol =3)
          
        } else {
          dk = hess.2(l = l.dk, p = p.dk, g = g.dk, s=s.dk, q = q1, zu = z1[i],
                      t1 = poisson.time.dk[i], t0 = t0.dk)
          ne = matrix(rep(0,9), nrow = 3, ncol =3)
        }
        mat.in[1:3, 1:3] = ne
        mat.in[4:6, 4:6] = dk
        mat[[i]] = mat.in
      } else {
        if(n.ne[i] == 1){
          ne = hess.1(l = l.ne, p = p.ne, g=g.ne, s = s.ne, q = q2, zu = z2[i], 
                      t1 = poisson.time.ne[i], t0 = t0.ne)
          dk = matrix(rep(0,9), nrow = 3, ncol = 3)
        } else {
          ne = hess.2(l = l.ne, p = p.ne, g=g.ne, s = s.ne, q = q2, zu = z2[i], 
                      t1 = poisson.time.ne[i], t0 = t0.ne)
          dk = matrix(rep(0,9), nrow = 3, ncol = 3)
        }
        mat.in[1:3, 1:3] = ne
        mat.in[4:6, 4:6] = dk
        mat[[i]] = mat.in
      }
    }
    o.mat = Reduce('+', mat)
    return(o.mat)
  }
  
  lik.mat = function(theta.dk, theta.ne, theta.dk.ml, theta.ne.ml){
    l.dk = theta.dk[1]; p.dk = theta.dk[2]; g.dk = theta.dk[3]; s.dk = theta.dk[4]
    l.ne = theta.ne[1]; p.ne= theta.ne[2]; g.ne = theta.ne[3]; s.ne = theta.ne[4]
    
    l.dk.ml = theta.dk.ml[1]; p.dk.ml = theta.dk.ml[2]
    g.dk.ml = theta.dk.ml[3]; s.dk.ml = theta.dk.ml[4]
    
    l.ne.ml = theta.ne.ml[1]; p.ne.ml = theta.ne.ml[2]
    g.ne.ml = theta.ne.ml[3]; s.ne.ml = theta.ne.ml[4]
    
    mat = list()
    
    for(i in 1:max(N1, N2)){
      if(i == 1){
        t0.ne = 0
        t0.dk = 0
      } else{
        t0.ne = poisson.time.ne[i-1]
        t0.dk = poisson.time.dk[i-1]
      }
      
      if(i <= min(N1,N2)){
        if(n.dk[i] == 1){
          lik.pl.dk = likdiff.1(l = l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
          lik.ml.dk = likdiff.1(l = l.dk.ml, p =  p.dk.ml, g=g.dk.ml, 
                                s=s.dk.ml, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
        } else {
          lik.pl.dk = likdiff.2(l = l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
          lik.ml.dk = likdiff.2(l = l.dk.ml, p = p.dk.ml, g = g.dk.ml, s = s.dk.ml,
                                q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
        }
        if(n.ne[i]==1){
          lik.pl.ne = likdiff.1(l = l.ne, p= p.ne, g = g.ne, s = s.ne, q = q2, 
                                zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
          lik.ml.ne = likdiff.1(l = l.ne.ml, p = p.ne.ml, g = g.ne.ml, s = s.ne.ml,
                                q = q2, zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
        } else{
          lik.pl.ne = likdiff.2(l = l.ne, p= p.ne, g = g.ne, s = s.ne, q = q2, 
                                zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
          lik.ml.ne = likdiff.2(l = l.ne.ml, p = p.ne.ml, g = g.ne.ml, s = s.ne.ml,
                                q = q2, zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
        }
        mat.pl = matrix(c(lik.pl.ne, lik.pl.dk),nrow = 6, ncol = 1) 
        mat.ml = matrix(c(lik.ml.ne, lik.ml.dk), nrow = 6, ncol = 1)
        mat[[i]] = mat.pl %*% t(mat.ml)
      } else if(N1>N2){
        if(n.dk[i] == 1){
          lik.pl.dk = likdiff.1(l = l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
          lik.ml.dk = likdiff.1(l = l.dk.ml, p = p.dk.ml, g = g.dk.ml, 
                                s = s.dk.ml, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
        } else {
          lik.pl.dk = likdiff.2(l = l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
          lik.ml.dk = likdiff.2(l = l.dk.ml, p = p.dk.ml, g = g.dk.ml,s = s.dk.ml,
                                q = q1, zu = z1[i],
                                t1 = poisson.time.dk[i], t0 = t0.dk)
        }
        mat.pl = matrix(c(0,0,0, lik.pl.dk), nrow = 6, ncol = 1)
        mat.ml = matrix(c(0,0,0, lik.ml.dk), nrow = 6, ncol = 1)
        mat[[i]] = mat.pl %*% t(mat.ml)
      } else {
        if(n.ne[i] == 1){
          lik.pl.ne = likdiff.1(l = l.ne, p= p.ne, g = g.ne, s = s.ne, q = q2, 
                                zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
          lik.ml.ne = likdiff.1(l = l.ne.ml, p = p.ne.ml, g = g.ne.ml, s = s.ne.ml,
                                q = q2, zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
        } else {
          lik.pl.ne = likdiff.2(l = l.ne, p= p.ne, g = g.ne, s = s.ne, q = q2, 
                                zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
          lik.ml.ne = likdiff.2(l = l.ne.ml, p = p.ne.ml, g = g.ne.ml, s = s.ne.ml,
                                q = q2, zu = z2[i], t1 = poisson.time.ne[i], t0 = t0.ne)
        }
        mat.pl = matrix(c(lik.pl.ne, 0,0,0), nrow = 6, ncol = 1)
        mat.ml = matrix(c(lik.ml.ne, 0,0,0), nrow = 6, ncol = 1)
        mat[[i]] = mat.pl %*% t(mat.ml)
      }
    }
    o.mat = Reduce('+', mat)
    return(o.mat)
  }
  
  ### Gradient for objective function in SLSQP algorithm ###
  
  grad.func=function(theta, pne, pdk){
    l.dk = theta[1]; p.dk = pdk; g.dk = theta[2]; s.dk = theta[3]
    l.ne = theta[4]; p.ne = pne; g.ne = theta[5]; s.ne = theta[6]
    
    vecs = list()
    for(i in 1:max(N1,N2)){
      if(i == 1){
        t0.dk = 0
        t0.ne = 0
      } else {
        t0.dk = poisson.time.dk[i-1]
        t0.ne = poisson.time.ne[i-1]
      }
      if(i <= min(N1,N2)){
        if(n.dk[i]==1){
          diff.dk = likdiff.1(l= l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                              t1 = poisson.time.dk[i], t0 = t0.dk)
        } else {
          diff.dk = likdiff.2(l= l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                              t1 = poisson.time.dk[i], t0 = t0.dk)
        }
        if(n.ne[i]==1){
          diff.ne = likdiff.1(l=l.ne, p=p.ne, g=g.ne, s=s.ne, q=q2, zu = z2[i],
                              t1 = poisson.time.ne[i], t0 = t0.ne)
        } else {
          diff.ne = likdiff.2(l=l.ne, p=p.ne, g=g.ne, s=s.ne, q=q2, zu = z2[i],
                              t1 = poisson.time.ne[i], t0 = t0.ne)
        }
      } else if(N1>N2){
        if(n.dk[i]==1){
          diff.dk = likdiff.1(l= l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                              t1 = poisson.time.dk[i], t0 = t0.dk)
        } else {
          diff.dk = likdiff.2(l= l.dk, p = p.dk, g = g.dk, s = s.dk, q = q1, zu = z1[i],
                              t1 = poisson.time.dk[i], t0 = t0.dk)
        }
        diff.ne = matrix(c(0,0,0), nrow = 3, ncol = 1)
      } else {
        if(n.ne[i]==1){
          diff.ne = likdiff.1(l=l.ne, p=p.ne, g=g.ne, s=s.ne, q=q2, zu = z2[i],
                              t1 = poisson.time.ne[i], t0 = t0.ne)
        } else {
          diff.ne = likdiff.2(l=l.ne, p=p.ne, g=g.ne, s=s.ne, q=q2, zu = z2[i],
                              t1 = poisson.time.ne[i], t0 = t0.ne)
        }
        diff.dk = matrix(c(0,0,0), nrow = 3, ncol = 1)
      }
      vecs[[i]]=matrix(c(diff.dk, diff.ne), nrow = 6, ncol = 1)
    }
    vec.tot = Reduce('+', vecs)
    return(c(-vec.tot[1], -vec.tot[2], -vec.tot[3], 
             -vec.tot[4], -vec.tot[5], -vec.tot[6]))
  }  
  
  
  ### Modified Profile Likelihood ###
  
  
  mpl.func.pi = function(param){
    pi.ne = param[1]
    pi.dk = param[2]
    
    start = c(ml.est.dk[c(1,3,4)], ml.est.ne[c(1,3,4)])
    
    if(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/ml.est.dk[4])^(1/ml.est.dk[3])>=1 | 
       is.nan(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/ml.est.dk[4])^(1/ml.est.dk[3]))){
      s = ml.est.dk[4]
      while(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/s)^(1/ml.est.dk[3])>=1 | 
            is.nan(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/s)^(1/ml.est.dk[3]))){
        s = s*1.5
      }
      start[3] = s
    } else {
      start = start
    }
    
    if(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/ml.est.ne[4])^(1/ml.est.ne[3])>=1 | 
       is.nan(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/ml.est.ne[4])^(1/ml.est.ne[3]))){
      s = ml.est.ne[4]
      while(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/s)^(1/ml.est.ne[3])>=1 | 
            is.nan(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/s)^(1/ml.est.ne[3]))){
        s = s*1.5
      }
      start[6] = s
    } else{
      start = start
    }
    
    likfunc = function(theta, pne, pdk){
      l.dk = theta[1]; p.dk = pdk; g.dk = theta[2]; s.dk = theta[3]
      l.ne = theta[4]; p.ne = pne; g.ne = theta[5]; s.ne = theta[6]
      
      a.dk = matrix(nrow = max(N1, N2), ncol=1); a.ne=matrix(nrow=max(N1,N2), ncol=1)
      b.dk = matrix(nrow = max(N1, N2), ncol=1); b.ne=matrix(nrow=max(N1,N2), ncol=1)
      c.dk = matrix(nrow = max(N1, N2), ncol=1); c.ne=matrix(nrow=max(N1,N2), ncol=1)
      
      for(i in 1:max(N1,N2)){
        if(i == 1){
          t0.dk = 0
          t0.ne = 0
        } else {
          t0.dk = poisson.time.dk[i-1]
          t0.ne = poisson.time.ne[i-1]
        }
        
        if(i <= min(N1,N2)){
          if(n.dk[i]==1){
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = log((1/s.dk)*(1+g.dk*z1[i]/s.dk)^((-1/g.dk)-1))
          } else{
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(1-p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = 0
          }
          if(n.ne[i]==1){
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = log((1/s.ne)*(1+g.ne*z2[i]/s.ne)^((-1/g.ne)-1))
          } else {
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(1-p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = 0
          }
        } else if(N1>N2){
          if(n.dk[i]==1){
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = log((1/s.dk)*(1+g.dk*z1[i]/s.dk)^((-1/g.dk)-1))
          } else{
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(1-p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = 0
          }
          a.ne[i] = 0
          b.ne[i] = 0
          c.ne[i] = 0
        } else {
          if(n.ne[i]==1){
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = log((1/s.ne)*(1+g.ne*z2[i]/s.ne)^((-1/g.ne)-1))
          } else{
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(1-p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = 0
          }
          a.dk[i] = 0
          b.dk[i] = 0
          c.dk[i] = 0
        }
      }
      lik = sum(a.ne)+sum(b.ne)+sum(c.ne)+sum(a.dk)+sum(b.dk)+sum(c.dk)
      return(-lik)
    }
    ## Setup functions for objective and constraints + gradients ##
    obj_func = function(x, p_ne, p_dk){
      val = likfunc(theta = x, pne = p_ne, pdk = p_dk)
      return(val)
    }
    
    obj_func_grad = function(x, p_ne, p_dk){
      a = grad.func(theta = x, pne = p_ne, pdk = p_dk)
      return(a)
    }
    
    constraint.func = function(x, p_ne, p_dk){
      l.dk = x[1]; p.dk = p_dk; g.dk = x[2]; s.dk = x[3]
      l.ne = x[4]; p.ne = p_ne; g.ne = x[5]; s.ne = x[6]
      
      a.dk.1 = -(1+g.dk*max(max(z1, na.rm = TRUE), q1)/s.dk) +10^(-5)
      a.dk.2 = p.dk*(1+g.dk*max(max(z1, na.rm = TRUE), q1)/s.dk)^(1/g.dk)-(1-10^(-6))
      
      a.ne.1 = -(1+g.ne*max(max(z2, na.rm = TRUE), q2)/s.ne) +10^(-5)
      a.ne.2 = p.ne*(1+g.ne*max(max(z2, na.rm = TRUE), q2)/s.ne)^(1/g.ne)-(1-10^(-6))
      
      condition = c(a.dk.1, a.dk.2, a.ne.1, a.ne.2)
      
      return(condition)
    }
    
    grad.constraint = function(x, p_ne, p_dk){
      l.dk = x[1]; p.dk = p_dk; g.dk = x[2]; s.dk = x[3]
      l.ne = x[4]; p.ne = p_ne; g.ne = x[5]; s.ne = x[6]
      
      condition.dk.1= c(0, -max(max(z1, na.rm = TRUE), q1)/s.dk, 
                        g.dk*max(max(z1, na.rm = TRUE), q1)/(s.dk^2),0,0,0)
      condition.dk.2= c(0, p.dk * (1 + g.dk * max(max(z1, na.rm = TRUE),q1) / s.dk) ^ (1 / g.dk) * (-1 / g.dk ^ 2 * log1p(g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk) + 1 / g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk / (1 + g.dk * max(max(z1, na.rm = TRUE), q1) / s.dk)),
                        -p.dk * (1 + g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk) ^ (0.1e1 / g.dk) *max(max(z1, na.rm = TRUE), q1) / s.dk ^ 2 / (1 + g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk), 0,0,0)
      
      
      condition.ne.1 = c(0,0,0,0, -max(max(z2, na.rm = TRUE),q2)/s.ne, 
                         g.ne*max(max(z2, na.rm = TRUE),q2)/(s.ne^2))
      condition.ne.2 = c(0,0,0,0, p.ne * (1 + g.ne * max(max(z2, na.rm = TRUE),q2) / s.ne) ^ (1 / g.ne) * (-1 / g.ne ^ 2 * log1p(g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne) + 1 / g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne / (1 + g.ne * max(max(z2, na.rm = TRUE), q2) / s.ne)),
                         -p.ne * (1 + g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne) ^ (0.1e1 / g.ne) *max(max(z2, na.rm = TRUE), q2) / s.ne ^ 2 / (1 + g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne))
      
      condition.final = rbind(condition.dk.1, condition.dk.2, condition.ne.1, condition.ne.2)
      return(condition.final)
    }
    
    opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-6, maxeval = 100)
    
    o = nloptr(x0 = start, p_ne = pi.ne, p_dk = pi.dk,
               eval_f = obj_func, eval_g_ineq = constraint.func, eval_grad_f = obj_func_grad,
               eval_jac_g_ineq = grad.constraint, opts = opts1,
               lb = c(10^(-10), -Inf, 10^(-10), 10^(-10), -Inf, 10^(-10)),
               ub = c(Inf, Inf, Inf, Inf, Inf, Inf)
    )
    theta.pl.dk = c(o$solution[1], pi.dk, o$solution[2], o$solution[3])
    theta.pl.ne = c(o$solution[4], pi.ne, o$solution[5], o$solution[6])
    
    likmat = lik.mat(theta.dk = theta.pl.dk, theta.ne = theta.pl.ne, 
                     theta.dk.ml = c(ml.est.dk[1], ml.pi.dk, ml.est.dk[3], ml.est.dk[4]),
                     theta.ne.ml = c(ml.est.ne[1], ml.pi.ne, ml.est.ne[3], ml.est.ne[4]))
    obsmat = obs.mat(theta.dk = theta.pl.dk, theta.ne = theta.pl.ne)
    
    mpl.lik = (-o$objective+(1/2)*log(abs(det(-obsmat)))
               -log(abs(det(likmat))))
    return(-mpl.lik)
  }
  
  local.opts = list("algorithm" = "NLOPT_GN_ISRES")
  opts1 = list("algorithm" = "NLOPT_LN_AUGLAG", "maxtime" = 120, "local_opts" = local.opts,
               "xtol_rel" = 10^(-15))
  
  optr = nloptr(x0 = c(ml.pi.ne, ml.pi.dk), eval_f = mpl.func.pi, 
                lb = c(0.7*ml.pi.ne, 0.7*ml.pi.dk), 
                ub = c(2*ml.pi.ne, 2*ml.pi.dk),
                opts = opts1)
  
  mpl.lik = optr$objective
  mpl.pars = optr$solution
  
  prof.lik.func = function(param){
    pi.ne = param[1]
    pi.dk = param[2]
    
    start = c(ml.est.dk[c(1,3,4)], ml.est.ne[c(1,3,4)])
    
    if(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/ml.est.dk[4])^(1/ml.est.dk[3])>=1 | 
       is.nan(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/ml.est.dk[4])^(1/ml.est.dk[3]))){
      s = ml.est.dk[4]
      while(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/s)^(1/ml.est.dk[3])>=1 | 
            is.nan(pi.dk*(1+ml.est.dk[3]*max(max(z1, na.rm = TRUE), q1)/s)^(1/ml.est.dk[3]))){
        s = s*1.5
      }
      start[3] = s
    } else {
      start = start
    }
    
    if(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/ml.est.ne[4])^(1/ml.est.ne[3])>=1 | 
       is.nan(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/ml.est.ne[4])^(1/ml.est.ne[3]))){
      s = ml.est.ne[4]
      while(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/s)^(1/ml.est.ne[3])>=1 | 
            is.nan(pi.ne*(1+ml.est.ne[3]*max(max(z2, na.rm = TRUE), q2)/s)^(1/ml.est.ne[3]))){
        s = s*1.5
      }
      start[6] = s
    } else{
      start = start
    }
    
    likfunc = function(theta, pne, pdk){
      l.dk = theta[1]; p.dk = pdk; g.dk = theta[2]; s.dk = theta[3]
      l.ne = theta[4]; p.ne = pne; g.ne = theta[5]; s.ne = theta[6]
      
      a.dk = matrix(nrow = max(N1, N2), ncol=1); a.ne=matrix(nrow=max(N1,N2), ncol=1)
      b.dk = matrix(nrow = max(N1, N2), ncol=1); b.ne=matrix(nrow=max(N1,N2), ncol=1)
      c.dk = matrix(nrow = max(N1, N2), ncol=1); c.ne=matrix(nrow=max(N1,N2), ncol=1)
      
      for(i in 1:max(N1,N2)){
        if(i == 1){
          t0.dk = 0
          t0.ne = 0
        } else {
          t0.dk = poisson.time.dk[i-1]
          t0.ne = poisson.time.ne[i-1]
        }
        
        if(i <= min(N1,N2)){
          if(n.dk[i]==1){
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = log((1/s.dk)*(1+g.dk*z1[i]/s.dk)^((-1/g.dk)-1))
          } else{
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(1-p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = 0
          }
          if(n.ne[i]==1){
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = log((1/s.ne)*(1+g.ne*z2[i]/s.ne)^((-1/g.ne)-1))
          } else {
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(1-p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = 0
          }
        } else if(N1>N2){
          if(n.dk[i]==1){
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = log((1/s.dk)*(1+g.dk*z1[i]/s.dk)^((-1/g.dk)-1))
          } else{
            a.dk[i] = -l.dk*(poisson.time.dk[i]-t0.dk)+log(l.dk*(poisson.time.dk[i]-t0.dk))
            b.dk[i] = log(1-p.dk*(1+g.dk*q1/s.dk)^(1/g.dk))
            c.dk[i] = 0
          }
          a.ne[i] = 0
          b.ne[i] = 0
          c.ne[i] = 0
        } else {
          if(n.ne[i]==1){
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = log((1/s.ne)*(1+g.ne*z2[i]/s.ne)^((-1/g.ne)-1))
          } else{
            a.ne[i] = -l.ne*(poisson.time.ne[i]-t0.ne)+log(l.ne*(poisson.time.ne[i]-t0.ne))
            b.ne[i] = log(1-p.ne*(1+g.ne*q2/s.ne)^(1/g.ne))
            c.ne[i] = 0
          }
          a.dk[i] = 0
          b.dk[i] = 0
          c.dk[i] = 0
        }
      }
      lik = sum(a.ne)+sum(b.ne)+sum(c.ne)+sum(a.dk)+sum(b.dk)+sum(c.dk)
      return(-lik)
    }
    
    obj_func = function(x, p_ne, p_dk){
      val = likfunc(theta = x, pne = p_ne, pdk = p_dk)
      return(val)
    }
    
    obj_func_grad = function(x, p_ne, p_dk){
      a = grad.func(theta = x, pne = p_ne, pdk = p_dk)
      return(a)
    }
    
    constraint.func = function(x, p_ne, p_dk){
      l.dk = x[1]; p.dk = p_dk; g.dk = x[2]; s.dk = x[3]
      l.ne = x[4]; p.ne = p_ne; g.ne = x[5]; s.ne = x[6]
      
      a.dk.1 = -(1+g.dk*max(max(z1, na.rm = TRUE), q1)/s.dk) +10^(-10)
      a.dk.2 = p.dk*(1+g.dk*max(max(z1, na.rm = TRUE), q1)/s.dk)^(1/g.dk)-(1-10^(-6))
      
      a.ne.1 = -(1+g.ne*max(max(z2, na.rm = TRUE), q2)/s.ne) +10^(-10)
      a.ne.2 = p.ne*(1+g.ne*max(max(z2, na.rm = TRUE), q2)/s.ne)^(1/g.ne)-(1-10^(-6))
      
      condition = c(a.dk.1, a.dk.2, a.ne.1, a.ne.2)
      
      return(condition)
    }
    
    grad.constraint = function(x, p_ne, p_dk){
      l.dk = x[1]; p.dk = p_dk; g.dk = x[2]; s.dk = x[3]
      l.ne = x[4]; p.ne = p_ne; g.ne = x[5]; s.ne = x[6]
      
      condition.dk.1= c(0, -max(max(z1, na.rm = TRUE), q1)/s.dk, 
                        g.dk*max(max(z1, na.rm = TRUE), q1)/(s.dk^2),0,0,0)
      condition.dk.2= c(0, p.dk * (1 + g.dk * max(max(z1, na.rm = TRUE),q1) / s.dk) ^ (1 / g.dk) * (-1 / g.dk ^ 2 * log1p(g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk) + 1 / g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk / (1 + g.dk * max(max(z1, na.rm = TRUE), q1) / s.dk)),
                        -p.dk * (1 + g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk) ^ (0.1e1 / g.dk) *max(max(z1, na.rm = TRUE), q1) / s.dk ^ 2 / (1 + g.dk *max(max(z1, na.rm = TRUE), q1) / s.dk), 0,0,0)
      
      
      condition.ne.1 = c(0,0,0,0, -max(max(z2, na.rm = TRUE),q2)/s.ne, 
                         g.ne*max(max(z2, na.rm = TRUE),q2)/(s.ne^2))
      condition.ne.2 = c(0,0,0,0, p.ne * (1 + g.ne * max(max(z2, na.rm = TRUE),q2) / s.ne) ^ (1 / g.ne) * (-1 / g.ne ^ 2 * log1p(g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne) + 1 / g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne / (1 + g.ne * max(max(z2, na.rm = TRUE), q2) / s.ne)),
                         -p.ne * (1 + g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne) ^ (0.1e1 / g.ne) *max(max(z2, na.rm = TRUE), q2) / s.ne ^ 2 / (1 + g.ne *max(max(z2, na.rm = TRUE), q2) / s.ne))
      
      condition.final = rbind(condition.dk.1, condition.dk.2, condition.ne.1, condition.ne.2)
      return(condition.final)
    }
    opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-6, maxeval = 1000)
    o = nloptr(x0 = start, p_ne = pi.ne, p_dk = pi.dk,
               eval_f = obj_func, eval_g_ineq = constraint.func, eval_grad_f = obj_func_grad,
               eval_jac_g_ineq = grad.constraint, opts = opts1,
               lb = c(10^(-10), -Inf, 10^(-10), 10^(-10), -Inf, 10^(-10)),
               ub = c(Inf, Inf, Inf, Inf, Inf, Inf)
    )
    return(o$objective)
  }
  
  bound.func = function(param1){
    val = 2*(-mpl.lik+prof.lik.func(param1)) - qchisq(p = 0.95, df = 1)
    return(val)
  }
  
  obj.func1 = function(param1){
    pi.ne = param1[1]
    pi.dk = param1[2]
    
    return(pi.ne-pi.dk)
  }
  
  obj.func2 = function(param1){
    pi.ne = param1[1]
    pi.dk = param1[2]
    
    return(pi.dk-pi.ne)
  }
  
  local.opts = list("algorithm" = "NLOPT_GN_ISRES")
  opts1 = list("algorithm" = "NLOPT_LN_AUGLAG", "maxtime" = 120, "local_opts" = local.opts,
               "xtol_rel" = 10^(-10))
  lower.opt = nloptr(x0 = mpl.pars, eval_f = obj.func1, eval_g_ineq = bound.func,
                     lb = 0.2*mpl.pars, ub = 3*mpl.pars, opts = opts1)
  
  upper.opt = nloptr(x0 = mpl.pars, eval_f = obj.func2, eval_g_ineq = bound.func,
                     lb = 0.2*mpl.pars, ub = 3*mpl.pars, opts = opts1)
  
  ests = c(lower.opt$objective, mpl.pars[1]-mpl.pars[2], -upper.opt$objective)
  return(ests)
}





