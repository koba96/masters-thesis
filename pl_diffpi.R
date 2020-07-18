#install.packages('plotrix')

library('nloptr')
library('plotrix')
library('extRemes')
library('tidyverse')

##### Functions for computing maximum likelihood estimates 
##### to initialize profile likelihood

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

##### Function to compute confidence intervals for the difference pi_ne - pi_dk using the profile likelihood #####

pl_diffpi = function(dataDK, dataNE, uDK, uNE, q, measure){
  
  measure1 = measure
  thresh1 = uDK
  thresh2 = uNE
  x1 = q
  
  N1 = length(dataDK[[measure1]])
  data1 = dataDK[dataDK[[measure1]] != -1, ]
  data.s = data1[data1[[measure1]]<thresh1, ]
  z1 = -data.s[[measure1]]+thresh1
  q1 = -x1+thresh1
  y1 = length(z1)
  
  poisson_time_dk = as.POSIXct(dataDK[["Date"]], tz = "", 
                               tryFormats = "%Y-%m-%d %H:%M:%OS")
  b00.dk = as.POSIXct(paste(str_sub(dataDK[['Date']][1], 1, 10), "00:00:01"),
                      tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  poisson.time.dk = as.numeric(difftime(poisson_time_dk, b00.dk, units = "hours"))
  t.final.dk = as.numeric(difftime(
    as.POSIXct(paste(str_sub(dataDK[['Date']][length(dataDK[,1])],1,10), "23:59:59"),
               tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS"), 
    b00.dk, units = "hours"))
  
  N2 = length(dataNE[[measure1]])
  data2 = dataNE[dataNE[[measure1]] != -1, ]
  data.s = data2[data2[[measure1]]<thresh2, ]
  z2 = -data.s[[measure1]]+thresh2
  q2 = -x1+thresh2
  y2 = length(z2)
  
  poisson_time_ne = as.POSIXct(dataNE[["Date"]], tz = "", 
                               tryFormats = "%Y-%m-%d %H:%M:%OS")
  b00.ne = as.POSIXct(paste(str_sub(dataNE[['Date']][1], 1, 10), "00:00:01"),
                      tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  poisson.time.ne = as.numeric(difftime(poisson_time_ne, b00.ne, units = "hours"))
  
  t.final.ne = as.numeric(difftime(
    as.POSIXct(paste(str_sub(dataNE[['Date']][length(dataNE[,1])],1,10), "23:59:59"),
               tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS"), 
    b00.ne, units = "hours"))
  
  ml.est.dk.all = pois.est(data = dataDK,  
                           thresh = uDK, measure = measure1, x = q1)$par
  ml.est.ne.all = pois.est(data = dataNE,  
                           thresh = uNE, measure = measure1, x = q2)$par
  
  ml.est.dk = ml.est.dk.all[c(1,3,4)]
  ml.est.ne = ml.est.ne.all[c(1,3,4)]
  
  proflik = function(pis, param = FALSE){
    pi.dk = pis[1]
    pi.ne = pis[2]
    start = c(ml.est.ne, ml.est.dk)
    
    if(start[2]< -0.5){
      start[2] = -0.5
    }
    
    if(start[5]< -0.5){
      start[5] = -0.5
    }
    
    if(pi.dk*(1+ml.est.dk[2]*max(max(z1), q1)/ml.est.dk[3])^(1/ml.est.dk[2])>=1 | 
       is.nan(pi.dk*(1+ml.est.dk[2]*max(max(z1), q1)/ml.est.dk[3])^(1/ml.est.dk[2]))){
      s = ml.est.dk[3]
      while(pi.dk*(1+ml.est.dk[2]*max(max(z1), q1)/s)^(1/ml.est.dk[2])>=1 | 
            is.nan(pi.dk*(1+ml.est.dk[2]*max(max(z1), q1)/s)^(1/ml.est.dk[2]))){
        s = s*1.5
      }
      start[6] = s
    } else {
      start = start
    }
    
    if(pi.ne*(1+ml.est.ne[2]*max(max(z2), q2)/ml.est.ne[3])^(1/ml.est.ne[2])>=1 | 
       is.nan(pi.ne*(1+ml.est.ne[2]*max(max(z2), q2)/ml.est.ne[3])^(1/ml.est.ne[2]))){
      s = ml.est.ne[3]
      while(pi.ne*(1+ml.est.ne[2]*max(max(z2), q2)/s)^(1/ml.est.ne[2])>=1 | 
            is.nan(pi.ne*(1+ml.est.ne[2]*max(max(z2), q2)/s)^(1/ml.est.ne[2]))){
        s = s*1.5
      }
      start[3] = s
    } else{
      start = start
    }
    
    
    likfunc = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      check1 = pi.ne*(1+g2*max(max(z2), q2)/s2)^(1/g2)
      check2 = pi.dk*(1+g1*max(max(z1), q1)/s1)^(1/g1)
      
      a.1 = -l1*t.final.dk + N1*log(l1*t.final.dk)
      b.1 = y1*(log(pi.dk)+(1/g1)*log(1+g1*q1/s1))
      c.1 = (N1-y1)*log(1-pi.dk*(1+g1*q1/s1)^(1/g1))
      d.1 = sum(log(devd(z1, scale = s1, shape = g1, threshold = 0, type = "GP")))
      
      a.2 = -l2*t.final.ne + N2*log(l2*t.final.ne)
      b.2 = y2*(log(pi.ne)+(1/g2)*log(1+g2*q2/s2))
      c.2 = (N2-y2)*log(1-pi.ne*(1+g2*q2/s2)^(1/g2))
      d.2 = sum(log(devd(z2, scale = s2, shape = g2, threshold = 0, type = "GP")))
      
      lik1 = a.1+b.1+c.1+d.1
      lik2 = a.2+b.2+c.2+d.2
      lik = lik1 + lik2
      return(-lik)
    }
    
    grad.func = function(par){
      n1 = y1
      n2 = y2
      p1 = pi.dk
      p2 = pi.ne
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      diff.l.dk = -t.final.dk + N1/l1
      diff.g.dk.1 = n1 * (-1 / g1 ^ 2 * log1p(g1 * q1 / s1) + 1 / g1 * q1 / s1 / (1 + g1 * q1 / s1)) - (N1 - n1) * p1 * (1 + g1 * q1 / s1) ^ (1 / g1) * (-1 / g1 ^ 2 * log1p(g1 * q1 / s1) + 1 / g1 * q1 / s1 / (1 + g1 * q1 / s1)) / (1 - p1 * (1 + g1 * q1 / s1) ^ (1 / g1))
      diff.g.dk.2 = sum(1 / g1 ^ 2 * log1p(1 / s1 * g1 * z1) + (-1 / g1 - 1) / s1 * z1 / (1 + 1 / s1 * g1 * z1))
      diff.s.dk.1 = -n1 * q1 / s1 ^ 2 / (1 + g1 * q1 / s1) + (N1 - n1) * p1 * (1 + g1 * q1 / s1) ^ (1 / g1) * q1 / s1 ^ 2 / (1 + g1 * q1 / s1) / (1 - p1 * (1 + g1 * q1 / s1) ^ (1 / g1))
      diff.s.dk.2 = sum((-1 / s1 ^ 2 * (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1) - 1 / s1 ^ 3 * (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1) * (-1 / g1 - 1) * g1 * z1 / (1 + 1 / s1 * g1 * z1)) * s1 / (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1))
      
      diff.l.ne = -t.final.ne + N2/l2
      diff.g.ne.1 = n2 * (-1 / g2 ^ 2 * log1p(g2 * q2 / s2) + 1 / g2 * q2 / s2 / (1 + g2 * q2 / s2)) - (N2 - n2) * p2 * (1 + g2 * q2 / s2) ^ (1 / g2) * (-1 / g2 ^ 2 * log1p(g2 * q2 / s2) + 1 / g2 * q2 / s2 / (1 + g2 * q2 / s2)) / (1 - p2 * (1 + g2 * q2 / s2) ^ (1 / g2))
      diff.g.ne.2 = sum(1 / g2 ^ 2 * log1p(g2 * z2 / s2) + (-1 / g2 - 1) * z2 / s2 / (1 + g2 * z2 / s2))
      diff.s.ne.1 = -n2 * q2 / s2 ^ 2 / (1 + g2 * q2 / s2) + (N2 - n2) * p2 * (1 + g2 * q2 / s2) ^ (1 / g2) * q2 / s2 ^ 2 / (1 + g2 * q2 / s2) / (1 - p2 * (1 + g2 * q2 / s2) ^ (1 / g2))
      diff.s.ne.2 = sum((-1 / s2 ^ 2 * (1 + g2 * z2 / s2) ^ (-1 / g2 - 1) - 1 / s2 ^ 3 * (1 + g2 * z2 / s2) ^ (-1 / g2 - 1) * (-1 / g2 - 1) * g2 * z2 / (1 + g2 * z2 / s2)) * s2 / (1 + g2 * z2 / s2) ^ (-1 / g2 - 1))
      
      vec = c(diff.l.ne, diff.g.ne.1+diff.g.ne.2, diff.s.ne.1+diff.s.ne.2,
              diff.l.dk, diff.g.dk.1+diff.g.dk.2, diff.s.dk.1+diff.s.dk.2)
      
      return(-vec)
    }
    
    constraint.func = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      ne.1 = -(1+g2*max(max(z2), q2)/s2)+10^(-4)
      ne.2 = pi.ne*(1+g2*max(max(z2), q2)/s2)^(1/g2) - 1
      
      dk.1 = -(1+g1*max(max(z1), q1)/s1)+10^(-4)
      dk.2 = pi.dk*(1+g1*max(max(z1), q1)/s1)^(1/g1) - 1
      return(c(ne.1, ne.2, dk.1, dk.2))
    }
    
    grad.constraint = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      a.ne.1 = c(0, -max(max(z2), q2)/s2, g2*max(max(z2), q2)/(s2^2), 0, 0 ,0)
      a.ne.2 = c(0, pi.ne * (1 + g2 * max(max(z2), q2) / s2) ^ (1 / g2) * (-1 / g2 ^ 2 * log1p(g2 * max(max(z2), q2) / s2) + 1 / g2 * max(max(z2), q2) / s2 / (1 + g2 * max(max(z2), q2) / s2)), -pi.ne * (1 + g2 * max(max(z2), q2) / s2) ^ (0.1e1 / g2) * max(max(z2), q2) / s2 ^ 2 / (1 + g2 * max(max(z2), q2) / s2), 0,0,0)
      
      a.dk.1 = c(0, 0, 0, 0, -max(max(z1), q1)/s1, g1*max(max(z1), q1)/(s1^2))
      a.dk.2 = c(0, 0, 0, 0, pi.dk * (1 + g1 * max(max(z1), q1) / s1) ^ (1 / g1) * (-1 / g1 ^ 2 * log1p(g1 * max(max(z1), q1) / s1) + 1 / g1 * max(max(z1), q1) / s1 / (1 + g1 * max(max(z1), q1) / s1)), -pi.dk * (1 + g1 * max(max(z1), q1) / s1) ^ (0.1e1 / g1) * max(max(z1), q1) / s1 ^ 2 / (1 + g1 * max(max(z1), q1) / s1))
      return(rbind(a.ne.1, a.ne.2, a.dk.1, a.dk.2))
    }
    local_opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-10, "maxeval" = 1000)
    opts1 = list("algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-10, "maxeval" = 1000,
                 "local_opts" = local_opts1)
    optr = nloptr(x0 = start, eval_f = likfunc, eval_grad_f = grad.func,
                  eval_g_ineq = constraint.func, eval_jac_g_ineq = grad.constraint,
                  opts = opts1, lb = c(10^(-20), -0.5, 10^(-20), 
                                       10^(-20), -0.5, 10^(-20)),
                  ub = c(Inf, Inf, Inf, Inf, Inf, Inf))
    p.lik = optr$objective
    if(param == FALSE){
      return(p.lik)
    } else {
      return(optr$solution)
    }
  }
  
  pi.ml.dk = ml.est.dk.all[2]*(1+ml.est.dk[2]*q1/ml.est.dk[3])^(-1/ml.est.dk[2])
  pi.ml.ne = ml.est.ne.all[2]*(1+ml.est.ne[2]*q2/ml.est.ne[3])^(-1/ml.est.ne[2])
  
  prof.est = auglag(x0 = c(pi.ml.dk, pi.ml.ne), fn = proflik, 
                    lower = c(10^(-20), 10^(-20)), 
                    upper = c(1, 1),
                    localsolver = "Cobyla", localtol = 10^(-10))  
  pi.prof = prof.est$par
  prof.lik = -prof.est$value
  
  par.prof = proflik(pis = pi.prof, param = TRUE)
  
  est = pi.prof[2]-pi.prof[1]
  
  obj.func1 = function(pis){
    pi.dk = pis[1]
    pi.ne = pis[2]
    
    val = pi.ne - pi.dk
    return(val)
  }
  
  obj.func2 = function(pis){
    pi.dk = pis[1]
    pi.ne = pis[2]
    
    val = pi.ne - pi.dk
    return(-val)
  }
  
  proflik2 = function(pis, param = FALSE){
    pi.dk = pis[1]
    pi.ne = pis[2]
    start = par.prof
    
    if(start[2]< -0.5){
      start[2] = -0.5
    }
    
    if(start[5]< -0.5){
      start[5] = -0.5
    }
    
    if(pi.dk*(1+par.prof[5]*max(max(z1), q1)/par.prof[6])^(1/par.prof[5])>=1 |
       is.nan(pi.dk*(1+par.prof[5]*max(max(z1), q1)/par.prof[6])^(1/par.prof[5]))){
      s = 1
      while(pi.dk*(1+par.prof[5]*max(max(z1), q1)/s)^(1/par.prof[5])>=1 | 
            is.nan(pi.dk*(1+par.prof[5]*max(max(z1), q1)/s)^(1/par.prof[5]))){
        s = s*1.5
      }
      start[6] = s
    } else {
      start = start
    }
    
    if(pi.ne*(1+par.prof[2]*max(max(z2), q2)/par.prof[3])^(1/par.prof[2])>=1 | 
       is.nan(pi.ne*(1+par.prof[2]*max(max(z2), q2)/par.prof[3])^(1/par.prof[2]))){
      s = 1
      while(pi.ne*(1+par.prof[2]*max(max(z2), q2)/s)^(1/par.prof[2])>=1 |
            is.nan(pi.ne*(1+par.prof[2]*max(max(z2), q2)/s)^(1/par.prof[2]))){
        s = s*1.5
      }
      start[3] = s
    } else{
      start = start
    }
    
    likfunc = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      a.1 = -l1*t.final.dk + N1*log(l1*t.final.dk)
      b.1 = y1*(log(pi.dk)+(1/g1)*log(1+g1*q1/s1))
      c.1 = (N1-y1)*log(1-pi.dk*(1+g1*q1/s1)^(1/g1))
      d.1 = sum(log(devd(z1, scale = s1, shape = g1, threshold = 0, type = "GP")))
      
      a.2 = -l2*t.final.ne + N2*log(l2*t.final.ne)
      b.2 = y2*(log(pi.ne)+(1/g2)*log(1+g2*q2/s2))
      c.2 = (N2-y2)*log(1-pi.ne*(1+g2*q2/s2)^(1/g2))
      d.2 = sum(log(devd(z2, scale = s2, shape = g2, threshold = 0, type = "GP")))
      
      lik1 = a.1+b.1+c.1+d.1
      lik2 = a.2+b.2+c.2+d.2
      lik = lik1 + lik2
      return(-lik)
    }
    
    grad.func = function(par){
      n1 = y1
      n2 = y2
      p1 = pi.dk
      p2 = pi.ne
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      diff.l.dk = -t.final.dk + N1/l1
      diff.g.dk.1 = n1 * (-1 / g1 ^ 2 * log1p(g1 * q1 / s1) + 1 / g1 * q1 / s1 / (1 + g1 * q1 / s1)) - (N1 - n1) * p1 * (1 + g1 * q1 / s1) ^ (1 / g1) * (-1 / g1 ^ 2 * log1p(g1 * q1 / s1) + 1 / g1 * q1 / s1 / (1 + g1 * q1 / s1)) / (1 - p1 * (1 + g1 * q1 / s1) ^ (1 / g1))
      diff.g.dk.2 = sum(1 / g1 ^ 2 * log1p(1 / s1 * g1 * z1) + (-1 / g1 - 1) / s1 * z1 / (1 + 1 / s1 * g1 * z1))
      diff.s.dk.1 = -n1 * q1 / s1 ^ 2 / (1 + g1 * q1 / s1) + (N1 - n1) * p1 * (1 + g1 * q1 / s1) ^ (1 / g1) * q1 / s1 ^ 2 / (1 + g1 * q1 / s1) / (1 - p1 * (1 + g1 * q1 / s1) ^ (1 / g1))
      diff.s.dk.2 = sum((-1 / s1 ^ 2 * (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1) - 1 / s1 ^ 3 * (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1) * (-1 / g1 - 1) * g1 * z1 / (1 + 1 / s1 * g1 * z1)) * s1 / (1 + 1 / s1 * g1 * z1) ^ (-1 / g1 - 1))
      
      diff.l.ne = -t.final.ne + N2/l2
      diff.g.ne.1 = n2 * (-1 / g2 ^ 2 * log1p(g2 * q2 / s2) + 1 / g2 * q2 / s2 / (1 + g2 * q2 / s2)) - (N2 - n2) * p2 * (1 + g2 * q2 / s2) ^ (1 / g2) * (-1 / g2 ^ 2 * log1p(g2 * q2 / s2) + 1 / g2 * q2 / s2 / (1 + g2 * q2 / s2)) / (1 - p2 * (1 + g2 * q2 / s2) ^ (1 / g2))
      diff.g.ne.2 = sum(1 / g2 ^ 2 * log1p(g2 * z2 / s2) + (-1 / g2 - 1) * z2 / s2 / (1 + g2 * z2 / s2))
      diff.s.ne.1 = -n2 * q2 / s2 ^ 2 / (1 + g2 * q2 / s2) + (N2 - n2) * p2 * (1 + g2 * q2 / s2) ^ (1 / g2) * q2 / s2 ^ 2 / (1 + g2 * q2 / s2) / (1 - p2 * (1 + g2 * q2 / s2) ^ (1 / g2))
      diff.s.ne.2 = sum((-1 / s2 ^ 2 * (1 + g2 * z2 / s2) ^ (-1 / g2 - 1) - 1 / s2 ^ 3 * (1 + g2 * z2 / s2) ^ (-1 / g2 - 1) * (-1 / g2 - 1) * g2 * z2 / (1 + g2 * z2 / s2)) * s2 / (1 + g2 * z2 / s2) ^ (-1 / g2 - 1))
      
      vec = c(diff.l.ne, diff.g.ne.1+diff.g.ne.2, diff.s.ne.1+diff.s.ne.2,
              diff.l.dk, diff.g.dk.1+diff.g.dk.2, diff.s.dk.1+diff.s.dk.2)
      
      return(-vec)
    }
    
    constraint.func = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      ne.1 = -(1+g2*max(max(z2), q2)/s2)+10^(-3)
      ne.2 = pi.ne*(1+g2*max(max(z2), q2)/s2)^(1/g2) - 1
      
      dk.1 = -(1+g1*max(max(z1), q1)/s1)+10^(-3)
      dk.2 = pi.dk*(1+g1*max(max(z1), q1)/s1)^(1/g1) - 1
      return(c(ne.1, ne.2, dk.1, dk.2))
    }
    
    grad.constraint = function(par){
      l2 = par[1]; g2 = par[2]; s2 = par[3]
      l1 = par[4]; g1 = par[5]; s1 = par[6]
      
      a.ne.1 = c(0, -max(max(z2), q2)/s2, g2*max(max(z2), q2)/(s2^2), 0, 0 ,0)
      a.ne.2 = c(0, pi.ne * (1 + g2 * max(max(z2), q2) / s2) ^ (1 / g2) * (-1 / g2 ^ 2 * log1p(g2 * max(max(z2), q2) / s2) + 1 / g2 * max(max(z2), q2) / s2 / (1 + g2 * max(max(z2), q2) / s2)), -pi.ne * (1 + g2 * max(max(z2), q2) / s2) ^ (0.1e1 / g2) * max(max(z2), q2) / s2 ^ 2 / (1 + g2 * max(max(z2), q2) / s2), 0,0,0)
      
      a.dk.1 = c(0, 0, 0, 0, -max(max(z1), q1)/s1, g1*max(max(z1), q1)/(s1^2))
      a.dk.2 = c(0, 0, 0, 0, pi.dk * (1 + g1 * max(max(z1), q1) / s1) ^ (1 / g1) * (-1 / g1 ^ 2 * log1p(g1 * max(max(z1), q1) / s1) + 1 / g1 * max(max(z1), q1) / s1 / (1 + g1 * max(max(z1), q1) / s1)), -pi.dk * (1 + g1 * max(max(z1), q1) / s1) ^ (0.1e1 / g1) * max(max(z1), q1) / s1 ^ 2 / (1 + g1 * max(max(z1), q1) / s1))
      return(rbind(a.ne.1, a.ne.2, a.dk.1, a.dk.2))
    }
    local_opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-10, "maxeval" = 1000)
    opts1 = list("algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-10, "maxeval" = 1000,
                 "local_opts" = local_opts1)
    optr = nloptr(x0 = c(start), eval_f = likfunc, eval_grad_f = grad.func,
                  eval_g_ineq = constraint.func, eval_jac_g_ineq = grad.constraint,
                  opts = opts1, lb = c(10^(-20), -0.5, 10^(-20), 
                                       10^(-20), -0.5, 10^(-20)),
                  ub = c(Inf, Inf, Inf, Inf, Inf, Inf))
    p.lik = optr$objective
    if(param == FALSE){
      return(p.lik)
    } else {
      return(optr$solution)
    }
  }
  
  constraint.func1 = function(pis){
    pi.dk = pis[1]
    pi.ne = pis[2]
    
    val = 2*(prof.lik + proflik2(c(pi.dk, pi.ne))) - qchisq(p = 0.95, df = 1)
    return(-val)
  }
  
  init = c(pi.prof)
  
  lb.pars = auglag(x0 = init, fn = obj.func1, hin = constraint.func1,
                   lower = c(10^(-20), 10^(-20)), 
                   upper = c(1-10^(-6), 1-10^(-6)),
                   localsolver = "COBYLA", localtol = 10^(-10))$par
  
  lb = lb.pars[2] - lb.pars[1]
  
  est = pi.prof[2]-pi.prof[1]
  
  ub.pars = auglag(x0 = init, fn = obj.func2, hin = constraint.func1,
                   lower = c(10^(-20), 10^(-20)), 
                   upper = c(1-10^(-6), 1-10^(-6)),
                   localsolver = "COBYLA", localtol = 10^(-10))$par
  ub = ub.pars[2] - ub.pars[1]
  
  ci = c(lb, est, ub)
  return(ci)
}




