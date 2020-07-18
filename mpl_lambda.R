library('extRemes')
library('tidyverse')
library('nloptr')

#### Functions to produce maximum likelihood estimates
#### to initialize the modified profile likelihood

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


#### Modified profile likelihood ######

mplfunc.alt = function(data, measure1, thresh1, x1){
  poiss_time = as.POSIXct(data[["Date"]], tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  b00 = as.POSIXct(paste(str_sub(data[['Date']][1], 1, 10), "00:00:01"),
                   tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  b0 = b00
  b1 = 0
  bend = as.POSIXct(paste(str_sub(data[['Date']][length(data[,1])],1,10), "23:59:59"),
                    tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  
  
  t.final = as.numeric(difftime(bend, b00, units = "hours"))
  
  poisson.t = as.numeric(difftime(poiss_time, b0, units = "hours"))
  nanz = which(is.na(poisson.t))
  
  for (num in nanz){
    poisson.t[nanz] = (poisson.t[nanz+1]+poisson.t[nanz-1])/2
  }
  
  for(i in 2:length(poisson.t)){
    if(poisson.t[i]==poisson.t[i-1]){
      poisson.t[i] = poisson.t[i-1]+1/60
    } else {
      poisson.t[i] = poisson.t[i]
    }
  }
  
  data = cbind(data, poisson.t)
  
  
  data1 = data[data[[measure1]] != -1, ]
  N = length(data[[measure1]])
  data2 = data1[data1[[measure1]]<thresh1, ]
  z = -data2[[measure1]]+thresh1
  y = length(z)
  q1 = -x1+thresh1
  poisson.time = data1[data1[[measure1]]<thresh1, ][["poisson.t"]]
  ml.est = pois.est(data = data1,  
                    thresh=thresh1, measure = measure1, x = x1)
  ml.par = ml.est$par[3:4]
  
  ml.lambda = ml.est$`Crash intensity`/(t.final-b1)
  unname(ml.lambda)
  unname(ml.par)
  
  grad_vec = function(theta, t1, t0, l, zu){
    g = theta[1]
    s = theta[2]
    z = zu
    q = q1
    cg3 = c(-l * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) * (t1 - t0) - 1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s) + 1 / g ^ 2 * log1p(g * z / s) + (-1 / g - 1) * z / s / (1 + g * z / s),l * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (t1 - t0) - q / s ^ 2 / (1 + g * q / s) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) * s / (1 + g * z / s) ^ (-1 / g - 1))
    return(cg3)
  }
  
  l_mat = function(theta, theta0, qu, zu, t1, t0){
    l = theta[1]; l0=theta0[1]
    g = theta[2]; g0=theta0[2]
    s = theta[3]; s0=theta0[3]
    q = qu
    z = zu
    cg=matrix(c((-l * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) * (t1 - t0) - 1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s) + 1 / g ^ 2 * log1p(g * z / s) + (-1 / g - 1) * z / s / (1 + g * z / s)) * (-l0 * (1 + g0 * q / s0) ^ (1 / g0) * (-1 / g0 ^ 2 * log1p(g0 * q / s0) + 1 / g0 * q / s0 / (1 + g0 * q / s0)) * (t1 - t0) - 1 / g0 ^ 2 * log1p(g0 * q / s0) + 1 / g0 * q / s0 / (1 + g0 * q / s0) + 1 / g0 ^ 2 * log1p(g0 * z / s0) + (-1 / g0 - 1) * z / s0 / (1 + g0 * z / s0)),(l * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (t1 - t0) - q / s ^ 2 / (1 + g * q / s) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) * s / (1 + g * z / s) ^ (-1 / g - 1)) * (-l0 * (1 + g0 * q / s0) ^ (1 / g0) * (-1 / g0 ^ 2 * log1p(g0 * q / s0) + 1 / g0 * q / s0 / (1 + g0 * q / s0)) * (t1 - t0) - 1 / g0 ^ 2 * log1p(g0 * q / s0) + 1 / g0 * q / s0 / (1 + g0 * q / s0) + 1 / g0 ^ 2 * log1p(g0 * z / s0) + (-1 / g0 - 1) * z / s0 / (1 + g0 * z / s0)),(-l * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) * (t1 - t0) - 1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s) + 1 / g ^ 2 * log1p(g * z / s) + (-1 / g - 1) * z / s / (1 + g * z / s)) * (l0 * (1 + g0 * q / s0) ^ (1 / g0) * q / s0 ^ 2 / (1 + g0 * q / s0) * (t1 - t0) - q / s0 ^ 2 / (1 + g0 * q / s0) + (-1 / s0 ^ 2 * (1 + g0 * z / s0) ^ (-1 / g0 - 1) - 1 / s0 ^ 3 * (1 + g0 * z / s0) ^ (-1 / g0 - 1) * (-1 / g0 - 1) * g0 * z / (1 + g0 * z / s0)) * s0 / (1 + g0 * z / s0) ^ (-1 / g0 - 1)),(l * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (t1 - t0) - q / s ^ 2 / (1 + g * q / s) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) * s / (1 + g * z / s) ^ (-1 / g - 1)) * (l0 * (1 + g0 * q / s0) ^ (1 / g0) * q / s0 ^ 2 / (1 + g0 * q / s0) * (t1 - t0) - q / s0 ^ 2 / (1 + g0 * q / s0) + (-1 / s0 ^ 2 * (1 + g0 * z / s0) ^ (-1 / g0 - 1) - 1 / s0 ^ 3 * (1 + g0 * z / s0) ^ (-1 / g0 - 1) * (-1 / g0 - 1) * g0 * z / (1 + g0 * z / s0)) * s0 / (1 + g0 * z / s0) ^ (-1 / g0 - 1))),nrow=2,ncol=2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      return(cg)
  }
  
  j_mat = function(theta, qu, zu, t1, t0){
    l = theta[1]
    g = theta[2]
    s = theta[3]
    q = qu
    z = zu
    cg1 = matrix(c(-l * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) ^ 2 * (t1 - t0) - l * (1 + g * q / s) ^ (1 / g) * (2 / g ^ 3 * log1p(g * q / s) - 2 / g ^ 2 * q / s / (1 + g * q / s) - 1 / g * q ^ 2 / s ^ 2 / (1 + g * q / s) ^ 2) * (t1 - t0) + 2 / g ^ 3 * log1p(g * q / s) - 2 / g ^ 2 * q / s / (1 + g * q / s) - 1 / g * q ^ 2 / s ^ 2 / (1 + g * q / s) ^ 2 - 2 / g ^ 3 * log1p(g * z / s) + 2 / g ^ 2 * z / s / (1 + g * z / s) - (-1 / g - 1) * z ^ 2 / s ^ 2 / (1 + g * z / s) ^ 2,l * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) * (t1 - t0) - l * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 * (t1 - t0) + q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 - 1 / g * z / s ^ 2 / (1 + g * z / s) - (-1 / g - 1) * z / s ^ 2 / (1 + g * z / s) + (-1 / g - 1) * z ^ 2 / s ^ 3 / (1 + g * z / s) ^ 2 * g,l * (1 + g * q / s) ^ (1 / g) * q / s ^ 2 / (1 + g * q / s) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) * (t1 - t0) - l * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 * (t1 - t0) + q ^ 2 / s ^ 3 / (1 + g * q / s) ^ 2 - 1 / g * z / s ^ 2 / (1 + g * z / s) - (-1 / g - 1) * z / s ^ 2 / (1 + g * z / s) + (-1 / g - 1) * z ^ 2 / s ^ 3 / (1 + g * z / s) ^ 2 * g,-l * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 * (t1 - t0) - 2 * l * (1 + g * q / s) ^ (1 / g) * q / s ^ 3 / (1 + g * q / s) * (t1 - t0) + l * (1 + g * q / s) ^ (1 / g) * q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 * (t1 - t0) * g + 2 * q / s ^ 3 / (1 + g * q / s) - q ^ 2 / s ^ 4 / (1 + g * q / s) ^ 2 * g + (2 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) + 4 / s ^ 4 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s) + 1 / s ^ 5 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) ^ 2 * g ^ 2 * z ^ 2 / (1 + g * z / s) ^ 2 - 1 / s ^ 5 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g ^ 2 * z ^ 2 / (1 + g * z / s) ^ 2) * s / (1 + g * z / s) ^ (-1 / g - 1) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) / (1 + g * z / s) ^ (-1 / g - 1) + (-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) / s / (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)),nrow=2,ncol=2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          return(cg1)
  }
  
  j_obs = function(theta1, q, times){
    js = list()
    for (i in 1:length(z)){
      if(i == 1){
        js[[i]] = j_mat(theta = theta1, zu = z[i], 
                        t1 = times[i], t0 = 0, qu = q)
      } else {
        js[[i]] = j_mat(theta = theta1, zu = z[i],
                        t1 = times[i], t0 = times[i-1], qu = q)
      }
    }
    sum.js = Reduce('+', js)
    return(sum.js)
  }
  
  i_obs = function(theta1, thetaml, q, times){
    ls = list()
    for (i in 1:length(z)){
      if(i == 1){
        ls[[i]] = l_mat(theta = theta1, theta0 = thetaml, zu = z[i], 
                        t1 = times[i], t0 = 0, qu = q)
      } else {
        ls[[i]] = l_mat(theta = theta1, theta0 = thetaml, zu = z[i],
                        t1 = times[i], t0 = times[i-1], qu = q)
      }
    }
    sum.i = Reduce('+', ls)
    return(sum.i)
  }
  
  grad.func = function(theta1, times, qu, lambda){
    grad = list()
    for (i in 1:length(z)){
      if(i == 1){
        grad[[i]] = grad_vec(theta = theta1, l = lambda, t1 = times[i], t0 = b1,
                             zu = z[i])
      } else {
        grad[[i]] = grad_vec(theta = theta1, l = lambda, t1 = times[i], t0 = times[i-1],
                             zu = z[i])
      }
    }
    grad_mat = Reduce('+', grad)
    return(grad_mat)
  }
  
  mpl.lik.func = function(l1){
    
    start = ml.par
    
    if(ml.par[1]< -0.5){
      start[1]=-0.5
    }
    
    s1 = ml.par[2]
    while((1+ml.par[1]*max(max(z), q1)/s1)<=0){
      s1 = s1*1.05
    }
    
    start[2]=s1
    
    loglikgp = function(theta, lambda = l1){
      l = lambda
      g = theta[1]
      s = theta[2]
      tim = poisson.time
      a = matrix(nrow = 1, ncol = length(z))
      b = matrix(nrow = 1, ncol = length(z))
      c = matrix(nrow = 1, ncol = length(z))
      for(i in 1:length(z)){
        if(i == 1){
          a[i] = -(l*(1+g*q1/s)^(1/g))*(tim[i]-b1)
          b[i] = log((l*(1+g*q1/s)^(1/g))*(tim[i]-b1))
          c[i] = log((1/s)*(1+g*z[i]/s)^((-1/g)-1))
        } else {
          a[i] = -(l*(1+g*q1/s)^(1/g))*(tim[i]-tim[i-1])
          b[i] = log((l*(1+g*q1/s)^(1/g))*(tim[i]-tim[i-1]))
          c[i] = log((1/s)*(1+g*z[i]/s)^((-1/g)-1))
        }
      }
      
      logL = sum(a)+sum(b)+sum(c)
      return(-logL)
    }
    
    eval_func = function(x, l1){
      g = x[1]
      s = x[2]
      return(loglikgp(theta = c(g,s), lambda = l1))
    }
    
    eval_ineq = function(x, l1){
      g = x[1]
      s = x[2]
      constr = -(1+g*max(max(z),q1)/s)
      return(constr)
    }
    
    eval_grad_func = function(x, l1){
      g = x[1]
      s = x[2]
      grad = grad.func(theta1 = c(g,s), times = poisson.time, qu = q1, lambda = l1)
      return(-grad)
    }
    
    eval_ineq_func_jac = function(x, l1){
      g = x[1]
      s = x[2]
      grad = c(-max(max(z),q1)/s, g*max(max(z),q1)/(s^2))
      return(grad)
    }
    
    opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-7, maxeval = 1000)
    
    o = nloptr(x0 = start, l1 = l1,
               eval_f= eval_func,eval_g_ineq = eval_ineq, opts = opts1,
               eval_grad_f = eval_grad_func, eval_jac_g_ineq = eval_ineq_func_jac,
               lb = c(-Inf, 0), ub = c(Inf, Inf))
    
    j_obs_pl = -j_obs(theta1 = c(l1, o$solution), q = q1, 
                      times = poisson.time)
    i_obs_ml_pl = i_obs(theta1 = c(l1, o$solution), thetaml = c(ml.lambda, ml.par),
                        q = q1, times = poisson.time)
    
    lik = (-loglikgp(theta = o$solution, lambda = l1)+
             (1/2)*log(abs(det(j_obs_pl)))-log(abs(det(i_obs_ml_pl))))  
    
    return(-lik)
  }
  
  prof.lik.func = function(l1){
    
    start = ml.par
    if(ml.par[1]< -0.5){
      start[1]= -0.5
    }
    s1 = ml.par[2]
    while((1+ml.par[1]*max(max(z), q1)/s1)<=0){
      s1 = s1*1.5
    }
    
    start[2]=s1
    
    if(l1 <= 0){
      lik = Inf
    } else {
      loglikgp = function(theta, lambda = l1){
        l = lambda
        g = theta[1]
        s = theta[2]
        tim = poisson.time
        a = matrix(nrow = 1, ncol = length(z))
        b = matrix(nrow = 1, ncol = length(z))
        c = matrix(nrow = 1, ncol = length(z))
        for(i in 1:length(z)){
          if(i == 1){
            a[i] = -(l*(1+g*q1/s)^(1/g))*(tim[i]-b1)
            b[i] = log((l*(1+g*q1/s)^(1/g))*(tim[i]-b1))
            c[i] = log((1/s)*(1+g*z[i]/s)^((-1/g)-1))
          } else {
            a[i] = -(l*(1+g*q1/s)^(1/g))*(tim[i]-tim[i-1])
            b[i] = log((l*(1+g*q1/s)^(1/g))*(tim[i]-tim[i-1]))
            c[i] = log((1/s)*(1+g*z[i]/s)^((-1/g)-1))
          }
        }
        
        logL = sum(a)+sum(b)+sum(c)
        return(-logL)
      }
      
      eval_func = function(x, l1){
        g = x[1]
        s = x[2]
        return(loglikgp(theta = c(g,s), lambda = l1))
      }
      
      eval_ineq = function(x, l1){
        g = x[1]
        s = x[2]
        constr = -(1+g*max(max(z),q1)/s)
        return(constr)
      }
      
      eval_grad_func = function(x, l1){
        g = x[1]
        s = x[2]
        grad = grad.func(theta1 = c(g,s), times = poisson.time, qu = q1, lambda = l1)
        return(-grad)
      }
      
      eval_ineq_func_jac = function(x, l1){
        g = x[1]
        s = x[2]
        grad = c(-max(max(z),q1)/s, g*max(max(z),q1)/(s^2))
        return(grad)
      }
      
      opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-7, maxeval = 1000)
      
      o = nloptr(x0 = start, l1 = l1,
                 eval_f= eval_func,eval_g_ineq = eval_ineq, opts = opts1,
                 eval_grad_f = eval_grad_func, eval_jac_g_ineq = eval_ineq_func_jac,
                 lb = c(-Inf, 0), ub = c(Inf, Inf))
      lik = o$objective
      
    }
    return(lik)
  }
  
  optr = optim(par = ml.lambda, fn = mpl.lik.func, lower = ml.lambda*0.9, 
               upper = ml.lambda*2, method = "Brent")
  
  mpl_est = optr$par
  mpl_lik = optr$value
  
  bound.func = function(lb){
    diff = 2*(-mpl_lik+prof.lik.func(lb))-qchisq(p = 0.95, df = 1)
    return(-diff)
  }
  
  mini.func = function(lb){
    return(lb)
  }
  
  maxi.func = function(lb){
    return(-lb)
  }
  
  lower = auglag(x0 = mpl_est/2, fn = mini.func, hin = bound.func, lower = 10^(-20),
                 localtol = 10^(-10), upper = 0.75*mpl_est)$par
  upper = auglag(x0 = mpl_est*1.5, fn = maxi.func, hin = bound.func, lower = mpl_est*1.25,
                 localtol = 10^(-10))$par
  
  return(c(t.final*lower, t.final*mpl_est, t.final*upper))
}





