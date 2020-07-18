library('nloptr')
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


##### Profile likelihood confidence intervals ######

proflikest.alt = function(data, thresh1, measure1, x1){
  data1 = data[data[[measure1]] != -1, ]
  N = length(data[[measure1]])
  data2 = data1[data1[[measure1]]<thresh1, ]
  z = -data2[[measure1]]+thresh1
  q = -x1+thresh1
  y = length(z)
  
  b00 = as.POSIXct(paste(str_sub(data[['Date']][1], 1, 10), "00:00:01"),
                   tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  
  bend = as.POSIXct(paste(str_sub(data[['Date']][length(data[,1])],1,10), "23:59:59"),
                    tz = "", tryFormats = "%Y-%m-%d %H:%M:%OS")
  
  t = as.numeric(difftime(bend, b00, units = "hours"))
  
  ml.est = pois.est(data = data1,
                    thresh=thresh1, measure = measure1, x = x1)
  ml.lambda = ml.est$`Crash intensity`/t
  ml.ci = ml.est$`CI Crash Intensity`/t
  ml.par = ml.est$par[2:4]
  unname(ml.lambda)
  unname(ml.par)
  
  opt.con.lik = function(l){
    if(l <=0){
      lik = Inf
    } else {
      loglikgp = function(theta){
        p = theta[1]
        g = theta[2]
        s = theta[3]
        a = (-l*t*(1+g*q/s)^(1/g))/(p)
        b = N*(log(l)+(1/g)*log(1+g*q/s)-log(p)+log(t))
        c = y*log(p)+(N-y)*log(1-p)
        d = sum(log((1/s)*(1+g*z/s)^(-1/g-1)))
        logL = a+b+c+d
        return(-logL)
      }
      
      grad.func = function(theta){
        n = y
        p = theta[1]
        g = theta[2]
        s = theta[3]
        diff.p = l / p ^ 2 * t * (1 + g * q / s) ^ (0.1e1 / g) - N / p + n / p - (N - n) / (1 - p)
        diff.g1 = -l / p * t * (1 + g * q / s) ^ (1 / g) * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s)) + N * (-1 / g ^ 2 * log1p(g * q / s) + 1 / g * q / s / (1 + g * q / s))
        diff.g2 = sum(1 / g ^ 2 * log1p(g * z / s) + (-1 / g - 1) * z / s / (1 + g * z / s))
        diff.s1 = l / p * t * (1 + g * q / s) ^ (0.1e1 / g) * q / s ^ 2 / (1 + g * q / s) - N * q / s ^ 2 / (1 + g * q / s)
        diff.s2 = sum((-1 / s ^ 2 * (1 + g * z / s) ^ (-1 / g - 1) - 1 / s ^ 3 * (1 + g * z / s) ^ (-1 / g - 1) * (-1 / g - 1) * g * z / (1 + g * z / s)) * s / (1 + g * z / s) ^ (-1 / g - 1))
        
        return(-c(diff.p, diff.g1 + diff.g2, diff.s1+diff.s2))
      }
      
      constraint.func = function(theta){
        p = theta[1]
        g = theta[2]
        s = theta[3]
        
        return(-c(1+g*max(max(z),q)/s-10^(-10)))
      }
      
      grad.constraint = function(theta){
        p = theta[1]
        g = theta[2]
        s = theta[3]
        
        a = -c(0, max(max(z),q)/s, -g*max(max(z),q)/(s^2))
        return(a)
      }
      
      opts1 = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-10, maxeval = 1000)
      o = nloptr(x0 = ml.par, eval_f = loglikgp, eval_grad_f = grad.func,
                 eval_g_ineq = constraint.func, eval_jac_g_ineq = grad.constraint,
                 opts = opts1, lb = c(0, -0.5, 10^(-50)), ub = c(1, Inf, Inf))
      lik = o$objective
    }
    return(lik)
  }
  
  est = optim(par = ml.lambda, fn = opt.con.lik, method = "Brent", upper = 2.5*ml.lambda,
              lower = 0.05*ml.lambda)
  
  prof.lik = -est$value
  l.prof = est$par
  unname(prof.lik)
  ci.bound = prof.lik - (1/2)*qchisq(p = 0.95, df = 1)
  
  bound.func = function(l1){
    val = 2*(prof.lik + opt.con.lik(l = l1))-qchisq(p = 0.95, df = 1)
    return(-val)
  }
  
  mini.func = function(l1){
    return(l1)
  }
  
  maxi.func = function(l1){
    return(-l1)
  }
  
  
  lowerc = auglag(x0 = 0.5*l.prof, fn = mini.func, hin = bound.func, lower = 0,
                  upper = l.prof - 0.001*l.prof, localtol = 10^(-5))$par
  upperc = auglag(x0 = 1.2*l.prof, fn = maxi.func, hin = bound.func, lower = 1.01*l.prof,
                  upper = 4*l.prof, localtol = 10^(-5))$par
  
  prof.sum = matrix(c(lowerc*t, l.prof*t, upperc*t), nrow = 1, ncol = 3)
  colnames(prof.sum) = c("lower", "estimate", "upper")
  
  return(prof.sum)

}


  
  
  
  