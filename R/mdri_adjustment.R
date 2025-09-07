# Copyright 2022 UNAIDS, World Health Organization, and individual contributors.
# Author: Eduard Grebe
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(tidyverse)

weibull_survival = function(t, alpha, beta = 4.5) {
  exp(-(t/alpha)^beta)
}

PR <- function(t, pr_parms) {
  1/(1 + exp(-(pr_parms[1] + pr_parms[2] * t + pr_parms[3] * t^2 + pr_parms[4] *
                 t^3)))
}

PprimeR = function(t, pr_parms, alpha, beta = 4.5) {
  PR(t, pr_parms) * weibull_survival(t, alpha, beta)
}

# For FRR estimation rho(t)
weibull_density = function(t, alpha, beta, bigT, max_t = 5000) {
  if (t < bigT) {
    return(NULL)
    } else {
      integ <- cubature::adaptIntegrate(
        weibull_survival,
        lowerLimit = bigT,
        upperLimit = max_t,
        alpha = alpha,
        beta = beta
      )$integral
    return(weibull_survival(t, alpha, beta) / integ)
    }
}

PR_weighted_frr <- function(t, pr_parms, alpha, beta, bigT, max_t = 5000) {
  if (t < bigT) {
    return(NULL)
  } else {
    PR_weighted = PR(t, pr_parms) * weibull_density(t, alpha, beta, bigT, max_t)
    return(PR_weighted)
  }
}

frr_untreated_weighted <- function(pr_parms, bigT, alpha, beta, max_t = 5000) {
    integ <- cubature::adaptIntegrate(
      PR_weighted_frr,
      lowerLimit = bigT,
      upperLimit = max_t,
      pr_parms = pr_parms,
      alpha = alpha,
      beta = beta,
      bigT = bigT
    )$integral
  return(integ)
}

# quantile_weibull = function(q = 0.5, alpha, beta = 4.5) {
#   (-log(q)/(alpha^beta))^(-beta)
# }

weibull_minus = function(t, alpha, beta, q = 0.5) {
  weibull_survival(t, alpha, beta) - q
}

weibull_quantile = function(q, alpha, beta, max_t = 4000) {
  if (q < 0 | q > 1) {stop("please provide a value of q in [0,1]")}
  t = uniroot(weibull_minus,
          interval = c(0,max_t),
          alpha = alpha,
          beta = beta,
          q = 1 - q)
  return(t$root)
}

find_scale = function(target_median = 365, beta = 2.5, max_t = 4000) {
  f_delta <- function(alpha) {
    realised_median = weibull_quantile(q = 0.5, alpha = alpha, beta = beta)
    delta = target_median - realised_median
    return(delta)
  }
  success <- FALSE
  max_alpha <- 100
  while (!success) {
    alpha <- try(uniroot(f_delta, lower = 1, upper = max_alpha)$root, silent = TRUE)
    if (class(alpha) == "try-error") {
        max_alpha <- max_alpha + 100 #* 2
        #print(paste0("Max alpha increased to:", max_alpha))
      } else {
        success <- TRUE
      }
  }
  return(alpha)
}

find_shape_quantile = function(target = 365.25, quantile = 0.25, alpha = 422, max_t = 20000) {
  f_delta <- function(beta) {
    realised = weibull_quantile(q = quantile, alpha = alpha, beta = beta)
    delta = target - realised
    return(delta)
  }
  success <- FALSE
  max_beta <- 50
  while (!success) {
    beta <- try(uniroot(f_delta, lower = 0.01, upper = max_beta)$root, silent = TRUE)
    if (class(beta) == "try-error") {
      max_beta <- max_beta + 10
    } else {
      success <- TRUE
    }
  }
  return(beta)
}

find_shape = function(target_90th = 730, alpha = 422, max_t = 4000, max_beta = 50) {
  #browser()
  f_delta <- function(beta) {
    realised_90th = weibull_quantile(q = 0.9, alpha = alpha, beta = beta)
    delta = target_90th - realised_90th
    return(delta)
  }
  beta = tryCatch({
    uniroot(f_delta, lower = 0.5, upper = max_beta)$root
  }, error = function(err) {
    return(NA)
  })



  return(beta)
}

plot_weibull = function(from = 0, to = 2000, alpha = 422, beta = 2.5, vertline_T = 730.5, vertline_median) {
  #browser()
  curve(exp(-(x/alpha)^beta), from = from, to = to, xlab = "t", ylab = "S(t)")
  abline(v = vertline_T, col = "red")
  abline(v = vertline_median, col = "blue")
}



