#' @title MCMC-ABC for an SIR epidemic model
#'
#' @description Performs the Markov chain Approximate Bayesian computation
#' algorithm for the SIR stochastic model.
#'
#' @param I.obs a positive integer valued vector of daily (weekly, etc.) case
#'     counts. Days with zero cases must be included.
#' @param s0 the total population size excluding the initial infectious
#'     individuals.
#' @param i0 the number of infectious individuals at the beginning of the
#'     epidemic.
#' @param max.infections the maximum number of infections to occur in a single
#'     epidemic. If the epidemic tends to generate considerably more infections
#'     than the observed epidemic, it stops immediately and the proposed
#'     particle in the ABC algorithm is rejected. Use especially for the ABC
#'     algorithm when the  population is large (> 10,000) to reduce the amount
#'     of time spent on generating a single epidemic. The value of maximal
#'     infections must be set manually regarding the data I.obs.
#' @param times a desired size of the Markov chain
#' @param burn.in the number of iterations to be run before the actual
#'     Markov chain begins, assuring that the chain converges to the
#'    stationary distribution.
#' @param tolerance a positive number specifying the tolerance for accepting
#'     proposed particles.
#' @param transf a function transforming the case count. Must be specified
#'     with a single parameter - a case counts vector. If other parameters are
#'     intended to be used in the \code{transf} function (e.g.the population
#'     size, or the number of initial infectious individuals), they must be
#'     defined globally.
#' @param kern a kernel function (non-negative , symmetric around zero, with
#'     maximum at zero). Must be specified without any additional parameters,
#'     i.e. Must be a function of only the value at which we wish to evaluate
#'     the kernel.
#' @param dist.metr the specification of the Mahalanobis distance metric between
#'     the observed and simulated summary statistic. Must be specified as an
#'     inverse of a positive definite matrix of the distance metric or the
#'     string "euclidean". In the latter case, standard euclidean distance is
#'     used. Use the matrix option only if the dimension of the summary
#'     statistics is known beforehands!
#' @param prior a named list specifying the prior distribution of the epidemic
#'     parameters 'lambda' and 'mu'. Both distributions are considered
#'     independent! Therefore the joint distribution is a product of the
#'     marginals. For the exact syntax see the Details section.
#' @param prior.params a named list specifying parameters of the prior
#'     distribution of 'lambda' and 'mu'. These prior parameters are passed to
#'     the functions specified in the parameter \code{prior}. For the exact
#'     syntax see the Details section.
#' @param Sigma a positive definite covariance matrix of the Markov proposal 
#'     distribution, which is assumed to be a bivariate Gaussian.
#'     #' @param other.mod An alternative function sampling from the model, which 
#'     returns a list with element of infection times \code{I} and logical value
#'     \code{stopped}, which is \code{TRUE} if the \code{max.infections} was
#'     reached. See examples of the \code{sir.ABC} function for more details.
#' @param other.mod.params A named lists of parameters supplied to the 
#'     \code{other.mod} function.
#' @details
#'
#' The parameter \code{prior} must be of the form:
#'
#' \code{list(lambda.samp = sampler of 'lambda', 
#'            mu.samp = sampler of 'mu', 
#'            lambda.dens = prior density of 'lambda',
#'            mu.dens = prior density of 'mu')}
#'
#' @returns a list consisting of five to seven elements:
#'     \itemize{
#'         \item \code{accept.parts} a matrix containing the accepted particles
#'         \item \code{accept.parts.adj} a matrix containing accepted particles
#'         adjusted by the linear regression
#'         \item \code{summary.stats} a list of summary statistics
#'     }
#'    If a kernel was used, there is also a vector \code{kernel.weights}. If a
#'    transformation is used, then the list of the non-transformed case counts
#'    is returned too as the \code{daily.cases.list}.
#'
#' @examples
#' epi.obs <- c(6, 0, 1, 4, 7, 5, 9, 14, 17, 11, 10, 13, 6, 9, 5, 2, 5, 7, 11, 9,
#'  12, 19, 13, 12, 13, 16, 19, 10, 19, 19, 7, 12, 11, 11, 9, 10, 10, 13, 18, 5,
#'  4, 8, 3, 1, 3, 3, 2, 2, 0)
#' 
#' # The variant without kernel, with the uniform prior distribution
#' # Unif(0.1, 2.5), Unif(0.1, 1.1) for 'lambda' and 'mu' respectively. No
#' # transformation is used
#' set.seed(78)
#' sir.MCMC.ABC(
#'   I, 
#'   s0 = s0,
#'   i0 = i0, times = 2 * 1e4,
#'   burn.in = 100,
#'   tolerance = 80,
#'   prior = c(lambda.samp = runif, mu.samp = runif, 
#'             lambda.dens = dunif, mu.dens = dunif), 
#'   prior.params = list(lambda = c(min = 0.01, max = 2),
#'                       mu = c(min = 0.01, max = 1.5))
#' )
#'
#' # The transformation of the case counts onto the vector with the total size
#' # of the epidemic as the first element and the total length of the epidemic
#' # as the second element.
#' # element
#' transf.total <- function (x) {
#'    non.zero <- x != 0
#'    return(c(sum(x[non.zero]), 
#'             ifelse(any(non.zero), max(which(non.zero)), 0)))
#' }
#' 
#' # The variant with the gaussian kernel, with the prior distribution
#' # Unif(0.1, 2.5) and Unif(0.1, 1.5) for 'lambda' and 'mu' respectively. The
#' # transformation 'transf.total' is used, where we give the double weight
#' # to the total size (the first element of the summary statistics) by the
#' # 'dist.metr' parameter.
#' set.seed(78)
#' sir.MCMC.ABC(
#'   I$obs.cases,
#'   s0 = s0, 
#'   i0 = i0, 
#'   times = 2 * 1e4, 
#'   max.init.times = 30, 
#'   burn.in = 100, 
#'   tolerance = 350,
#'   kern = dnorm,
#'   dist.metr = diag(c(2, 1)),
#'   prior = c(lambda.samp = runif, mu.samp = runif, 
#'             lambda.dens = dunif, mu.dens = dunif),
#'   prior.params = list(lambda = c(min = 0.1, max = 2.5), 
#'                       mu = c(min = 0.1, max = 1.1))
#'  )

sir.MCMC.ABC <- function (
    I.obs, s0, i0, max.infections = s0,
    times = 100, max.init.times = 30, burn.in = floor(0.1 * times), tolerance = 1e2,
    transf = NULL, kern = NULL, dist.metr = "euclidean",
    prior = list(lambda.samp = runif, mu.samp = runif, lambda.dens = NULL,
                 mu.dens = NULL),
    Sigma = diag(c(1, 1)),
    prior.params = list(...),
    other.mod = NULL,
    other.mod.params = list(...)
    ) {
  
  # Functions for sampling a candidate particles and calculating density of the
  # Markov proposal distribution (a Gaussian) ==============================================
  
  # Creates a function for sampling from the bivariate normal distribution
  rbinorm <- function (mu1, sigma1, mu2, sigma2, rho)
  {
    x1 <- rnorm(1, mu1, sigma1)
    x2 <- rnorm(1, mu2 + (sigma2 / sigma1) * rho * (x1 - mu1),
                sqrt((1 - rho ^ 2) * sigma2 ^ 2))
    return(c(x1, x2))
  }
  
  # Input parameters check =====================================================
  
  # Stores the information, whether the user  wants to use a kernel function,
  # because it is tested later in multiple conditions
  use.kern.flag <- !is.null(kern)
  
  int.param.check <- c(s0, i0, times, max.infections, max.init.times, burn.in)
  if (any(int.param.check <= 0| int.param.check %% 1 != 0)) {
    stop("Parameters 's0', 'i0', 'max.infections', 'times',
         'max.init.times' and 'burn.in' must be positive integer values.")
  }
  if (!is.function(prior$lambda.samp) | !is.function(prior$mu.samp)) {
    stop("Parameter 'prior' must be a properly named list.")
  }
  if (!is.function(kern) & use.kern.flag) {
    stop("Parameter 'kern' must be an appropriate functions.")
  }
  if (!is.numeric(I.obs) | sum(I.obs %% 1, na.rm = TRUE) != 0) {
    stop("Parameter 'I.obs' must be an integer valued vector.")
  }
  # If transformation of case counts not specified, we define the transformation
  # function as an identity.
  use.transf.flag <- !is.null(transf)
  if (use.transf.flag) {
    if (!is.function(transf)) {
      stop("Parameter 'transf' must be a function of the case count.")
    }
  } else {
    transf <- function (x) {return(x)}
  }
  
  # If function measuring the distance of summary statistics not specified,
  # we define it as an euclidean distance. If it is specified, we define it as
  # a bilinear form of the observed and simulated summary statistic. Attention!
  # It is not checked, whether the matrix is positive definite or not.
  if (is.matrix(dist.metr)) {
    dist.fun <- function (obs.trans, samp.trans) {
      return(sqrt(t(obs.trans - samp.trans) %*% dist.metr %*%
                    (obs.trans - samp.trans)))
    }
  } else {
    if (dist.metr != "euclidean") {
      warning("The function computing distance between summary statistics must be
             specified by a positive definite matrix of the distance metric
             between the observed and simulated summary statistic or as a string
             'euclidean'. By default the Euclidean metric was used.")
    }
    dist.fun <- function (obs.trans, samp.trans) {
      return(sqrt(sum((obs.trans - samp.trans) ^ 2)))
    }
  }
  
  # Creates a kenel function based on the specifications in the 'kern'. If no
  # kernel function is specified, a constant function returning 1 is created.
  if (use.kern.flag) {
    
    # Normalises the input kernel function
    kernel.normalised <- function (x) {
      if (abs(x) > tolerance) {
        return(0)
      } else {
        return(kern(x / tolerance) / kern(0))
      }
    }
  } else {
    kernel.normalised <- function (x) {
      return(1)
    }
  }
  
  # Attention! Not all parameters checked properly, especially 'prior.params'.
  # Injection attack possible.
  
  # Preparation of general variables of the algorithm ==========================
  
  # Setting up the particular function, which we want to like to generate the 
  # epidemic. Either the 'gener.sir' function from the EpidemicABC package or
  # any other function with outputs compatible to the algorithm implementation.
  if (!is.null(other.mod)) {
    gener.epi <- other.mod
    model.params <- other.mod.params
  } else {
    gener.epi <- gener.sir
    model.params <- list(s0 = s0, i0 = i0, max.infections = max.infections)
  }
  
  # Finds out the length of the observed epidemic and transforms it
  max.day.obs <- length(I.obs) - 1
  obs.trans <- transf(I.obs)
  
  # Extracts parameters 'sigma1', 'sigma2', and 'rho' from 'Sigma'
  sigma1 <- sqrt(Sigma[1, 1])
  sigma2 <- sqrt(Sigma[2, 2])
  rho <- Sigma[1, 2] / (sigma1 * sigma2)
  
  # Converts the 'I.obs' parameter into data frame for better manipulation
  max.day.obs <- length(I.obs) - 1
  daily.cases.obs <- data.frame(
    day = 0:max.day.obs,
    obs.cases = I.obs
  ) # data frame type allows for comparing with simulated epidemics
  
  # Allocates the output parameters
  accept.parts <- matrix(data = NA, nrow = times, ncol = 2,
                         dimnames = list(NULL, c("lambda", "mu")))
  summary.stats <- vector(mode = "list", length = times)
  kernel.weights <- rep(NA, times)
  # If there is a transformation, summary statistics and the case counts must be
  # stored separately, because their dimension may differ due to the
  # transformation. If there is no transformation, the list of summary
  # statistics contains only vectors of daily case counts, whereas the list of
  # daily case counts contains data frames with additional column specifying the
  # day.
  daily.cases.list <- summary.stats
  
  # Allocates a local variable to store proposed particles lambda and mu
  parts <- c(NA, NA)
  
  # Flag inicating, whether to store the previous particle in the Markov chain
  store.prev <- TRUE
  
  # Prepares the argument lists for evaluating the prior density function
  lambda.prior.dens.args <- as.list(c(x = NA, prior.params$lambda))
  mu.prior.dens.args <- as.list(c(x = NA, prior.params$mu))
  
  # Initialisation of the algorithm ============================================
  distance <- Inf
  init.counter <- 0
  cont.init <- TRUE
  
  # Generates the initial particle of the MCMC algorithm. Maximum number of
  # iteration is 'max.init.times', then an error is displayed.
  while (cont.init & init.counter < max.init.times) {
    # Proposes initial particles lambda and mu
    # lambda stored at the first position, mu at the second position
    parts[1] <- do.call(prior$lambda.samp, args = as.list(c(1, prior.params$lambda)))
    parts[2] <- do.call(prior$mu.samp, args = as.list(c(1, prior.params$mu)))
    
    # Samples from the model and processes its output. Continuous infection
    # times are converted into case counts.
    epi.samp <- do.call(gener.epi, args = c(as.list(parts), model.params))
    I.samp <- epi.samp$I[!is.na(epi.samp)]
    daily.cases.samp <- as.data.frame(table(floor(I.samp), dnn = list("day")),
                                      responseName = "samp.cases")
    daily.cases.samp$day <- as.numeric(levels(daily.cases.samp$day))
    max.day.samp <- daily.cases.samp$day[nrow(daily.cases.samp)]
    # data.frame in this context is faster  than data.table
    
    # Makes a one-sided outer join of the data frame with the sampled case
    # counts and a data frame with only column "day" in order to add missing 
    # days with zero counts.
    if (max.day.samp != nrow(daily.cases.samp) - 1) {
      daily.cases.samp <- merge(data.frame(day = 0:max.day.samp), 
                                daily.cases.samp, by = "day", all.x = TRUE)
      daily.cases.samp$samp.cases[is.na(daily.cases.samp$samp.cases)] <- 0
    }
    
    # If the sampled epidemic lasted for more days than the observed one, the
    # observed epidemic is expanded by these days with zero cases.
    if (max.day.samp > max.day.obs) {
      I.obs <- c(I.obs, rep(0, max.day.samp - max.day.obs))
      max.day.obs <- max.day.samp
      # If the transformation function calculates one element of the vector on 
      # the basis of multiple other elements of the case count, the simplest 
      # way is to recalculate the the transformation of the whole vector of 
      # observed case counts
      obs.trans <- transf(I.obs)
      
    } else if (max.day.samp < max.day.obs) {
      daily.cases.samp <- rbind(
        daily.cases.samp,
        data.frame(day = (max.day.samp + 1):max.day.obs, samp.cases = 0)
      )
    }
    
    # Transforms the case count.
    samp.trans <- transf(daily.cases.samp$samp.cases)
    
    # Computes the distance
    distance <- dist.fun(obs.trans, samp.trans)
    print(distance)  # detete before submitting
    
    # Accepts the proposed particle pair with probability 
    # 'K_{tolerance}(distance)'
    threshold <- kernel.normalised(distance)
    U <- runif(1)
    if (U <= threshold) {
      cont.init <- FALSE
    }
    
    # Increase the counter of iterations needed for initialisation
    init.counter <- init.counter + 1
  }
  
  # If none of the proposed initial particles generated an epidemic close
  # enough to the observed one, tolerance might be too low, or it is just
  # the seed generating not close enough epidemics.
  if (init.counter == max.init.times & cont.init) {
    stop("Initialisation of the algorithm failed. Try changing the seed,
    increasing or the parameter 'max.init.times', or decreasing the tolerance.")
  }
  
  # Stores the initial particle together with summary statistics, daily case
  # counts and kernel weight into the pre-allocated matrices/lists. Kernel
  # weights of the initial particle are set as 1.
  accept.parts[1, ] <- parts
  kernel.weights[1] <- threshold
  summary.stats[[1]] <- samp.trans
  daily.cases.list[[1]] <- daily.cases.samp
  
  # Compute the prior probability of the initial particle. Using the assumption
  # that prior distributions of mu and lambda are independent.
  lambda.prior.dens.args[[1]] <- parts[1]
  mu.prior.dens.args[[1]] <- parts[2]
  prior.prob <- do.call(prior$lambda.dens, args = lambda.prior.dens.args) *
    do.call(prior$mu.dens, args = mu.prior.dens.args)
  
  
  # The burn-in period =========================================================
  
  for (total.counter in 1:burn.in) {
    # Proposes candidate particles lambda and mu
    # lambda stored at the first position, mu at the second position
    parts.cand <- rbinorm(parts[1], sigma1, parts[2], sigma2, rho)
    
    # Compute the prior probability of the candidate particle
    lambda.prior.dens.args[[1]] <- parts.cand[1]
    mu.prior.dens.args[[1]] <- parts.cand[2]
    prior.prob.cand <- do.call(prior$lambda.dens,
                               args = lambda.prior.dens.args) *
      do.call(prior$mu.dens, args = mu.prior.dens.args)
    
    # If the candidate particle is inside the support of the prior density,
    # the ABC and MCMC step of the algorithm is executed, i.e. an epidemic is
    # generated, the distance is checked and the candidate particle is (with
    # corresponding probability) either added into the chain or not. If it is
    # outside of the support, we keep the current particle (stay in the current
    # state of the chain) and proceed with the next iteration.
    if (prior.prob.cand > 0) {
      # Samples from the model and processes its output. Continuous infection
      # times are converted into case counts.
      epi.samp <- do.call(gener.epi, args = c(as.list(parts), model.params))
      
      if (!epi.samp$stopped) {
        # If there were less infections than specified in the 'max.infections',
        # the output of the epidemic model is processed. Continuous infection
        # times are converted into case counts.
        I.samp <- epi.samp$I[!is.na(epi.samp)]
        daily.cases.samp <- as.data.frame(table(floor(I.samp), dnn = list("day")),
                                          responseName = "samp.cases")
        daily.cases.samp$day <- as.numeric(levels(daily.cases.samp$day))
        max.day.samp <- daily.cases.samp$day[nrow(daily.cases.samp)]
        # data.frame in this context is faster  than data.table
        
        # Makes a one-sided outer join of the data frame with the sampled case
        # counts and a data frame with only column "day" in order to add missing 
        # days with zero counts.
        if (max.day.samp != nrow(daily.cases.samp) - 1) {
          daily.cases.samp <- merge(data.frame(day = 0:max.day.samp), 
                                    daily.cases.samp, by = "day", all.x = TRUE)
          daily.cases.samp$samp.cases[is.na(daily.cases.samp$samp.cases)] <- 0
        }
        
        # If the sampled epidemic lasted for more days than the observed one, the
        # observed epidemic is expanded by these days with zero cases.
        if (max.day.samp > max.day.obs) {
          I.obs <- c(I.obs, rep(0, max.day.samp - max.day.obs))
          max.day.obs <- max.day.samp
          # If the transformation function calculates one element of the vector on 
          # the basis of multiple other elements of the case count, the simplest 
          # way is to recalculate the the transformation of the whole vector of 
          # observed case counts
          obs.trans <- transf(I.obs)
          
        } else if (max.day.samp < max.day.obs) {
          daily.cases.samp <- rbind(
            daily.cases.samp,
            data.frame(day = (max.day.samp + 1):max.day.obs, samp.cases = 0)
          )
        }
        
        # Transforms the case count.
        samp.trans <- transf(daily.cases.samp$samp.cases)
        
        # Computes the distance
        distance.cand <- dist.fun(obs.trans, samp.trans)
        
        # Accepts or rejects proposed lambda and mu according to the algorithm
        # variant
        if (distance.cand < tolerance) {
          # Computes the probability of accepting the candidate particle. The
          # Markov proposal density does not need to be computed, since we use
          # the Gaussian distribution only, which is symmetrical.
          treshold <- min(
            1,
            (prior.prob.cand * kernel.normalised(distance.cand)) /
              (prior.prob * kernel.normalised(distance))
          )
          
          # Depending on the result of accepting, stores either the candidate particle
          # or the previous particle together with summary statistics and daily case
          # counts into the pre-allocated matrices/lists.
          U <- runif(1)
          if (U < treshold) {
            
            # Stores the candidate particle, which is now part of the Markov chain
            parts <- parts.cand
            accept.parts[1, ] <- parts
            summary.stats[[1]] <- samp.trans
            kernel.weights[1] <- kernel.normalised(distance.cand)
            daily.cases.list[[1]] <- daily.cases.samp
            # Prior probability of the candidate particle is now prior probability of
            # the Markov chain particle. The same with the distance.
            prior.prob <- prior.prob.cand
            distance <- distance.cand
          }
        }
      }
    }
    
    print(total.counter) # Shows, where the algorithm is. Delete before submitting.
  }
  
  # Main cycle of the algorithm ================================================
  
  for (total.counter in 2:times) {
    # Proposes candidate particles lambda and mu
    # lambda stored at the first position, mu at the second position
    parts.cand <- rbinorm(parts[1], sigma1, parts[2], sigma2, rho)
    
    # Compute the prior probability of the candidate particle
    lambda.prior.dens.args[[1]] <- parts.cand[1]
    mu.prior.dens.args[[1]] <- parts.cand[2]
    prior.prob.cand <- do.call(prior$lambda.dens,
                               args = lambda.prior.dens.args) *
      do.call(prior$mu.dens, args = mu.prior.dens.args)
    
    # If the candidate particle is inside the support of the prior density,
    # the ABC and MCMC step of the algorithm is executed, i.e. an epidemic is
    # generated, the distance is checked and the candidate particle is (with
    # corresponding probability) either added into the chain or not. If it is
    # outside of the support, we keep the current particle (stay in the current
    # state of the chain) and proceed with the next iteration.
    if (prior.prob.cand > 0) {
      # Samples from the model and processes its output. Continuous infection
      # times are converted into case counts.
      epi.samp <- do.call(gener.epi, args = c(as.list(parts), model.params))
      
      if (!epi.samp$stopped) {
        # If there were less infections than specified in the 'max.infections',
        # the output of the epidemic model is processed. Continuous infection 
        # times are converted into case counts.
        I.samp <- epi.samp$I[!is.na(epi.samp)]
        daily.cases.samp <- as.data.frame(
          table(floor(epi.samp$I), dnn = list("day"), useNA = "no"),
          responseName = "samp.cases"
        )
        daily.cases.samp$day <- as.numeric(levels(daily.cases.samp$day))
        max.day.samp <- daily.cases.samp$day[nrow(daily.cases.samp)]
        # data.frame in this context is faster than data.table
        
        # Makes a one-sided outer join of the data frame with the sampled case
        # counts and a data frame with only column "day" in order to add missing 
        # days with zero counts.
        if (max.day.samp != nrow(daily.cases.samp) - 1) {
          daily.cases.samp <- merge(data.frame(day = 0:max.day.samp), 
                                    daily.cases.samp, by = "day", all.x = TRUE)
          daily.cases.samp$samp.cases[is.na(daily.cases.samp$samp.cases)] <- 0
        }
        
        # If the sampled epidemic lasted for more days than the observed one, the
        # observed epidemic is expanded by these days with zero cases.
        if (max.day.samp > max.day.obs) {
          I.obs <- c(I.obs, rep(0, max.day.samp - max.day.obs))
          max.day.obs <- max.day.samp
          # If the transformation function calculates one element of the vector on 
          # the basis of multiple other elements of the case count, the simplest 
          # way is to recalculate the the transformation of the whole vector of 
          # observed case counts
          obs.trans <- transf(I.obs)
        } else if (max.day.samp < max.day.obs) {
          daily.cases.samp <- rbind(
            daily.cases.samp,
            data.frame(day = (max.day.samp + 1):max.day.obs, samp.cases = 0)
          )
        }
        
        # Transforms the case count.
        samp.trans <- transf(daily.cases.samp$samp.cases)
        
        # Computes the distance
        distance.cand <- dist.fun(obs.trans, samp.trans)
        
        # Accepts or rejects proposed lambda and mu according to the algorithm
        # variant
        if (distance.cand < tolerance) {
          # Computes the probability of accepting the candidate particle. The
          # Markov proposal density does not need to be computed, since we use
          # the Gaussian distribution only, which is symmetrical.
          treshold <- min(
            1,
            (prior.prob.cand * kernel.normalised(distance.cand)) /
              (prior.prob * kernel.normalised(distance))
          )
          
          # Depending on the result of accepting, stores either the candidate particle
          # or the previous particle together with summary statistics and daily case
          # counts into the pre-allocated matrices/lists.
          U <- runif(1)
          if (U < treshold) {
            
            # Stores the candidate particle, which is now part of the Markov chain
            parts <- parts.cand
            accept.parts[total.counter, ] <- parts
            summary.stats[[total.counter]] <- samp.trans
            kernel.weights[total.counter] <- kernel.normalised(distance.cand)
            daily.cases.list[[total.counter]] <- daily.cases.samp
            # Prior probability of the candidate particle is now prior probability of
            # the Markov chain particle. The same with the distance.
            prior.prob <- prior.prob.cand
            distance <- distance.cand
            
            # Previous particle will not be stored
            store.prev <- FALSE
          }
        }
        
        
      } else {
        # We end up here, if the epidemic generated had more cases than
        # specified in 'max.infect'. In that case, the candidate particle is
        # rejected and the previous particle is accepted.
        store.prev <- TRUE
      }
    }
    
    if (store.prev) {
      # Stores the previous particle
      accept.parts[total.counter, ] <- parts
      summary.stats[[total.counter]] <- summary.stats[[total.counter - 1]]
      daily.cases.list[[total.counter]] <- daily.cases.list[[total.counter - 1]]
    }
    
    # Sets the indicator of storing previous particle back to TRUE for the next
    # iteration
    store.prev <- TRUE
    print(total.counter) # Shows, where the algorithm is. Delete before submitting.
  }
  
  # Linear regression adjustment ===============================================
  
  # Adds the days, when there were zero cases into all data frames with daily
  # cases
  padded <- lapply(
    daily.cases.list,
    function (x) {
      len <- nrow(x)
      if (len < max.day.obs + 1){
        return(rbind(x, data.frame(day = len:max.day.obs, samp.cases = 0)))
      }
      return(x)
    }
  )
  
  # Creates the design matrix with (s_i - s_{obs}) as the i-th row.
  padded.transf <- lapply(
    padded,
    function (x) {return(transf(x$samp.cases) - obs.trans)}
  )
  covariates <- matrix(unlist(padded.transf), ncol = length(padded.transf[[1]]),
                       byrow = TRUE)
  
  # Tries to fit a linear model. There might be too little particles to fit
  # the model. If a kernel function was used, the weighted least squares are
  # used.
  if (use.kern.flag) {
    adj.model <- try(lm(accept.parts ~ covariates, weights = kernel.weights),
                     silent = TRUE)  # Possible to use just lm.fit()
  } else {
    adj.model <- try(lm(accept.parts ~ covariates), silent = TRUE)
    # Possible to use just lm.fit()
  }
  
  # If the fitting was successful, the columns of the design matrix might be
  # linearly dependent, typically couple of the last columns, which are mostly
  # padded by zeros.
  if (class(adj.model)[1] != "try-error") {
    independent <- !is.na(adj.model$coefficients[, 1])
    beta.hat <- adj.model$coefficients[independent, ]
    if (!is.null(dim(beta.hat))) {
      beta.hat <- as.matrix(beta.hat[-1, ])
      # If the summary statistic is one-dimensional, beta.hat is a column instead
      # a row and we need to transpose it.
      if (ncol(beta.hat) == 1) {
        beta.hat <- t(beta.hat)
      }
      accept.parts.adj <- accept.parts - covariates[, independent[-1], drop = FALSE] %*% beta.hat
    } else {
      class(adj.model)[1] <- "try-error"
    }
  }
  
  # Returning values ===========================================================
  
  # Prepares the return value list
  return.values <- list(
    accept.parts = accept.parts,
    summary.stats = summary.stats
  )
  
  if (class(adj.model)[1] != "try-error") {
    return.values$accept.parts.adj <- accept.parts.adj
  } else {
    warning("There was not enough accepted particles to fit a linear regression
            model.")
  }
  
  # If the user requires usage of kernel, we add the kernel weights into the
  # list of return values
  if (use.kern.flag) {
    return.values$kernel.weights <- kernel.weights
  }
  
  # If the user specifies some transformation, we additionally return the
  # non-transformed case counts
  if (use.transf.flag) {
    return.values$samp.counts <- daily.cases.list
  }
  
  return(return.values)
}