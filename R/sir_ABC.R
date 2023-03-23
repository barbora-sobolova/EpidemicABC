#' @title ABC for an SIR epidemic model
#'
#' @description Performs the classic Approximate Bayesian computation algorithm
#' for the SIR stochastic model.
#'
#' @param I.obs a positive integer valued vector of daily (weekly, etc.) case
#'     counts. Days with zero cases must be included.
#' @param n.pop the total population size excluding the initial infectious
#'     individuals.
#' @param m the number of infectious individuals at the beginning of the
#'     epidemic.
#' @param max.infections the maximum number of infections to occur in a single
#'     epidemic. If the epidemic tends to generate considerably more infections
#'     than the observed epidemic, it stops immediately and the proposed
#'     particle in the ABC algorithm is rejected. Use especially for the ABC
#'     algorithm when the  population is large (> 10,000) to reduce the amount
#'     of time spent on generating a single epidemic. The value of maximal
#'     infections must be set manually regarding the data I.obs.
#' @param times a desired size of the sample from the approximate posterior
#'     distribution.
#' @param max.times a maximal number of iterations.
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
#' @param importance a named list specifying the importance distribution of
#'     'lambda' and 'mu'. For the exact syntax see the Details section.
#' @param prior.params a named list specifying parameters of the prior
#'     distribution of 'lambda' and 'mu'. These prior parameters are passed to
#'     the functions specified in the parameter \code{prior}. For the exact
#'     syntax see the Details section.
#' @param importance.params a named list specifying parameters of the
#'     importance distribution of 'lambda' and 'mu'. These prior parameters are
#'     passed to the functions specified in the parameter \code{importance}.
#'     For the exact syntax see the Details section.
#' @param other.mod An alternative function sampling from the model, which 
#'     returns a list with element of infection times \code{I} and logical value
#'     \code{stopped}, which is \code{TRUE} if the \code{max.infections} was
#'     reached. See examples for more details.
#' @param other.mod.params A named lists of parameters supplied to the 
#'     \code{other.mod} function.
#'
#' @details
#'
#' The parameter \code{prior} must be either of the form:
#'
#' \code{list(lambda.samp = prior distribution of 'lambda',
#'            mu.samp = prior distribution of 'mu')}
#'
#' for the Rejection sampling variant, or:
#'
#' \code{list(lambda.samp = prior sampler of 'lambda',
#'            mu.samp = prior sampler of 'mu',
#'            lambda.dens = prior density of 'lambda',
#'            mu.dens = prior density of 'mu')}
#'
#' for the Importance sampling variant. In case of the Importance sampling, the
#' sampler and the density functions accept the same parameters passed from
#' \code{prior.params}.
#'
#' The parameter \code{importance} must be of the form:
#'
#' \code{list(lambda.dens = importance density of 'lambda',
#'       mu.dens = importance density of 'mu')}
#'
#' The sampler and the density functions accept the same parameters passed from
#' \code{importance.params}.
#'
#' The parameter \code{prior.params} must be of the form:
#'
#' \code{list(lambda = c(params for the prior sampler and the density of lambda),
#'            mu = c(params for the prior sampler and the density of mu))}
#'
#'
#' The parameter \code{importance.params} must be of the form:
#'
#' \code{list(lambda = c(params for the importance sampler and the density of lambda),
#'            mu = c(params for the importance sampler and the density of mu))}
#'
#' @returns a list consisting of five to seven elements:
#'     \itemize{
#'         \item \code{accept.parts} a matrix containing the accepted particles
#'         \item \code{accept.parts.adj} a matrix containing accepted particles
#'         adjusted by the linear regression
#'         \item \code{summary.stats} a list of summary statistics
#'         \item \code{n.iteration} the total number of iterations
#'         \item \code{n.accepted} the number of iterations, which resulted in
#'         accepting the proposed particle
#'     }
#'    In case of the importance sampling, there is an additional vector of the
#'    importance weights \code{IS.weights}. If a kernel was used, there is also
#'    a vector \code{kernel.weights}. If a transformation is used, then the list
#'    of the non-transformed case counts is returned too as the
#'    \code{daily.cases.list}.
#' @examples
#' 
#' epi.obs <- c(6, 0, 1, 4, 7, 5, 9, 14, 17, 11, 10, 13, 6, 9, 5, 2, 5, 7, 11, 9,
#'  12, 19, 13, 12, 13, 16, 19, 10, 19, 19, 7, 12, 11, 11, 9, 10, 10, 13, 18, 5,
#'   4, 8, 3, 1, 3, 3, 2, 2, 0)
#'   
#' # The transformation of the case counts onto the estimates of the offspring
#' # mean m.hat and the offspring variance sigma2.hat of the branching process
#' # approximating the general stochastic epidemic.The estimation uses case
#' # counts up to the highest point of the epidemic. If this highest point was
#' # the first day, the epidemic tends to cease very early, which does not
#' # qualitatively correspond to a longer epidemic, which we often observe (and
#' # which is in this case simulated using the seed 79). However the offspring
#' # mean and variance may be within the tolerance region and the proposed
#' # particle may be accepted. So the summary statistic is set as a sufficiently
#' # large number, which leads to the rejection.
#' i0 <- 5
#' transf <- function (x) {
#'
#'   N <- which.max(x)
#'   if (N == 1) {
#'   return(c(1e7, 1e7))
#'   }
#'   x.incr <- c(x[1] - i0, x[2:N])
#'   x.shifted <- c(i0, x.incr[-N])
#'
#'   m.hat <- sum(x.incr) / sum(x.shifted)
#'   sigma2.hat <- mean(x.shifted * (ifelse(x.shifted == 0,
#'                                          1, x.incr / x.shifted) - m.hat) ^ 2)
#'   return(c(m.hat, sigma2.hat))
#' }
#'
#' # The Importance sampling variant without kernel, with the uniform prior
#' # distribution Unif(0.1, 2.5), Unif(0.1, 1.5) for 'lambda' and 'mu'
#' # respectively and the importance distribution Gamma(2, 2) for both
#' # 'lambda' and 'mu'. The transformation to the branching process
#' # characteristics is used, where we put 10-times more weight onto the 
#' # offspring mean via the 'dist.metr' parameter.
#' set.seed(75)
#' sir.ABC(
#'   epi.obs, 
#'   s0 = s0, 
#'   i0 = i0,
#'   times = 1e2, 
#'   transf = transf,
#'   dist.metr = diag(c(10, 1)),
#'   max.times = 2e3,
#'   tolerance = 0.8, 
#'   kern = NULL,
#'   prior = c(lambda.dens = dunif, mu.dens = dunif),
#'   prior.params = list(lambda = c(min = 0.1, max = 2.5), 
#'                       mu = c(min = 0.1, max = 1.1)), 
#'   importance = list(lambda.samp = rgamma, mu.samp = rgamma,
#'   lambda.dens = dgamma, mu.dens = dgamma),
#'   importance.params = list(lambda = c(shape = 2, rate = 2),
#'                            mu = c(shape = 2, rate = 2))
#' )
#'
#' # The Rejection sampling variant with the gaussian kernel, with the prior
#' # distribution Unif(0.1, 2.5) and Unif(0.1, 1.1) for 'lambda' and 'mu'
#' # respectively. No transformation is used.
#' set.seed(76)
#' sir.ABC(
#'   epi.obs,
#'   s0 = 1000,
#'   i0 = 5,
#'   max.infections = 500,
#'   times = 1e2,
#'   max.times = 2e3,
#'   tolerance = 70,
#'   kern = dnorm,
#'   prior = c(lambda.samp = runif, mu.samp = runif),
#'   prior.params = list(lambda = c(min = 0.1, max = 2.5),
#'                       mu = c(min = 0.1, max = 1.1))
#' )
#' 
#' # The Rejection sampling variant with the gaussian kernel, with the prior
#' # distribution Unif(0.1, 2.5) and Unif(0.1, 1.1) for 'lambda' and 'mu'
#' # respectively. An alternative function from sampling from the model is used
#' # namely the 'ssa.adaptivetau' from the 'adaptivetau' library. No 
#' # transformation is used.
#' 
#' # Setting up the alternative sampling function
#' library(adaptivetau)
#' 
#' transitions <- list(
#'   c(S = -1, I = +1),  # infection
#'   c(I = -1, R = +1)  # recovery
#'   ) 
#'   
#' rates <- function (x, params, t) {
#'   s0 <- x["S"] + x["I"] + x["R"]  # total population size
#'   c(params$lambda * x["S"] * x["I"] / s0, # rate of infection
#'   params$mu * x["I"])  # rate of recovery
#'   }
#'     
#'  epi.len <- length(epi.obs) + 20
#'     
#'  model <- function (lambda, mu, init.values, transitions, rateFunc,
#'   max.duration) {
#'   epi <- as.data.frame(
#'   ssa.adaptivetau(init.values = init.values, transitions = transitions,
#'                   rateFunc = rateFunc, tf = max.duration,
#'                   params = list(lambda = lambda, mu = mu))
#'                   )
#'  return(list(I = c(rep(0, i0), epi$time[c(0, -diff(epi$S)) > 0]),
#'   stopped = FALSE))  # We ignore the max.infection parameter in this example.
#' }
#' 
#' # Algorithm itself
#' sir.ABC(
#'   epi.obs, 
#'   s0 = s0, 
#'   i0 = i0, 
#'   times = 1e2, 
#'   max.times = 2e3, 
#'   tolerance = 70, 
#'   kern = dnorm,
#'   prior = c(lambda.samp = runif, mu.samp = runif),
#'   prior.params = list(lambda = c(min = 0.1, max = 2.5), 
#'                       mu = c(min = 0.1, max = 1.1)
#'                       ),
#'   other.mod = model,
#'   other.mod.params = list(transitions = transitions, rateFunc = rates,
#'                           max.duration = epi.len,
#'                           init.values = c(S = s0, I = i0, R = 0))
#' )

sir.ABC <- function (
    I.obs, s0, i0, max.infections = s0, times = 100, max.times = 1e4, 
    tolerance = 1e2, transf = NULL, kern = NULL, dist.metr = "euclidean",
    prior = list(lambda.samp = runif, mu.samp = runif, lambda.dens = NULL,
                 mu.dens = NULL),  
    importance = list(lambda.samp = NULL, mu.samp = NULL, lambda.dens = NULL,
                      mu.dens = NULL), 
    prior.params = list(...), 
    importance.params = list(...),
    other.mod = NULL,
    other.mod.params = list(...)) {
  
  # Functions for different algorithm variants =================================
  K.IS.ABC <- function() {
    # Algorithm variant performing importance sampling with a kernel function
    
    threshold <- kernel.normalised(distance)
    # Accepts the proposed particle pair with probability 
    # 'K_{tolerance}(distance)'
    U <- runif(1)
    if (U <= threshold) {
      # Updates the counter of accepted particles and saves them together
      # with summary statistics, daily case counts, kernel weights and
      # importance sampling weights into the pre-allocated matrices/lists. The 
      # global assignment operator '<<-' is used to avoid copying the 
      # parameters, which would be otherwise passed to the function
      accept.counter <<- accept.counter + 1
      kernel.weights[accept.counter] <<- threshold
      accept.parts[accept.counter, ] <<- parts
      summary.stats[[accept.counter]] <<- samp.trans
      daily.cases.list[[accept.counter]] <<- daily.cases.samp
      IS.weights[accept.counter] <<- single.IS.weight
    }
  }
  
  K.RE.ABC <- function() {
    # Algorithm variant performing rejection sampling with a kernel function
    threshold <- kernel.normalised(distance)
    
    # Accepts the proposed particle pair with probability 
    # 'K_{tolerance}(distance)'
    U <- runif(1)
    if (U <= threshold) {
      # Updates the counter of accepted particles and saves them together
      # with summary statistics, daily case counts and kernel weights into the 
      # pre-allocated matrices/lists. The global assignment operator '<<-' is 
      # used to avoid copying the parameters, which would be otherwise passed to
      # the function
      accept.counter <<- accept.counter + 1
      kernel.weights[accept.counter] <<- threshold
      accept.parts[accept.counter, ] <<- parts
      summary.stats[[accept.counter]] <<- samp.trans
      daily.cases.list[[accept.counter]] <<- daily.cases.samp
    }
  }
  
  IS.ABC <- function() {
    # Algorithm variant performing importance sampling
    
    # Updates the counter of accepted particles and saves them together
    # with summary statistics, daily case counts and importance sampling weights
    # into the pre-allocated matrices/lists. The global assignment operator 
    # '<<-' is used to avoid copying the parameters, which would be otherwise
    # passed to the function
    accept.counter <<- accept.counter + 1
    accept.parts[accept.counter, ] <<- parts
    summary.stats[[accept.counter]] <<- samp.trans
    daily.cases.list[[accept.counter]] <<- daily.cases.samp
    IS.weights[accept.counter] <<- single.IS.weight
  }
  
  RE.ABC <- function() {
    # Algorithm variant performing rejection sampling with
    
    # Updates the counter of accepted particles and saves them together
    # with summary statistics and daily case counts into the pre-allocated 
    # matrices/lists. The global assignment operator '<<-' is used to avoid 
    # copying the parameters, which would be otherwise passed to the function
    accept.counter <<- accept.counter + 1
    accept.parts[accept.counter, ] <<- parts
    summary.stats[[accept.counter]] <<- samp.trans
    daily.cases.list[[accept.counter]] <<- daily.cases.samp
  }
  
  # Input parameters check =====================================================
  
  # Stores the information, whether the user  wants to use a kernel function or
  # importance sampling, because it is tested later in multiple conditions
  use.kern.flag <- !is.null(kern)
  use.importance.flag <- !is.null(unlist(importance))
  
  int.param.check <- c(s0, i0, times, max.times)
  if (any(int.param.check <= 0 | int.param.check %% 1 != 0)) {
    stop("Parameters 's0', 'i0', 'times' and 'max.times' must be positive 
         integer values.")
  }
  if ((!is.function(prior$lambda.samp) || !is.function(prior$mu.samp)) &&
      !use.importance.flag) {
    stop("Parameter 'prior' must be a properly named list of funtions.")
  }
  if (!is.function(kern) && use.kern.flag) {
    stop("Parameter 'kern' must be an appropriate functions.")
  }
  if ((!is.function(importance$lambda.samp) || !is.function(importance$mu.samp))
      && use.importance.flag) {
    stop("Parameters 'importance$lambda' and 'importance$mu' must be functions.")
  }
  if (!is.numeric(I.obs) || sum(I.obs %% 1, na.rm = TRUE) != 0) {
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
  
  # Attention! Not all parameters checked properly, especially 'prior.params'
  # and 'importance.params'. Injection attack possible.
  
  # Algorithm choice ===========================================================
  
  # Decides which type of an ABC algorithm to use based on the specifications 
  # in the 'kern' and 'importance' parameters.
  if (use.importance.flag) {
    
    # Sets the sampling distribution as the importance distribution.
    lambda.sampler <- importance$lambda.samp
    mu.sampler <- importance$mu.samp
    lambda.sampler.args <- as.list(c(n = 1, importance.params$lambda))
    mu.sampler.args <- as.list(c(n = 1, importance.params$mu))
    
    # Prepares the argument lists for evaluating the prior and importance
    # density functions
    lambda.prior.dens.args <- as.list(c(x = NA, prior.params$lambda))
    mu.prior.dens.args <- as.list(c(x = NA, prior.params$mu))
    lambda.importance.dens.args <- as.list(c(x = NA, importance.params$lambda))
    mu.importance.dens.args <- as.list(c(x = NA, importance.params$mu))
    
    # Allocates the output parameter 'weights'
    IS.weights <- rep(NA, times)
    
    if (use.kern.flag) {
      # This branch is the algorithm using a kernel function and the importance 
      # sampling
      
      # Allocates a vector of kernel weights
      kernel.weights <- rep(NA, times)
      # Normalises the input kernel function
      kernel.normalised <- function (x) {
        if (abs(x) > tolerance) {
          return(0)
        } else {
          return(kern(x / tolerance) / kern(0))
        }
      }
      
      # Sets the algorithm variant as Kernel-Importance-Sampling-ABC
      alg.variant <- K.IS.ABC
    } else {
      # This branch is the algorithm using the importance sampling only
      
      # Sets the algorithm variant as Importance-Sampling-ABC
      alg.variant <- IS.ABC
    }
    
  } else {
    
    # Sets the sampling distribution as the prior distribution.
    lambda.sampler <- prior$lambda
    lambda.sampler.args <- as.list(c(n = 1, prior.params$lambda))
    mu.sampler <- prior$mu
    mu.sampler.args <- as.list(c(n = 1, prior.params$mu))
    
    if (use.kern.flag) {
      # This branch is the algorithm using a kernel function and the rejection 
      # sampling
      
      # Allocates a vector of kernel weights
      kernel.weights <- rep(NA, times)
      # Normalises the input kernel function
      kernel.normalised <- function (x) {
        if (abs(x) > tolerance) {
          return(0)
        } else {
          return(kern(x / tolerance) / kern(0))
        }
      }
      
      # Sets the algorithm variant as Kernel-Rejection-Sampling-ABC
      alg.variant <- K.RE.ABC
    } else {
      # This branch is the algorithm using the rejection sampling only
      
      # Sets the algorithm variant as Rejection-Sampling-ABC
      alg.variant <- RE.ABC
    }
  }
  
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
  
  # Sets the variable, which indicates, whether to sample an epidemic from the
  # model. It stays always TRUE for the rejection sampling variants. For the 
  # importance sampling, it is set as False every time, there is proposed
  # with 0 prior probability.
  do.simulate <- TRUE
  
  # Finds out the length of the observed epidemic and transforms it
  max.day.obs <- length(I.obs) - 1
  obs.trans <- transf(I.obs)
  
  # Allocates the output parameters
  accept.counter <- 0
  total.counter <- 0
  accept.parts <- matrix(data = NA, nrow = times, ncol = 2,
                         dimnames = list(NULL, c("lambda", "mu")))
  summary.stats <- vector(mode = "list", length = times)
  # If there is a transformation, summary statistics and the case counts must be 
  # stored separately, because their dimension may differ due to the 
  # transformation. If there is no transformation, the list of summary 
  # statistics contains only vectors of daily case counts, whereas the list of 
  # daily case counts contains data frames with additional column specifying the 
  # day.
  daily.cases.list <- summary.stats
  
  
  # Allocates a local variable to store proposed particles lambda and mu
  parts <- c(NA, NA)
  
  # Main cycle of the algorithm ================================================
  
  while (accept.counter < times & total.counter < max.times) {
    
    # Proposing particles lambda and mu
    # lambda stored at the first position, mu at the second position
    parts[1] <- do.call(lambda.sampler, args = lambda.sampler.args) 
    parts[2] <- do.call(mu.sampler, args = mu.sampler.args)
    
    # If the importance weights of the proposed particle are 0 (i.e. the prior 
    # distribution is 0, we do not need to generate a sample epidemic)
    if (use.importance.flag) {
      
      # Completes the argument lists for evaluation of the prior and the 
      # importance density
      lambda.prior.dens.args[[1]] <- parts[1]
      mu.prior.dens.args[[1]] <- parts[2]
      lambda.importance.dens.args[[1]] <- parts[1]
      mu.importance.dens.args[[1]] <- parts[2]
      
      single.IS.weight <- 
        do.call(prior$lambda.dens, args = lambda.prior.dens.args) /
        do.call(importance$lambda.dens, args = lambda.importance.dens.args) *
        do.call(prior$mu.dens, args = mu.prior.dens.args) /
        do.call(importance$mu.dens, args = mu.importance.dens.args)
      # zero becomes FALSE, non-zero number becomes TRUE
      do.simulate <- as.logical(single.IS.weight)  
    }
    
    if (do.simulate) {
      # Samples from the model. If there were more infections, than specified in 
      # the 'max.infections', the whole sample is discarded and new iteration 
      # starts.
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
        distance <- dist.fun(obs.trans, samp.trans) # Delete before submitting.
        print(distance)
        
        # Accepts or rejects proposed lambda and mu according to the algorithm 
        # variant
        if (distance < tolerance) {
          alg.variant()
        }
      }
      
      print(total.counter) # Shows, where the algorithm is. Delete before submitting.
      total.counter <- total.counter + 1
    }
  }
  
  # Tests whether there was at least one accepted particle
  if (accept.counter == 0) {
    stop("Number of accepted particles is 0. Consider increasing tolerance
         or maximal number of iteration.")
  }
  
  # Linear regression adjustment ===============================================
  
  # Crops the NA elements of all variables needed for the linear model in case 
  # that the number of accepted particles is less than required.
  accept.parts <- accept.parts[1:accept.counter, ]
  daily.cases.list <- daily.cases.list[1:accept.counter]
  summary.stats <- summary.stats[1:accept.counter]
  
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
    kernel.weights <- kernel.weights[1:accept.counter]
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
      warning("It was impossible to perform the linear regression adjustment, 
              due to the lack of variability of the summary statistics.")
      class(adj.model)[1] <- "try-error"
    }
  }
  
  # Returning values ===========================================================
  
  # Prepares the return value list
  return.values <- list(
    accept.parts = accept.parts,
    summary.stats = summary.stats, 
    n.iterations = total.counter,
    n.accepted = accept.counter
  )
  
  if (class(adj.model)[1] != "try-error") {
    return.values$accept.parts.adj <- accept.parts.adj
  } else {
    warning("There was not enough accepted particles to fit a linear regression 
            model.")
  }
  
  # If a user requires importance sampling, we add the importance weights
  # into the list of return values
  if (use.importance.flag) {
    return.values$IS.weights <- IS.weights[1:accept.counter]
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