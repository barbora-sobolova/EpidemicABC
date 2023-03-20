#' @title ABC for an SVEIRD epidemic model
#' 
#' @description Performs the classic Approximate Bayesian computation algorithm
#' for the SVEIRD stochastic model.
#' 
#' @param I.obs A positive integer valued vector of daily case counts. Days
#' with zero cases must be included.
#' @param s0 A positive integer specifying the total population size
#' excluding the initial infectious individuals.
#' @param i0 A positive integer specifying the number of infectious individuals.
#' @param v0 A positive integer specifying the number of vaccinated individuals.
#' @param r0 A positive integer specifying the number of recovered individuals.
#' @param max.infections A positive integer specifying the maximum number of 
#' infections to occur in a single epidemic. If the epidemic tends to generate
#' more infections, it stops immediately and the proposed particle in the ABC 
#' algorithm is rejected. Use especially for the ABC algorithm when the 
#' population is large to reduce the amount of time spent on generating a
#' single epidemic. The value of maximal infections must be set manually 
#' regarding the data I.obs. 
#' @param times A desired size of the sample from the approximate posterior
#' distribution
#' @param max.times A maximal number of iterations
#' @param tolerance A positive number specifying the tolerance for accepting 
#' proposed particles.
#' @param transf A named list of functions transforming the case count. All 
#' functions must be specified with one parameters - a case counts vector. The
#' susceptible population size and the number of initial infectious individuals 
#' can be used too, since these variables exist within the 'sir.ABC' function. 
#' Their names are 's0' and 'i0' respectively.
#' @param kern A kernel non-negative function, symmetric
#' around zero, with maximum at zero. Must be specified without any additional 
#' parameters, i.e. Must be a function of only the value at which we wish to 
#' evaluate the kernel.
#' @param dist.metr Specification of the Mahalanobis distance metric between 
#' the observed and simulated summary statistic. Must be specified as an 
#' inverse of a positive definite matrix of the distance metric or the string 
#' "euclidean". In the latter case, standard euclidean distance is used. 
#' Use the matrix option only if the dimension of the summary statistics is 
#' known beforehands!
#' @param prior A named list specifying the prior distribution of 'lambda' and 
#' 'mu'. Both distributions are considered independent! Therefore the joint
#' distribution is a product of the marginals. Must be either of the form:
#' list(lambda.samp = prior distribution of 'lambda',
#' mu.samp = prior distribution of 'mu', nu.samp = prior distribution of 'nu')
#' for the Rejection sampling variant, or of the form: 
#' list(lambda.dens = prior density of 'lambda', 
#' mu.dens = prior density of 'mu', nu.dens = prior density of 'nu')
#' for the Importance sampling variant. In case of the Importance sampling, 
#' the sampler and the density functions accept the same parameters passed 
#' from 'prior.params'.
#' @param importance A named list specifying the importance distribution of 
#' 'lambda', 'mu' and 'nu'. Must be of the form: 
#' list(lambda.samp = importance sampler of 'lambda',
#' mu.samp = importance sampler of 'mu', lambda.dens = importance density of 'lambda',
#' mu.dens = importance density of 'mu'). The sampler and the density functions
#' accept the same parameters passed from 'prior.params'.
#' @param prior.params A named list specifying parameters of the prior 
#' distribution of 'lambda' and 'mu'. These prior parameters are passed to 
#' the functions specified in the parameter 'prior'. Must be of the form:
#' list(lambda = c(params for the prior sampler and density of lambda),
#' mu = c(params for the prior sampler and density of mu)).
#' @param importance.params A named list specifying parameters of the 
#' importance distribution of 'lambda' and 'mu'. These prior parameters are 
#' passed to the functions specified in the parameter 'importance'. Must be of 
#' the form: 
#' list(lambda = c(params for the importance sampler and density of lambda),
#' mu = c(params for the importance sampler and density of mu)).
#' 
#' @returns A list consisting of five to seven elements: a numerical matrix 
#' containing accepted particles, a numerical matrix containing accepted
#' particles adjusted by the linear regression, a list of summary statistics, 
#' total number of iterations and number of iterations, which resulted in 
#' accepting the proposed particle. In case of importance sampling, there is an 
#' additional numerical matrix of importance weights. If a transformation
#' is used, then list of non-transformed case counts is returned too.
#' @examples
#' 
#' transf <- function (x) {
#'  log(x + 1)
#'}
#' gse <- gener.sir(n.pop = 100, m = 5, lambda = 1.5, mu = 0.5)
#' I <- as.vector(table(floor(gse$I)))
#' 
#' # The Importance sampling variant without kernel, with the prior 
#' # distribution
#' Unif(0.5, 2), Unif(0.1, 1.1) for lambda and mu respectively and the 
#' importance distribution Gamma(1.5, 1) for both lambda and mu.
#' sir.ABC(
#' I, n.pop = 100, m = 5, times = 1e2, max.times = 1e3, tolerance = 50, 
#' prior = c(lambda.samp = runif, mu.samp = runif, 
#'           lambda.dens = dunif, mu.dens = dunif),
#' prior.params = list(lambda = c(min = 0.5, max = 2), 
#'                     mu = c(min = 0.1, max = 1.1)), 
#' importance = list(lambda.samp = rgamma, mu.samp = rgamma,
#'                   lambda.dens = dgamma, mu.dens = dgamma),
#' importance.params = list(lambda = c(shape = 1.5, rate = 1),
#'                          mu = c(shape = 1.5, rate = 1))
#' )
#' 
#' The Rejection sampling variant with the standardised normal kernel, with 
#' the prior distribution Unif(0.5, 2), Unif(0.1, 1.1) for lambda and mu 
#' respectively and the importance distribution Gamma(1.5, 1) for both lambda 
#' and mu. A transformation 'transf' is used.
#' sir.ABC(
#' I, n.pop = 100, m = 5, times = 1e2, max.times = 1e3, tolerance = 3, 
#' kern = dunif, transf = transf,
#' prior = c(lambda.samp = runif, mu.samp = runif),
#' prior.params = list(lambda = c(min = 0.5, max = 2), 
#'                     mu = c(min = 0.1, max = 1.1))
#' )

sveird.ABC <- function (
    I.obs, s0, i0, v0 = 0, r0 = 0, max.infections = 1e4, times = 100, 
    max.times = 1e4, max.duration = 100, tolerance = 1e2, 
    transf = list(I = NULL, VI = NULL, RI = NULL, RVI = NULL, V = NULL, R = NULL, D = NULL), 
    kern = NULL, V.obs = NULL, R.obs = NULL, D.obs = NULL, VI.obs = NULL,
    RI.obs = NULL, RVI.obs = NULL,
    compart.weights = c(I = 1, VI = NA, RI = NA, RVI = NA, V = NA, R = NA, D = NA),
    dist.metr = list(I = "euclidean", VI = NULL, RI = NULL, RVI = NULL, V = NULL,
                     R = NULL, D = NULL),
    prior = list(lambda.samp = runif, mu.samp = runif, nu.samp = runif,
                 psi.samp = runif, phi.samp = runif, kappa.samp = runif,
                 omega.samp = runif, lambda.dens = NULL, mu.dens = NULL,
                 nu.dens = NULL, psi.dens = NULL, phi.dens = NULL,
                 kappa.dens = NULL, omega.dens = NULL),
    importance = list(lambda.samp = NULL, mu.samp = NULL, nu.samp = NULL,
                      psi.samp = NULL, phi.samp = NULL, kappa.samp = NULL,
                      omega.samp = NULL, lambda.dens = NULL, mu.dens = NULL, 
                      nu.dens = NULL, psi.dens = NULL, phi.dens = NULL,
                      kappa.dens = NULL, omega.dens = NULL), 
    prior.params = list(...), 
    importance.params = list(...)) {
  
  # Functions for different algorithm variants =================================
  K.IS.ABC <- function() {
    # Algorithm variant performing importance sampling with a kernel function
    
    threshold <- kernel.normalised(distance.total)
    # Accepts the proposed particle pair with probability 
    # 'K_{tolerance}(distance)'
    U <- runif(1)
    if (U <= threshold) {
      
      # Updates the counter of accepted particles and saves them together
      # with summary statistics, daily case counts, kernel weights and
      # importance sampling weights into the pre-allocated matrices/lists. The 
      # global assignment operator '<<-' is used to avoid copying the 
      # parameters, which would be otherwise passed to the function.
      accept.counter <<- accept.counter + 1
      kernel.weights[accept.counter] <<- threshold
      accept.parts[accept.counter, ] <<- parts
      IS.weights[accept.counter] <<- single.IS.weight
      summary.stats$I[[accept.counter]] <<- results.I$samp.trans
      daily.list$I[[accept.counter]] <<- results.I$daily.samp
      if (!is.null(VI.obs)) {
        summary.stats$VI[[accept.counter]] <<- results.VI$samp.trans
        daily.list$VI[[accept.counter]] <<- results.VI$daily.samp
      }
      if (!is.null(RI.obs)) {
        summary.stats$RI[[accept.counter]] <<- results.RI$samp.trans
        daily.list$RI[[accept.counter]] <<- results.RI$daily.samp
      }
      if (!is.null(RVI.obs)) {
        summary.stats$RVI[[accept.counter]] <<- results.RVI$samp.trans
        daily.list$RVI[[accept.counter]] <<- results.RVI$daily.samp
      }
      if (!is.null(R.obs)) {
        summary.stats$R[[accept.counter]] <<- results.R$samp.trans
        daily.list$R[[accept.counter]] <<- results.R$daily.samp
      }
      if (!is.null(V.obs)) {
        summary.stats$V[[accept.counter]] <<- results.V$samp.trans
        daily.list$V[[accept.counter]] <<- results.V$daily.samp
      }
      if (!is.null(D.obs)) {
        summary.stats$D[[accept.counter]] <<- results.D$samp.trans
        daily.list$D[[accept.counter]] <<- results.D$daily.samp
      }
    }
  }
  
  K.RE.ABC <- function() {
    # Algorithm variant performing rejection sampling with a kernel function
    threshold <- kernel.normalised(distance.total)
    
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
      summary.stats$I[[accept.counter]] <<- results.I$samp.trans
      daily.list$I[[accept.counter]] <<- results.I$daily.samp
      if (!is.null(VI.obs)) {
        summary.stats$VI[[accept.counter]] <<- results.VI$samp.trans
        daily.list$VI[[accept.counter]] <<- results.VI$daily.samp
      }
      if (!is.null(RI.obs)) {
        summary.stats$RI[[accept.counter]] <<- results.RI$samp.trans
        daily.list$RI[[accept.counter]] <<- results.RI$daily.samp
      }
      if (!is.null(RVI.obs)) {
        summary.stats$RVI[[accept.counter]] <<- results.RVI$samp.trans
        daily.list$RVI[[accept.counter]] <<- results.RVI$daily.samp
      }
      if (!is.null(R.obs)) {
        summary.stats$R[[accept.counter]] <<- results.R$samp.trans
        daily.list$R[[accept.counter]] <<- results.R$daily.samp
      }
      if (!is.null(V.obs)) {
        summary.stats$V[[accept.counter]] <<- results.V$samp.trans
        daily.list$V[[accept.counter]] <<- results.V$daily.samp
      }
      if (!is.null(D.obs)) {
        summary.stats$D[[accept.counter]] <<- results.D$samp.trans
        daily.list$D[[accept.counter]] <<- results.D$daily.samp
      }
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
    IS.weights[accept.counter] <<- single.IS.weight
    summary.stats$I[[accept.counter]] <<- results.I$samp.trans
    daily.list$I[[accept.counter]] <<- results.I$daily.samp
    if (!is.null(VI.obs)) {
      summary.stats$VI[[accept.counter]] <<- results.VI$samp.trans
      daily.list$VI[[accept.counter]] <<- results.VI$daily.samp
    }
    if (!is.null(RI.obs)) {
      summary.stats$RI[[accept.counter]] <<- results.RI$samp.trans
      daily.list$RI[[accept.counter]] <<- results.RI$daily.samp
    }
    if (!is.null(RVI.obs)) {
      summary.stats$RVI[[accept.counter]] <<- results.RVI$samp.trans
      daily.list$RVI[[accept.counter]] <<- results.RVI$daily.samp
    }
    if (!is.null(R.obs)) {
      summary.stats$R[[accept.counter]] <<- results.R$samp.trans
      daily.list$R[[accept.counter]] <<- results.R$daily.samp
    }
    if (!is.null(V.obs)) {
      summary.stats$V[[accept.counter]] <<- results.V$samp.trans
      daily.list$V[[accept.counter]] <<- results.V$daily.samp
    }
    if (!is.null(D.obs)) {
      summary.stats$D[[accept.counter]] <<- results.D$samp.trans
      daily.list$D[[accept.counter]] <<- results.D$daily.samp
    }
  }
  
  RE.ABC <- function() {
    # Algorithm variant performing rejection sampling with
    
    # Updates the counter of accepted particles and saves them together
    # with summary statistics and daily case counts into the pre-allocated 
    # matrices/lists. The global assignment operator '<<-' is used to avoid 
    # copying the parameters, which would be otherwise passed to the function
    accept.counter <<- accept.counter + 1
    accept.parts[accept.counter, ] <<- parts
    summary.stats$I[[accept.counter]] <<- results.I$samp.trans
    daily.list$I[[accept.counter]] <<- results.I$daily.samp
    if (!is.null(VI.obs)) {
      summary.stats$VI[[accept.counter]] <<- results.VI$samp.trans
      daily.list$VI[[accept.counter]] <<- results.VI$daily.samp
    }
    if (!is.null(RI.obs)) {
      summary.stats$RI[[accept.counter]] <<- results.RI$samp.trans
      daily.list$RI[[accept.counter]] <<- results.RI$daily.samp
    }
    if (!is.null(RVI.obs)) {
      summary.stats$RVI[[accept.counter]] <<- results.RVI$samp.trans
      daily.list$RVI[[accept.counter]] <<- results.RVI$daily.samp
    }
    if (!is.null(R.obs)) {
      summary.stats$R[[accept.counter]] <<- results.R$samp.trans
      daily.list$R[[accept.counter]] <<- results.R$daily.samp
    }
    if (!is.null(V.obs)) {
      summary.stats$V[[accept.counter]] <<- results.V$samp.trans
      daily.list$V[[accept.counter]] <<- results.V$daily.samp
    }
    if (!is.null(D.obs)) {
      summary.stats$D[[accept.counter]] <<- results.D$samp.trans
      daily.list$D[[accept.counter]] <<- results.D$daily.samp
    }
  }
  
  # Definitions of methods of data processing ==================================
  processing <- function (compartment, max.day.obs.inner, obs.vector,
                                 obs.trans.inner, transf.inner, dist.fun.inner) {
    compartment.times <- unlist(compartment)
    if (length(compartment.times) != 0) {
      daily.samp <- as.data.frame(table(floor(compartment.times), dnn = list("day")), 
                                  responseName = "samp")
      daily.samp$day <- as.numeric(levels(daily.samp$day))
      max.day.samp <- daily.samp$day[nrow(daily.samp)]
      
      # Makes a one-sided outer join of the data frame with the sampled case
      # counts and a data frame with only column "day" in order to add missing 
      # days with zero counts.
      if (max.day.samp != nrow(daily.samp) - 1) {
        daily.samp <- merge(data.frame(day = 0:max.day.samp), daily.samp, 
                            by = "day", all.x = TRUE)
        daily.samp$samp[is.na(daily.samp$samp)] <- 0
      }
      
      # If the sampled epidemic lasted for more days than the observed one, the
      # observed epidemic is expanded by these days with zero cases.
      if (max.day.samp > max.day.obs.inner) {
        obs.vector <- c(obs.vector, rep(0, max.day.samp - max.day.obs.inner))
        max.day.obs.inner <- max.day.samp
        # If the transformation function calculates one element of the vector on 
        # the basis of multiple other elements of the case count, the simplest 
        # way is to recalculate the the transformation of the whole vector of 
        # observed case counts
        obs.trans.inner <- transf.inner(obs.vector)
        
      } else if (max.day.samp < max.day.obs.inner) {
        daily.samp <- rbind(
          daily.samp,
          data.frame(day = (max.day.samp + 1):max.day.obs.inner, samp = 0)
        )
      }
 #     cat("Max day samp:", max.day.samp, "\n")
    } else {
      daily.samp <- data.frame(day = 0:max.day.obs.inner, samp.cases = 0)
#      cat("Max day samp:", 0, "\n")
    }
    
    # Transforms the case count.
    samp.trans <- transf.inner(daily.samp$samp)
    
    # Computes the distance.
    distance <- dist.fun.inner(obs.trans.inner, samp.trans)
    
    

    # cat("Max day obs:", max.day.obs.inner, "\n")
    # cat("Daily samp:\n")
    # print(daily.samp)
    # cat("Daily obs:\n")
    # print(obs.vector)
    # 
    # cat("Samp trans:", length(samp.trans), "\n")
    # cat("Obs trans:", length(obs.trans.inner), "\n")
    # cat("Distance:", distance, "\n")
    
    return(list(samp.trans = samp.trans, obs.trans = obs.trans.inner, 
                distance = distance, max.day.obs = max.day.obs.inner, 
                obs.vector = obs.vector, daily.samp = daily.samp))
  }
  
  # Input parameters check =====================================================
  
  # Stores the information, whether the user  wants to use a kernel function or
  # importance sampling, because it is tested later in multiple conditions
  use.kern.flag <- !is.null(kern)
  use.importance.flag <- !is.null(unlist(importance))
  
  int.param.check <- c(s0, i0, v0, r0, times, max.times)
  if (any(int.param.check < 0 | int.param.check %% 1 != 0)) {
    stop("Parameters 's0', 'i0', 'times' and 'max.times' 
          (and possibly 'v0', 's0') must be positive integer values.")
  }
  if (!use.importance.flag &&
      (!is.function(prior$lambda.samp) || !is.function(prior$mu.samp) ||
       !is.function(prior$nu.samp) || !is.function(prior$psi.samp) ||
       !is.function(prior$phi.samp) || !is.function(prior$kappa.samp) ||
       !is.function(prior$omega.samp))) {
    stop("Parameter 'prior' must be a properly named list.")
  }
  if (!is.function(kern) & use.kern.flag) {
    stop("Parameter 'kern' must be an appropriate functions.")
  }
  if (use.importance.flag && 
      (!is.function(importance$lambda.samp) || 
       !is.function(importance$mu.samp) || !is.function(importance$nu.samp) ||
       !is.function(importance$psi.samp) || !is.function(importance$phi.samp) ||
       !is.function(importance$kappa.samp) || 
       !is.function(importance$omega.samp))) {
    stop("Parameters 'importance$lambda' and 'importance$mu' must be functions.")
  }
  if (!is.numeric(I.obs) || sum(I.obs %% 1, na.rm = TRUE) != 0) {
    stop("Parameter 'I.obs' must be an integer valued vector.")
  }
  
  if (!is.list(transf)) {
    stop("Parameter 'transf' must be a named list containing appropriate 
         functions or NULL values.")
  }
  
  if (any(!(names(transf) %in% c("I", "VI", "RI", "RVI", "V", "R", "D")))) {
    stop("Names of the list 'trans' must be from the set 
         {'I', 'VI', 'RI', 'RVI', 'V', 'R', 'D'}")
  }
  
  if (!is.null(RVI.obs) && is.null(RI.obs) && is.null(VI.obs)) {
    warning("Reinfections with respect to vaccination ('RVI') supplied, but 
            either the counts of all reinfections('RI'), or the counts of all 
            infections after vaccination ('VI'), not supplied. 'RVI' will be ignored")
    RVI.obs <- NULL
    compart.weights["RVI"] <- NA
    transf$RVI <- NULL
    dist.metr$RVI <- NULL
  }
  
  # If transformation of case counts not specified, we define the transformation
  # function as an identity.
  
  if (is.null(transf$I)) {
    transf$I <- function (x) {return(x)}
  }
  
  # We possibly transform other counts. If their transformation is 
  # not specified, we use the identity function
  if (is.null(R.obs) && !is.null(transf$R)) {
    warning("Observed recoveries not supplied. The transformation for the
                recoveries will be ignored.")
  } else if (!is.null(R.obs) && is.null(transf$R)) {
    transf$R <- function (x) {return(x)}
  }
  if (is.null(V.obs) && !is.null(transf$V)) {
    warning("Observed vaccinations not supplied. The transformation for the
                vaccinations will be ignored.")
  } else if (!is.null(V.obs) && is.null(transf$V)) {
    transf$V <- function (x) {return(x)}
  }
  if (is.null(D.obs) && !is.null(transf$D)) {
    warning("Observed deaths not supplied. The transformation for the
                deaths will be ignored.")
  } else if (!is.null(D.obs) && is.null(transf$D)) {
    transf$D <- function (x) {return(x)}
  }
  if (is.null(VI.obs) && !is.null(transf$VI)) {
    warning("Observed infections with respect to the vaccination not supplied. 
    The transformation for the infections after vaccination will be ignored.")
  } else if (!is.null(VI.obs) && is.null(transf$VI)) {
    transf$VI <- function (x) {return(x)}
  }
  if (is.null(RI.obs) && !is.null(transf$RI)) {
    warning("Observed infections with respect to the recoveries not supplied. 
    The transformation for the infections after recovery will be ignored.")
  } else if (!is.null(RI.obs) && is.null(transf$RI)) {
    transf$RI <- function (x) {return(x)}
  }
  if (is.null(RVI.obs) && !is.null(transf$RVI)) {
    warning("Observed infections with respect to the recoveries not supplied. 
    The transformation for the infections after recovery will be ignored.")
  } else if (!is.null(RVI.obs) && is.null(transf$RVI)) {
    transf$RVI <- function (x) {return(x)}
  }
  
  # If function measuring the distance of summary statistics not specified,
  # we define it as an euclidean distance. If it is specified, we define it as
  # a bilinear form of the observed and simulated summary statistic. Attention!
  # It is not checked, whether the matrix is positive definite!
  
  if (!is.list(dist.metr)) {
    stop("Parameter 'dist.metr' must be a named list containing appropriate 
         matrices, the 'euclidean' strings, or NULL values.")
  }
  
  if (any(!(names(dist.metr) %in% c("I", "VI", "RI", "RVI", "V", "R", "D")))) {
    stop("Names of the list 'dist.metr' must be from the set 
         {'I', 'VI', 'RI', 'RVI', 'V', 'R', 'D'}")
  }
  
  if (is.matrix(dist.metr$I)) {
    if (any(dist.metr$I != t(dist.metr$I))) {
      stop("The matrix in the element 'I' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      I = function (x, y) {sqrt(t(x - y) %*% dist.metr$I %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$I) && dist.metr$I != "euclidean") {
      warning("Element 'I' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- list(I = function (x, y) {sqrt(sum((x - y) ^ 2))})
  }
  
  if (is.matrix(dist.metr$V)) {
    if (any(dist.metr$V != t(dist.metr$V))) {
      stop("The matrix in the element 'V' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      V = function (x, y) {sqrt(t(x - y) %*% dist.metr$V %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$V) && dist.metr$V != "euclidean") {
      warning("Element 'V' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(V = function (x, y) {sqrt(sum((x - y) ^ 2))}))
  }
  
  if (is.matrix(dist.metr$R)) {
    if (any(dist.metr$R != t(dist.metr$R))) {
      stop("The matrix in the element 'R' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      R = function (x, y) {sqrt(t(x - y) %*% dist.metr$R %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$R) && dist.metr$R != "euclidean") {
      warning("Element 'R' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(R = function (x, y) {sqrt(sum((x - y) ^ 2))}))
  }
  
  if (is.matrix(dist.metr$D)) {
    if (any(dist.metr$D != t(dist.metr$D))) {
      stop("The matrix in the element 'D' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      D = function (x, y) {sqrt(t(x - y) %*% dist.metr$D %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$D) && dist.metr$D != "euclidean") {
      warning("Element 'D' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(D = function (x, y) {sqrt(sum((x - y) ^ 2))}))
  }
  
  if (is.matrix(dist.metr$VI)) {
    if (any(dist.metr$VI != t(dist.metr$VI))) {
      stop("The matrix in the element 'VI' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      VI = function (x, y) {sqrt(t(x - y) %*% dist.metr$VI %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$VI) && dist.metr$VI != "euclidean") {
      warning("Element 'VI' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(VI = function (x, y) {sqrt(sum((x - y) ^ 2))}))
  }
  if (is.matrix(dist.metr$RI)) {
    if (any(dist.metr$RI != t(dist.metr$RI))) {
      stop("The matrix in the element 'RI' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      RI = function (x, y) {sqrt(t(x - y) %*% dist.metr$RI %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$RI) && dist.metr$RI != "euclidean") {
      warning("Element 'RI' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(RI = function (x, y) {sqrt(sum((x - y) ^ 2))}))
  }
  if (is.matrix(dist.metr$RVI)) {
    if (any(dist.metr$RVI != t(dist.metr$RVI))) {
      stop("The matrix in the element 'RVI' of the 'dist.metr' list is not symmetrical.")
    }
    
    dist.fun <- list(
      RVI = function (x, y) {sqrt(t(x - y) %*% dist.metr$RVI %*% (x - y))}
      )
  } else {
    if (!is.null(dist.metr$RVI) && dist.metr$RVI != "euclidean") {
      warning("Element 'RVI' of the 'dist.metr' list must be a symmetrical matrix,
      or the 'euclidean' string. The euclidean distance will be used." )
    }
    dist.fun <- c(dist.fun, list(RVI = function (x, y) {sqrt(sum((x - y) ^ 2))}))
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
    nu.sampler <- importance$nu.samp
    psi.sampler <- importance$psi.samp
    phi.sampler <- importance$phi.samp
    kappa.sampler <- importance$kappa.samp
    omega.sampler <- importance$omega.samp
    lambda.sampler.args <- as.list(c(n = 1, importance.params$lambda))
    mu.sampler.args <- as.list(c(n = 1, importance.params$mu))
    nu.sampler.args <- as.list(c(n = 1, importance.params$nu))
    psi.sampler.args <- as.list(c(n = 1, importance.params$psi))
    phi.sampler.args <- as.list(c(n = 1, importance.params$phi))
    kappa.sampler.args <- as.list(c(n = 1, importance.params$kappa))
    omega.sampler.args <- as.list(c(n = 1, importance.params$omega))
    
    # Prepares the argument lists for evaluating the prior and importance
    # density functions
    lambda.prior.dens.args <- as.list(c(x = NA, prior.params$lambda))
    mu.prior.dens.args <- as.list(c(x = NA, prior.params$mu))
    nu.prior.dens.args <- as.list(c(x = NA, prior.params$nu))
    psi.prior.dens.args <- as.list(c(x = NA, prior.params$psi))
    phi.prior.dens.args <- as.list(c(x = NA, prior.params$phi))
    kappa.prior.dens.args <- as.list(c(x = NA, prior.params$kappa))
    omega.prior.dens.args <- as.list(c(x = NA, prior.params$omega))
    lambda.importance.dens.args <- as.list(c(x = NA, importance.params$lambda))
    mu.importance.dens.args <- as.list(c(x = NA, importance.params$mu))
    nu.importance.dens.args <- as.list(c(x = NA, importance.params$nu))
    psi.importance.dens.args <- as.list(c(x = NA, importance.params$psi))
    phi.importance.dens.args <- as.list(c(x = NA, importance.params$phi))
    kappa.importance.dens.args <- as.list(c(x = NA, importance.params$kappa))
    omega.importance.dens.args <- as.list(c(x = NA, importance.params$omega))
    
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
    nu.sampler <- prior$nu
    nu.sampler.args <- as.list(c(n = 1, prior.params$nu))
    psi.sampler <- prior$psi
    psi.sampler.args <- as.list(c(n = 1, prior.params$psi))
    phi.sampler <- prior$phi
    phi.sampler.args <- as.list(c(n = 1, prior.params$phi))
    kappa.sampler <- prior$kappa
    kappa.sampler.args <- as.list(c(n = 1, prior.params$kappa))
    omega.sampler <- prior$omega
    omega.sampler.args <- as.list(c(n = 1, prior.params$omega))
    
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
  
  # Sets the variable, which indicates, whether to sample an epidemic from the
  # model. It stays always TRUE for the rejection sampling variants. For the 
  # importance sampling, it is set as False every time, there is proposed
  # with 0 prior probability.
  do.simulate <- TRUE
  
  # Finds out the length of the observed epidemic and transforms it
  max.day.obs <- c(I = length(I.obs) - 1)
  obs.trans <- list(I = transf(I.obs))
  
  # Allocates the output parameters
  accept.counter <- 0
  total.counter <- 0
  accept.parts <- matrix(data = NA, nrow = times, ncol = 7,
                         dimnames = list(NULL, c("lambda", "mu", "nu", "psi",
                                                 "phi", "kappa", "omega")))
  # If there is a transformation, summary statistics and the case counts must be 
  # stored separately, because their dimension may differ due to the 
  # transformation. If there is no transformation, the list of summary 
  # statistics contains only vectors of daily case counts, whereas the list of 
  # daily case counts contains data frames with additional column specifying the 
  # day.
  summary.stats <- list(I = vector(mode = "list", length = times))
  
  if (!is.null(VI.obs)) {
    if (is.na(compart.weights["VI"])) {
      warning("Weight of distance for the infections after the vaccination
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["VI"] <- 0
    }
    summary.stats <-c(summary.stats, list(VI = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, VI = length(VI.obs) - 1)
    obs.trans <- c(obs.trans, list(VI = transf$VI(VI.obs)))
  }
  if (!is.null(RI.obs)) {
    if (is.na(compart.weights["RI"])) {
      warning("Weight of distance for the infections after the recovery
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["RI"] <- 0
    }
    summary.stats <- c(summary.stats, list(RI = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, RI = length(RI.obs) - 1)
    obs.trans <- c(obs.trans, list(RI = transf$RI(RI.obs)))
  }
  if (!is.null(RVI.obs)) {
    if (is.na(compart.weights["RVI"])) {
      warning("Weight of distance for the infections after the recovery
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["RVI"] <- 0
    }
    summary.stats <- c(summary.stats, list(RVI = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, RVI = length(RVI.obs) - 1)
    obs.trans <- c(obs.trans, list(RVI = transf$RVI(RVI.obs)))
  }
  if (!is.null(V.obs)) {
    if (is.na(compart.weights["V"])) {
      warning("Weight of distance for the V compartment 
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["V"] <- 0
    }
    summary.stats <- c(summary.stats, list(V = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, V = length(V.obs) - 1)
    obs.trans <- c(obs.trans, list(V = transf$V(V.obs)))
  }
  if (!is.null(R.obs)) {
    if (is.na(compart.weights["R"])) {
      warning("Weight of distance for the R compartment 
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["R"] <- 0
    }
    summary.stats <- c(summary.stats, list(R = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, R = length(R.obs) - 1)
    obs.trans <- c(obs.trans, list(R = transf$R(R.obs)))
  }
  if (!is.null(D.obs)) {
    if (is.na(compart.weights["D"])) {
      warning("Weight of distance for the D compartment 
              (variable 'compart.weights') not specified. It will be set as 0.")
      compart.weights["D"] <- 0
    }
    summary.stats <- c(summary.stats, list(D = summary.stats[[1]]))
    max.day.obs <- c(max.day.obs, D = length(D.obs) - 1)
    obs.trans <- c(obs.trans, list(D = transf$D(D.obs)))
  }
  
  daily.list <- summary.stats

  # Ensures that the compartment weights are in the right order
  compart.weights <- compart.weights[c("I", "VI", "RI", "RVI", "V", "R", "D")]
  compart.weights <- compart.weights[!is.na(compart.weights)]
  
  # Allocates a local variable to store proposed particles lambda and mu
  parts <- rep(NA, 7)
  
  # Main cycle of the algorithm ================================================
  
  while (accept.counter < times & total.counter < max.times) {
    # Proposing particles lambda and mu
    # lambda stored at the first position, mu at the second position
    parts[1] <- do.call(lambda.sampler, args = lambda.sampler.args) 
    parts[2] <- do.call(mu.sampler, args = mu.sampler.args)
    parts[3] <- do.call(nu.sampler, args = nu.sampler.args)
    parts[4] <- do.call(psi.sampler, args = psi.sampler.args)
    parts[5] <- do.call(phi.sampler, args = phi.sampler.args)
    parts[6] <- do.call(kappa.sampler, args = kappa.sampler.args)
    parts[7] <- do.call(omega.sampler, args = omega.sampler.args)
    
    # If the importance weights of the proposed particle are 0 (i.e. the prior 
    # distribution is 0, we do not need to generate a sample epidemic)
    if (use.importance.flag) {
      
      # Completes the argument lists for evaluation of the prior and the 
      # importance density
      lambda.prior.dens.args[[1]] <- parts[1]
      mu.prior.dens.args[[1]] <- parts[2]
      nu.prior.dens.args[[1]] <- parts[3]
      psi.prior.dens.args[[1]] <- parts[4]
      phi.prior.dens.args[[1]] <- parts[5]
      kappa.prior.dens.args[[1]] <- parts[6]
      omega.prior.dens.args[[1]] <- parts[7]
      lambda.importance.dens.args[[1]] <- parts[1]
      mu.importance.dens.args[[1]] <- parts[2]
      nu.importance.dens.args[[1]] <- parts[3]
      psi.importance.dens.args[[1]] <- parts[4]
      phi.importance.dens.args[[1]] <- parts[5]
      kappa.importance.dens.args[[1]] <- parts[6]
      omega.importance.dens.args[[1]] <- parts[7]
      
      single.IS.weight <- 
        do.call(prior$lambda.dens, args = lambda.prior.dens.args) /
        do.call(importance$lambda.dens, args = lambda.importance.dens.args) *
        do.call(prior$mu.dens, args = mu.prior.dens.args) /
        do.call(importance$mu.dens, args = mu.importance.dens.args) *
        do.call(prior$nu.dens, args = nu.prior.dens.args) /
        do.call(importance$nu.dens, args = nu.importance.dens.args) *
        do.call(prior$psi.dens, args = psi.prior.dens.args) /
        do.call(importance$psi.dens, args = psi.importance.dens.args) *
        do.call(prior$phi.dens, args = phi.prior.dens.args) /
        do.call(importance$phi.dens, args = phi.importance.dens.args) *
        do.call(prior$kappa.dens, args = kappa.prior.dens.args) /
        do.call(importance$kappa.dens, args = kappa.importance.dens.args) *
        do.call(prior$omega.dens, args = omega.prior.dens.args) /
        do.call(importance$omega.dens, args = omega.importance.dens.args)
      do.simulate <- as.logical(single.IS.weight)
    }
    
    # Samples from the model. If there were more infections, than specified in 
    # the 'max.infections', the whole sample is discarded and new iteration 
    # starts.
    if (do.simulate) {
      epi.samp <- gener.sveird(
        s0 = s0, i0 = i0, v0 = v0, r0 = r0, max.duration = max.duration,
        lambda = parts[1], mu = parts[2], nu = parts[3], psi = parts[4], 
        phi = parts[5], kappa = parts[6], omega = parts[7], 
        max.infections = max.infections
      )
      
      if (!epi.samp$stopped) {
        # If there were less infections than specified in the 'max.infections',
        # the output of the epidemic model is processed. Continuous infection 
        # times are converted into case counts.
        
        results.I <- processing(
          compartment = epi.samp$I, max.day.obs.inner =  max.day.obs["I"], 
          obs.vector =  I.obs, obs.trans.inner = obs.trans[["I"]],
          transf.inner = transf$I, dist.fun.inner = dist.fun$I
        )
        obs.trans$I <- results.I$obs.trans
        max.day.obs["I"] <- results.I$max.day.obs
        I.obs <- results.I$obs.vector
        distance.total <- c(I = results.I$distance)
        if (!is.null(VI.obs)) {
          
          # We extract only those infection times, which happened after the
          # vaccination.
          prev.vac.times <- unlist(
            lapply(epi.samp$V, function (x) {ifelse(length(x) == 0, Inf, x[length(x)])}))
          compart.after <- mapply(function (x, times) {x[x > times]},
                                  epi.samp$I, prev.vac.times)
          # if (length(compart.after) != length(prev.vac.times)) {
          #   print(length(epi.samp$V))
          #   print(length(lapply(epi.samp$V, function (x) {ifelse(length(x) == 0, Inf, x[length(x)])})))
          #   print(prev.vac.times)
          #   cat("Length of comp (VI):", length(compart.after), "\n")
          #   cat("Length of times:", length(prev.vac.times), "\n")
          # }
          


          results.VI <- processing(
            compartment = compart.after, max.day.obs.inner =  max.day.obs["VI"],
            obs.vector =  VI.obs, obs.trans.inner = obs.trans[["VI"]],
            transf.inner = transf$VI, dist.fun.inner = dist.fun$VI
          )
          obs.trans$VI <- results.VI$obs.trans
          max.day.obs["VI"] <- results.VI$max.day.obs
          VI.obs <- results.VI$obs.vector
          distance.total <- c(distance.total, VI = results.VI$distance)
        }
        if (!is.null(RI.obs)) {
          compart.after <- lapply(epi.samp$I, function (x, times) {x[-length(x)]})

          # We extract only the reinfections (second and later infections).
          results.RI <- processing(
            compartment = compart.after, max.day.obs.inner =  max.day.obs["RI"],
            obs.vector =  RI.obs, obs.trans.inner = obs.trans[["RI"]],
            transf.inner = transf$RI, dist.fun.inner = dist.fun$RI
          )
          obs.trans$RI <- results.RI$obs.trans
          max.day.obs["RI"] <- results.RI$max.day.obs
          RI.obs <- results.RI$obs.vector
          distance.total <- c(distance.total, RI = results.RI$distance)
        }
        if (!is.null(RVI.obs)) {

          # In 'compart.after', there is still the list of reinfection times and
          # in the 'prev.vac.times', there is the border, after which the
          # individuals are vaccinated.
          compart.after <- mapply(function (x, times) {x[x > times]},
                                  compart.after, times = prev.vac.times)
          if (length(compart.after) != length(prev.vac.times)) {
            cat("Length of comp (RVI):", length(compart.after), "\n")
            cat("Length of times:", length(prev.vac.times), "\n")
          }
          
          results.RVI <- processing(
            compartment = compart.after, max.day.obs.inner =  max.day.obs["RVI"],
            obs.vector =  RVI.obs, obs.trans.inner = obs.trans[["RVI"]],
            transf.inner = transf$RVI, dist.fun.inner = dist.fun$RVI
          )
          obs.trans$RVI <- results.RVI$obs.trans
          max.day.obs["RVI"] <- results.RVI$max.day.obs
          RVI.obs <- results.RVI$obs.vector
          distance.total <- c(distance.total, RVI = results.RVI$distance)
        }
        if (!is.null(V.obs)) {
          results.V <- processing(
            compartment = epi.samp$V, max.day.obs.inner =  max.day.obs["V"], 
            obs.vector =  V.obs, obs.trans.inner = obs.trans[["V"]],
            transf.inner = transf$V, dist.fun.inner = dist.fun$V
            )
          obs.trans$V <- results.V$obs.trans
          max.day.obs["V"] <- results.V$max.day.obs
          V.obs <- results.V$obs.vector
          distance.total <- c(distance.total, V = results.V$distance)
        }
        if (!is.null(R.obs)) {
          results.R <- processing(
            compartment = epi.samp$R, max.day.obs.inner =  max.day.obs["R"], 
            obs.vector =  R.obs, obs.trans.inner = obs.trans[["R"]],
            transf.inner = transf$R, dist.fun.inner = dist.fun$R
            )
          obs.trans$R <- results.R$obs.trans
          max.day.obs["R"] <- results.R$max.day.obs
          R.obs <- results.R$obs.vector
          distance.total <- c(distance.total, R = results.R$distance)
        }
        if (!is.null(D.obs)) {
          if (all(is.na(epi.samp$D))) {  # Can get rid of this condition by
                                         # rewriting the 'gener.sveird' function
                                         # to have NULL for the empty individuals
                                         # instead of NA.
            epi.samp$D <- list()
          }
          results.D <- processing(
            compartment = epi.samp$D, max.day.obs.inner =  max.day.obs["D"], 
            obs.vector =  D.obs, obs.trans.inner = obs.trans[["D"]],
            transf.inner = transf$D, dist.fun.inner = dist.fun$D
            )
          obs.trans$D <- results.D$obs.trans
          max.day.obs["D"] <- results.D$max.day.obs
          D.obs <- results.D$obs.vector
          distance.total <- c(distance.total, D = results.D$distance)
        }
        
        distance.total <- sum(distance.total * compart.weights)
        print(distance.total) # delete before submitting
        
        # Accepts or rejects proposed lambda and mu according to the algorithm 
        # variant
        if (distance.total < tolerance) {
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

  # Returning values ===========================================================
  
  # Prepares the return value list
  return.values <- list(
    accept.parts = accept.parts,
    summary.stats = summary.stats, 
    n.iterations = total.counter,
    n.accepted = accept.counter
  )
  
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
  return.values$samp.counts <- daily.list
  
  return(return.values)
}