#' @title Generates a course of an SEIR epidemic.
#'
#' @description Generates a course of an SEIR epidemic. Returns vector of
#' infection, exposition and removal times.
#' 
#' @import base
#' @import stats
#'
#' @param n.pop the total population size excluding the initial infectious
#'     individuals.
#' @param m the number of initial infectious individuals at the beginning of the
#'     epidemic.
#' @param lambda the rate of the Poisson process, which generates the points
#'     when the infection can occur.
#' @param mu the rate of the exponential distribution of the infectious periods.
#' @param nu the rate of the exponential distribution of the latent periods.
#' @param max.infections how many infections can occur before the epidemic
#'     stops. Use especially for the ABC algorithm  when the population is large
#'     to reduce the amount of time spent on generating a single epidemic, which
#'     is too different from the observed one anyway. The value of maximal
#'     infections must be set manually regarding the data.
#'
#' @returns A list of three numerical vectors \code{E}, \code{I} and \code{R}
#'     and a logical value \code{stopped}. Initial infectious individuals are
#'     included in the first \code{m} positions. Hence all three vectors are of
#'     length \code{m} + \code{n.pop}.
#'     \itemize{
#'     \item \code{I} the infection times for each individual
#'     \item \code{E} the exposition times for each individual
#'     \item \code{R} the recovery times for each individual
#'     \item \code{stopped} logical value, if the epidemic was stopped after
#'     exceeding the \code{max.infections} number
#'     }
#'
#' @examples
#' gener.seir(n.pop = 1000, m = 5, lambda = 0.6, mu = 0.3, nu = 1)

gener.seir <- function (lambda, mu, nu, n.pop = 100, m = 1,
                        max.infections = n.pop) {
  
  # Generates the infectious periods of the initial infectious individuals.
  infect.per <- rexp(m, mu)
  
  # Allocates the infection and exposition times, which are 0 for every 
  # initially infectious individual and NA for everyone susceptible.
  infect.time <- rep(c(0, NA), times = c(m, n.pop))
  exposed.time <- infect.time
  
  # Allocates the recovery times, which are 0 + infectious period for
  # every initially infectious individual and NA for everyone susceptible.
  recov.time <- c(infect.per, rep(NA, n.pop))
  
  # Generates the list of the Poisson process points for the initial infectious 
  # individuals. In the following algorithm, always the earliest point will be
  # selected as the infection time.
  contact.times <- lapply(infect.per, poiss.proc.points, lambda = lambda)
  
  # Starts measuring the time of events and sets the counter of infections.
  act.time <- 0
  infections <- 0
  
  # Keeps track of the latent and active cases. At the beginning, there are no
  # latent cases and active cases are those initially infectious.
  # Active cases are those, who are infectious and can potentially infect
  # someone else, i.e. their element in the list 'contact.times' is non-empty.
  still.active <- unlist(lapply(contact.times, length)) != 0
  contact.times <- contact.times[still.active]
  
  # The epidemic continues if there are at least one active or latent case, the 
  # maximum number of infections does not exceed the 'max.infections' parameter
  # and at least one individual has not been in the Infectious compartment yet.
  while (length(contact.times) != 0 && infections <= max.infections) {
    
    # The individual, who has the earliest event of the Poisson process
    # potentially infects another individual.
    
    carrier <- which.min(lapply(contact.times, FUN = choose.first))
    
    # Updates the actual time.
    act.time <- contact.times[[carrier]][1]
    # Chooses a new infected. If still susceptible, the infection was
    # successful, new points of the Poisson process are generated and the
    # times of the exposition, infectiousness and recovery are set.
    
    new.infected <- sample(1:(n.pop + m), size = 1)
    if (is.na(infect.time[new.infected])) {
      infections <- infections + 1
      exposed.time[new.infected] <- act.time
      single.latent.per <- rexp(1, nu)
      single.infect.per <- rexp(1, mu)
      infect.time[new.infected] <- act.time + single.latent.per
      recov.time[new.infected] <-
        infect.time[new.infected] + single.infect.per
      # New Poisson process points are generated. If there is none of them,
      # the empty vector is ignored.
      new.jumps <- poiss.proc.points(lambda, single.infect.per)
      if (length(new.jumps) != 0) {
        contact.times<- c(contact.times, 
                          list(new.jumps + act.time + single.latent.per))
      }
    }
    
    # The 'used' point of the Poisson process is removed. If no contact
    # remains, the corresponding individual is removed from the active cases.
    contact.times[[carrier]] <- contact.times[[carrier]][-1]
    if (length(contact.times[[carrier]]) == 0) {
      contact.times <- contact.times[-carrier]
    }
  }
  return(list(E = exposed.time, I = infect.time, R = recov.time,
              stopped = infections > max.infections))
}
