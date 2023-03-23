#' @title Generate an SIR epidemic
#'
#' @description Generates a course of an SIR epidemic. Returns vector of
#' infection and removal times.
#'
#' @import base
#' @import stats
#'
#' @param s0 the total population size excluding the initial infectious
#'     individuals.
#' @param m the number of initial infectious individuals at the beginning of the
#'     epidemic.
#' @param lambda the rate of the Poisson process, which generates the points
#'     when the infection can occur.
#' @param mu the rate of the exponential distribution of the infectious periods.
#' @param max.infections how many infections can occur before the epidemic
#'     stops. Use especially for the ABC algorithm  when the population is large
#'     to reduce the amount of time spent on generating a single epidemic, which
#'     is too different from the observed one anyway. The value of maximal
#'     infections must be set manually regarding the data.
#'
#' @returns A list of two numerical vectors \code{I} and \code{R} and a logical
#'     value \code{stopped}. Initial infectious individuals are included in the
#'     first \code{i0} positions. Hence both vectors are of length
#'     \code{i0} + \code{s0}.
#'     \itemize{
#'     \item \code{I} the infection times for each individual
#'     \item \code{R} the recovery times for each individual
#'     \item \code{stopped} logical value, if the epidemic was stopped after
#'     exceeding the \code{max.infections} number
#'     }
#'
#' @examples
#' set.seed(80)
#' gener.sir(s0 = 1000, i0 = 5, lambda = 0.6, mu = 0.5)

gener.sir <- function (lambda, mu, s0 = 100, i0 = 1, max.infections = s0) {

  # Generates the infectious periods
  infect.per <- rexp(i0, mu)
  
  # Allocates the infection times, which are 0 for every initially infectious
  # individual and NA for everyone susceptible.
  infect.time <- rep(c(0, NA), times = c(i0, s0))
  
  # Allocates the recovery times, which are 0 + infectious period for
  # every initially infectious individual and NA for everyone susceptible.
  recov.time <- c(infect.per, rep(NA, s0))
  
  # Generates the list of the Poisson process points for the initial infectious 
  # individuals. In the following algorithm, always the earliest point will be
  # selected as the infection time.
  contact.times <- lapply(infect.per, poiss.proc.points, lambda = lambda)
  
  # Starts measuring the time of events and sets the counter of infections.
  act.time <- 0
  infections <- 0
  
  # Some of the individuals might have not made any contact. Therefore we remove
  # them.
  still.active <- unlist(lapply(contact.times, length)) != 0
  contact.times <- contact.times[still.active]
  
  # The epidemic continues if there are at least one active case, the maximum 
  # number of infections does not exceed the 'max.infections' parameter and
  # at least one individual has not been in the Infectious compartment yet.
  while (length(contact.times) != 0 && infections <= max.infections) {
    
    # The individual, who has the earliest event of the Poisson process
    # potentially infects another individual.
    
    carrier <- which.min(lapply(contact.times, FUN = choose.first))
    
    # Updates the actual time.
    act.time <- contact.times[[carrier]][1]
    
    # Chooses a new infected. If still susceptible, the infection was successful
    # and new points of the Poisson process are generated.
    new.infected <- sample(1:(s0 + i0), size = 1)
    if (is.na(infect.time[new.infected])) {
      infections <- infections + 1
      single.infect.per <- rexp(1, mu)
      infect.time[new.infected] <- act.time
      recov.time[new.infected] <- act.time + single.infect.per
      new.jumps <- poiss.proc.points(lambda, single.infect.per)
      
      # The newly infected individual is added into the active cases only if
      # he makes at least one contact. Otherwise, he is ignored.
      if (length(new.jumps) != 0) {
        contact.times <- c(
          contact.times,
          list(new.jumps + act.time)
        )
      }
    }
    
    # The 'used' point of the Poisson process is removed. If no contact remains,
    # the corresponding individual is removed from the active cases.
    contact.times[[carrier]] <- contact.times[[carrier]][-1]
    if (length(contact.times[[carrier]]) == 0) {
      contact.times <- contact.times[-carrier]
    }
  }
  
  return(list(I = infect.time, R = recov.time,
              stopped = infections > max.infections))
}