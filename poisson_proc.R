choose.first <- function (x) {x[1]}

poiss.proc.points <- function (lambda, duration) {
  #' @title Generate jump times of a Poisson process
  #'
  #' @description The function generates a sequence of jump times of a Poisson
  #' process of rate \code{lambda}.
  #'
  #' @import base
  #'
  #' @param lambda the rate of the Poisson process.
  #' @param duration a numeric value specifying the maximal duration of the
  #'   process.
  #'
  #' @returns a numeric vector of jump times of the Poisson process. If no event
  #'     occurred until the \code{duration}, a vector of length 0 is returned.
  #'
  #' @examples poiss.proc.points(lambda = 1, duration = 2.5)
  
  # Allocates the vector of Poisson process jumps times.
  poiss.proc <- c()
  
  # Generates the exponential times of the Poisson process
  while (sum(poiss.proc) < duration) {
    poiss.proc <- c(poiss.proc,
                    rexp(1, rate = lambda))
  }
  
  # Returns the vector without the last element, since the last event occurred
  # after the maximal duration
  return(cumsum(poiss.proc[-length(poiss.proc)]))
}