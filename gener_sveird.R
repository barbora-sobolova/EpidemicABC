#' @title Generates a course of an SVEIRD epidemic.
#' 
#' @description Generates a course of an SVEIRD epidemic. Returns lists of
#' infection, exposition, removal, death, vaccination and renewed 
#' susceptibility times.
#'
#' @param lambda the rate of the Poisson process, which generates the points
#'     when the infection can occur.
#' @param mu the rate of the exponential distribution of the time from
#'     the start of the infectiousness to the time of recovery.
#' @param nu the rate of the exponential distribution of the latent periods.
#' @param psi the rate of the exponential distribution of the time before
#'     the first vaccine dose.
#' @param omega the rate of the exponential distribution of the time, after
#'     which the immunity gained by the recovery from the disease fades out
#'     and the individual becomes susceptible again.
#' @param kappa the rate of the exponential distribution of the time from
#'     the start of the infectiousness to the time of death.
#' @param max.duration the maximum duration of the epidemic, which must be set
#'     manually, since in this epidemic model does not guarantee that the
#'     algorithm ends in a reasonable time.
#' @param s0 the initial susceptible population size.
#' @param i0 the number of initial infectious individuals.
#' @param v0 the number of initial vaccinated individuals.
#' @param r0 the number of initial recovered individuals.
#' @param max.infections how many infections can occur before the epidemic
#'     stops. Use especially for the ABC algorithm  when the population is large
#'     to reduce the amount of time spent on generating a single epidemic, which
#'     is too different from the observed one anyway. The value of maximal
#'     infections must be set manually regarding the data.
#' @param all.comps indicator, whether the outcome should include entry times
#'     of all compartments. If \code{FALSE}, only the 'E', 'I', 'R' and 'D'
#'     compartments are included.
#'
#' @returns A list of three lists \code{I}, \code{V}, \code{R}, one vector 
#'     \code{D}, one numeric value \code{duration} one character vector and a 
#'     logical value \code{stopped}. If \code{all.comps = TRUE}, additional 
#'     lists \code{S} and \code{V} are included. Initial infectious individuals
#'     are included in the first \code{i0} positions of the lists/vector. Hence 
#'     all three/five lists and both vectors are of length 
#'     \code{i0} + \code{s0}. The elements in the six lists are in the 
#'     decreasing order for the sake of convenience of the algorithm.
#'     \itemize{
#'     \item \code{S} the renewed susceptibility times for each individual
#'     \item \code{E} the exposition times for each individual
#'     \item \code{I} the infection times for each individual
#'     \item \code{V} the vaccination times for each individual
#'     \item \code{R} the recovery times for each individual
#'     \item \code{D} the death times for each individual
#'     \item \code{stopped} logical value, if the epidemic was stopped after
#'     exceeding the \code{max.infections} number
#'     \item \code{duration} the duration of the epidemic, which is less or 
#'     equal to the \code{max.duration} number.
#'     }
#'
#' @examples
#' gener.sveird(lambda = 2, mu = 0.5, nu = 0.25, psi = 0.2, omega = 0.05,
#'              phi = 0.1, max.duration = 100, kappa = 0.1, s0 = 1000,
#'              i0 = 1, v0 = 5, r0 = 10, all.comps = TRUE)

gener.sveird <- function (lambda, mu, nu, psi, phi, omega, kappa, max.duration, 
                          s0 = 100, i0 = 1, v0 = 0, r0 = 0, 
                          max.infections = 1e4, all.comps = TRUE) {
  
  # Functions resolving the contacts and infection =============================
  
  # These functions resolve the process of contact and passing on the infection
  # between an infectious individual and an individual in different states.
  # They are implemented in the way that they work with variables, which are
  # global in the environment corresponding to the gener.sveird() function and
  # they also use the global assignment operator '<<-'. A better implementation,
  # when the global variables would be passed to the functions as the formal
  # parameters, and the benchmarking, whether this implementation is faster,
  # is left as future work.
  
  infect.fun <- function () {
    # A successful infection occurs. The newly infected individual is assigned 
    # the "E" state, the number of infections increases. The newly infected
    # individual is assigned the exposed time infection time and the 
    # recovery/death time. New points of the Poisson process are generated.
    
    state[new.infected] <<- "E"
    infections <<- infections + 1
    single.latent.per <- rexp(1, nu)
    
    # Recovery or death may occur. We choose the earlier event.
    single.infect.per.toR <- rexp(1, mu)
    single.infect.per.toD <- rexp(1, kappa)
    
    
    if (single.infect.per.toR < single.infect.per.toD) {
      # Recovery first
      recov.time[[new.infected]] <<- c(
        act.time + single.latent.per + single.infect.per.toR,
        recov.time[[new.infected]]
      )
      removal.time[new.infected] <<- recov.time[[new.infected]][1]
      new.jumps <- poiss.proc.points(lambda, single.infect.per.toR)
    } else {
      # Death first
      death.time[new.infected] <<- 
        act.time + single.latent.per + single.infect.per.toD
      removal.time[new.infected] <<- death.time[new.infected]
      new.jumps <- poiss.proc.points(lambda, single.infect.per.toD)
    }
    
    # If the newly generated Poisson process points are not an empty vector, we
    # add the newly infected individual into the 'latent.cases' set, otherwise
    # we ignore him.
    if (length(new.jumps) != 0) {
      contact.times.latent <<- c(
        contact.times.latent,
        list(new.jumps + act.time +
               single.latent.per)
      )
      latent.cases <<- c(latent.cases, new.infected)
    }
    
    exposed.time[[new.infected]] <<- c(act.time, exposed.time[[new.infected]])
    infect.time[[new.infected]] <<- c(act.time + single.latent.per, 
                                      infect.time[[new.infected]])
  }
  
  Scontact <- function () {
    # The contact was made between an infectious and a susceptible individual.
    # The susceptible individual might have already been vaccinated. In the 
    # following, we generate the candidate vaccination time and compare it 
    # with the possible infection time.
    
    # The candidate first dose time is the last time entering the susceptible
    # state (which is 0 for those, susceptible from the beginning) plus the
    # exponentially distributed time with parameter 'psi'.
    vac.time.cand <- 
      ifelse(is.null(suscep.time[[new.infected]][1]), 
             0, suscep.time[[new.infected]][1]) + rexp(1, psi)
    
    if (contact.times[[active.min]][1] < vac.time.cand) {
      # The infection occurred earlier than the vaccine dose. The candidate 
      # vaccination time is ignored and the individual is infected
      # (enters the latent period), which corresponds to the real-world
      # situation, when an infection would usually result in postponing the 
      # vaccination.
      
      # Updates the actual time.
      act.time <<- contact.times[[active.min]][1]
      
      infect.fun()
    } else {
      # The vaccine dose was given before the infection. The susceptible
      # individual transits to the V state and we resolve the contact by the
      # Vcontact() function, which then might lead back to the Scontact().
      # Eventually, the infection either occurs or is prevented and the
      # recursion ends. An infinite (or very long) loop should not occur in a 
      # meaningful setting. Otherwise, the epidemic is in some way degenerate 
      # and the input parameters should be reconsidered (e.g. the cycle of 
      # vaccination and loosing the acquired immunity is unproportionally faster
      # than the "susceptibility -> infection -> recovery -> susceptibility" 
      # cycle).
      
      # Updates the actual time and state of the individual.
      vac.time[[new.infected]] <<- c(vac.time.cand, vac.time[[new.infected]])
      act.time <<- vac.time.cand
      state[new.infected] <<- "V"
      Vcontact()
    }
  }
  
  Vcontact <- function () {
    # The contact was made between an infectious and a vaccinated individual.
    # The vaccinated individual might have already lost the immunity. In the 
    # following, we generate the candidate renewed susceptibility time and 
    # compare it with the possible infection time.
    
    # The candidate second dose time is the last time entering the first dose
    # state plus the exponentially distributed time with parameter 'phi1'.
    suscep.time.cand <- vac.time[[new.infected]][1] + rexp(1, phi)
    
    if (contact.times[[active.min]][1] < suscep.time.cand) {
      # The infection occurred earlier than the loss of immunity, then the
      # individual was protected by the vaccine and nothing happens.
      
      # Updates the actual time.
      act.time <<- contact.times[[active.min]][1]
    } else {
      # The immunity was lost before the infection. The vaccinated  individual 
      # transits to the S state and we resolve the contact by the Scontact() 
      # function, which then might lead back to the Vcontact(). Eventually, the 
      # infection either occurs or is prevented and the recursion ends.
      # An infinite (or very long) loop should not occur in a meaningful 
      # setting. Otherwise, the epidemic is in some way degenerate and the 
      # input parameters should be reconsidered (e.g. the cycle of vaccination
      # and loosing the acquired immunity is unproportionally faster than the 
      # "susceptibility -> infection -> recovery -> susceptibility" cycle).
      
      # Updates the actual time and state of the individual.
      suscep.time[[new.infected]] <<- c(suscep.time.cand, 
                                        suscep.time[[new.infected]])
      act.time <<- suscep.time.cand
      state[new.infected] <<- "S"
      Scontact()
    }      
  }
  
  Icontact <- function () {
    # The contact was made between two infectious individuals. The contacted one
    # might have already recovered/died. In the  following, we check the state.
    
    if (removal.time[new.infected] < contact.times[[active.min]][1]) {
      # The individual has already recovered/died.
      
      if (is.na(death.time[new.infected])) {
        # The individual has recovered. The contact is resolved by the 
        # Rcontact() function.
        act.time <- recov.time[[new.infected]][1]
        state[new.infected] <<- "R"
        Rcontact()
      } else {
        # The individual has died. Nothing happens
        state[new.infected] <<- "D"
        act.time <<- contact.times[[active.min]][1]
      }
    } else {
      # The individual is still infectious, therefore nothing happens
      act.time <<- contact.times[[active.min]][1]
    }
  }
  
  Rcontact <- function () {
    # The contact was made between an infectious and a recovered individual. The
    # gained immunity might have already worn out and the recovered individual 
    # might be susceptible again. In the following, we generate the candidate
    # renewed susceptibility time and compare it with the possible infection 
    # time.
    
    # The candidate renewed susceptibility time is the last time entering the 
    # Recovered state plus the exponentially distributed time with parameter 
    # 'omega'.
    
    if (is.na(next.trans.time[new.infected])) {
      next.trans.time[new.infected] <<- 
        recov.time[[new.infected]][1] + rexp(1, omega)
    }
    
    
    if (contact.times[[active.min]][1] < next.trans.time[new.infected]) {
      # The infection occurred earlier than the loss of immunity. Therefore 
      # nothing happens
      
      # Updates the actual time.
      act.time <<- contact.times[[active.min]][1]
      
    } else {
      # The individual became susceptible again before the infection.
      
      # Updates the actual time and state of the individual.
      suscep.time[[new.infected]] <<- c(next.trans.time[new.infected], 
                                        suscep.time[[new.infected]])
      act.time <<- next.trans.time[new.infected]
      next.trans.time[new.infected] <<- NA
      state[new.infected] <<- "S"
      Scontact()
    }
  }
  
  EDcontact <- function () {
    # The contact was made between an infectious and an exposed individual, or 
    # between an infectious and a dead individual. Nobody can be infected again,
    # therefore nothing happens.
    
    act.time <<- contact.times[[active.min]][1]
  }
  
  # Allocation of variables ====================================================
  n.pop <- s0 + i0 + v0 + r0
  
  # Keeps track of the current state of each individual. The possible states
  # are "S", "E", "I" , "V", "R", "D".
  state <- c(rep("I", i0), rep("S", s0), rep("V", v0), rep("R", r0))
  
  # Stores the exponential time of the R -> S transition.
  next.trans.time <- rep(NA, n.pop)
  
  # Allocates the times of the possible renewed susceptibility.
  suscep.time <- vector(mode = "list", length = n.pop)
  
  # Allocates the infection and exposition times, which are lists due to the
  # possibility of a repeated infection. At the beginning, all elements of the
  # lists are empty (NULL) except for the initially infectious individuals
  # whose times of exposition and infectiousness are set as 0.
  infect.time <- suscep.time
  infect.time[1:i0] <- as.list(rep(0, i0))
  exposed.time <- infect.time
  
  # Allocates the times of vaccination and possible recovery times. All are 
  # list, since the events can occur more than once. We do not allocate a new 
  # list, rather we reuse already existing empty list 'suscep.time'. If 'v0' 
  # or 'r0' are non-zero, we have to fill the corresponding positions by zero.
  vac.time <- suscep.time
  if (v0 > 0) {
    vac.time[s0 + i0 + 1:v0] <- as.list(rep(0, v0))
  }
  recov.time <- suscep.time
  if (r0 > 0) {
    recov.time[1:r0 + n.pop - r0] <- as.list(rep(0, r0))
  }
  
  # Allocates the death times, which is not a list but only a simple vector,
  # because unlike the other states, death can happen only once. The 
  # 'removal.time' vector stores for each individual the time of completing the
  # infectious period, i.e. either the death or recovery times. This variable is
  # not returned by the function, it only helps us in the process of deciding,
  # which individual is still infectious and updating the active cases.
  death.time <- next.trans.time
  removal.time <- death.time
  
  # Decides the outcome of the initial infectious individuals after the
  # completion of their infectious period. After it, they can become either
  # recovered or die.
  infect.per.toR <- rexp(i0, mu)
  infect.per.toD <- rexp(i0, kappa)
  death.first <- infect.per.toR < infect.per.toD
  
  death.time[1:i0][death.first] <- infect.per.toD[death.first]
  recov.time[1:i0][!death.first] <- as.list(infect.per.toR[!death.first])
  
  infect.per <- ifelse(death.first, infect.per.toD, infect.per.toR)
  removal.time[1:i0] <- infect.per
  removal.time[i0 + s0 + v0 + 1:r0] <- 0
  
  # Generates the Poisson process points for the initial infectious individuals.
  contact.times <- lapply(infect.per, poiss.proc.points, lambda = lambda)
  
  # May delete before submitting
  rm(infect.per.toD, infect.per.toR, death.first)
  
  # Allocates the list of the Poisson process points for the exposed 
  # individuals. These points are stored separately from 'jump times' to avoid
  # searching through them. 
  contact.times.latent <- list()
  
  # Starts measuring the time of events and sets the counter of infections.
  act.time <- 0
  infections <- 0
  
  # Keeps track of the latent cases. At the beginning, there are no latent
  # cases. Here active cases are those, who are infectious and can potentially 
  # infect someone else, i.e. their element in the list 'contact.times' is 
  # non-empty.
  
  still.active <- unlist(lapply(contact.times, length)) != 0
  contact.times <- contact.times[still.active]
  latent.cases <- c()
  
  # Keeps track, whether there are any active/latent cases.
  any.active <- length(unlist(contact.times)) != 0
  any.latent <- FALSE
  
  # The main cycle of the algorithm ============================================
  
  # The epidemic continues if there are at least one active or latent case,
  # the maximum number of infections does not exceed the 'max.infections'
  # parameter and the duration of the epidemic does not exceed the 'max.duration'
  # parameter.
  while ((any.active || any.latent) != 0 && infections <= max.infections
         && act.time < max.duration) {
    
    # Finds out, which event occurred first - a new infection, or a transition
    # from the Exposed state to the Infectious state. The conditions depend on
    # whether there are any latent and/or active cases.
    if (any.active && any.latent) {
      
      # Finds the earliest time of potential infection.
      active.min <- which.min(lapply(contact.times, FUN = choose.first))
      
      # Finds the earliest time of possible transition from the Exposed, to 
      # the Infectious compartment.
      latent.min <- which.min(lapply(infect.time[latent.cases], choose.first))
      if (contact.times[[active.min]][1] < infect.time[[latent.cases[latent.min]]][1]) {
        infection <- TRUE
      } else {
        infection <- FALSE
      }
      
    } else{
      if (any.active) {
        infection <- TRUE
        # Finds the earliest time of potential infection.
        active.min <- which.min(lapply(contact.times, FUN = choose.first))
      } else {
        infection <- FALSE
        # Finds the earliest time of possible transition from the Exposed, to 
        # the Infectious compartment
        latent.min <- which.min(lapply(infect.time[latent.cases], choose.first))
      }
    }
    
    if (infection) {
      # Infection occurred first.
      new.infected <- sample(1:n.pop, size = 1)
      
      # The success of the infection depends on the state, in which the
      # potentially newly infected individual is.
      
      switch(state[new.infected], "S" = Scontact(), "V" = Vcontact(),
             "E" = EDcontact(), "R" = Rcontact(), "I" = Icontact(),
             "D" = EDcontact())
      
      # The 'used' point of the Poisson process is removed. If it was the last
      # point (contact), the individual is removed from the 'contact.times'
      # list.
      contact.times[[active.min]] <- contact.times[[active.min]][-1]
      if (length(contact.times[[active.min]]) == 0) {
        contact.times <- contact.times[-active.min]
      }
      
    } else {
      # Transition from the Exposed state to the Infectious state occurred
      # first.
      
      state[latent.cases[latent.min]] <- "I"
      act.time <- infect.time[[latent.cases[latent.min]]][1]
    }
    
    # Updates the latent cases, which are those, who are exposed and have not
    # become infectious yet. The Poisson process points of the newly infectious
    # individuals are added to the 'contact.times' list.
    still.latent <- 
      unlist(lapply(infect.time[latent.cases], choose.first)) > act.time
    
    contact.times <- c(contact.times, contact.times.latent[!still.latent])
    contact.times.latent <- contact.times.latent[still.latent]
    
    state[latent.cases[!still.latent]] <-  "I"
    latent.cases <- latent.cases[still.latent]
    
    # Checks if there are any latent and/or active cases. If both FALSE, the
    # cycle ends.
    any.active <- length(contact.times) != 0
    any.latent <- length(latent.cases) != 0
  }
  
  # Evaluation of possible transfers other than 'S -> E' and 'E -> I', which 
  # occurred before the end of the epidemic ====================================
  
  # The evaluation is done only if the user is interested in compartments 'S'
  # and 'V' too.
  epi.end <- min(act.time, max.duration)
  rm(act.time)
  
  if (all.comps) {
    
    # Finds the particles, which transferred from 'I' to 'R' before the end.
    in.state <- which(state == "I")
    
    moves <- 
      removal.time[in.state] < epi.end & !is.na(death.time[in.state])
    
    # Evaluates possible 'R -> S' transitions. For each particle in state 'R'
    # we generate the candidate 'R -> S' transition times. If they are earlier
    # than the epidemic end, we append them to the 'suscep.time' list.
    in.state <- c(in.state[moves], which(state == "R"))
    if (length(in.state) != 0) {
      suscep.time.cand <- removal.time[in.state] + rexp(length(in.state), 
                                                        omega)
      moves <- suscep.time.cand < epi.end
      if (any(moves)) {
        # For 'moves' all equal FALSE, we obtain error message.
        suscep.time[in.state[moves]] <- mapply(
          c, as.list(suscep.time.cand[moves]),
          suscep.time[in.state[moves]],
          SIMPLIFY = FALSE)
      }
    }
    
    # Evaluates possible 'S -> V' transitions. For each particle in state 'S'
    # we generate the candidate 'S -> V' transition times. If they are earlier
    # than the epidemic end, we append them to the 'vac.time' list.
    in.state <- c(in.state[moves], which(state == "S"))
    len <- length(in.state)
    # For vectors of length 0 we obtain errors.
    if (len != 0) {
      # Some vaccination times might be non-existent, so we perceive them as 0.
      vac.time.cand <- 
        unlist(
          lapply(
            suscep.time[in.state], 
            function (x) {
              ifelse(length(x) == 0, 0, x[1])
            })) +
        rexp(len, psi)
      moves <- vac.time.cand < epi.end
    } else {
      moves <- rep(FALSE, len)
    }
    
    if (any(moves)) {
      # For 'moves' all equal FALSE, we obtain error message.
      vac.time[in.state[moves]] <- mapply(
        c, as.list(vac.time.cand[moves]),
        vac.time[in.state[moves]],
        SIMPLIFY = FALSE
      )
    }
    
    # Evaluates possible 'V -> S' transitions. For each particle in state 'V'
    # we generate the candidate 'V -> S' transition times. If they are earlier
    # than the epidemic end, we append them to the 'suscep.time' list.
    in.state <- c(which(state == "V"), in.state[moves])
    # For vectors of length 0 we obtain errors.
    len <- length(in.state)
    if (len != 0) {
      # Some times of renewed susceptibility might be non-existent, so we perceive
      # them as 0.
      suscep.time.cand <- 
        unlist(
          lapply(
            vac.time[in.state], 
            function (x) {
              ifelse(length(x) == 0, 0, x[1])
            })) + rexp(len, phi)
      moves <- suscep.time.cand < epi.end
    } else {
      moves <- rep(FALSE, len)
    }
    
    if (any(moves)) {
      # For 'moves' all equal FALSE, we obtain error message.
      suscep.time[in.state[moves]] <- mapply(
        c, as.list(suscep.time.cand[moves]),
        suscep.time[in.state[moves]],
        SIMPLIFY = FALSE
      )
    }
    
    # All particles that moved from 'V' to 'S' might further move in the cycle
    # 'S -> V -> S'. The cycle ends, when no particle can move anymore.
    while (any(moves)) {
      
      # The 'S -> V' direction
      in.state <- in.state[moves]
      len <- length(in.state)
      if (len != 0) {
        # Some vaccination times might be still non-existent, so we perceive 
        # them as 0.
        vac.time.cand <-
          unlist(
            lapply(
              suscep.time[in.state],
              function (x) {
                ifelse(length(x) == 0, 0, x[1])
              })) +
          rexp(len, psi)
        moves <- vac.time.cand < epi.end
      } else {
        moves <- rep(FALSE, len)
      }
      
      if (any(moves)) {
        # For 'moves' all equal FALSE, we obtain error message.
        vac.time[in.state[moves]] <- mapply(
          c, as.list(vac.time.cand[moves]),
          vac.time[in.state[moves]],
          SIMPLIFY = FALSE
        )
      }
      
      # The 'V -> S' direction
      in.state <- in.state[moves]
      len <- length(in.state)
      if (len != 0) {
        # Some times of renewed susceptibility might be non-existent, so we 
        # perceive them as 0.
        suscep.time.cand <-
          unlist(
            lapply(
              vac.time[in.state],
              function (x) {
                ifelse(length(x) == 0, 0, x[1])
              })) + rexp(len, phi)
        moves <- suscep.time.cand < epi.end
      } else {
        moves <- rep(FALSE, len)
      }
      
      if (any(moves)) {
        # For 'moves' all equal FALSE, we obtain error message.
        suscep.time[in.state[moves]] <- mapply(
          c, as.list(suscep.time.cand[moves]),
          suscep.time[in.state[moves]],
          SIMPLIFY = FALSE)
      }
    }
    rm(in.state, moves)
  }
  
  # Removal of events after the end of the epidemic ============================
  
  # There can be at most one individual exposed after the epidemic's end.
  # If this is the case, the label of such individual is still stored in the 
  # 'new.infected' variable. But in case of zero infections, we must be aware
  # that the variable 'new.infected' does not exist.
  if (infections > 0 && !is.null(exposed.time[[new.infected]][1]) && 
      exposed.time[[new.infected]][1] > max.duration) {
    exposed.time[[new.infected]] <- exposed.time[[new.infected]][-1]
    # This individual is also the only one, which can have two infection times 
    # after the end of the epidemic. So we remove one infection time. The 
    # potential second one is removed in the next code chunk. The same applies 
    # for the recovery times.
    infect.time[[new.infected]] <- infect.time[[new.infected]][-1]
    if (is.na(death.time[new.infected])) {
      recov.time[[new.infected]] <- recov.time[[new.infected]][-1]
    }
    
    # There can be more than two vaccination/susceptibility times for this
    # individual, so we remove them all here.
    if (all.comps) {
      suscep.time[[new.infected]] <-
        suscep.time[[new.infected]][suscep.time[[new.infected]] < max.duration]
      vac.time[[new.infected]] <-
        vac.time[[new.infected]][vac.time[[new.infected]] < max.duration] 
    }
  }
  
  # Removes the infection times after the epidemic end.
  after <- which(unlist(lapply(infect.time, function (x) {
    ifelse(length(x) == 0, 0, x[1])
  })) > epi.end)
  infect.time[after] <- lapply(infect.time[after], function (x) {x[-1]})
  
  # Removes the recovery and death times after the epidemic end.
  after <- which(removal.time > epi.end)
  recovered <- is.na(death.time[after])
  recov.time[after[recovered]] <- 
    lapply(recov.time[after[recovered]], function (x) {x[-1]})
  death.time[after[!recovered]] <- NA
  
  # Returning values ===========================================================
  if (all.comps) {
    return(list(S = suscep.time, V = vac.time, E = exposed.time,
                I = infect.time, R = recov.time, D = death.time, state = state,
                duration = epi.end, stopped = infections > max.infections))
  } else {
    return(list(E = exposed.time, I = infect.time, R = recov.time, 
                D = death.time, state = state, duration = epi.end, 
                stopped = infections > max.infections))
  }
}