seir_ode <- function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  inc<-Y[5]
  N <- sum(Y[1:4])
  
  beta<-par[1]
  sigma<-par[2]
  gamma<-par[3]
  
  dYdt<-vector(length=4)
  dYdt[1]= -beta*I*S/N 
  dYdt[2]= beta*I*S/N - sigma*E
  dYdt[3]= sigma*E - gamma*I
  dYdt[4]= gamma*I
  dYdt[5] = beta*I*S/N
  
  return(list(dYdt))
}

simulate_seir_process <- function(n_indivs, pars, times){
  # Set parameter values
  R0 <- pars["R0"]
  sigma<-pars["sigma"];
  gamma<-pars["gamma"];
  beta <- R0*gamma
  N <- n_indivs
  I0 <- pars["I0"]
  
  init <- c(N-I0,0,I0,0,0)
  t<-times
  par<-c(beta,sigma,gamma)
  # Solve system using lsoda
  sol<-lsoda(init,t,seir_ode,par)
  # Plot solution
  sol <- as.data.frame(sol)
  colnames(sol) <- c("time","S","E","I","R","cumulative_incidence")
  incidence <- diff(c(0,sol$cumulative_incidence/N))
  prevalence <- (sol$E + sol$I)/N
  sol_melted <- reshape2::melt(sol, id.vars="time")
  sol_melted$value <- sol_melted$value/N
  
  p <- ggplot(sol_melted) + 
    geom_line(aes(x=time,y=value,col=variable)) + 
    ylab("Per capita") + 
    xlab("Date") +
    theme_bw()
  p_inc <- ggplot(data.frame(x=t,y=incidence,y1=prevalence)) + geom_line(aes(x=x,y=y),col="red") +
    geom_line(aes(x=x,y=y1),col="blue") +
    ylab("Per capita incidence (red) and prevalence (blue)") + 
    xlab("Date") +
    theme_bw()
  
  return(list(plot=p, incidence_plot=p_inc, incidence=incidence, 
              seir_outputs=sol,prevalence=prevalence,seir_outputs_melted=sol_melted,
              overall_prob_infection=sum(incidence)))  
}

simulate_infection_times <- function(n, p_infected, incidence){
  scaled_incidence <- incidence/sum(incidence)
  are_infected <- numeric(n)
  infection_times <- numeric(n)
  for(i in 1:n){
    infection <- rbinom(1,1, p_infected)
    are_infected[i] <- infection
    if(infection == 1){
      t_inf <- sample(1:length(incidence), 1, prob=scaled_incidence)
      infection_times[i] <- t_inf
    } else {
      infection_times[i] <- -1
    }
  }
  return(infection_times)
}

## Viral kinetics function
model_func_tinf <- function(ts, tinf, tp, viral_peak, wane){
  y <- numeric(length(ts))
  growth <- viral_peak/tp
  y[ts < tinf] <- 0
  y[ts <= tp + tinf & ts > tinf] <- growth*(ts[ts <= tp + tinf & ts > tinf]-tinf)
  y[ts > tp + tinf] <- viral_peak + wane*tp - wane*(ts[ts > tp + tinf] - tinf)
  y
}

model_trajectory <- function (pars, times, logSigma = FALSE) 
{
  mu <- pars["mu"]
  dp <- pars["dp"]
  tp <- pars["tp"]
  ts <- pars["ts"]
  m <- pars["m"]
  y0 <- pars["y0"]
  eff_y0 <- pars["eff_y0"]
  t_i <- pars["t_i"]
  x <- pars["x"]
  c <- pars["c"]
  tau <- pars["tau"]
  order <- pars["order"]
  primed <- pars["primed"]
  mod <- pars["mod"]
  lower_bound <- pars["lower_bound"]
  if (logSigma) {
    sigma <- exp(pars["sigma"])
    beta <- exp(pars["beta"])
    y0_mod <- exp(pars["y0_mod"])
  }
  else {
    sigma <- pars["sigma"]
    beta <- pars["beta"]
    y0_mod <- pars["y0_mod"]
  }
  boost_limit <- pars["boost_limit"]
  cr <- mu - sigma * x
  prime_cr <- c - beta * x
  if (cr < 0) 
    cr <- 0
  if (prime_cr < 0) 
    prime_cr <- 0
  mu <- cr + prime_cr * primed
  titre_dependent <- 1
  if (y0_mod >= -999) {
    if (y0 >= boost_limit) {
      titre_dependent <- (1 - y0_mod * boost_limit)
    }
    else {
      titre_dependent <- (1 - y0_mod * y0)
    }
  }
  if (titre_dependent < 0) 
    titre_dependent <- 0
  seniority <- 1 - tau * (order - 1)
  if (seniority < 0) 
    seniority <- 0
  mu <- seniority * titre_dependent * mod * mu
  if (mu < 0) 
    mu = 0
  y <- numeric(length(times))
  i <- 1
  for (t in times) {
    tmp <- 0
    if (t <= t_i) 
      tmp = 0
    else if (t > t_i & t <= (t_i + tp)) 
      tmp = (mu/tp) * (t - t_i)
    else if (t > (tp + t_i) & t <= (ts + t_i + tp)) 
      tmp = ((-(dp * mu)/ts) * (t) + ((mu * dp)/ts) * (t_i + 
                                                         tp) + mu)
    else tmp = (-m * (t) + m * (t_i + tp + ts) + (1 - dp) * 
                  mu)
    y[i] <- tmp + eff_y0
    if (y[i] < lower_bound) 
      y[i] = lower_bound
    i <- i + 1
  }
  return(y)
}
