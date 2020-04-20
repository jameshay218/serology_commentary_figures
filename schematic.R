library(tidyverse)
library(ggplot2)
library(deSolve)
library(patchwork)
library(ggpubr)
source("schematic_functions.R")

set.seed(2)

## Entire popuation size (this is MA in 2018 from US Census Bureau)
population_n <- 6900000

## Sample size
n <- 1000

## Epidemic growth rate
growth_rate <- 0.1

## Duration of epidemic so far in days
times <- seq(0,365,by=1)

## Simulate the epidemic process
#epidemic_process <- simulate_epidemic_process(population_n,0.1,times)
#epidemic_process$plot
## R0 of 2.5, infectious period of 5 days, incubation period of 5.1 days, I0 of 10 days
seir_pars <- c("R0"=2.4,"gamma"=1/5,"sigma"=1/5.1,"I0"=10)
epidemic_process <- simulate_seir_process(population_n,seir_pars,times)
epidemic_process$plot

## Simulate infection times
infection_times <- simulate_infection_times(n, epidemic_process$overall_prob_infection, 
                                            epidemic_process$incidence)

seir_dat <- epidemic_process$seir_outputs
beta <- seir_pars["R0"]*seir_pars["gamma"]
seir_dat$N <- rowSums(seir_dat[,c("S","E","I","R")])
seir_dat$FOI <- beta*seir_dat$I/seir_dat$N
seir_dat[,c("S","E","I","R")] <- seir_dat[,c("S","E","I","R")]/seir_dat$N
seir_dat <- seir_dat[,c("time","S","E","I","R","FOI")]

ab_pars <- c("lower_bound"=0,"S"=1,"EA"=0,"MAX_TITRE"=13,
          "mu"=8,"tp"=14,"dp"=0.5,"ts"=20,"m"=0.01,"beta"=0.6, "c"=4,
          "sigma"=1,"y0_mod"=-10000,"boost_limit"=0,"tau"=0.05,
          "order"=1,"primed"=0,"mod"=1,
          "x"=0,"t_i"=10,"y0"=0,"eff_y0"=0)

## Vector of times to solve the model over
times <- seq(0,365,by=0.1)

all_trajectories <- matrix(0, nrow=length(infection_times),ncol=length(times))

for(i in seq_along(infection_times)){
  if(infection_times[i] > 0){
    ab_pars["t_i"] <- infection_times[i]
    all_trajectories[i,] <- model_trajectory(ab_pars,times)
    all_trajectories[i,(ab_pars["t_i"]*10-5):(ab_pars["t_i"]*10)] <- NA
  }
}

traj_melted <- reshape2::melt(all_trajectories[1:100,])
colnames(traj_melted) <- c("individual","time","Ab titer")

seir_dat_melt <- reshape2::melt(seir_dat, id.vars="time")
colnames(seir_dat_melt)[2] <- "Compartment"
seir_key <- c("S"="Susceptible","E"="Exposed","I"="Infected","R"="Removed","FOI"="Force of infection")
seir_dat_melt$Compartment <- seir_key[seir_dat_melt$Compartment]
seir_dat_melt$Compartment <- factor(seir_dat_melt$Compartment,levels=seir_key)

pA <- ggplot(seir_dat_melt) + 
  geom_line(aes(x=time,y=value,col=Compartment)) +
  scale_color_manual(values=c("black","#403891ff","6b4596ff","#f9a242ff","#5DC863FF")) +
  #scale_color_viridis_d(option="B") +
  guides(color=guide_legend(keyheight=0.2,default.unit="inch",ncol=2)) +
  theme_pubr() +
  theme(legend.position=c(0.7,0.45),
        panel.grid=element_blank(),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.box.background = element_blank(),
        legend.background =  element_rect(color=NA,fill=NA)) +
  xlab("") +
  ylab("Per capita") +
  scale_x_continuous(expand=c(0,0),labels=seq(0,365,by=50),breaks=seq(0,365,by=50)) +
  labs(tag="A")
pB <- ggplot(traj_melted) + 
  geom_tile(aes(x=time, y=individual, fill=`Ab titer`),col=NA) + 
  scale_x_continuous(expand=c(0,0),breaks=seq(0,3650,by=500),labels=seq(0,365,by=50)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_viridis_c(option="B",na.value="white") + 
  theme_minimal() +
  xlab("Time (day)") +
  ylab("Individual") +
  theme(legend.position=c(0,1),
        legend.justification = c(0,1),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8),
        legend.box.background = element_rect(color="white",fill="white"))+
  labs(tag="B")
left_plot <- pA / pB + plot_layout(heights=c(1,2))

## Pretend cross sectional studies etc
## pD
sample_1 <- all_trajectories[1:200,900]
sample_1[sample_1 < 2] <- 0
sample_1[sample_1 >= 2] <- 1
sample_1[is.na(sample_1)] <- 0
sample_2 <- all_trajectories[1:400,1100]
sample_2[sample_2 < 2] <- 0
sample_2[sample_2 >= 2] <- 1
sample_2[is.na(sample_2)] <- 0

confint_1 <- prop.test(sum(sample_1),length(sample_1), conf.level=0.95,correct=FALSE)$conf.int[1:2]
confint_2 <- prop.test(sum(sample_2),length(sample_2), conf.level=0.95,correct=FALSE)$conf.int[1:2]

dat1 <- data.frame("Proportion"=sum(sample_1/length(sample_1)), "Lower"=confint_1[1], "Upper"=confint_1[2],"Study"="Time 1")
dat2 <- data.frame("Proportion"=sum(sample_2/length(sample_2)), "Lower"=confint_2[1], "Upper"=confint_2[2],"Study"="Time 2")
dat <- rbind(dat1, dat2)
dat$Study <- as.factor(dat$Study)
dat_join <- data.frame(y=c(sum(sample_1/length(sample_1)),sum(sample_2/length(sample_2))),Study=c("Time 1","Time 2"))
dat_join$Study <- as.factor(dat$Study)
pDi <- ggplot(dat) + 
  geom_point(aes(x=Study,y=Proportion),size=1) + 
  geom_errorbar(aes(x=Study,ymin=Lower,ymax=Upper),size=0.5,width=0.1) + 
  geom_line(data=dat_join, aes(x=Study,y=y,group=1),linetype="dashed") +
  ylab("Prevalence of seropositivity") +
  theme_pubr() +
  ggtitle("Cross-sectional studies") +
  theme(
    plot.title=element_text(size=8,face="bold",hjust=0.5),
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8)) +
  scale_x_discrete(expand=c(0.2,0.2))+
  xlab("") +
  scale_y_continuous(limits=c(0,0.5))+ 
  labs(tag="C")
pDii <- ggplot(dat) + 
  geom_bar(aes(x=Study,y=Proportion),stat="identity",width=0.5,fill="black") + 
  scale_y_continuous(limits=c(0,0.5))  +
  theme_pubr() +
  theme(
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8)) +
  scale_x_discrete(expand=c(0.2,0.2))+
  ylab("Proportion seropositive") +
  xlab("")
pDii
pD <- pDi/pDii 


## Pretend cross sectional studies etc
## pD
titres_1 <- all_trajectories[1:50,900]
titres_1[is.na(titres_1)] <- 0
titres_1 <- titres_1 + rnorm(length(titres_1), 0, 0.25)

titres_2 <- all_trajectories[1:50,1100]
titres_2[is.na(titres_2)] <- 0
titres_2 <- titres_2 + rnorm(length(titres_2), 0, 0.25)

dat_titre <- data.frame(t1=titres_1, t2=titres_2,individual=1:length(titres_1))
dat_titre <- reshape2::melt(dat_titre,id.vars="individual")
key <- c("t1"="Visit 1", "t2" = "Visit 2")
dat_titre$variable <- key[dat_titre$variable]


dat1 <- data.frame("Proportion"=sum(sample_1/length(sample_1)), "Lower"=confint_1[1], "Upper"=confint_1[2],"Study"="Visit 1")
dat2 <- data.frame("Proportion"=sum(sample_2/length(sample_2)), "Lower"=confint_2[1], "Upper"=confint_2[2],"Study"="Visit 2")
dat <- rbind(dat1, dat2)
dat$Study <- as.factor(dat$Study)
dat_join <- data.frame(y=c(sum(sample_1/length(sample_1)),sum(sample_2/length(sample_2))),Study=c("Visit 1","Visit 2"))
dat_join$Study <- as.factor(dat$Study)

pEi <- ggplot(dat) + 
  geom_point(aes(x=Study,y=Proportion),size=1) + 
  geom_errorbar(aes(x=Study,ymin=Lower,ymax=Upper),size=0.5,width=0.1) + 
  geom_line(data=dat_join, aes(x=Study,y=y,group=1),linetype="dashed") +
  ylab("Prevalence of seropositivity") +
  ggtitle("Follow-up individuals") +
  theme_pubr() +
  scale_x_discrete(expand=c(0.2,0.2))+
  xlab("") +
  theme(
    plot.title=element_text(size=8,hjust=0.5,face="bold"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8)) +
  scale_y_continuous(limits=c(0,0.5))+ 
  labs(tag="D")
pEii <- ggplot(dat_titre) + 
  geom_point(aes(x=variable,y=value,group=individual),shape=4,size=1)+
  geom_line(aes(x=variable,y=value,group=individual),size=0.15) +
  geom_hline(yintercept=2,linetype="dotted") +
  scale_x_discrete(expand=c(0.2,0.2))+
  xlab("") +
  theme_pubr() +
  theme(
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8)) +
  ylab("Antibody titer") 
pE <- pEi/pEii

top_right <- pD | pE

(left_plot | top_right) + plot_layout(widths=c(1,1.5))

## Solve the model using both the R and Cpp implementations 
ab_pars["t_i"] <- 14
ab_pars["m"] <- 0.01
ab_pars["ts"] <- 14
times <- seq(0,75,by=0.1)
titre_trajectory_cpp <- model_trajectory_cpp(ab_pars,times)
y <- smooth.spline(titre_trajectory_cpp,spar=0.4)$y
y[y < 0] <- 0
trajectory_dat <- data.frame(y=y,time=times,Quantity="Antibody titer")

viral_load <- model_func_tinf(times, 5, 5, 7, 0.5)
y1 <- smooth.spline(viral_load,spar=0.4)$y
y1[y1 < 0] <- 0
viral_dat <- data.frame(y=y1, time=times, Quantity="LOG10 viral load")
all_dats <- rbind(trajectory_dat, viral_dat)

pFi <- ggplot(all_dats) + geom_line(aes(x=time, y=y,col=Quantity)) + 
  geom_hline(yintercept = 2,linetype="dotted") + 
  ylab("Antibody titer (green)") +
  xlab("Time since infection") +
  scale_y_continuous(sec.axis=sec_axis(~.,name="LOG10 viral load/ml (blue)")) +
  scale_color_manual(values=c("#5DC863FF", "#3B528BFF"))+
  ggtitle("Follow-up sampling of infected individuals") +
  theme_pubr() +
  theme(legend.position=c(0.8,0.8),
        plot.title=element_text(hjust=0.5,size=8,face="bold"),
    legend.text=element_text(size=7),
    legend.title=element_text(size=7),
    axis.text=element_text(size=8),
    axis.title=element_text(size=8)) +
  labs(tag="E")
pFi

right <- (top_right / pFi) + plot_layout(heights=c(2,1))

pdf("tmp2.pdf",height=7,width=8)
(left_plot | right) + plot_layout(widths=c(1,1))
dev.off()

png("tmp2.png",height=7,width=8,res=600,units="in")
(left_plot | right) + plot_layout(widths=c(1,1))
dev.off()

