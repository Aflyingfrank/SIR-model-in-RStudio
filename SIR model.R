library(deSolve)
library(ggplot2)
library(dplyr)

#initial values for SIR model
t <- 300       #time for model (unit in days)
b <- 0.184   #value for beta (contact rate)
g <- 0.1   #value for gamma (recovery rate)
P <- 100000   #value for assumed population

#input initial values into the model
SIR.model <- function(t, b, g)
  require(deSolve)
  initial_values <- c(S=1-1e-6, I=1e-6, R=0) #initial population for three classes
  parameters <- c(beta=b, gamma=g)
  time <- seq(0, t, by=t/(5*length(1:t)))

#define the differential equations
equation <- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
      dS <- -beta*S*I
      dI <- beta*S*I-gamma*I
      dR <- gamma*I
      return(list(c(dS, dI, dR)))
      })
}

#calculate output for proportion of S,I,R
output <- ode(y=initial_values, times=time, func=equation, parms=parameters)
output.df <- as.data.frame(output)

#calculate output for the I population
output.I <- output.df[, c("time", "I")]
output.I <- output.I %>% mutate(I.P=I*P)

#set theme for the graph
require(ggplot2)
mytheme <- theme_bw() +
  theme(text=element_text(colour="black")) +
  theme(panel.grid = element_line(colour = "white")) +
  theme(panel.background = element_rect(fill = "lightgrey"))
theme_set(mytheme)

#set title and subtitle for the graph
title <- bquote("SIR Model for Influenza")
subtit <- bquote(list(beta==.(parameters[1]), ~gamma==.(parameters[2])))

#plot the graph for proportion of S,I,R
graph_SIR <- ggplot(output.df, aes(x=time))+
  ggtitle("")+
  geom_line(aes(y=S, colour="Susceptible"))+
  geom_line(aes(y=I, colour="Infected"))+
  geom_line(aes(y=R, colour="Recovered"))+
  ylab(label="Proportion (%)")+
  xlab(label="Time (days)")+
  theme(plot.title = element_text(size=16, hjust = 0.5))+
  theme(plot.subtitle = element_text(size=12, hjust = 0.5))+
  theme(legend.justification=c(4.5,0), legend.position=c(0.95,0.4))+
  theme(legend.title=element_text(size=12, face="bold", hjust=0.5),
        legend.background = element_rect(fill='#FFFFFF', size=0.5, linetype="solid"),
        legend.text=element_text(size=10),
        legend.key=element_rect(colour="#FFFFFF", fill='lightgrey', size=0.25, linetype="solid"))+
  scale_colour_manual("Classes", breaks=c("Susceptible", "Infected", "Recovered"), values=c("blue", "red", "darkgreen"))
print(graph_SIR)

#plot the graph for the I population
graph_infection <- ggplot(output.I, aes(x=time, y=I.P))+
  ggtitle("")+
  geom_line(color="red")+
  ylab(label="Population")+
  xlab(label="Time (days)")+
  theme(plot.title = element_text(size=16, hjust = 0.5))
print(graph_infection)
