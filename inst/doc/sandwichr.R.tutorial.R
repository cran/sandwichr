## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, message=FALSE, fig.cap = " ", fig.path='figs/')

## -----------------------------------------------------------------------------
# Install the sandwichr package
# install.packages("sandwichr")

## -----------------------------------------------------------------------------
# Import the sandwichr package and other packages
library("sandwichr")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ape)

## -----------------------------------------------------------------------------
# Input data from shapefiles
sim.sampling.name <- system.file("extdata", "sim.sampling.shp", 
                                package="sandwichr")
sim.ssh.name <- system.file("extdata", "sim.ssh.shp", 
                           package="sandwichr")
sim.reporting.name <- system.file("extdata", "sim.reporting.shp", 
                                 package="sandwichr")

sim.data <- load.data.shp(sampling.file=sim.sampling.name, 
                      ssh.file=sim.ssh.name,
                      reporting.file=sim.reporting.name)

# Sampling
head(sim.data[[1]])
class(sim.data[[1]])
attributes(sim.data[[1]])

# Stratification
head(sim.data[[2]])
class(sim.data[[2]])
attributes(sim.data[[2]])

# Reporting
head(sim.data[[3]])
class(sim.data[[3]])
attributes(sim.data[[3]])

## -----------------------------------------------------------------------------
sim.dists <- as.matrix(dist(cbind(sim.data[[1]]$Lon, sim.data[[1]]$Lat)))
sim.dists.inv <- 1/sim.dists
diag(sim.dists.inv) <- 0
Moran.I(sim.data[[1]]$Value, sim.dists.inv)

## -----------------------------------------------------------------------------
# Prepare the stratification layer(s) for evaluation
sim.join <- ssh.data.shp(object=sim.data[[1]], ssh.lyr=sim.data[[2]], ssh.id="X")
head(sim.join)

## -----------------------------------------------------------------------------
# Calculate the geographical detector q-statistic
ssh.test(object=sim.join, y="Value", x=c("X"), test="factor")

## -----------------------------------------------------------------------------
p = ggerrorplot(sim.join, x = "X", y = "Value", 
            desc_stat = "mean_sd", color = "black",
            add = "violin", add.params = list(color = "darkgray")
            )

p + theme(axis.title.x = element_blank())

sim.join %>%                                        
  group_by(X) %>%                         
  summarise_at(vars(Value),              
               list(name = mean))

sim.join %>%
  group_by(X) %>%
  summarise(mean=mean(Value), sd=sd(Value), n=n())

## -----------------------------------------------------------------------------
# Perform the SSH based spatial prediction
sim.sw <- sandwich.model(object=sim.data, sampling.attr="Value", type="shp")
head(sim.sw$object)
summary(sim.sw)

## ----fig.align="center", fig.width=8, fig.height=3----------------------------
# Plot the estimated mean values and standard errors
ggplot2::autoplot(object=sim.sw)

## -----------------------------------------------------------------------------
# Perform k-fold cross validation
set.seed(0)
sim.cv <- sandwich.cv(object=sim.data, sampling.attr="Value", k=5, type="shp", ssh.id.col="X")
sim.cv

## -----------------------------------------------------------------------------
# Input data from text files
bc.sampling_ssh.name <- system.file("extdata", "bc_sampling_ssh.csv", 
                                package="sandwichr")
bc.reporting_ssh.name <- system.file("extdata", "bc_reporting_ssh.csv", 
                                 package="sandwichr")

bc.data <- load.data.txt(sampling_ssh.file=bc.sampling_ssh.name, 
                         reporting_ssh.file=bc.reporting_ssh.name)

# Sampling-stratification
head(bc.data[[1]])    
class(bc.data[[1]])

# Reporting-stratification
head(bc.data[[2]])    
class(bc.data[[2]])

## -----------------------------------------------------------------------------
bc.dists <- as.matrix(dist(cbind(bc.data[[1]]$X, bc.data[[1]]$Y)))
bc.dists.inv <- 1/bc.dists
diag(bc.dists.inv) <- 0
Moran.I(bc.data[[1]]$Incidence, bc.dists.inv)

## -----------------------------------------------------------------------------
# Prepare the stratification layer for evaluation
bc.join <- ssh.data.txt(object=bc.data)
head(bc.join)

## -----------------------------------------------------------------------------
# Calculate the geographical detector q-statistic
ssh.test(object=bc.join, y="Incidence", x="SSHID", test="factor", type="txt")

## -----------------------------------------------------------------------------
p = ggerrorplot(bc.data[[1]], x = "SSHID", y = "Incidence", 
            desc_stat = "mean_sd", color = "black",
            add = "violin", add.params = list(color = "darkgray")
            )

p + scale_x_discrete(labels=c("1" = "Urban", "2" = "Rural")) + 
  theme(axis.title.x = element_blank()) + labs(y="Breast Cancer Incidence (Rate per 100,000)")

bc.data[[1]] %>%                                        
  group_by(SSHID) %>%                         
  summarise_at(vars(Incidence),              
               list(name = mean))

bc.data[[1]] %>%
  group_by(SSHID) %>%
  summarise(mean=mean(Incidence), sd=sd(Incidence), n=n())

## -----------------------------------------------------------------------------
# Perform the SSH based spatial prediction
bc.sw <- sandwich.model(object=bc.data, sampling.attr="Incidence", type="txt", 
                        ssh.id.col="SSHID", ssh.weights=list(c(1,2), c("W1","W2")))
head(bc.sw$object)
summary(bc.sw)

## -----------------------------------------------------------------------------
# Perform k-fold cross validation
set.seed(0)
bc.cv <- sandwich.cv(object=bc.data, sampling.attr="Incidence", k=5, type="txt", 
                     ssh.id.col="SSHID", reporting.id.col="GBCODE", 
                     ssh.weights=list(c(1,2), c("W1","W2")))
bc.cv

