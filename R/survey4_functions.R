########################################
# Survey 4 Functions
########################################

# packages
library(ggplot2)
library(foreign)
library(xtable)
library(reshape)
library(haven)
library(plyr)
library(foreign)
library(systemfit)
library(lfe)

########################################
# Data Processing Functions

# adding up columns with NAs
psum <- function(..., na.rm=T) { 
  x <- list(...)
  rowSums(matrix(unlist(x), ncol=length(x)), na.rm=T)
} 

# function to relabel
relab <- function(old.var, old.labs, new.labs){
  map <- setNames(object = old.labs, nm = new.labs)
  return(as.numeric(names(map[unlist(eval(old.var))])))
}

# function to standardize
vscale <- function(var, vig=d$V){
  var <- eval(var)
  vig <- eval(vig)
  sv <- c()
  for (i in 1:length(var)){
    v <- vig[i]
    sv[i] <- (var[i] - mean(var[vig==v], na.rm = T))/
      sd(var[vig==v], na.rm=T)
  }
  return(sv)
}

# function to organize multiple versions of the same question
multiq <- function(v1, v2, v3, v4, v5, v6){
  res <- rep(NA, length(v1))
  res[!is.na(eval(v1))] <- eval(v1)[!is.na(eval(v1))]
  res[!is.na(eval(v2))] <- eval(v2)[!is.na(eval(v2))]
  res[!is.na(eval(v3))] <- eval(v3)[!is.na(eval(v3))]
  res[!is.na(eval(v4))] <- eval(v4)[!is.na(eval(v4))]
  res[!is.na(eval(v5))] <- eval(v5)[!is.na(eval(v5))]
  res[!is.na(eval(v6))] <- eval(v6)[!is.na(eval(v6))]
  return(res)
}

# function to sort out the labels for regions
region.lab <- function(region1, region2){
  t <- ifelse(eval(region1)==1 |
                eval(region2)==1, 1, 0)
  t[is.na(t)] <- 0
  return(t)
}
#####################################
# Analysis Functions

# robust SE function
robustse <- function(model){ 
  # OLS with robust SE
  # written by Baobao Zhang
  require("sandwich")
  require("lmtest")
  model$newse<-vcovHC(model,type = "HC3")
  coef <- coeftest(model,model$newse)
  return(as.matrix(coef))
}

# function to perform OLS regression and extract coef.
myreg <- function(Y, Z=d$Z, V=d$V, row=2, print.out=TRUE){
  Y <- as.numeric(eval(Y))
  Z <- as.numeric(eval(Z))
  V <- as.numeric(eval(V))
  res <- matrix(NA, max(V), 2)
  for (i in 1:max(V)){
    regout <- robustse(lm(Y[V==i]~Z[V==i]))
    res[i,] <- regout[row,1:2]
    if (print.out) {print(regout)}
  }
  return(res)
}

# function to perform 2SLS regression
manualiv <- function (instrument,treatment,outcome){
  # written by Baobao Zhang
  first <- lm(treatment ~ instrument) # first stage
  treatment.hat <- first$fitted
  second <- lm(outcome ~ treatment.hat)
  X.adj <- cbind(rep(1, length(treatment)), treatment)
  second$residuals <- outcome - X.adj %*% coef(second)
  return(robustse(second))
}

# function to perform 2SLS regression and extract the coefficients
myreg.iv <- function(Y, D=d$D, Z=d$Z, V=d$V, row=2, print.out=TRUE){
  Y <- as.numeric(eval(Y))
  Z <- as.numeric(eval(Z))
  V <- as.numeric(eval(V))
  D <- as.numeric(eval(D))
  res <- matrix(NA, max(V), 2)
  for (i in 1:max(V)){
    regout <- manualiv(instrument = Z[V==i], 
                       treatment = D[V==i], outcome = Y[V==i])
    res[i,] <- regout[row,1:2]
    if (print.out) {print(regout)}
  }
  return(res)
}

#########################################
# Graphics Functions

# theme for graphics
theme_bb <- function (base_size = 12, base_family = "",where="bottom") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position=where, 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), strip.background = element_rect(fill = "grey80", 
                                                                                                                colour = "grey50"))
}

# Plot density function
myden <- function(V=d$V, Z=d$Z, Y, v, glab=c("Non-democracy", "Democracy"),
                  color1="orange1", color1b="orange3",
                  color2="steelblue4", ylab="Density",
                  alpha=0.75, adjust=1.5, mytheme=theme_bb(),
                  x.limits = NULL, x.breaks = NULL, x.labels = NULL,
                  title, xlab){
  V <- as.numeric(eval(V))
  Z <- eval(Z)
  Y <- eval(Y)
  data <- data.frame(Z[V==v], Y[V==v])
  names(data) <- c("x", "y")
  if (is.null(x.breaks) & is.null(x.labels) & is.null(x.limits)){
    ggplot(data, aes(x=y, fill=x)) + 
      geom_density(alpha=alpha,adjust=adjust) +
      geom_vline(xintercept = mean(data$y[data$x=="Non-democracy"],
                                   na.rm=T), 
                 linetype="longdash", colour=color1b, size=1)+
      geom_vline(xintercept = mean(data$y[data$x=="Democracy"],
                                   na.rm=T), 
                 linetype="longdash", colour=color2, size=1) + ylab(ylab)+
      xlab(xlab) + mytheme + ggtitle(title) +
      scale_fill_manual(name="Regime Type", 
                        values=c(color1, color2))
  } else{
    ggplot(data, aes(x=y, fill=x)) + 
    geom_density(alpha=alpha,adjust=adjust) +
      geom_vline(xintercept = mean(data$y[data$x=="Non-democracy"],
                                   na.rm=T), 
                 linetype="longdash", colour=color1b, size=1)+
      geom_vline(xintercept = mean(data$y[data$x=="Democracy"],
                                   na.rm=T), 
                 linetype="longdash", colour=color2, size=1) + ylab(ylab)+
      scale_x_continuous(limits=x.limits, breaks=x.breaks,
                         labels=x.labels)+
      xlab(xlab) + mytheme + ggtitle(title) +
      scale_fill_manual(name="Regime Type", 
                        values=c(color1, color2))
  }
}

# Plot histograms
myhist <- function(X, Y, glab=c("Non-democracy", "Democracy"), 
                   color1="orange1", color2="steelblue4",
                   outline="black", mytheme=theme_bb(),
                   title, xlab, ylab){
  data <- na.omit(data.frame(eval(X), eval(Y)))
  data$x <- factor(data[,1], labels=glab)
  data$y <- factor(data[,2], labels=c("No", "Yes"))
  ggplot(data=data, aes(x=y, fill=x))+
    geom_histogram(position="dodge", color=outline)+
    xlab(xlab)+ylab(ylab)+
    mytheme+
    ggtitle(title)+
    scale_fill_manual(values=c(color1, color2),
                      name="Regime Type")
}

# ggplot2 themes


theme_bb <- function (base_size = 12, base_family = "",where="bottom") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position=where, 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), strip.background = element_rect(fill = "grey80", 
                                                                                                                                                                                                                                                                                                                                    colour = "grey50"))
}


theme_bare <- function (base_size = 12, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position="none", 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), strip.background = element_rect(fill = "grey80", 
                                                                                                                colour = "grey50"))
}

theme_bb3 <- function (base_size = 12, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position="right", 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), strip.background = element_rect(fill = "grey80", 
                                                                                                                colour = "grey50"))
}


theme_bb4 <- function (base_size = 12, base_family = "", height, width) 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position="right", 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), 
          strip.background = element_rect(fill = "grey80", colour = "grey50"),
          legend.key.height = unit(height, "cm"),
          legend.key.width = unit(width, "cm"))
}

theme_bbtop <- function (base_size = 12, base_family = "", height, width) 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position="top", 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), 
          strip.background = element_rect(fill = "grey80", colour = "grey50"))
}

theme_nolegend <- function (base_size = 12, base_family = "", height, width) 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          legend.position="none", 
          axis.ticks = element_line(colour = "black"), 
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey90", size = 0.2), 
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80",  colour = "grey50"), 
          strip.background = element_rect(fill = "grey80", colour = "grey50"))
}

