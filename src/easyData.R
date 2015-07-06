#---------------------------------------------------------------------------
#
# Functions for making very simple data sets for testing lpSVM material.
#
# Copyright (C) Sean Holden 2014-15.
#
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Make some easy data using two normals.
#---------------------------------------------------------------------------
easyData <- function(m1, m2) {
	 
    require("MASS")

    y <- c(rep(1, m1),
           rep(-1, m2))

    x1 <- mvrnorm(m1, c(1.5,3), matrix(c(1,0.2,0.2,5), nrow=2, ncol=2))
    x2 <- mvrnorm(m2, c(-2,-2), matrix(c(1,0,0,1),     nrow=2, ncol=2))
    x <- rbind(x1,x2)

    p <- sample(1:(m1+m2))
    y <- y[p]
    x <- x[p,]

    list(X=x, y=y)
}

#---------------------------------------------------------------------------
# Make some easy data using parity.
#---------------------------------------------------------------------------
easyDataP <- function(m1, m2) {
	 
    require("MASS")

    y <- c(rep(1, 2*m1),
           rep(-1, 2*m2))

    x1 <- mvrnorm(m1, c(1.5,3),   matrix(c(1,0.2,0.2,5), nrow=2, ncol=2))
    x2 <- mvrnorm(m1, c(-2,-2), matrix(c(1,0,0,1), nrow=2, ncol=2))
    x3 <- mvrnorm(m2, c(1,-1),   matrix(c(2,0.2,0.1,3), nrow=2, ncol=2))
    x4 <- mvrnorm(m2, c(-2.5,1.5), matrix(c(1,0,0,0.75), nrow=2, ncol=2))
    x <- rbind(x1,x2,x3,x4)

    p <- sample(1:(2*(m1+m2)))
    y <- y[p]
    x <- x[p,]

    list(X=x, y=y)
}

#---------------------------------------------------------------------------
# Display the data produced by easyData.
#---------------------------------------------------------------------------
plotEasyData <- function(data) {

    y <- data$y
    x1 <- data$X[y==1,]
    x2 <- data$X[y==-1,]
  
    range1 <- range(data$X[,1])
    range2 <- range(data$X[,2])

    plot(range1, range2, type="n")
    points(x1, col="red", pch=21)
    points(x2, col="blue", pch=22)
}

#---------------------------------------------------------------------------
# Make 3-class data.
#---------------------------------------------------------------------------
easyData3 <- function(m1,m2,m3) {
	require("MASS")

    y <- c(rep(1, m1), rep(2, m2), rep(3, m3))

    x1 <- mvrnorm(m1, c(1.5,3),   matrix(c(1,0.2,0.2,5), nrow=2, ncol=2))
    x2 <- mvrnorm(m2, c(-2,-2), matrix(c(1,0,0,1), nrow=2, ncol=2))
    x3 <- mvrnorm(m3, c(1,-1),   matrix(c(2,0.2,0.1,3), nrow=2, ncol=2))
    x <- rbind(x1,x2,x3)

    p <- sample(1:(m1+m2+m3))
    y <- y[p]
    x <- x[p,]

    list(X=x, y=y)
}

#---------------------------------------------------------------------------
# Plot 3-class data.
#---------------------------------------------------------------------------
plotData3 <- function(data) {
    y <- data$y
    x1 <- data$X[y==1,]
    x2 <- data$X[y==2,]
    x3 <- data$X[y==3,]
    range1 <- range(data$X[,1])
    range2 <- range(data$X[,2])

    plot(range1, range2, type="n")
    points(x1, col="red", pch=21)
    points(x2, col="blue", pch=22)
    points(x3, col="green", pch=23)
}

#---------------------------------------------------------------------------
# Make 4-class data.
#---------------------------------------------------------------------------
easyData4 <- function(m1,m2,m3,m4) {
	require("MASS")

    y <- c(rep(1, m1), rep(2, m2), rep(3, m3), rep(4,m4))

    x1 <- mvrnorm(m1, c(1.5,3),   matrix(c(1,0.2,0.2,3), nrow=2, ncol=2))
    x2 <- mvrnorm(m2, c(-2,-2), matrix(c(1,0,0,1), nrow=2, ncol=2))
    x3 <- mvrnorm(m3, c(1,-1),   matrix(c(1.5,0.2,0.1,0.5), nrow=2, ncol=2))
    x4 <- mvrnorm(m4, c(-2.5,1.5), matrix(c(1,0,0,0.75), nrow=2, ncol=2))
    x <- rbind(x1,x2,x3,x4)

    p <- sample(1:(m1+m2+m3+m4))
    y <- y[p]
    x <- x[p,]

    list(X=x, y=y)
}

#---------------------------------------------------------------------------
# Plot 4-class data.
#---------------------------------------------------------------------------
plotData4 <- function(data) {
    y <- data$y
    x1 <- data$X[y==1,]
    x2 <- data$X[y==2,]
    x3 <- data$X[y==3,]
    x4 <- data$X[y==4,]
    range1 <- range(data$X[,1])
    range2 <- range(data$X[,2])

    plot(range1, range2, type="n")
    points(x1, col="red", pch=21)
    points(x2, col="blue", pch=22)
    points(x3, col="green", pch=23)
    points(x4, col="purple", pch=24)
}

