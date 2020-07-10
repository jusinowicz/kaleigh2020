#=============================================================================
# Example code for fitting an exponential model of the general form 
#			y = a*exp (b*x)
# in R by fitting a linear model (LM) to log-transformed data. 
#
# 1. Sample data is generated from a simple exponential model with added 
#	 noise and some points removed. 
# 2. Data are log-transformed
# 3. A LM is fit to the log-transformed data.
# 4. Check out the fit.
# 5. Simulate data from the fitted model and graph it all together. 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(lubridate)

#=============================================================================
#1. Simulate some data. Use an exponential of the general form y = a*exp (b*x). 
# This could represent e.g. Kleiber's law when b = -1/kT, and a = B0*M^(3/4)
#=============================================================================

as = 10
bs = 0.5
xs = seq(0.1,10,0.1)
gn = rnorm(length(xs) ) #Add Gaussian noise
ys = as*exp(bs*xs+gn ) #Case A: noise added inside the exp()

ggplot( data.frame(xs=xs,ys=ys), aes(x=xs,y=ys)) +geom_point()

#Note: The variance scales with the value of X: as X increases, the variance 
#also increases. This is a classic sign of an exponential process and indicates 
#that the data should be log-transformed. Notice the difference with: 
gn2 = rnorm(length(xs),1,20 ) #Add Gaussian noise with larger variance
ys2 = as*exp(bs*xs )+gn2 #Case B: noise added outside the exp()
ggplot( data.frame(xs=xs,ys=ys2), aes(x=xs,y=ys2)) +geom_point()
#Just something to keep in mind in the future when you're thinking about how 
#to fit data. Does your variance look more like case A, or case B? 

#=============================================================================
#2. Log-transform data, i.e. take the log of the response, y: 
#	Taking the log of y = a*exp(b*x) shows how this becomes linear: 
#    log(y) = log(a*exp(b*x)) [ use that log (a*b) = log(a)+log(b) ]
#	 		= log(a) + log(exp(b*x))  [use that log (exp(x) ) = x]
#			= log(a) + b*x 
#
#  There are two important things to note about this relationship: 
#  		The intercept is in term of log(a), so we will ned to transform it back 
#		with exp(). 
#		The coefficient b is in terms of the original scale, x, so b will 
#		already be on the correct scale. 
#  These relationships will be useful when working backwards from the fitted LM
#=============================================================================
lys = log(ys)

#=============================================================================
#3. Fit a linear model to the transformed data
#=============================================================================
model1 = lm ( lys~xs )

#=============================================================================
#4. Check out the fit 
#=============================================================================
summary (model1)
coef(model1) #Just the coefficients 

#=============================================================================
# 5. Simulate data from the fitted model and graph it all together.
#=============================================================================
#	 Based on the relationships in (2), the fitted intercept stored in 
#	 coef(model1)[1] is the log(a): Taking exp( )
af = exp(coef(model1)[1] )

#	 Based on the relationships in (2), the fitted coefficient b stored in 
#	 coef(model1)[2] is b
bf = coef(model1)[2]

############################################################
#There are two ways to fit a model which are equivalent in this case. 
############################################################
#First, use the fitted coefficients and the model structure:  
yp1 = af*exp(bf*xs) 

#Second, use predict. It is worth learning this approach at some point since it 
#is MUCH easier to use with more complicated model structures (e.g. GLMMs). 
#Predict needs a data.frame with all variables, with same names as supplied to
#the model in lm( ) 
newdata = data.frame(xs = xs) 
#Predicted data (i.e. the fitted curve over the specified interval)
yp2 = predict(model1, newdata=newdata, type="response")
#Note: This will be on the linear scale! Just do exp(yp2) to get it on the 
#original scale

#Version 1, on the original scale 
ggplot() + geom_point( data = data.frame(xs=xs,ys=ys), 
			mapping = aes(x=xs,y=ys, color= "Data")  ) + #Data
		   geom_line( data = data.frame(xs=xs,yp=yp1), 
		   	mapping = aes(x=xs, y=yp1, color = "Fit")  ) 

#Version 2, The linear log-log plot:
ggplot() + geom_point( data = data.frame(lxs=xs,lys=lys), 
			mapping = aes(x=xs,y=lys, color= "Data")  ) + #Data
		   geom_line( data = data.frame(xs=xs,yp=yp2), 
		   	mapping = aes(x=xs, y=yp2, color = "Fit")  ) 
#Note how the points are about equally distant above and below the fitted line in this
#graph. If your points are generally closer at one end of the line then the model structure 
#might not be correct. 


