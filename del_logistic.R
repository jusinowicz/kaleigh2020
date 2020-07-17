#=============================================================================
#Implementation of a logistic growth function with a time delay
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(deSolve)

#=============================================================================
# Define the population dynamics through the following function
#=============================================================================
del_log= function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{		
		if (times < tau){ 
			lag1 = 1
		}else{
			lag1 = lagvalue(times-tau)
		}
			# Birth - Death = Logistic growth - consumption
			dN1 = N1 * ( gN1 *(1 - lag1/K_N1 ) )  

	  	list( c(dN1) )
		})	

}  

#=============================================================================
# Set values of the population parameters and time lag tau 
#=============================================================================
gN1 =  1.2#Growth of resource N1
K_N1 = 10  #Carrying capacity of resource N1
tau = 1 #Time lag of density dependent feedback. 

parms = list(
			gN1 = gN1, K_N1 = K_N1, tau=tau
		 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 100
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
minit =  (c(N1 = 1))

del_out = dede(y=minit, times=times, func=del_log, parms=parms, atol = 1e-9)

#=============================================================================
# Plot
#=============================================================================
plot( del_out [,2], t="l", ylab = "Time", xlab = "Population density")


