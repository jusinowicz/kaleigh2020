#=============================================================================
# R code to implement and explore an ODE for resource-consumer dynamics 
# between algae and photosynthetic diazotrophs. The model assumes that two
# Nitrogen pools are present: 1 in the water column that only algae consume, 
# and 1 that is internal to nitrogen-fixing diazotrophs. This second source
# can become available to algae when diazotrophs die. 
# 
# Eventually the demographic rates will become temperature dependent. 
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(deSolve)


#=============================================================================
# Define the population dynamics through the following function
# This model includes two consumers: A(lgae) and D(iazotrophs)
# This model includes two resources: N1 and N2
# A consumes N1 and N2 in the form of dead D
# D consumes N2 only from an internal fixation process
# Competition for light is introduced as a per-capita shading factor (essentially
# an interspecific competition coefficient)
#=============================================================================

ad_mod = function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{	#####Consumer dynamics
		dA = A * ( rA *(cA1*N1 +cA2*D*muD)  - aij*D- muA ) #Algae: consumption of N2 is mortality of Diazotroph
		dD = D * ( rD *(cD2*N2) - aji*A - muD )			 #Diazotroph

		####Resource dynamics 
		dN1 = N1 * ( gN1 *(1 - N1/K_N1 ) - cA1*A )  #Logistic growth - consumption
		dN2 = N2 * ( gN2 *(1 - N2/K_N2 ) - cA2*D ) #How does Nitrogen "grow" in availability for D? 
	  	list( c(dA,dD,dN1,dN2) )
		})	

}  

#=============================================================================
# Set values of the population parameters
#=============================================================================
rA = 2 
cA1 = 0.5
cA2 = 0.5
muA = 0.5 
rD =  1
cD2 = 1 
muD = 0.8 
gN1 =  100
K_N1 = 100  
gN2 = 100
K_N2 = 100 

parms = list(
			rA = rA, cA1 = cA1, cA2 = cA2, aij = 0.8, muA = muA, 
			rD = rD, cD2 = cD2, muD = muD, aji = 0.8, 
			gN1 = gN1, K_N1 = K_N1,  
			gN2 = gN2, K_N2 = K_N2 
		 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 100
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
minit =  c(A =1, D=1, N1 = K_N1, N2=K_N2 )
out_temp = ode(y=minit, times=times, func=ad_mod, parms=parms)
