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
			#Algae: consumption of N2 is mortality of Diazotroph, competition
			#for light via aij*D
			dA = A * ( rA *(cA1*N1 +cA2*D*muD)  - aij*D- muA ) 
			#Diazotroph: Competition for light via aji*A
			dD = D * ( rD *(cD2*N2) - aji*A - muD )			 

			####Resource dynamics 
			# Birth - Death = Logistic growth - consumption
			dN1 = N1 * ( gN1 *(1 - N1/K_N1 ) - cA1*A )  
			
			# This is an internal process of D. 
			# How does Nitrogen "grow" in availability for D? 
			dN2 = N2 * ( gN2 *(1 - N2/K_N2 ) - cA2*D ) 

	  	list( c(dA,dD,dN1,dN2) )
		})	

}  

#=============================================================================
# Set values of the population parameters
#=============================================================================
rA = 2 #Algae intrinsic growth
cA1 = 0.5 #Algal consumption of N1
cA2 = 0.5 #Algal consumption of N2
aij = 0.8 #Light competition from D
muA = 0.5 # Mortality of A

rD =  1 #Diazotroph intrinsic growth
cD2 = 1 #D consumption of N2
aji = 0.8 #Light competition from A
muD = 0.8 #D mortalityt

gN1 =  100 #Growth of resource N1
K_N1 = 100  #Carrying capacity of resource N1
gN2 = 100 #Growth of resource N2 
K_N2 = 100 #Carrying capacity of resource N2

parms = list(
			rA = rA, cA1 = cA1, cA2 = cA2, aij = aij, muA = muA, 
			rD = rD, cD2 = cD2, muD = muD, aji = aji, 
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
ad_out = ode(y=minit, times=times, func=ad_mod, parms=parms)

#=============================================================================
# Plot
#=============================================================================

plot( ad_out [,2], t="l", ylab = "Time", xlab = "Population density")
cuse = c("red","blue","green")
for ( n in 3:5) {
	lines (ad_out[,n], col = cuse[n] )
}


