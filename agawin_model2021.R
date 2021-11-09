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
# This model includes one consumers: Pn and Pf 
# 	Pn is a consumer of Nitrogen and light. 
#	Pf is a consumer of Nitrogen and light, but also fixes N. 
# This model includes one resource: R 
# 	R can be added at a constant rate. 
#	R is produced by the mortality of both Pn and Pf
#
# The growth rates of both species are determined by light intensity. 
# Light intensity is a function of depth, and decreases as densities of both
# species increase. The extra set of equations in the model definition 
# calculate a new average growth rate each time step based on these factors.  
#=============================================================================

agawin_mod = function(times,sp,parms){
	with( as.list(c(parms, sp )), #args = parms (list of parameters), sp (matrix of outputs (N) -- add a line each time step)
		{	
			##### Determine species-specific average growth rates as a function
			##### of resource availibility and light intensity (eqs. 7, 9-10)
			Iout = Iin *(exp( zM*(-Kbg - kN*Pn - kF*Pf) ) ) #eq 7 at zM #Iout = average amount of light as a function of depth
			       #k = light attenuation, P is abundance, for species N and F (F is fixer) 
		
			#Average growth rate of Pn (eq 9)
			gn_ave = gmn * (R/(Mn+R))*( ( log(Hn+Iin) - log(Hn+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			#Average growth rate of Pf (eq 10)
			gf_ave = (gmfn2 *Mf + gmfn03*R)/(Mf+R)*( ( log(Hf+Iin) - log(Hf+Iout) ) / 
			(log(Iin)-log(Iout)) )

			#Average growth rates for each nitrogen source, separated
			gfno3_ave = gmfn2 * (R/(Mf+R))*( ( log(Hf+Iin) - log(Hf+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			gfn2_ave = gmfn03 * (R/(Mf+R))*( ( log(Hf+Iin) - log(Hf+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			##### Consumer dynamics 
			#Algae (Eq 1) : 
			dPn = Pn * ( gn_ave - mn) #dPn calculates change at next time point, so uses all variables from the last time step; Jacob has factored out PN
			#N fixer (Eq 2): 
			dPf = Pf * ( gf_ave - mf )			 
			
			####Resource (N) (Eq 3): 
			dR = D*(Rin - R) - Qn*gn_ave*Pn - Pf*gfno3_ave*Qf + ef*Pf*gfn2_ave*Qf

	  	list( c(dPn,dPf,dR) )
		})	

}  


#=============================================================================
# Set values of the population parameters
#=============================================================================
###From table 3
###Non-fixer (there are two of these, chlorella for now)
mn = 0.014 #mortality

#These terms are all for the growth rate gn_ave. 
gmn = 0.06 #Max gr
Mn = 0.308 #Half saturation constants 
Hn = 80

###Fixer
mf =0.014 #mortality

#These terms are all for the growth rates gfno3 and gfn2. There are two sets of
#these for high and low light: 
#Use the high light treatment for now.
gmfn03 = .084
gmfn2 =  0.025
#Half saturation constants 
#All of Hs and Ms are set to be equal by the authors. 
Hf = 70
Mf = 1
# Hfn2 = 70
# Hfno3 = 70
# Mfn2 = 1
# Mfno3 = 1 

###Resource
Qn = 0.079#cell quota of N
Qf = 0.092#cell quota of N

D = 0.014 #dilution rate in experiment. 
Iin = 40 #Light in, HL
Rin = 8 #These are the treatment levels of N input. 0, 0.1, 0.5, 8. 

###Parametes for light intensity calculation (eqs 7-10)
#Growth rate is spatially (i.e. depth) dependent. The average growth rates (eqs 4 - 6) for 
#each species are light dependent. Light intensity is dependent on population 
#densities, depth, and some physical parameters as described by equation 7.
# 
#The integral in eq. 8 could be implemented directly in this code. However, 
#the authors show that this integral can be solved analytically which leads to 
#equations 9-10. These equations only require the endpoints of light intensity: 
#i.e., the light at depth 0, and at the mixing depth z_M, at the bottom of the 
#container.
Kbg = 4.75 #4.75 #Background turbidity. In a few cases this is 11
zM = 0.05 #mixing depth
kF = 4.86 #light attenutation constants. 
kN = 5.68 #This is determined by eq 13 for Synechococcus
ef = 8 #

#Put these parameters all together in a list to pass to deSolve with the model
#definition. 
parms = list(
			mn = mn, gmn = gmn, Mn = Mn, Hn = Hn, 
			mf = mf, gmfn03 = gmfn03, gmfn2 = gmfn2, Hf = Hf, Mf = Mf,  
			Qn = Qn, Qf = Qf, D= D, Iin = Iin, Rin = Rin, Kbg = Kbg, zM = zM, kF = kF,
			kN = kN, ef=ef
		 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 3000
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
#Use ICs to help set the scenario. E.g., when Pn or Pf =0, 
#this is a monoculture experiment. Otherwise, use ICs from those reported in 
#Table 1
minit =  c(Pn = 0.04,Pf = 0.16, R = Rin)

#Run the model
agawin_out = ode(y=minit, times=times, func=agawin_mod, parms=parms)

#=============================================================================
# Plot
#=============================================================================
#Black is Pn (Cyanothece) 
#Red is Pf (Chlorella)
#Blue is R (N)

plot( agawin_out [,2], t="l", ylab = "Time", xlab = "Population density", 
	ylim = c(0, max(agawin_out[,2:4]) ) )
cuse = c("black","black","red","blue","green")
for ( n in 3:4) {
	lines (agawin_out[,n], col = cuse[n] )
}


