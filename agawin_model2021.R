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
	with( as.list(c(parms, sp )), #args = parms (list of parameters), sp (matrix of outputs (Pf, Pn, R) -- add a line of values at each time step)
		{	
			##### Determine species-specific average growth rates as a function
			##### of resource availibility and light intensity (eqs. 7, 9-10)
			Iout = Iin *(exp( zM*(-Kbg - kN*Pn - kF*Pf) ) ) #eq 7 at zM #Iout = average amount of light as a function of depth
			       #k = light attenuation, P is abundance, for species N and F (F is fixer) 
		
			#Average growth rate of Pn (eq 9)
			gn_ave = gmn * (R/(Mn+R))*( ( log(Hn+Iin) - log(Hn+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			#Average growth rate of Pf (eq 10)
			# gf_ave = (gmfn2 *Mf + gmfno3*R)/(Mf+R)*( ( log(Hf+Iin) - log(Hf+Iout) ) / 
			# (log(Iin)-log(Iout)) )

			#Average growth rates for each nitrogen source, separated
			gfno3_ave =  gmfno3 *  (R/(Mfno3+R))*( ( log(Hfno3+Iin) - log(Hfno3+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			gfn2_ave = gmfn2 * (R/(Mfn2+R))*( ( log(Hfn2+Iin) - log(Hfn2+Iout) ) / 
			(log(Iin)-log(Iout) ) )

			gf_ave = gfn2_ave + gfno3_ave

			##### Consumer dynamics 
			#Algae (Eq 1) : 
			dPn = Pn * ( gn_ave - mn) #dPn calculates change at next time point, so uses all variables from the last time step; Jacob has factored out PN; this is stored as the new population size at this time step
			#N fixer (Eq 2): 
			dPf = Pf * ( gf_ave - mf )			 
			
			####Resource (N) (Eq 3): 
			dR = D*(Rin - R) - Qn*gn_ave*Pn - Pf*gfno3_ave*Qf + ef*Pf*gfn2_ave*Qf 


	  	list( c(dPn,dPf,dR) ) #also store these outputs in the matrix; these are the "dynamic variables", which they calculate as the growth rates they are and then back calculate N at current time step from

	  	list( c(dPn,dPf,dR, Iout=Iout, gn_ave=gn_ave, gfno3_ave  = gfno3_ave,gfn2_ave =gfn2_ave ) )

		})	

}  


#=============================================================================
# Set values of the population parameters (high light scenario; i.e. no synechococcus)
#=============================================================================
###From table 3
###Non-fixer (there are two of these, chlorella for now)
mn = 0.014 #mortality

#These terms are all for the growth rate gn_ave. 
gmn = 0.06 #Max gr
Mn = 0.308 #Half saturation constants 
Hn = 80

###Fixer
mf = 0.014 #mortality

#These terms are all for the growth rates gfno3 and gfn2. There are two sets of
#these for high and low light: 
#Use the high light treatment for now.
gmfno3 = .084
gmfn2 =  0.025
#Half saturation constants 
#Don't assume that Hs and Ms are set to be equal, as they are by the authors. 
Hfno3 = 70
Hfn2 = 70
Mfno3 = 1
Mfn2 = 0.5

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
#container. these values in text of table 1
Kbg = 4.75 #4.75 #Background turbidity. In a few cases this is 11
zM = 0.05 #mixing depth
kF = 4.86 #light attenutation constants. 
kN = 5.68 #This is determined by eq 13 for Synechococcus
ef = 8 #epsilon

#Put these parameters all together in a list to pass to deSolve with the model
#definition. 
parms = list(
			mn = mn, gmn = gmn, Mn = Mn, Hn = Hn, 
			mf = mf, gmfn03 = gmfn03, gmfn2 = gmfn2, Hfno3  = Hfno3 , Mfno3  = Mfno3 ,  
			Hfn2  = Hfn2, Mfn2  = Mfn2 , 
			Qn = Qn, Qf = Qf, D = D, Iin = Iin, Rin = Rin, Kbg = Kbg, zM = zM, kF = kF,
			kN = kN, ef = ef
		 )

#=============================================================================
# Run the model with initial conditions (IC)
#=============================================================================
tend = 3000 #t final
delta1 = 0.1 #size of time step
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times) 
#Use ICs to help set the scenario. E.g., when Pn or Pf =0, 
#this is a monoculture experiment. Otherwise, use ICs from those reported in 
#Table 1

minit =  c(Pn = 10,Pf = 0.16, R = Rin) #model initial values

Pn_ic =0.04
Pf_ic =0.16
Iout = Iin-1
gn_ic = gmn * (Rin/(Mn+Rin))*( ( log(Hn+Iin) - log(Hn+Iout) ) / 
			(log(Iin)-log(Iout) ) )
gfno3_ic = gmfno3 *  (Rin/(Mfno3+Rin)	)*( ( log(Hfno3+Iin) - log(Hfno3+Iout) ) / 
			(log(Iin)-log(Iout) ) )
gfn2_ic = gmfn2 * (Rin/(Mfn2+Rin))*( ( log(Hfn2+Iin) - log(Hfn2+Iout) ) / 
			(log(Iin)-log(Iout) ) )

minit =  c(Pn = Pn_ic ,Pf = Pf_ic, R = Rin, Iout=Iout, gn_ave=gn_ic, 
				gfno3_ave  = gfno3_ic,gfn2_ave = gfn2_ic )


#Run the model
agawin_out = ode(y = minit, times = times, func = agawin_mod, parms = parms)

#=============================================================================
# Plot
#=============================================================================
#Black is Pn (Chorella) 
#Red is Pf (Cyanothece)
#Blue is R (N)

plot( agawin_out [,2], t="l", ylab = "pop density", xlab = "time", 
	ylim = c(0, max(agawin_out[,2:4]) ) )
cuse = c("black","black","red","blue","green")
for ( n in 3:4) {
	lines (agawin_out[,n], col = cuse[n] )
}

#get it read to plot with ggplot
agawin_out <- as.data.frame(agawin_out)

#plot
agawin_out %>% 
  ggplot() +
  geom_point(aes(x = time, y = Pn), colour = "green") +
  # geom_point(aes(x = time, y = gfno3_ave), colour = "blue") +
  # geom_point(aes(x = time, y = gfn2_ave), colour = "green") +
  geom_point(aes(x = time, y = Pf), colour = "orange") + #total growth rate of fixer
  geom_point(aes(x = time, y = R), colour = "black")
  