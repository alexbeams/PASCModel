rm(list=ls()) #clears the working directory to start from a clean slate each time


# This code is similar to pasc.R but includes a highly susceptible sub-population to
# examine the effects of heterogeneous susceptibility to infection

# Information we know which places constraints on parameters in model:
# vaccine efficacy against PASC, = 1-(1-epsilon)*(1-eta), 
# eta=efficacy against infection, epsilon=efficacy against PASC conditional on breakthrough infection
# Q <- 1/50 #fraction of all infections that have resulted in long COVID
# I_active <- .02 #fraction of population with active COVID infection

# functions:

# 1. getout: produces a dataframe called "out" containing the ODE solution. It includes
#    a state uexp which solves for PASC assuming exponential duration, as well as Sigma,
#    the per-capita PASC acquisition rate for the integral equation. It needs all of the
#    parameters fed in, except for the latency and precision of the integral eqns. 
# 2. getu: acts on out to produce PASC prevalence for different distributions. It needs Sigma, 
#    as well as the mean duration of PASC. You can tune the precison 
#    of the PASC distributions independent of the ODE model for transmission
#    on the assumption that PASC doesn't affect transmission

getout <- function(tau,eta,epsilon,Q,I_active,gamma=1,delta=1/52,alpha=1/52,s=0.1){

	# This code solves the heterogeneous SIP ODE model
 	
	require(deSolve)	#load in the ODE solver package

	#run the simulation for 100 weeks to eliminate transients - what time should the variant emerge?
	Tinv <- 1000

	##############
	# Parameters #		#description of parameter, [units of parameter]
	##############


	# Model parameters
	VE <- 1-(1-eta)*(1-epsilon)

	pars <- list(
	beta    = 0,            # depends on eta;  transmission rate, [1/population density/week]
	deltabeta_low=0.0,          # change in transmisison rate, [1/population density/week]
	deltabeta_high=0.0,          # change in transmisison rate, [1/population density/week]
	s = s,                # selection coefficient
	sw = 0,
	gamma   = gamma,            # recovery rate, [1/week]
	delta   = delta,        # rate of adaptive immune loss, [1/week]
	deltaeta = 0.0,	        # change in immune efficacy, [dimensionless]
	alpha   = alpha,         # per-capita vaccination rate, [1/week]
	tau     = tau,           #mean duration of long Covid, [week]
	phi     = 0,         #probability of developing long Covid from S0, [dimensionless]
	epsilon  = 0,           #efficacy against PASC conditional on breakthrough infection (in immunized people), [dimensionless]
	L	= 10,		#mean delay between infection and onset of PASC, [weeks]
	mu      = 0/52/75)      #per-capita birth/death rate, [1/week], don't have the formula from Overleaf with this part yet

	# calculate beta to satisfy calibration target I*, assuming gamma, delta, 
	# alpha are known, and we want to supply various values of eta.
	# However, there is a subpopulation of size p_risk with elevated susceptibility rho 
	# choose rho to accomodate I_active, delta, alpha, gamma, eta, p_risk, and beta
	# should be the case that when rho=1, we get the same beta value as before

	# this is how we chose beta in the homogeneous case:
	getBeta_old = function(eta){with(pars,{ 
		a = (eta-1)*I_active^2 - eta*I_active+I_active
		b = (eta-1)*(I_active-1)*alpha - (delta+gamma)*I_active + delta
		c = -gamma*(alpha+delta)

		return( max( (-b + sqrt(b^2-4*a*c))/2/a,
		(-b - sqrt(b^2-4*a*c))/2/a))
		
	})}
	
	# this is assuming I_active = 0.02 as in figure 3. Need to pull the more complete formula from Maple if we want to use that
	getzeta = function(beta,eta,sigma){ 
	zeta= (3.773698829*10^(-14)*sqrt(((1.56249999900001*10^30 + 4.569760001*10^29*beta^4 + 3.51520000000000*10^30*beta^3 + 8.44999999700001*10^30*beta^2 + 6.49999999700001*10^30*beta)*eta^2 + (-9.139520001*10^29*beta^4 - 7.03040000000000*10^30*beta^3 - 1.85900000000000*10^31*beta^2 - 1.95000000000000*10^31*beta - 6.25000000000000*10^30)*eta + 4.569760000*10^29*beta^4 + 3.51520000110308*10^30*beta^3 + 1.01400000008924*10^31*beta^2 + 1.30000000122887*10^31*beta + 6.25000000783408*10^30)*sigma^2 + ((-5.167344000*10^28*beta^4 - 1.15934000000000*10^29*beta^3 - 9.55500000000000*10^28*beta^2 - 3.06250000000000*10^28*beta)*eta^3 + (-7.578771201*10^29*beta^4 - 3.91167400000000*10^30*beta^3 - 7.97225000000000*10^30*beta^2 - 4.75625000000000*10^30*beta - 6.25000000100000*10^28)*eta^2 + (1.670774560*10^30*beta^4 + 8.11740800000000*10^30*beta^3 + 1.25528000000000*10^31*beta^2 + 1.06250000000000*10^30*beta - 9.50000000000000*10^30)*eta - 8.612240001*10^29*beta^4 - 4.08980000002945*10^30*beta^3 - 4.48500000050159*10^30*beta^2 + 3.62499999962937*10^30*beta + 6.25000000093175*10^30)*sigma + 1.623076000*10^26*(beta + 0.9615384615)^2*beta^2*eta^4 + (4.966612560*10^28*beta^4 + 9.55308900000000*10^28*beta^3 + 4.65745000000000*10^28*beta^2 + 6.12500000100000*10^26*beta)*eta^3 + (4.84068499900000*10^29*beta^2 + 1.88075000000000*10^29*beta + 3.059498261*10^29*beta^4 + 5.94231820000000*10^29*beta^3 + 6.25000000100000*10^26)*eta^2 + (7.32662000100000*10^29*beta^2 + 1.35415000000000*10^30*beta - 7.615472593*10^29*beta^4 - 1.42262484000000*10^30*beta^3 + 1.92499999900000*10^29)*eta - 1.43750000027229*10^30*beta - 1.26187499952968*10^30*beta^2 + 7.32550000287585*10^29*beta^3 + 4.057690000*10^29*beta^4 + 1.56249999907113*10^30) + (-0.4807692308*eta^2 + (-25.51020408*sigma + 26.4423077)*eta + 25.51020408*sigma - 25.96153847)*beta^2 + (-0.4622781066*eta^2 + 25.39699311*eta - 23.43466973)*beta + (47.17123535*sigma - 0.9434247071)*eta - 94.34247071*sigma + 50.94493418)/(beta*((eta^2 - 2.*eta + 0.9999999998)*beta^2 + (0.9615384614*eta^2 + (-51.02040816*sigma - 1.86420722)*eta + 51.02040816*sigma + 0.9026687591)*beta + (-98.11616953*sigma + 1.96232339)*eta + 98.11616953*sigma - 1.96232339))

		return(zeta) 		
	}

	# looks to be working: getzeta(getBeta_old(.55),.55,.2) = 1, as it must
	# now, beta in the getzeta call describes the transmission rate experienced by low-risk individuals, and beta*zeta is for the high risk individuals
	# to capture the effect of heterogeneity for a given eta, choose smaller values of beta than getBeta(eta) - the high-risk people will have transmission rate beta*zeta
	# for eta=0.55, decr=0.65 produces an approximately 20-80 situation: 80% of new infections come from 20% of the population

	# let's say the low-risk people have a lower transmission rate,
	# smaller than the average by a factor decr:

	# transmission rate for low-risk group:
	getBeta_low <- function(eta,decr) getBeta_old(eta)* decr	

	# transmission rate for high-risk group:
	getBeta_high <- function(eta,decr) getBeta_low(eta,decr)*getzeta(getBeta_low(eta,decr),eta,0.2)

	#calculate S and P in terms of alpha, gamma, delta, eta, beta, zeta, sigma, and I_active

	getS_old <- function(eta){
		with(pars, {
			S = 1 + ((alpha-gamma)*I_active - alpha)/(getBeta(eta)*I_active*(1-eta)+alpha+delta)  - I_active
			return(S)
		})

	}

	getP_old <- function(eta){
		with(pars, {
			P =  -((alpha-gamma)*I_active - alpha)/(getBeta(eta)*I_active*(1-eta)+alpha+delta) 
			return(P)
		})
	}


	# this looks like it's working correctly. In these function calls, 
	# replacing beta=getBeta_low(eta) with getBeta_old(eta), 
	# we see getS_old(.55)=getS_het(.55,getzeta(getBeta_old(.55),.55,.2))
	# just as they must
	getS_het <- function(eta,decr){
		with(pars,{ 
			beta=getBeta_low(eta,decr)
			zeta=getzeta(getBeta_low(eta,decr),eta,0.2)
			S = -(beta*(zeta + 1)*(-1 + eta)*I_active^2 + (-(zeta + 1)*(-1 + eta)*beta + alpha*eta - alpha - delta - gamma)*I_active - alpha*eta + alpha + delta)*delta/(I_active^2*zeta*(-1 + eta)*beta^2 + I_active*(alpha + delta)*(zeta + 1)*(-1 + eta)*beta + (alpha + delta)*(alpha*eta - alpha - delta))
			return(S)
		})

	}

	getP_het <- function(eta,decr){
		with(pars,{ 
			beta=getBeta_low(eta,decr)
			zeta=getzeta(getBeta_low(eta,decr),eta,0.2)
			P = (-beta^2*zeta*(-1 + eta)*I_active^3 + beta*((-zeta - 1)*alpha + beta*zeta)*(-1 + eta)*I_active^2 + ((1 - eta)*alpha^2 + ((zeta + 1)*(-1 + eta)*beta + delta)*alpha - delta*gamma)*I_active + alpha*((-1 + eta)*alpha - delta))/(I_active^2*zeta*(-1 + eta)*beta^2 + I_active*(alpha + delta)*(zeta + 1)*(-1 + eta)*beta + (alpha + delta)*((-1 + eta)*alpha - delta))
			return(P)
		})

	}
	

	# notice that getS_het(.55) > getS_old(.55), and getP_het(.55) < getP_old(.55). Heterogeneity translates to
	# lower levels of immunity/higher levels of susceptibility in the population because infections are targeted toward the at-risk group
	# the incidence is concentrated among a subset of the population
	# increasing the variance in susceptibility makes this effect even more pronounced

	# we have our overall states S, P calibrated, but need to get the relative number of high- and low-risk individuals right

	getx_het <- function(eta,decr){
		with(pars,{ 
			beta=getBeta_low(eta,decr)
			zeta=getzeta(getBeta_low(eta,decr),eta,0.2)
			x = (-I_active*(-1 + eta)*(I_active - 1)*beta^2 + ((-alpha*eta + alpha + delta + gamma)*I_active + alpha*eta - alpha - delta)*beta + gamma*(alpha + delta))/(beta*(zeta - 1)*(I_active*(zeta + 1)*(-1 + eta)*(I_active - 1)*beta + (alpha*eta - alpha - delta - gamma)*I_active - alpha*eta + alpha + delta))
			return(x)
		})
	}


	gety_het <- function(eta,decr){
		with(pars,{ 
			beta=getBeta_low(eta,decr)
			zeta=getzeta(getBeta_low(eta,decr),eta,0.2)
			sigma=0.2
			numery = beta^3*eta*sigma*zeta^2*I_active^2 + alpha*beta^2*eta*sigma*zeta^2*I_active - beta^3*eta*sigma*zeta*I_active^2 - beta^3*eta*zeta*I_active^3 - beta^3*sigma*zeta^2*I_active^2 + beta^2*delta*eta*sigma*zeta^2*I_active - alpha*beta^2*eta*zeta*I_active^2 - alpha*beta^2*sigma*zeta^2*I_active + beta^3*eta*zeta*I_active^2 + beta^3*sigma*zeta*I_active^2 + beta^3*zeta*I_active^3 - beta^2*delta*sigma*zeta^2*I_active + alpha^2*beta*eta*sigma*zeta - alpha*beta^2*eta*sigma*I_active + alpha*beta^2*eta*zeta*I_active - alpha*beta^2*eta*I_active^2 + alpha*beta^2*zeta*I_active^2 + alpha*beta*delta*eta*sigma*zeta - beta^3*zeta*I_active^2 - beta^2*delta*eta*sigma*I_active - beta^2*delta*eta*I_active^2 + beta^2*delta*zeta*I_active^2 + beta^2*gamma*zeta*I_active^2 - alpha^2*beta*eta*sigma - alpha^2*beta*eta*I_active - alpha^2*beta*sigma*zeta + alpha*beta^2*eta*I_active + alpha*beta^2*sigma*I_active - alpha*beta^2*zeta*I_active + alpha*beta^2*I_active^2 - alpha*beta*delta*eta*sigma - alpha*beta*delta*eta*I_active - 2*alpha*beta*delta*sigma*zeta + alpha*beta*gamma*zeta*I_active + beta^2*delta*eta*I_active + beta^2*delta*sigma*I_active - beta^2*delta*zeta*I_active + beta^2*delta*I_active^2 - beta*delta^2*sigma*zeta + beta*delta*gamma*zeta*I_active + alpha^2*beta*eta + alpha^2*beta*sigma + alpha^2*beta*I_active - alpha*beta^2*I_active + alpha*beta*delta*eta + 2*alpha*beta*delta*sigma + 2*alpha*beta*delta*I_active + alpha*beta*gamma*I_active - beta^2*delta*I_active + beta*delta^2*sigma + beta*delta^2*I_active + beta*delta*gamma*I_active - alpha^2*beta + alpha^2*gamma - 2*alpha*beta*delta + 2*alpha*delta*gamma - beta*delta^2 + delta^2*gamma 
			denomy = (beta^2*eta*zeta*I_active^2 + alpha*beta*eta*zeta*I_active - beta^2*zeta*I_active^2 + beta*delta*eta*zeta*I_active + alpha*beta*eta*I_active - alpha*beta*zeta*I_active + beta*delta*eta*I_active - beta*delta*zeta*I_active + alpha^2*eta - alpha*beta*I_active + alpha*delta*eta - beta*delta*I_active - alpha^2 - 2*alpha*delta - delta^2)*beta*(zeta - 1)*I_active

			y = numery/denomy
			return(y)
		})
	}
	# now, we have our model calibrated with underlying heterogeneity

	getSL <- function(eta,decr) (1-getx_het(eta,decr))*getS_het(eta,decr)
	getSH <- function(eta,decr) getx_het(eta,decr)*getS_het(eta,decr)
	getIL <- function(eta,decr) (1-gety_het(eta,decr))*I_active
	getIH <- function(eta,decr) gety_het(eta,decr)*I_active
	getPL <- function(eta,decr) 0.8 - getSL(eta,decr) - getIL(eta,decr)
	getPH <- function(eta,decr) 0.2 - getSH(eta,decr) - getIH(eta,decr)
	# IH and IL are 0.2-SH-PH and 0.8-SL-PL, respectively
	
	#calibrate phi for the chosen value of eta, along with the background values of VE and fr (Q)
	# Q is the ratio of PASC incidence to total incidence, and is assumed known
	# some algebra:
	# Q = (phi*(beta_low*SL + beta_high*SH)*I + phi*(1-epsilon)*(1-eta)*I*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH)*I + (1-eta)*I*(beta_low*PL+beta_high*PH))
	# Q = (phi*(beta_low*SL + beta_high*SH) + phi*(1-epsilon)*(1-eta)*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH) + (1-eta)*(beta_low*PL+beta_high*PH))
	# Q = phi* ((beta_low*SL + beta_high*SH) + (1-epsilon)*(1-eta)*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH) + (1-eta)*(beta_low*PL+beta_high*PH))
	# Q = phi* (beta_s(x)*S + (1-epsilon)*(1-eta)*beta_p(y)*P) / (beta_s(x)*S+(1-eta)*beta_p(y)*P) 

	# phi = Q * (beta_s(x)*S+(1-eta)*beta_p(z)*P) / (beta_s(x)*S + (1-epsilon)*(1-eta)*beta_p(z)*P)


	getPhi <- function(eta,decr,onlyhighriskgetpasc){
		numerSterm = getBeta_low(eta,decr)*getSL(eta,decr)*onlyhighriskgetpasc + getBeta_high(eta,decr)*getSH(eta,decr)
		numerPterm = (1-eta)*( getBeta_low(eta,decr)*getPL(eta,decr)*onlyhighriskgetpasc + getBeta_high(eta,decr)*getPH(eta,decr) )		
		denomSterm = getBeta_low(eta,decr)*getSL(eta,decr)*1 + getBeta_high(eta,decr)*getSH(eta,decr)

		denomPterm = (1-epsilon)*(1-eta)*( getBeta_low(eta,decr)*getPL(eta,decr)*1 + getBeta_high(eta,decr)*getPH(eta,decr) )		

		Phi =  Q * (numerSterm+numerPterm)/(denomSterm+denomPterm)
		return(Phi)
	}


	#############
	# timesteps #
	#############

	# how long to run (weeks)
	tEnd <-2000		
	# timestep (weeks)
	DT <- 1

	##########################
	# differential equations #
	##########################

	F <- function(Time, state, Pars) {
		with(as.list(c(state, Pars)), {

			beta_low_fun = function(Time) beta_low + deltabeta_low * ( plogis(sw*(Time-Tinv)) )
			beta_high_fun = function(Time) beta_high + deltabeta_high * ( plogis(sw*(Time-Tinv)) )
			Eta  = function(Time) eta - deltaeta * ( plogis(sw*(Time-Tinv))) 
	
			dSL	<- -beta_low_fun(Time)*SL*(IL+IH) + delta*PL + mu*0.8 - mu*SL - alpha*SL
			dIL	<-  beta_low_fun(Time)*SL*(IL+IH) + (1-Eta(Time))*beta_low_fun(Time)*PL*(IL+IH) - gamma*IL  - mu*IL
			dPL	<-  gamma*IL  - (1-Eta(Time))*beta_low_fun(Time)*PL*(IL+IH) - delta*PL  - mu*PL  + alpha*SL
			dSH	<- -beta_high_fun(Time)*SH*(IL+IH) + delta*PH + mu*0.2 - mu*SH - alpha*SH
			dIH	<-  beta_high_fun(Time)*SH*(IL+IH) + (1-Eta(Time))*beta_high_fun(Time)*PH*(IL+IH) - gamma*IH  - mu*IH
			dPH	<-  gamma*IH  - (1-Eta(Time))*beta_high_fun(Time)*PH*(IL+IH) - delta*PH  - mu*PH  + alpha*SH

			Sigma	<- phi*(beta_low*SL + beta_high*SH+(1-Eta(Time))*(1-epsilon)*(beta_low*PL+beta_high*PH))*(IL+IH)

			duexp	<- Sigma*(1-uexp)-uexp/tau

		return(list(c(dSL,dIL,dPL,dSH,dIH,dPH,duexp),c(Sigma)))
		})
	}


	
	#########################
	# set up and run solver #
	#########################

	#use the parameters calibrated to I_active, fr, VE, and the chosen eta 

	pars$eta = eta
	pars$beta_low = getBeta_low(eta,0.65)
	pars$beta_high = getBeta_high(eta,0.65)
	pars$epsilon = epsilon
	pars$phi = getPhi(eta,0.65,1) #zero in the third argument if only the high risk ppl get pasc, 1 if everyone can

	#first, run the system to steady state before introducing the variant

	Times <- seq(0, tEnd, by = DT)		#define timesteps
	y0 <- c(SL=getSL(Eta,0.65),
		IL=getIL(Eta,0.65),
		PL=getPL(Eta,0.65),
		SH=getSH(Eta,0.65),
		IH=getIH(Eta,0.65),
		PH=getPH(Eta,0.65),
		uexp=0)
	varnms = names(y0)
	outprelim <- ode(y0,Times,F,pars)

	y0 <- tail(outprelim,n=1)[2:8]
	y0 <- as.numeric(y0)
	names(y0) <- c(varnms)

	# this deltabeta business is deprecated now; each beta goes up by a factor of 1+s	
	# solve the system under 2 scenarios: all of the selective benefit in beta, and all of it in eta:
	deltabeta_low <- with(pars, s*pars$beta_low )
	deltabeta_high <- with(pars, s*pars$beta_high )

	#calculate the initial rate of increase of p = I1/(I1+I2) for the logistic equation
	sw <- s * as.numeric( with(pars, beta_low*(y0['SL']+(1-eta)*y0['PL'])+beta_high*(y0['SH']+(1-eta)*y0['PH']) ) )

	#update sw in pars 
	pars$sw <- sw

	#create 2 parameter sets: one for immune escape evolution, one for heightened transmission
	pars.transmission <- pars
	pars.transmission$deltabeta_low <- deltabeta_low
	pars.transmission$deltabeta_high <- deltabeta_high


	#pars.immescape <- pars
	#pars.immescape$deltaeta = deltaeta

	# simulate the system with the new parameter values to simulate variant emergence
	out <- as.data.frame(ode(y0, Times, F, pars.transmission )) 		
	#out.immescape    <- as.data.frame(ode(y0, Times, F, pars.immescape    ))

	#Beta = function(Time,spd,beta,deltabeta) beta + deltabeta * ( plogis(spd*(Time-200))  - plogis(spd*(Time-210)))
	#out$beta = sapply(out$time, function(x) Beta(x,pars$spd,pars$beta,pars$deltabeta)  )

	names(out)[9] <- 'Sigma'
	
	#the logistic approximation for variant turnover doesn't work great here. It's off a little at the start, but
	# gives the idea
	return(out[])

}

getu <- function(tau,alphatau,l,alphal,mod,out){

	################
	# PASC latency #
	################

	# define the PDF for the distributed delay between infection and long COVID symptom emergence
	h <- function(u) dgamma(u, shape=alphal, scale=l/alphal)

	#################
	# PASC duration #
	#################

	# define the survival probability for PASC duration:

	# Durations to plug into the pdfs for plotting the Sequelae duration:
	uvals = seq(0,250, length=1000)

	if(mod=='exp'){
		#Exponential Distribution:
		s <- function(u) 1-pexp(u,rate=1/tau)
		dens <- function(u) dexp(u,rate=1/tau)
		densvals <- sapply(uvals,dens)
	}
	
	if(mod=='gam'){	
		#Gamma Distribution:
		alpha.tau = alphatau  #No. of 'layers' in the Exponential Filter
		sigma = tau/alpha.tau
		s <- function(u) 1-pgamma(u, shape=alpha.tau, scale=sigma)
		dens <- function(u) dgamma(u,shape=alpha.tau,scale=sigma)
		densvals = sapply(uvals,dens)
	}

	if(mod=='weibull'){	
		#Weibull Distribution:
		alpha.tau = alphatau
		sigma = tau/gamma(1+1/alpha.tau)
		s <- function(u) 1-pweibull(u, shape=alpha.tau, scale=sigma)	#survival probability
		dens <- function(u) dweibull(u,shape=alpha.tau, scale=sigma)	#PDF
		densvals = sapply(uvals,dens)				#PDF values evaluated at uvals, for displaying the distribution
	}

	if(! mod %in% c('exp','gam','weibull')){print('Error! choose exp, gam, or weibull, or edit the code')}
			
	#####################
	# Integral Equation #
	#####################


	# The PASC prevalence is the solution to the 
	# Fredholm Integral Equation of the 2nd Kind: 
	# u(t) = f(t) - \int_0^\infty K(t,s) u(t-s) ds

	#We solve by iterating Equation 2 from the notes.

	# Calculate f(t), the inhomogeneous part of the FIE:

	Sigma = out$Sigma
	Times = out$time

	mu = 0
	#mort = 1-pexp(Times,rate=mu)
	#f = convolve(Sigma[-1],rev(s(Times[-1])*rev(mort[-1])*rev(diff(Times))),type='open')[1:length(Times)] 
	f = convolve(Sigma[-1],rev(s(Times[-1])*rev(diff(Times))),type='open')[1:length(Times)] 

	get.u = function(u){

		# include a latency distributions
		#ntgrnd = convolve((1-u[-1])*Sigma[-1], rev(h(Times[-1]))*rev(diff(Times)), type='open')[1:length(Times)]
		
		# don't include a latency distribution
		ntgrnd = (1-u[-1])*Sigma[-1]
		#u = f- convolve( lambda[-1]*u[-1], rev(s(Times[-1]))*rev(mort[-1])*rev(diff(Times)) , type='open' )[1:length(Times)]
		
		u = convolve(ntgrnd[-1], rev(s(Times[-1]))*rev(diff(Times)), type='open')[1:length(Times)]

		return(u)
	}

	# Solve the FIE with iteration. Make the initial guess u=f, and iterate away:
	u=rnorm(f,f,.01)
	u.infty = rep(Sigma[1]/(Sigma[1] + 1/tau), length(f) )
	for(i in 1:50) u = get.u(u)
	return(u)

}


# run the simulations that make the plots in the results section of the manuscript


##############################################
## VE_T, eta, or epsilon figure simulations ##
##############################################

# run different simulations to generate the 3-part figure
#
##1. run simulations for combined VE =0.55
VE = 0.55
etavals = seq(0,VE, length=30)
epsvals = sapply(etavals, function(x) 1-(1-VE)/(1-x)  )
parms = cbind(etavals,epsvals)
outlist = list()
f <- function(x) getout(52,as.numeric(x[1]),as.numeric(x[2]),1/50,.02,alpha=1/52)
outlist.ve <- apply(parms,1,f)

#2. run simulations with epsilon = 0.55
etavals = seq(0,0.99,length=30)
f <- function(x) getout(52,as.numeric(x), 0.55, 1/50, 0.02)
outlist.eps <- lapply(etavals,  f)

#3. run simulations with eta = 0.55
epsvals = seq(0,0.99,length=30)
f <- function(x) getout(52,0.55,as.numeric(x),1/50, 0.02)
outlist.eta <- lapply(epsvals,  f)




# write the funtions that generate the plots from the sims
# plot u~time; feed in a list of 'out' with different
# parameter values to see the prediction spread

getuplot <- function(outlist,ylims,mainlab){
        out = outlist[[1]]
        out$time = out$time-950
        uold = out$uexp
        xold = out$time

        grays = gray((length(outlist):0)/length(outlist),1)

        plot(uold~xold,type='l',lwd=3,
                        xlab='Time [weeks]', ylab='PASC Prevalence',
                        main=mainlab,
                        xlim=c(0,200),ylim=ylims,las=1,lty='dotted')

        for(i in 2:length(outlist)){
                outnew = outlist[[i]]
                outnew$time = outnew$time-950
                unew = outnew$u
                unew = unew
                xnew = outnew$time

                polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=grays[i])

                uold=unew
                xold=xnew
        }

        lines(uold~xold,lwd=3)
}

#like getuplot, but superimposes on a previous plot instead
# of drawing a new plotting window

require(scales)

adduplot <- function(outlist,col='gray'){
        out = outlist[[1]]
        out$time = out$time-950
        uold = out$uexp
        xold = out$time

        grays = alpha(col, alpha=rev(c((length(outlist):0)/length(outlist),1)))

        lines(uold~xold,type='l',lwd=1,lty='dotted')

        for(i in 2:length(outlist)){
                outnew = outlist[[i]]
                outnew$time = outnew$time-950
                unew = outnew$u
                unew = unew
                xnew = outnew$time

                polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=grays[i])

                uold=unew
                xold=xnew
        }

        lines(uold~xold,lwd=1)
}

# like getuplot, but just draws the individual trajectories, and not the shaded
# gradients between them

getulines <- function(outlist,ylims=c(0,.1),maintitle='Title goes here!'){

        outlist[[1]]$time = outlist[[1]]$time - 950
        plot(uexp~time,outlist[[1]],xlab='Time [weeks]',
                ylim=ylims, main=maintitle,type='l',
                ylab='PASC Prevalence')

        if(length(outlist)>1){
                for(i in 2:length(outlist)){
                        outlist[[i]]$time = outlist[[i]]$time-950
                        lines(uexp~time,outlist[[i]])

                }
        }

}

getuplot2 <- function(outlist,ylims,mainlab,axno=T,col='gray',axesvar=F){
        out = outlist[[1]]
        out$time = out$time-950
        uold = out$uexp
        xold = out$time

        grays = alpha(col,alpha=rev(c((length(outlist):0)/length(outlist),1)))

        plot(uold~xold,type='l',lwd=1,
                        xlab='', ylab='',
                        main=mainlab,
                        xlim=c(0,200),ylim=ylims,las=1,lty='dotted',axes=axesvar)


        for(i in 2:length(outlist)){
                outnew = outlist[[i]]
                outnew$time = outnew$time-950
                unew = outnew$u
                unew = unew
                xnew = outnew$time

                polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=grays[i])

                uold=unew
                xold=xnew
        }

        lines(uold~xold,lwd=1)
}


# generate the figures:

getfig3_het <- function(save=F){
        if(save==T){jpeg(file='pasc_preds_exp_het_only_get_pasc.jpeg', width = 2880, height = 960,
             pointsize = 36, quality = 100, bg = "white")
        }

        par(mfrow=c(1,3))
        getuplot2(outlist.ve,c(0,0.1),bquote('Outcomes for '~VE==0.55),col='black')
        text(x=150,y=0.025,bquote(epsilon==0.55~', '~eta==0))
        text(x=100,y=0.07,bquote(epsilon==0~', '~eta==0.55))
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
        mtext('Time [weeks]',line=-32,outer=T)
        mtext('PASC Prevalence',las=3,line=-1.5,side=2,outer=T)

        getuplot2(outlist.eps,c(0,0.1),bquote('Outcomes for '~epsilon==0.55),col='black')
        text(x=150,y=0.015,bquote(eta==1~', '~VE==1))
        text(x=100,y=0.07,bquote(eta==0~', '~VE==0.55))
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)

        getuplot2(outlist.eta,c(0,0.1),bquote('Outcomes for '~eta==0.55),col='black')
        text(x=150,y=0.02,bquote(epsilon==1~', '~VE==1))
        text(x=100,y=0.04,bquote(epsilon==0~', '~VE==0.55))
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)

        if(save==T){dev.off()}

}



getleggrad <- function(xleftdiff,xrightdiff,ytopdiff,ybottomdiff){

        # Gradient parameters for legend
        n <- 100  # Number of gradient steps
        eta_values <- seq(0, 0.55, length.out = n)
        colors <- gray(seq(0, 1, length.out = n))  # Gradient from white to black
        # Positioning for the legend
        xleft <- par("usr")[1]+xleftdiff  # Left boundary of the plot
        xright <- xleft + xrightdiff  # Width of the legend
        ytop <- par("usr")[4]+ytopdiff  # Top boundary of the plot
        ybottom <- ytop +ybottomdiff  # Height of the legend

        # Draw the gradient
        for (i in 1:n) {
          rect(xleft = xleft, ybottom = ybottom + (i-1)*(ytop-ybottom)/n,
               xright = xright, ytop = ybottom + i*(ytop-ybottom)/n,
               col = colors[i], border = NA)
        }


}



getfig3a <- function(save=F){
        if(save==T){pdf(file='pasc_preds_exp_het_a.pdf', width = 5, height = 5)
        }

        getuplot2(outlist.ve,c(0,0.1),bquote('Outcomes for '~VE[T]==0.55),col='black')
        text(x=150,y=0.027,bquote(eta==0.55~', '~epsilon==0))
        text(x=120,y=0.07,bquote(eta==0~', '~epsilon==0.55),srt=25)
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
        mtext('Time [weeks]',line=-28,outer=T)
        mtext('PASC Prevalence',las=3,line=-1.5,side=2,outer=T)

        getleggrad(20,20,-.025,-.025)
        # Add text labels for eta values
        text(x = par('usr')[1]+25, y = .085, labels = expression(eta == 0), adj = 0)
        text(x = par('usr')[1]+20, y = .05, labels = expression(eta == 0.55), adj = 0)
        if(save==T){dev.off()}
}


getfig3b <- function(save=F){
        if(save==T){pdf(file='pasc_preds_exp_het_b.pdf', width = 5, height = 5)
        }

        par(mfrow=c(1,1))
        getuplot2(outlist.eps,c(0,0.1),bquote('Outcomes for '~epsilon==0.55),col='black')
        text(x=150,y=0.015,bquote(eta==1~', '~VE[T]==1))
        text(x=110,y=0.065,bquote(eta==0~', '~VE[T]==0.55),srt=35)
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
        mtext('Time [weeks]',line=-24,outer=T)
        getleggrad(20,20,-.025,-.025)
        text(x = par('usr')[1]+23, y = .085, labels = expression(eta == 0), adj = 0)
        text(x = par('usr')[1]+23, y = .05, labels = expression(eta == 1), adj = 0)
        if(save==T){dev.off()}
}

getfig3c <- function(save=F){
        if(save==T){pdf(file='pasc_preds_exp_het_c.pdf', width = 5, height = 5)
        }

        par(mfrow=c(1,1))
        getuplot2(outlist.eta,c(0,0.1),bquote('Outcomes for '~eta==0.55),col='black')
        text(x=150,y=0.021,bquote(epsilon==1~', '~VE[T]==1))
        text(x=100,y=0.035,bquote(epsilon==0~', '~VE[T]==0.55),srt=5)
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
        getleggrad(20,20,-.025,-.025)
        text(x = par('usr')[1]+25, y = .085, labels = expression(epsilon == 0), adj = 0)
        text(x = par('usr')[1]+25, y = .05, labels = expression(epsilon == 1), adj = 0)

        if(save==T){dev.off()}

}


