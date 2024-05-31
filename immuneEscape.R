rm(list=ls()) #clears the working directory to start from a clean slate each time


# Information we know which places constraints on parameters in model:
#VE <- .7 #vaccine efficacy against PASC, = 1-(1-epsilon)*(1-eta), eta=efficacy against infection, epsilon=efficacy against PASC conditional on breakthrough infection
#fr <- 1/50 #fraction of all infections that have resulted in long COVID
#I_active <- .01 #fraction of population with active COVID infection



getOut <- function(Tau,a.duration,L,a.latent,VE,fr,I_active){

	# This function outputs the SIP ODE solution, and PASC prevalence
	# ODE is either SIP or immune escape model
	# PASC prevalence obtained as solution of a Fredholm Integral Eqn
	# assumes that PASC and infection dynamics are independent
	# nice for handling general duration of sequelae

	require(deSolve)	#load in the ODE solver package
	require(boot)		#load in inv.logit

	#run the simulation for 100 weeks to eliminate transients - what time should the variant emerge?
	Tinv <- 1000

	#selection coefficient
	sel <- .1

	##############
	# Parameters #		#description of parameter, [units of parameter]
	##############

	# Model parameters
	pars <- list(
	beta    = 0,            # depends on eta;  transmission rate, [1/population density/week]
	deltabeta=0.0,          # change in transmisison rate, [1/population density/week]
	sw = .4,                # selection coefficient
	gamma   = 1,            # recovery rate, [1/week]
	delta   = 1/30,        # rate of adaptive immune loss, [1/week]
	#eta     = .95,         # immune efficacy against infection, [dimensionless]
	deltaeta = 0.0,	        # change in immune efficacy, [dimensionless]
	alpha   = 1/50,         # per-capita vaccination rate, [1/week]
	tau     = Tau,           #mean duration of long Covid, [week]
	phi     = 0,         #probability of developing long Covid from S0, [dimensionless]
	epsilon  = 0,           #efficacy against PASC conditional on breakthrough infection (in immunized people), [dimensionless]
	L	= L,		#mean delay between infection and onset of sequelae, [weeks]
	mu      = 0/52/75)      #per-capita birth/death rate, [1/week], don't have the formula from Overleaf with this part yet

	# fix I_active, and calculate beta, assuming gamma, delta, alpha are known, and we want to supply various values of eta
	# choose beta to accomodate I_active, delta, alpha, gamma, eta:

	beta = with(pars, delta*gamma/(delta - I_active*(gamma+delta) ) )
	pars$beta = beta

	#calculate S and R

	S = with(pars, gamma/beta)

	R = with(pars, (1-gamma/beta)*gamma/(gamma+delta) )

	# choose eta based on the value of the selection coefficient

	eta = with(pars, 1-sw * S/R)

	#calibrate phi for the chosen value of eta, along with the background values of VE and fr

	Phi = fr * (S+(1-eta)*R)/(S + (1-VE)*R)



	################
	# PASC latency #
	################

	# define the PDF for the distributed delay between infection and long COVID symptom emergence
	alpha1 = a.latent
	h <- function(u) dgamma(u, shape=alpha1, scale=pars$L/alpha1)

	#####################
	# sequelae duration #
	#####################

	# define the survival probability for sequelae duration:

	# Durations to plug into the pdfs for plotting the Sequelae duration:
	tvals = seq(0,1000, length=3000)

	#Comment out the ones you don't want to use:

	#Exponential Distribution:
	s <- function(u) 1-pexp(u,rate=1/pars$tau)
	dens <- function(u) dexp(u,rate=1/pars$tau)
	densvals <- sapply(tvals,dens)

	dens.exp<-dens
	densvals.exp<-densvals

	#Gamma Distribution (uncomment the following lines to use):
	alpha = a.duration  #No. of 'layers' in the Exponential Filter
	sigma = pars$tau/alpha
	s <- function(u) 1-pgamma(u, shape=alpha, scale=sigma)
	dens <- function(u) dgamma(u,shape=alpha,scale=sigma)
	densvals = sapply(tvals,dens)

	#Weibull Distribution:
	#a = a.duration
	#sigma = pars$tau/gamma(1+1/a)
	#s <- function(u) 1-pweibull(u, shape=a, scale=sigma)	#survival probability
	#dens <- function(u) dweibull(u,shape=a,scale=sigma)	#PDF
	#densvals = sapply(tvals,dens)				#PDF values evaluated at uvals, for displaying the distribution

	latvals = sapply(tvals,h)

	#############
	# timesteps #
	#############

	tEnd <-2000		
	DT <- 1/2			#temporal resolution for solver; can be rather coarse

	##########################
	# differential equations #
	##########################


	# Immune Escape model #

	G <- function(Time, state, Pars) {
		with(as.list(c(state, Pars)), {


			dS	<- -beta*S*I1 - beta*S*I2 + delta*(R1+R2)
			
			dI1	<-  beta*S*I1 - gamma*I1
			dR1	<-  gamma*I1  - (1-eta)*beta*I2*R1 - delta*R1

			dI2	<- beta*S*I2 + (1-eta)*beta*R1*I2 - gamma*I2
			dR2	<-  gamma*I2 - delta*R2

			du_exp <- phi*beta*(S*(I1+I2)+(1-eta)*(1-epsilon)*R1*I2) * (1-u_exp) - u_exp/tau

		return(list(c(dS,dI1,dR1,dI2,dR2,du_exp)))
		})
	}




	#########################
	# set up and run solver #
	#########################

	#want to set up the system so it is properly calibrated to I_active, fr, VE, and the chosen eta 

	#update pars
	pars$eta = eta
	pars$beta = beta
	pars$epsilon = (VE-eta)/(1-eta)
	pars$phi = Phi

	#first, run the system to steady state before introducing the variant

	Times <- seq(0, tEnd, by = DT)		#define timesteps
	outprelim <- ode(c(S=S,I1=I_active,R1=R,I2=0,R2=0,u_exp=0),Times,G,pars)

	y0 <- tail(outprelim,n=1)[-1]
	y0 <- as.numeric(y0)
	names(y0) <- c('S','I1','R1','I2','R2','u_exp')

	y0['I2'] = 1e-8
	y0['I1'] = y0['I1']-1e-8	

	outnew <- as.data.frame(ode(y0, Times + DT + tail(Times,n=1), G, pars )) 		#solve the system! place solution in a dataframe called "out" (short for output)

	out <- rbind(outprelim,outnew)

	out$I <- out$I1 + out$I2
	out$Sigma <- with(as.list(c(pars,out)),phi*beta*I*(S+(1-eta)*(1-epsilon)*R1*I2/I) )

	getu <- function(out){

	# influx into the long COVID compartment, obtained from the ODE solution. Gets
	# fed into the Fredholm Integral Equation for the long COVID prevalence
	Sigma = out$Sigma
	Times = out$time
		
	#####################
	# Integral Equation #
	#####################


	# The sequelae prevalence is the solution to the 
	# Fredholm Integral Equation of the 2nd Kind: 
	# u(t) = f(t) - \int_0^\infty K(t,s) u(t-s) ds

	#We solve by iterating Equation 2 from the notes.

	# Calculate f(t), the inhomogeneous part of the FIE:

	mort = 1-pexp(Times,rate=pars$mu)
	f = convolve(Sigma[-1],rev(s(Times[-1])*rev(mort[-1])*rev(diff(Times))),type='open')[1:length(Times)] 

	#Define th:

	get.u = function(u){

		# the first ntgrnd uses a latent distribution, h; the second one does not
		ntgrnd = convolve((1-u[-1])*Sigma[-1], rev(h(Times[-1]))*rev(diff(Times)), type='open')[1:length(Times)]
		
		#ntgrnd = (1-u[-1])*Sigma[-1]
		#u = f- convolve( lambda[-1]*u[-1], rev(s(Times[-1]))*rev(mort[-1])*rev(diff(Times)) , type='open' )[1:length(Times)]
		u = convolve(ntgrnd[-1], rev(s(Times[-1]))*rev(diff(Times)), type='open')[1:length(Times)]

		return(u)
	}

	# Solve the FIE with iteration. Make the initial guess L=f, and iterate away:
	u=rnorm(f,f,.01)
	u.infty = rep(Sigma[1]/(Sigma[1] + 1/pars$tau), length(f) )
	for(i in 1:100) u = get.u(u)
	return(u)

	}

	u <- getu(out)

	out$u <- u

	#out <- tail(out,n=1500)

	return(list(out=out,densvals=densvals,densvals.exp=densvals.exp,latvals=latvals,tvals=tvals))

}


crud = getOut(50,6,10,10,.55,1/200,.02)
densvals=crud$densvals
densvals.exp=crud$densvals.exp
latvals=crud$latvals
tvals=crud$tvals
out=crud$out
out = out[out$time > 2000 & out$time < 2200,]
out$time = out$time - 2000

#creates Fig. 7 in the manuscript

getPlot <- function(save=F){

	if(save==T){jpeg(file='immescape.jpeg', width = 960, height = 960,
	     pointsize = 24, quality = 100, bg = "white")}

	par(mfrow=c(2,2),mar=c(5,5,4,2))
	plot(u ~ time,out, 
		type='l',xlab='Time [weeks]',ylab='',main='PASC',
		ylim=c(0,.01),las=1)
	title(ylab='Prevalence', line=3.5,cex.lab=1.2)
	lines( u_exp ~ time,out,lty='dotted')

	plot((I1+I2)~time,out,type='l',
		ylab='',xlab='Time [weeks]',main='Active Infections',
		ylim=c(0,.1),las=1)
	abline(h=0,lty='dotted')
	title(ylab='Prevalence', line=3.5,cex.lab=1.2)

	plot(tvals,densvals,type='l',xlim=c(0,100),
		main='PASC duration',xlab='Duration [weeks]',ylab='',las=1)
	lines(tvals,densvals.exp,lty='dotted')
	title(ylab='Prob. density', line=3.5,cex.lab=1.2)

	plot(tvals,latvals,type='l',xlim=c(0,30),
		main='PASC latency', xlab='Latency [weeks]',
		ylab='',las=1)
	title(ylab='Prob. density', line=3.5,cex.lab=1.2)


	if(save==T){dev.off()}
}

