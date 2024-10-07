rm(list=ls()) #clears the working directory to start from a clean slate each time


# This code is similar to pasc.R but includes a separate vaccination state
# to examine how differences between vaccine-induced and infection-derived
# immunity affect results 
 
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
	deltabeta=0.0,          # change in transmisison rate, [1/population density/week]
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
	mu      = 0/52/75,      #per-capita birth/death rate, [1/week], don't have the formula from Overleaf with this part yet
	theta	= .5,		# efficacy of vaccine relative to avg immune efficacy, [dimensionless]
	kappa = 2)		# rate of vaccine immunity loss relative to avg rate of infection-derived immunity loss, [dimensionless]


	# calculate beta to satisfy calibration target I*, assuming gamma, delta, 
	# alpha are known, and we want to supply various values of eta.
	# However, there is a subpopulation of size p_risk with elevated susceptibility rho 
	# choose rho to accomodate I_active, delta, alpha, gamma, eta, p_risk, and beta
	# should be the case that when rho=1, we get the same beta value as before

	# this is how we chose beta in the homogeneous case:
	getBeta = function(eta){with(pars,{ 
		a = (eta-1)*I_active^2 - eta*I_active+I_active
		b = (eta-1)*(I_active-1)*alpha - (delta+gamma)*I_active + delta
		c = -gamma*(alpha+delta)

		return( max( (-b + sqrt(b^2-4*a*c))/2/a,
		(-b - sqrt(b^2-4*a*c))/2/a))
		
	})}
	

	getS <- function(eta){
		with(pars, {
			S = 1 + ((alpha-gamma)*I_active - alpha)/(getBeta(eta)*I_active*(1-eta)+alpha+delta)  - I_active
			return(S)
		})

	}

	getP <- function(eta){
		with(pars, {
			P =  -((alpha-gamma)*I_active - alpha)/(getBeta(eta)*I_active*(1-eta)+alpha+delta) 
			return(P)
		})
	}

	# equilibrium proportion P(vaccinated | immune)
	# kappa is is delta_V/delta_P, and theta is eta_V/eta_P
	getx <- function(eta,kappa,theta){
		with(pars,{
			S = getS(eta)
			I=I_active
			P=getP(eta)
			beta=getBeta(eta)
			x = alpha*S/(-I*P*beta*eta*theta + I*P*beta*eta + P*delta*kappa + I*gamma - P*delta + S*alpha)
			return(x)
		})
	}

	# we have our overall states S, P calibrated, but need to get the relative number of high- and low-risk individuals right

	#calibrate phi for the chosen value of eta, along with the background values of VE and fr (Q)
	# Q is the ratio of PASC incidence to total incidence, and is assumed known
	# some algebra:
	# Q = (phi*(beta_low*SL + beta_high*SH)*I + phi*(1-epsilon)*(1-eta)*I*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH)*I + (1-eta)*I*(beta_low*PL+beta_high*PH))
	# Q = (phi*(beta_low*SL + beta_high*SH) + phi*(1-epsilon)*(1-eta)*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH) + (1-eta)*(beta_low*PL+beta_high*PH))
	# Q = phi* ((beta_low*SL + beta_high*SH) + (1-epsilon)*(1-eta)*(beta_low*PL+beta_high*PH)) / ((beta_low*SL + beta_high*SH) + (1-eta)*(beta_low*PL+beta_high*PH))
	# Q = phi* (beta_s(x)*S + (1-epsilon)*(1-eta)*beta_p(y)*P) / (beta_s(x)*S+(1-eta)*beta_p(y)*P) 

	# phi = Q * (beta_s(x)*S+(1-eta)*beta_p(z)*P) / (beta_s(x)*S + (1-epsilon)*(1-eta)*beta_p(z)*P)


	getPhi <- function(eta){
                Phi =  Q * (getS(eta)+(1-eta)*getP(eta))/(getS(eta) + (1-VE)*getP(eta))
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

			beta_fun = function(Time) beta + deltabeta * ( plogis(sw*(Time-Tinv)) )
			eta_fun  = etaV*x + etaP*(1-x)
			delta_fun = deltaV*x + deltaP*(1-x)
	
			dS	<- -beta_fun(Time)*S*I + delta_fun*P + mu - mu*S - alpha*S
			dI	<-  beta_fun(Time)*S*I + (1-eta_fun)*beta_fun(Time)*P*I - gamma*I  - mu*I
			dP	<-  gamma*I - (1-eta_fun)*beta_fun(Time)*P*I - delta_fun*P  - mu*P  + alpha*S
			dx	<- (((eta_fun*beta_fun(Time)*(theta - 1)*I - delta*(kappa - 1))*P - gamma*I - alpha*S)*x + alpha*S)/P
	
			Sigma	<- phi*(beta_fun(Time)*S +(1-eta_fun)*(1-epsilon)*(beta_fun(Time)*P))*I

			duexp	<- Sigma*(1-uexp)-uexp/tau

		return(list(c(dS,dI,dP,dx,duexp),c(Sigma,beta_fun(Time),eta_fun,delta_fun)))
		})
	}


	
	#########################
	# set up and run solver #
	#########################

	#use the parameters calibrated to I_active, fr, VE, and the chosen eta 

	pars$eta = eta
	pars$etaV = with(pars, theta*eta)
	pars$etaP = with(pars, eta*(theta*getx(eta,kappa,theta)-1)/(getx(eta,kappa,theta)-1))
	pars$deltaV = with(pars, kappa*delta)
	pars$deltaP = with(pars, delta*(kappa*getx(eta,kappa,theta)-1)/(getx(eta,kappa,theta)-1))

	pars$beta = getBeta(eta)
	pars$epsilon = epsilon
	pars$phi = getPhi(eta) #zero in the third argument if only the high risk ppl get pasc, 1 if everyone can

	#first, run the system to steady state before introducing the variant

	Times <- seq(0, tEnd, by = DT)		#define timesteps
	y0 <- c(S=getS(pars$eta),
		I=I_active,
		P=getP(pars$eta),
		x=getx(Eta,pars$kappa,pars$theta),
		uexp=0)
	varnms = names(y0)
	outprelim <- ode(y0,Times,F,pars)

	y0 <- tail(outprelim,n=1)[2:6]
	y0 <- as.numeric(y0)
	names(y0) <- c(varnms)

	# this deltabeta business is deprecated now; each beta goes up by a factor of 1+s	
	# solve the system under 2 scenarios: all of the selective benefit in beta, and all of it in eta:
	deltabeta <- with(pars, s*pars$beta )

	#calculate the initial rate of increase of p = I1/(I1+I2) for the logistic equation
	sw <- s * as.numeric( with(pars, beta*(y0['S']+(1-eta)*y0['P'])) ) 

	#update sw in pars 
	pars$sw <- sw

	#create 2 parameter sets: one for immune escape evolution, one for heightened transmission
	pars.transmission <- pars
	pars.transmission$deltabeta <- deltabeta


	#pars.immescape <- pars
	#pars.immescape$deltaeta = deltaeta

	# simulate the system with the new parameter values to simulate variant emergence
	out <- as.data.frame(ode(y0, Times, F, pars.transmission )) 		
	#out.immescape    <- as.data.frame(ode(y0, Times, F, pars.immescape    ))

	#Beta = function(Time,spd,beta,deltabeta) beta + deltabeta * ( plogis(spd*(Time-200))  - plogis(spd*(Time-210)))
	#out$beta = sapply(out$time, function(x) Beta(x,pars$spd,pars$beta,pars$deltabeta)  )

	names(out)[7] <- 'Sigma'
	names(out)[8] <- 'beta'
	names(out)[9] <- 'eta'	
	names(out)[10] <- 'delta'
	return(out)

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

                polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=NA) #grays[i])

                uold=unew
                xold=xnew
        }

        lines(uold~xold,lwd=1)
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

# generate the figures:

getfig3_vax <- function(save=F){
        if(save==T){jpeg(file='pasc_preds_exp_vax.jpeg', width = 2880, height = 960,
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


getfig3a <- function(save=F){
        if(save==T){jpeg(file='pasc_preds_exp_vax_a.jpeg', width = 5, height = 5, units='in', res=700)
        }

        getuplot2(outlist.ve,c(0,0.1),bquote('Outcomes for '~VE[T]==0.55),col='black')
        text(x=150,y=0.024,bquote(eta==0.55~', '~epsilon==0))
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
        if(save==T){jpeg(file='pasc_preds_exp_vax_b.jpeg', width = 5, height = 5, units='in', res=700)
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
        if(save==T){jpeg(file='pasc_preds_exp_vax_c.jpeg', width = 5, height = 5, units='in', res=700)
        }

        par(mfrow=c(1,1))
        getuplot2(outlist.eta,c(0,0.1),bquote('Outcomes for '~eta==0.55),col='black')
        text(x=150,y=0.02,bquote(epsilon==1~', '~VE[T]==1))
        text(x=100,y=0.03,bquote(epsilon==0~', '~VE[T]==0.55),srt=5)
        axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)
        axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
        getleggrad(20,20,-.025,-.025)
        text(x = par('usr')[1]+25, y = .085, labels = expression(epsilon == 0), adj = 0)
        text(x = par('usr')[1]+25, y = .05, labels = expression(epsilon == 1), adj = 0)

        if(save==T){dev.off()}

}


plot(eta~time,outlist.ve[[1]],type='l',ylim=c(0,1),
	xlab='Time [weeks]', ylab=expression(eta))
for(i in 2:30) lines(eta~time,outlist.ve[[i]])


plot(eta~time,outlist.eta[[1]],type='l',ylim=c(0,1),
	xlab='Time [weeks]', ylab=expression(eta))
for(i in 2:30) lines(eta~time,outlist.eta[[i]])

plot(eta~time,outlist.eps[[1]],type='l',ylim=c(0,1),
	xlab='Time [weeks]', ylab=expression(eta))
for(i in 2:30) lines(eta~time,outlist.eps[[i]])



plot(delta~time,outlist.ve[[1]],type='l',ylim=c(1/104,2/52),
	xlab='Time [weeks]', ylab=expression(delta))
for(i in 2:30) lines(delta~time,outlist.ve[[i]])


plot(delta~time,outlist.eta[[1]],type='l',ylim=c(1/104,2/52),
	xlab='Time [weeks]', ylab=expression(delta))
for(i in 2:30) lines(delta~time,outlist.eta[[i]])

plot(delta~time,outlist.eps[[1]],type='l',ylim=c(1/104,2/52),
	xlab='Time [weeks]', ylab=expression(delta))
for(i in 2:30) lines(delta~time,outlist.eps[[i]])

