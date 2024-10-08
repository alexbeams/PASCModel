rm(list=ls()) #clears the working directory to start from a clean slate each time


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

	# This code solves the SIP ODE model
 	
	require(deSolve)	#load in the ODE solver package
	require(boot)		#load in inv.logit

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
	mu      = 0/52/75)      #per-capita birth/death rate, [1/week], don't have the formula from Overleaf with this part yet

	# fix I_active, and calculate beta, assuming gamma, delta, 
	# alpha are known, and we want to supply various values of eta
	# choose beta to accomodate I_active, delta, alpha, gamma, eta:

	getBeta = function(eta){with(pars,{ 
		a = (eta-1)*I_active^2 - eta*I_active+I_active
		b = (eta-1)*(I_active-1)*alpha - (delta+gamma)*I_active + delta
		c = -gamma*(alpha+delta)

		return( max( (-b + sqrt(b^2-4*a*c))/2/a,
		(-b - sqrt(b^2-4*a*c))/2/a))
		
	})}

	#calculate S and P in terms of alpha, gamma, delta, eta, beta(eta), and I_active

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

	#calibrate phi for the chosen value of eta, along with the background values of VE and fr (Q)

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

			Beta = function(Time) beta + deltabeta * ( inv.logit(sw*(Time-Tinv)) )
			Eta  = function(Time) eta - deltaeta * ( inv.logit(sw*(Time-Tinv))) 
	
			dS	<- -Beta(Time)*S*I + delta*P + mu - mu*S - alpha*S
			dI	<-  Beta(Time)*S*I + (1-Eta(Time))*Beta(Time)*P*I - gamma*I  - mu*I
			dP	<-  gamma*I  - (1-Eta(Time))*Beta(Time)*I*P - delta*P  - mu*P  + alpha*S

			Sigma	<- Beta(Time)*phi*(S+(1-Eta(Time))*(1-epsilon)*P)*I

			duexp	<- Sigma*(1-uexp)-uexp/tau

		return(list(c(dS,dI,dP,duexp),c(Sigma)))
		})
	}


	
	#########################
	# set up and run solver #
	#########################

	#use the parameters calibrated to I_active, fr, VE, and the chosen eta 

	pars$eta = eta
	pars$beta = getBeta(eta)
	pars$epsilon = epsilon
	pars$phi = getPhi(eta)

	#first, run the system to steady state before introducing the variant

	Times <- seq(0, tEnd, by = DT)		#define timesteps
	outprelim <- ode(c(S=getS(Eta),I=I_active,P=getP(Eta),uexp=0),Times,F,pars)

	y0 <- tail(outprelim,n=1)[2:5]
	y0 <- as.numeric(y0)
	names(y0) <- c('S','I','P','uexp')

	# solve the system under 2 scenarios: all of the selective benefit in beta, and all of it in eta:
	deltabeta <- with(pars, s*pars$beta )
	#deltaeta  <- with(pars, s*(1-eta + y0['S']/y0['P']))

	#calculate the initial rate of increase of p = I1/(I1+I2) for the logistic equation
	sw <- s * as.numeric( with(pars, beta*(y0["S"] + (1-eta)*y0["P"]) ) )

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

	#Beta = function(Time,spd,beta,deltabeta) beta + deltabeta * ( inv.logit(spd*(Time-200))  - inv.logit(spd*(Time-210)))
	#out$beta = sapply(out$time, function(x) Beta(x,pars$spd,pars$beta,pars$deltabeta)  )

	names(out)[6] <- 'Sigma'
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

#1. run simulations for combined VE =0.55
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

#########################################
## different gammas figure simulations ##
#########################################

# run the VE model for a couple of different gammas to generate the next results figure

f <- function(x) getout(52,as.numeric(x[1]),as.numeric(x[2]),1/50,.02,gamma=.5,alpha=1/52)
outlist.ve.gam.5 <- apply(parms,1,f)
f <- function(x) getout(52,as.numeric(x[1]),as.numeric(x[2]),1/50,.02,gamma=2,alpha=1/52)
outlist.ve.gam2 <- apply(parms,1,f)

####################################
## changing I* figure simulations ##
####################################

# run the model for different values of I_active 
Ivals <- seq(0.001, 0.1, length=30)
VE = 0.55
etavals = seq(0,VE, length=10)
epsvals = sapply(etavals, function(x) 1-(1-VE)/(1-x)  )
parms = cbind(etavals,epsvals)

ivalparms = expand.grid(Ivals, c(.5,1,2))
ivallist <- list()
for(i in 1:dim(ivalparms)[1]) ivallist[[i]] <- t(apply(parms, 1, function(x) unlist(tail(getout(52,as.numeric(x[1]), as.numeric(x[2]), 1/50,ivalparms[i,1],alpha=1/52, gamma=ivalparms[i,2]),n=1))  ))

colnames(ivalparms) <- c('Iactive','gamma')

crud=ivallist

for(i in 1:length(ivallist)) crud[[i]] <- cbind(crud[[i]], ivalparms[i,])

cruddat = do.call(rbind, crud)

dat = cruddat
dat$eta = etavals 

######################################
## PASC decrease figure simulations ##
######################################

out.reduct = getout(52,0.65,1,1/100, .09, gamma=2,delta=1/52,alpha=1/52)
out.reduct$time = out.reduct$time - 900
out.reduct = out.reduct[out.reduct$time > 0 & out.reduct$time < 200,]

######################################
## immune escape figure simulations ##
######################################

# generated from immescape.R

#####################################
## last results figure simulations ##
#####################################


# run the VE model for different combinations of gammas, q's, and I_active's to generate the last figure 
# run on a stencil of parameter values, p/m 50% of the central parm value used
# this will help visualize how the 3 pieces of information interact to affect PASC 

# run simulations for combined VE =0.55
VE = 0.55
etavals = seq(0,VE, length=10)
epsvals = sapply(etavals, function(x) 1-(1-VE)/(1-x)  )
#parms = cbind(etavals,epsvals)

# make the stencil of parameter values
I_active.0 = c(.01, .02, .03)
g.0 = c(.5, 1, 1.5)
q.0 = I_active.0 #conveniently chose 1/50 for both 
stncl <- expand.grid(I_active.0,g.0,q.0)
colnames(stncl) <- c('I_active','g','q')
stncl <- stncl[,c('q','I_active','g')]
# make some indices for where these will go in the plotting window
stncind <- stncl
stncind[,1] <- stncind[,1] - .02
stncind[,2] <- stncind[,2] - .02
stncind[,3] <- stncind[,3] - 1
crud=stncind

stncind[crud < 0] = 1
stncind[crud == 0] = 2
stncind[crud > 0] = 3
rm(crud)

g <- function(x) apply(cbind(etavals,epsvals),1,function(y){ 
	x=as.numeric(x)
	y=as.numeric(y)	
	return(getout(52,y[1],y[2],x[1],x[2],x[3]  ))
	}
)

crud = list()
for(i in 1:27) crud[[i]] <- g(stncl[i,])

outlist.stncl <-crud 

###################################################################
## write some functions to generate figures from the simulations ##
###################################################################


getfiga <- function(x){

	if(x==1) for(i in 1:9) getuplot(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(gamma==~.(stncl[i,3])~', '~I==.(stncl[i,2])~', '~q==.(stncl[i,1])) )
	if(x==2) for(i in 10:18) getuplot(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(gamma==~.(stncl[i,3])~', '~I==.(stncl[i,2])~', '~q==.(stncl[i,1])) )
	if(x==3) for(i in 19:27) getuplot(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(gamma==~.(stncl[i,3])~', '~I==.(stncl[i,2])~', '~q==.(stncl[i,1])) )

}

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

		polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=NA ) #grays[i])

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

		polygon(c(xold, rev(xnew)), c(uold, rev(unew)), col=grays[i],border=NA )#grays[i])

		uold=unew
		xold=xnew
	}

	#lines(uold~xold,lwd=1)
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

getfig3 <- function(save=F){
	if(save==T){jpeg(file='pasc_preds_exp.jpeg', width = 14, height = 4,
	     units='in', res = 800)
	}

	par(mfrow=c(1,3))
	getuplot2(outlist.ve,c(0,0.1),bquote('Outcomes for '~VE[T]==0.55),col='black')
	text(x=150,y=0.028,bquote(eta==0.55~', '~epsilon==0))
	text(x=120,y=0.07,bquote(eta==0~', '~epsilon==0.55),srt=25)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)		
	axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
	mtext('Time [weeks]',line=-28,outer=T)
	mtext('PASC Prevalence',las=3,line=-1.5,side=2,outer=T)

	getleggrad(20,20,-.025,-.025)
	# Add text labels for eta values
	text(x = par('usr')[1]+25, y = .085, labels = expression(eta == 0), adj = 0)
	text(x = par('usr')[1]+20, y = .05, labels = expression(eta == 0.55), adj = 0)
	
	getuplot2(outlist.eps,c(0,0.1),bquote('Outcomes for '~epsilon==0.55),col='black')
	text(x=150,y=0.015,bquote(eta==1~', '~VE[T]==1))
	text(x=110,y=0.065,bquote(eta==0~', '~VE[T]==0.55),srt=35)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)		
	axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
	getleggrad(20,20,-.025,-.025)
	text(x = par('usr')[1]+23, y = .085, labels = expression(eta == 0), adj = 0)
	text(x = par('usr')[1]+23, y = .05, labels = expression(eta == 1), adj = 0)
	
	getuplot2(outlist.eta,c(0,0.1),bquote('Outcomes for '~eta==0.55),col='black')
	text(x=150,y=0.023,bquote(epsilon==1~', '~VE[T]==1))
	text(x=100,y=0.035,bquote(epsilon==0~', '~VE[T]==0.55),srt=5)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)		
	axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
	getleggrad(20,20,-.025,-.025)
	text(x = par('usr')[1]+25, y = .085, labels = expression(epsilon == 0), adj = 0)
	text(x = par('usr')[1]+25, y = .05, labels = expression(epsilon == 1), adj = 0)

	if(save==T){dev.off()}

}

getfig3a <- function(save=F){
	if(save==T){jpeg(file='pasc_preds_exp_a.jpeg', width = 5, height = 5, units='in', res=700)
	}

	getuplot2(outlist.ve,c(0,0.1),bquote('Outcomes for '~VE[T]==0.55),col='black')
	text(x=150,y=0.028,bquote(eta==0.55~', '~epsilon==0))
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
	if(save==T){jpeg(file='pasc_preds_exp_b.jpeg', width = 5, height = 5, units='in', res=700)
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
	if(save==T){jpeg(file='pasc_preds_exp_c.jpeg', width = 5, height = 5, units='in', res=700)
	}

	par(mfrow=c(1,1))
	getuplot2(outlist.eta,c(0,0.1),bquote('Outcomes for '~eta==0.55),col='black')
	text(x=150,y=0.023,bquote(epsilon==1~', '~VE[T]==1))
	text(x=100,y=0.035,bquote(epsilon==0~', '~VE[T]==0.55),srt=5)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)		
	axis(2,at=c(.0,.025,.05,.075,.1),labels=c(.0,.025,.05,.075,.1),las=2,pos=0)
	getleggrad(20,20,-.025,-.025)
	text(x = par('usr')[1]+25, y = .085, labels = expression(epsilon == 0), adj = 0)
	text(x = par('usr')[1]+25, y = .05, labels = expression(epsilon == 1), adj = 0)

	if(save==T){dev.off()}

}

getleggradcol <- function(xleftdiff,xrightdiff,ytopdiff,ybottomdiff,col){

	n <- 100 #no. of gradient steps
	colors = alpha(col,alpha=c(n:0)/n)

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


getfig4 <- function(save=F){
	if(save==T){jpeg(file='pasc_gamma_colors.jpeg', width = 4, height = 4, units='in', res=700)
	}
	
	par(mfrow=c(1,1))
	col1 = hcl.colors(n=3,palette='zissou 1')[1]
	col2 = hcl.colors(n=3,palette='zissou 1')[3]

	getuplot2(outlist.ve.gam2, c(0,.2),'' ,col=col1,axesvar=F)
	adduplot(outlist.ve.gam.5,col=col2)
	#text(x=150,y=0.005,bquote(1/gamma==2))
	#text(x=100,y=0.15,bquote(1/gamma==0.5))
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	mtext('Time [weeks]',line=-18,outer=T)
	mtext('PASC Prevalence',las=3,line=-1.5,side=2,outer=T)
	#legend(x=10,y=.20,legend=c(bquote(gamma^{-1}~'=0.5 weeks'),bquote(gamma^{-1}~'=2 weeks')),
	#	col=c(col1,col2),lty=1,lwd=4)
	getleggradcol(20,10,-.05,-.05,col1)
	getleggradcol(30,10,-.05,-.05,col2)
	text(x = par('usr')[1]+18, y = .1, labels = expression(eta == 0.55), adj = 0,cex=.5)
	text(x = par('usr')[1]+22, y = .17, labels = expression(eta == 0), adj = 0,cex=.5)
	text(x=150,y=0.1,labels = expression(gamma^{-1}==0.5~' weeks'),cex=0.5,col='black')
	text(x=160,y=0.026,labels = expression(gamma^{-1}==2~' weeks'),cex=0.5,col='white')


	if(save==T){dev.off()}
}


getfig5color <- function(save=F){
        if(save==T){jpeg(file='pasc_Iactive_colors.jpeg', width = 6, height = 6, units='in', res=700)
                }
        Iactive = Ivals
        #manually make a color plot showing the scenario endpoints as Iactive changes
        gamval=2
        uold = dat[dat$gamma==gamval & dat$eta==etavals[1],'uexp']
        cols=hcl.colors(n=3,palette='zissou 1')
        grays = alpha(cols[1], alpha=rev(c((10:1)/10,1)))

        plot(Iactive,uold,ylim=c(0,0.3),type='l',lty='dotted', ylab='PASC Prevalence after Variant Replacement',
                xlab='Prevalence of Active Infections before Variant Replacement',las=1)
        for(i in 2:10){
                unew = dat[dat$gamma==gamval & dat$eta==etavals[i],'uexp']
                polygon(c(Iactive, rev(Iactive)), c(uold, rev(unew)), col=grays[i],border=grays[i])
                uold=unew
        }
        #manually make a color plot showing the scenario endpoints as Iactive changes
        gamval=1
        uold = dat[dat$gamma==gamval & dat$eta==etavals[1],'uexp']
        cols=hcl.colors(n=3,palette='zissou 1')
        grays = alpha('black', alpha=rev(c((10:1)/10,1)))

        lines(Iactive,uold,ylim=c(0,0.3),lty='dotted')
        for(i in 2:10){
                unew = dat[dat$gamma==gamval & dat$eta==etavals[i],'uexp']
                polygon(c(Iactive, rev(Iactive)), c(uold, rev(unew)), col=grays[i],border=grays[i])
                uold=unew
        }
        #manually make a color plot showing the scenario endpoints as Iactive changes
        gamval=0.5
        uold = dat[dat$gamma==gamval & dat$eta==etavals[1],'uexp']
        cols=hcl.colors(n=3,palette='zissou 1')
        grays = alpha(cols[3], alpha=rev(c((10:1)/10,1)))

        lines(Iactive,uold,ylim=c(0,0.3),lty='dotted')
        for(i in 2:10){
                unew = dat[dat$gamma==gamval & dat$eta==etavals[i],'uexp']
                polygon(c(Iactive, rev(Iactive)), c(uold, rev(unew)), col=grays[i],border=grays[i])
                uold=unew
        }
#       legend('bottomright',legend=c(bquote(gamma^{-1}~'=0.5 weeks'),
#                               bquote(gamma^{-1}~'=1 week'),
#                               bquote(gamma^{-1}~'=2 weeks')),
#               col=c(cols[1],'black',cols[3]),lty=1,lwd=3)
        n <- 30 #no. of gradient steps
        colors = alpha(hcl.colors(3,'zissou 1')[1],alpha=c(n:0)/n)

        # Positioning for the legend
        xleft <- 0.005  # Left boundary of the plot
        xright <- 0.007  # Width of the legend
        ytop <- .27 # Top boundary of the plot
        ybottom <- .2  # Height of the legend

        # Draw the gradient
        for (i in 1:n) {
          rect(xleft = xleft, ybottom = ybottom + (i-1)*(ytop-ybottom)/n,
               xright = xright, ytop = ybottom + i*(ytop-ybottom)/n,
               col = colors[i], border = NA)
        }
        colors = alpha('black',alpha=c(n:0)/n)

        # Positioning for the legend
        xleft <- 0.007  # Left boundary of the plot
        xright <- 0.009  # Width of the legend
        ytop <- .27 # Top boundary of the plot
        ybottom <- .2  # Height of the legend

        # Draw the gradient
        for (i in 1:n) {
          rect(xleft = xleft, ybottom = ybottom + (i-1)*(ytop-ybottom)/n,
               xright = xright, ytop = ybottom + i*(ytop-ybottom)/n,
               col = colors[i], border = NA)
        }
        colors = alpha(hcl.colors(3,'zissou 1')[3],alpha=c(n:0)/n)

        # Positioning for the legend
        xleft <- 0.009  # Left boundary of the plot
        xright <- 0.011  # Width of the legend
        ytop <- .27 # Top boundary of the plot
        ybottom <- .2  # Height of the legend

        # Draw the gradient
        for (i in 1:n) {
          rect(xleft = xleft, ybottom = ybottom + (i-1)*(ytop-ybottom)/n,
               xright = xright, ytop = ybottom + i*(ytop-ybottom)/n,
               col = colors[i], border = NA)
}
        text(x = 0.006, y = .28, labels = expression(eta == 0), adj = 0,cex=.7)
        text(x = 0.004, y = .19, labels = expression(eta == 0.55), adj = 0,cex=.7)
        text(x=0.038,y=0.042,labels=expression(gamma^{-1}==2~' weeks'),adj=0,cex=1,col='white',srt=8)
        text(x=0.025,y=0.07,labels=expression(gamma^{-1}==1~' weeks'),adj=0,cex=1,col='white',srt=15)
        text(x=0.01,y=0.12,labels=expression(gamma^{-1}==0.5~' weeks'),adj=0,cex=1,col='black',srt=12)
        if(save==T){dev.off()}
}

#getfig5color(save=T)



getfig6 <- function(save=F){
	if(save==T){jpeg(file='pasc_decline.jpeg', width = 10, height = 5, units='in', res=700)}
	par(mfrow=c(1,2))
	plot(uexp~time,out.reduct,type='l',las=1,ylim=c(0,0.2),
		ylab='PASC Prevalence',xlab='Time [weeks]',lwd=3)
	lines(S~time,out.reduct,lty='dotted',lwd=3)
	lines(I~time,out.reduct,lty='dashed',lwd=3)
	legend('topleft',legend=c('PASC','Susceptible','Infected'),lty=c(1,3,2),lwd=3)
	plot(P~time,out.reduct,type='l',ylim=c(0,1),las=1,lwd=3,xlab='Time [weeks]')

	if(save==T){dev.off()}
}


getfigstencilfig_bw <- function(x,save=F,filnm='appendix_fig_1.jpeg'){
	if(save==T){jpeg(file=filnm, width = 960, height = 960,
	     pointsize = 24, quality = 100, bg = "white")
	}

	par(mfrow=c(3,3),mar=c(1,1,2,1),oma=c(5,5,1,1))	
	if(x==1){inds=1:9; ylabinds=c(1,4,7);xlabinds=c(7,8,9)}
	if(x==2){inds=10:18;ylabinds=c(10,13,16);xlabinds=c(16,17,18)}
	if(x==3){inds=19:27;ylabinds=c(19,22,25);xlabinds=c(25,26,27)}
	
	for(i in inds){ 
		getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
		if(i %in% ylabinds){axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2)} 
		if(i %in% xlabinds){axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1)} 
	

	}
	mtext('Time [weeks]',line=-48,outer=T)
	mtext('Proportion',las=3,line=3,side=2,outer=T)
	if(save==T){dev.off()}
}

getstencilfig <- function(){
	# manually run getuplot2 and adduplot to superimpose the 3 figures in getfiga on top of each other:
	par(mfrow=c(3,3),mar=c(1,1,2,1),oma=c(5,5,1,1))	

	i=19	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[10]])
	adduplot(outlist.stncl[[1]])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2)

	i=20	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[11]])
	adduplot(outlist.stncl[[2]])	

	i=21	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[12]])
	adduplot(outlist.stncl[[3]])	

	i=22	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[13]])
	adduplot(outlist.stncl[[4]])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2)

	i=23	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[14]])
	adduplot(outlist.stncl[[5]])	

	i=24	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[15]])
	adduplot(outlist.stncl[[6]])	


	i=25	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[16]])
	adduplot(outlist.stncl[[7]])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1)

	i=26	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[17]])
	adduplot(outlist.stncl[[8]])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1)

	i=27	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~I==.(stncl[i,2])))  
	adduplot(outlist.stncl[[18]])
	adduplot(outlist.stncl[[9]])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1)
}

getstencilfig2 <- function(save=F){
	if(save==T){jpeg(file='figure8_colors.jpeg', width = 960, height = 960,
		     pointsize = 24, quality = 100, bg = "white") }

	# try doing these in a different order, so that I* is changing within each window
	par(mfrow=c(3,3),mar=c(1,1,2,1),oma=c(5,5,1,1))	

	i=3	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[2]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[1]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	legend('topright',legend=c(bquote(I~'= 0.01'),bquote(I~'= 0.02'),bquote(I~'= 0.03')),col=rev(hcl.colors(3,palette='zissou 1')),
		lty=1,lwd=3)

	i=6	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[5]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[4]],col=hcl.colors(n=3,palette='zissou 1')[3])	

	i=9	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[8]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[7]],col=hcl.colors(n=3,palette='zissou 1')[3])	

	i=12	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[11]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[10]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)

	i=15	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[14]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[13]],col=hcl.colors(n=3,palette='zissou 1')[3])	

	i=18	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[17]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[16]],col=hcl.colors(n=3,palette='zissou 1')[3])	


	i=21	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[20]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[19]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)

	i=24	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[23]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[22]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)

	i=27	
	getuplot2(outlist.stncl[[i]],ylims=c(0,0.2),main=bquote(1/gamma==~.(round(1/stncl[i,3],2))~' wk, '~q==.(stncl[i,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[26]],col=hcl.colors(n=3,palette='zissou 1')[2])
	adduplot(outlist.stncl[[25]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)	
	mtext('Time [weeks]',line=-48,outer=T)
	mtext('PASC Prevalence',las=3,line=3,side=2,outer=T)

	if(save==T){dev.off()}

}


	
getstencilfig3 <- function(save=F){
	if(save==T){jpeg(file='figure8_colors.jpeg', width = 8, height = 8, units='in', res=700)}

	# try doing these in a different order, so that gamma is changing within each window and keep the color scheme identical with the earlier figures
	par(mfrow=c(3,3),mar=c(1,1,2,1),oma=c(5,5,1,1))	

	getuplot2(outlist.stncl[[7]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[1,2],2))~' '~italic(F)==.(stncl[7,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[4]],col='black')
	adduplot(outlist.stncl[[1]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	legend('top',legend=c(bquote(gamma^{-1}~'= 0.67 wks'),
				bquote(gamma^{-1}~'= 1 wk'),
				bquote(gamma^{-1}~'= 2 wks')),
			col=c(hcl.colors(3,palette='zissou 1')[1],'black',hcl.colors(3,palette='zissou 1')[3]),
		lty=1,lwd=3)

	getuplot2(outlist.stncl[[8]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[8,2],2))~' '~italic(F)==.(stncl[8,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[5]],col='black')
	adduplot(outlist.stncl[[2]],col=hcl.colors(n=3,palette='zissou 1')[3])	

	getuplot2(outlist.stncl[[9]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[9,2],2))~' '~italic(F)==.(stncl[9,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[6]],col='black')
	adduplot(outlist.stncl[[3]],col=hcl.colors(n=3,palette='zissou 1')[3])	

	getuplot2(outlist.stncl[[16]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[16,2],2))~' '~italic(F)==.(stncl[16,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[13]],col='black')
	adduplot(outlist.stncl[[10]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	getleggradcol(20,10,-.05,-.1,hcl.colors(3,'zissou 1')[1])
	getleggradcol(30,10,-.05,-.1,'black')
	getleggradcol(40,10,-.05,-.1,hcl.colors(3,'zissou 1')[3])
	text(x = par('usr')[1]+18, y = .05, labels = expression(eta == 0.55), adj = 0,cex=1)
	text(x = par('usr')[1]+22, y = .17, labels = expression(eta == 0), adj = 0,cex=1)

	getuplot2(outlist.stncl[[17]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[17,2],2))~' '~italic(F)==.(stncl[17,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[14]],col='black')
	adduplot(outlist.stncl[[11]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	text(x=140,y=0.08,labels=expression(gamma^{-1}==0.67~' wks'))
	text(x=110,y=0.04,labels=expression(gamma^{-1}==1~' wk'),col='white')
	text(x=150,y=0.025,labels=expression(gamma^{-1}==2~' wks'),col='black')

	getuplot2(outlist.stncl[[18]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[18,2],2))~' '~italic(F)==.(stncl[18,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[15]],col='black')
	adduplot(outlist.stncl[[12]],col=hcl.colors(n=3,palette='zissou 1')[3])	


	getuplot2(outlist.stncl[[25]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[25,2],2))~' '~italic(F)==.(stncl[25,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[22]],col='black')
	adduplot(outlist.stncl[[19]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(2,at=c(0,.05,.10,.15,.20),labels=c(0,.05,.10,.15,.20),las=2,pos=0)
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)

	getuplot2(outlist.stncl[[26]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[26,2],2))~' '~italic(F)==.(stncl[26,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[23]],col='black')
	adduplot(outlist.stncl[[20]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)

	getuplot2(outlist.stncl[[27]],ylims=c(0,0.2),main=bquote(I==~.(round(stncl[27,2],2))~' '~italic(F)==.(stncl[27,1])),col=hcl.colors(n=3,palette='zissou 1')[1])  
	adduplot(outlist.stncl[[24]],col='black')
	adduplot(outlist.stncl[[21]],col=hcl.colors(n=3,palette='zissou 1')[3])	
	axis(1,at=c(0,50,100,150,200),labels=c(0,50,100,150,200),las=1,pos=0)	
	mtext('Time [weeks]',line=-58,outer=T)
	mtext('PASC Prevalence',las=3,line=3,side=2,outer=T)

	if(save==T){dev.off()}

}
#getstencilfig3(save=T)
	
