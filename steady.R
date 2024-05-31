rm(list=ls())

# make plots for the steady state results section

#graph beta as a function of eta, for fixed I
# J represnts I^*, b/c Maple reserves I for sqrt(-1)

getBeta <- function(eta,J,alpha,gamma,delta){
	beta= (-J*alpha*eta + J*alpha + J*delta + gamma*J + alpha*eta - alpha - delta + sqrt(J^2*alpha^2*eta^2 - 2*J^2*alpha^2*eta - 2*J^2*alpha*delta*eta + 2*J^2*alpha*eta*gamma + 4*J^2*delta*eta*gamma - 2*J*alpha^2*eta^2 + J^2*alpha^2 + 2*J^2*alpha*delta - 2*J^2*alpha*gamma + J^2*delta^2 - 2*J^2*delta*gamma + J^2*gamma^2 + 4*J*alpha^2*eta + 4*J*alpha*delta*eta - 2*J*alpha*eta*gamma - 4*J*delta*eta*gamma + alpha^2*eta^2 - 2*J*alpha^2 - 4*J*alpha*delta + 2*J*alpha*gamma - 2*J*delta^2 + 2*J*delta*gamma - 2*alpha^2*eta - 2*alpha*delta*eta + alpha^2 + 2*alpha*delta + delta^2))/(2*J*(J*eta - J - eta + 1))
	return(beta)
}

getS <- function(eta,J,alpha,gamma,delta){
	beta = getBeta(eta,J,alpha,gamma,delta)
	S = -delta*(J - 1)/(J*beta + alpha + delta)
	return(S)

}

getP <- function(eta,J,alpha,gamma,delta){
	beta = getBeta(eta,J,alpha,gamma,delta)
	P = -(J - 1)*(J*beta + alpha)/(J*beta + alpha + delta)
	return(P)

}


getSFrac <- function(eta,J,alpha,gamma,delta){

	SFrac = getS(eta,J,alpha,gamma,delta)/(getS(eta,J,alpha,gamma,delta)+getP(eta,J,alpha,gamma,delta))
	return(SFrac)
}


getSIncFrac <- function(eta,J,alpha,gamma,delta){

	SIncFrac = getS(eta,J,alpha,gamma,delta)/(getS(eta,J,alpha,gamma,delta)+(1-eta)*getP(eta,J,alpha,gamma,delta))
	return(SIncFrac)
}

f = 1/20
getPhi <- function(eta,J,alpha,gamma,delta){
	Phi = ( getS(eta,J,alpha,gamma,delta) + (1-eta)*getP(eta,J,alpha,gamma,delta) )/( getS(eta,J,alpha,gamma,delta) + (1-VE_T)*getP(eta,J,alpha,gamma,delta) )
	Phi = f*Phi
	return(Phi)

}

#make the range of eta values

VE_T = 0.8

etavals <- seq(0,VE_T,length=200)

J = .02
alpha = 1/50
gamma = 1
delta = 1/50

betavals <- sapply(etavals, function(x) getBeta(x,J,alpha,gamma,delta)  )
phivals <- sapply(etavals, function(x) getPhi(x,J,alpha,gamma,delta) )
SFracvals <- sapply(etavals, function(x) getSFrac(x,J,alpha,gamma,delta)  )
SIncFracvals <- sapply(etavals, function(x) getSIncFrac(x,J,alpha,gamma,delta)  )


jpeg(file='steady_state_heterogeneity.jpeg',quality=100,height=960,width=960,pointsize=24,bg='white')
par(mfrow=c(2,2),mar=c(5,6,4,1))

plot(etavals,betavals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('Transmission rate ('~beta~')'),lwd=2)
plot(etavals,phivals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('Baseline PASC probability ('~phi~')'),lwd=2)
plot(etavals,SFracvals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('P(Susceptible | Uninfected)'),lwd=2)
plot(etavals,SIncFracvals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote(atop('Proportion of New Infections ','from Susceptibles')),lwd=2)

dev.off()


# make a plot with a few different values of gamma

J = .02
alpha = 1/50
gamma = .5
delta = 1/50

betavals.5 <- sapply(etavals, function(x) getBeta(x,J,alpha,gamma,delta)  )
phivals.5<- sapply(etavals, function(x) getPhi(x,J,alpha,gamma,delta) )
SFracvals.5  <- sapply(etavals, function(x) getSFrac(x,J,alpha,gamma,delta)  )
SIncFracvals.5 <- sapply(etavals, function(x) getSIncFrac(x,J,alpha,gamma,delta)  )

gamma = 2
betavals2 <- sapply(etavals, function(x) getBeta(x,J,alpha,gamma,delta)  )
phivals2 <- sapply(etavals, function(x) getPhi(x,J,alpha,gamma,delta) )
SFracvals2  <- sapply(etavals, function(x) getSFrac(x,J,alpha,gamma,delta)  )
SIncFracvals2 <- sapply(etavals, function(x) getSIncFrac(x,J,alpha,gamma,delta)  )

cols = hcl.colors(n=2,palette='zissou 1')

jpeg(file='steady_state_heterogeneity_colors.jpeg',quality=100,height=960,width=960,pointsize=24,bg='white')
par(mfrow=c(2,2),mar=c(5,6,4,1))

plot(etavals,betavals/1,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('Reproduction number ('~R[0]~')'),lwd=2,
	xlim=c(0,VE_T),ylim=c(0,5),las=1)
lines(etavals,betavals.5/.5,col=cols[2],lwd=2)
lines(etavals,betavals2/2,col=cols[1],lwd=2)
legend('topleft',legend=c(bquote(gamma^{-1}~'=0.5 weeks'),bquote(gamma^{-1}~'=1 week'),bquote(gamma^{-1}~'=2 weeks')),
	col=c(cols[1],'black',cols[2]),lwd=2)

plot(etavals,phivals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('Baseline PASC probability ('~phi~')'),lwd=2,
	xlim=c(0,VE_T),ylim=c(0,0.15),las=1)
lines(etavals,phivals.5,col=cols[2],lwd=2)
lines(etavals,phivals2,col=cols[1],lwd=2)


plot(etavals,SFracvals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote('P(Susceptible | Uninfected)'),lwd=2,
	xlim=c(0,VE_T),ylim=c(0.1,0.4),las=1)
lines(etavals,SFracvals.5,col=cols[2],lwd=2)
lines(etavals,SFracvals2,col=cols[1],lwd=2)


plot(etavals,SIncFracvals,type='l',xlab=bquote('Primary efficacy ('~eta~')'),
	ylab=bquote(atop('Proportion of New Infections ','from Susceptibles')),lwd=2,
	xlim=c(0,VE_T),ylim=c(0.25,0.7),las=1)
lines(etavals,SIncFracvals.5,col=cols[2],lwd=2)
lines(etavals,SIncFracvals2,col=cols[1],lwd=2)

dev.off()
