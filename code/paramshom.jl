################################################################################
# paramshom.jl
#
# Description:
#   Parameter file for the 2-type, 60-overlapping-generations model.
#   Contains global constants, arrays, and initial conditions read
#   by mainIM.jl. Adjust these values to modify the simulation's
#   economic environment (e.g., preferences, production, or transitions).
#
# Usage:
#   1) This file is included by "mainhom.jl".
#   2) Make any changes to the constants or arrays as needed, then re-run.
#
################################################################################

const hh=60 # number of generations
nn=2 # types per generation
const retire=46
hs=ones(hh+1)*0.96
#    hs[1:5].=0.9  # age dependent discounting is possible but not calibration in the paper
#    hs[6:10].=0.93
const beta=hs
const sigma=2 # coefficient of relative risk aversion
const nstates=4
const maxiter = 300 # max number of iterations 
const maxtime=200  # Simulation steps in each iteration
const maxsim= 10
const alpha=0.33 # capital share
const bond=[1,1,0.85,1] # payoffs of risky bond, 0.85 at s=3
const clevern=1
const frac=0.5 # fraction of sophisticated traders; corresponds to \Psi, e.g., Psi = 0, 0.5, 0.9
hs=zeros(2*hh)
hs[1:hh].=frac
hs[hh+1:2*hh].=1-frac
const hfrac=hs

# determine life-cycle labor profile
hen=zeros(hh+1)
for i=1:retire
 	hen[i]=1/100.0*exp(4.47+0.033*i-0.0006*i^2)
end
hen[retire]=hen[retire-1]/2
for i=retire+1:hh+1
	hen[i]=0.0
end
const en=hen[:]/(sum(hen[i] for i=1:hh+1))

# tfp shocks, depreciation shocks and probabilities
const tfp=[1,1.15,1,1.15]*6.0
const delta=[0.1, 0.1, 0.2, 0.2]           
const prob=[0.4, 0.4, 0.1, 0.1]

# steady state capital
const ststk=((1/beta[1]-(1.0-delta[1]))/(alpha*tfp[1]))^(1.0/(alpha-1.0))

const eps=0.00001 # very small number
global runflag::Int = 0


hs=zeros(Int64,clevern)
hs[1:clevern].=2
# hs[4:clevern].=1
const dim=hs
const maxdim=2
const tdeg=8
global coeff=zeros(hh,tdeg)
global lincoeff=zeros(hh,nn,2)
global bounds=zeros(2,clevern,maxdim)
global ah=zeros(hh,nn,maxsim*maxtime) # x-points for regression bond
global yh=zeros(hh,nn,maxsim*maxtime)
global proj=zeros(hh,nn*hh,maxdim)

# find starting coefficients
hpc=ones(hh)
for h=1:hh
    hpc[h]=1-1/(1+sum(prod(beta[i] for i=h:k) for k=h:hh))
end
# global lincoeff[:,2]=hpc[:]
const stapc=hpc
global oldk=zeros(hh+1,nn,maxsim)
global oldb=zeros(hh+1,nn,maxsim)

# Potentially placeholders for reading data from files (commented)
# global coeff=raread("coeffIM24T.dat")
# global xstar=raread("xstarIM24T.dat")
# ...
