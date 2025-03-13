using  LinearAlgebra
using NLsolve, Plots, RawArray
#  2 types with different forecasting abilities, 60 overlapping generations 

include("./params.jl") # is assumed to be in the same folder.
  
function cheby(x,pdeg)
	# evaluate Chebychev polynomial of degree pdeg at x
    res=ones(pdeg)
    res[1]=1.0
    res[2]=x
    for i=3:pdeg
       res[i] = 2.0*x*res[i-1]-res[i-2]
    end
    return res[pdeg]
end
function polyn(x)
	# evaluate polynomial function at 2D x
	res=zeros(tdeg)
	res[1]=1
	res[2]=cheby(x[1],2)
	res[3]=cheby(x[1],3)
	res[4]=cheby(x[1],4)
	res[5]=cheby(x[2],2)
	res[6]=cheby(x[2],3)
	res[7]=cheby(x[2],4)
	res[8]=x[1]*x[2]
return res
end
function dpolyn(x)
	# evaluate the Jacobian of the above function
	res=zeros(tdeg,2)	
	res[:,1]=[0,1,4*x[1], 12*x[1]^2-3,0,0,0,x[2]]
	res[:,2]=[0,0,0,0,1,4*x[2],12*x[2]^2-3,x[1]]
return res
end



function mut(x) # marginal utilty
	if x > eps
    	return 1/(x^sigma)
	else
		return eps^(-sigma)-sigma*(x-eps)*eps^(-sigma-1)
	end
end

function mut2(x) # derivative of marginal utility
	if x > eps
    	return -sigma/(x^(sigma+1))
	else
		return -sigma/(eps^(sigma+1))
	end
end


function mutinv(x) 
	if x > eps/100
    	return x^(-1/sigma)
	else
		return eps^(-1/sigma)-1/sigma*(x-eps)*eps^(-1/sigma-1)
	end
end


function ut(x)
    	return max(x,eps/10000)^(1-sigma)/(1-sigma)
end

function utinv(x)
	return (x*(1-sigma))^(1/(1-sigma))
end

function evalp(x)	
	# in x(hh+1,nn,nstates) cash at hand,  for Mr hh+1 irrelevant, Mr 1 only relevant in proj
	# out res[hh+1,nn,nstates], savings, zero for Mr hh+1, irrelevant for Mr 1 
res=zeros(hh+1,nn,nstates) 
if runflag == 0 
		for h=2:hh
			for n=1:nn
				res[h,n,:] .= stapc[h]*x[h,n,:]
			end
		end
else
    for h=2:hh
		for n=1:nn
			for s=1:nstates
				res[h,n,s]=lincoeff[h,n,1]+lincoeff[h,n,2]*x[h,n,s]
			end
			if h<clevern && n==1
				hx=ones(nstates,dim[h]) 
				for s=1:nstates
					y=ones(2*hh)
					y[1:hh]=x[1:hh,1,s]
					y[hh+1:2*hh]=x[1:hh,2,s]
					hx[s,1:dim[h]]=proj[h,:,1:dim[h]]'*y
					for j=1:2
						m= bounds[1,h,j]
						l= bounds[2,h,j]
						hx[s,j]=(hx[s,j]-m)/l
					end
					kstar = polyn(hx[s,:])
					res[h,n,s] += kstar'*coeff[h,:]
				end	
			end
		end
	end		
end
res[hh+1,:,:].=0.0
return res
end

function hff(x,h)
	# x[1:2*hh] contains cash at hand of 1-hh, not hh+1, ffh is help variable, denote agent
	res=zeros(2*hh)
	hx=ones(2)
	hx= proj[h,:,1:2]'*x
	for j=1:2
		m= bounds[1,h,j]
		l= bounds[2,h,j]
		hx[j]=(hx[j]-m)/l
	end
	deri=(dpolyn(hx))'*coeff[h,:]
	deri[1]=deri[1]/bounds[2,h,1]
	deri[2]=deri[2]/bounds[2,h,2]
	res=proj[h,:,1:2]*deri
	return res
#  return z[:,1]'
end

function devalp(x)	
	# in x(hh+1,nn,nstates), h=1 must be in there, hh+1 irr
	# out res[2*hh,2*hh,nstates] h=1 is ignored. h=hh+1 is all zero! Mr h's cah is h+1, since Mr 1 is 2 years old! 
res=zeros(2*hh,2*hh,nstates) 
if runflag == 0 
		for h=1:hh-1
			res[h,h+1,:] .= stapc[h+1]
			res[hh+h,hh+h+1,:] .= stapc[h+1]
		end
else
    for h=1:hh-1
		for s=1:nstates
			res[h,h+1,s]=lincoeff[h+1,1,2]
			res[hh+h,hh+h+1,s]=lincoeff[h+1,2,2]
		end
		if h+1<clevern
			for s=1:nstates
				help=zeros(2*hh)
				help[1:hh]=x[1:hh,1,s]  
				help[hh+1:2*hh]=x[1:hh,2,s]
				deri=hff(help,h+1)
			    res[h,:,s] += deri
			end
#				fun=hff(help)
#				for i=1:2*hh
#					help[i] += eps
#					zet=hff(help)
#			    	res[h,i,s] += (fun[1]-zet[1])/eps
#					help[i] += -eps
#				end
#			end
		end
	end
end
res[hh,:,:] .= 0.0
return res
end

function mpl(k,s)
	return tfp[s]*(1.0-alpha)*max(eps,k)^alpha
end

function mplk(k,s)
	return tfp[s]*(1.0-alpha)*alpha*max(eps,k)^(alpha-1)
end

function mpk(k,s)
	return tfp[s]*alpha*max(eps,k)^(alpha-1)+1-delta[s]
end

function mpkk(k,s)
	return tfp[s]*alpha*(alpha-1)*max(eps,k)^(alpha-2)
end


function f!(F,x,lolk)  # first order conditions
# for nn=2
# vars are: cap1, cap2, bond, bond2, q
# equations are: focc1, focc2, focb1, focb2, mc
	hlam=ones(hh)
	hlam[1:hh-1] .= lolk[2*hh+1]
	nextk=sum(x[h]*hfrac[h] for h=1:2*hh)
	q=x[4*hh+1]
	b=zeros(hh,nn)
	b[1:hh,1]=x[2*hh+1:2*hh+hh]
	b[1:hh,2]=x[3*hh+1:4*hh]
	hv=ones(hh,nn)
  	hv[:,1]=lolk[1:hh]-x[1:hh]-q*b[1:hh,1]
	hv[:,2]=lolk[hh+1:2*hh]-x[hh+1:2*hh]-q*b[1:hh,2]
	for h=1:hh
		F[h]=-mut(hv[h,1])
		F[hh+h]=-mut(hv[h,2])
		F[2*hh+h]=-hlam[h]*q*mut(hv[h,1])-(1-hlam[h])*b[h,1]
		F[3*hh+h]=-hlam[h]*q*mut(hv[h,2])-(1-hlam[h])*b[h,2]
	end
	F[4*hh+1]=frac*sum(b[h,1] for h=1:hh)+(1-frac)*sum(b[h,2] for h=1:hh)
	help=zeros(hh+1,2,nstates)
	hx=zeros(hh+1,2,nstates)
	for s=1:nstates
		r=mpk(nextk,s)
		wage=mpl(nextk,s)	
		for h=1:hh
			hx[h+1,1,s]=x[h]*r+en[h+1]*wage+b[h,1]*bond[s]
			hx[h+1,2,s]=x[hh+h]*r+en[h+1]*wage+b[h,2]*bond[s]
		end
		hx[1,:,s] .= en[1]*wage
	end
	help[:,:,:]=evalp(hx)	
	for s=1:nstates
		r=mpk(nextk,s)
		for h=1:hh
			mu1=mut(hx[h+1,1,s]-help[h+1,1,s])
			mu2=mut(hx[h+1,2,s]-help[h+1,2,s])
			F[h] += beta[h]*prob[s]*r*mu1
			F[hh+h] +=  beta[h]*prob[s]*r*mu2
			F[2*hh+h] +=  hlam[h]*beta[h]*prob[s]*bond[s]*mu1
			F[3*hh+h] +=  hlam[h]*beta[h]*prob[s]*bond[s]*mu2
       	end
    end
end

function df!(J,x,lolk)  # analytic Jacobian of first order conditions
	J .= 0
	hlam=ones(hh)
	hlam[1:hh-1] .= lolk[2*hh+1]
	nextk=sum(x[h]*hfrac[h] for h=1:2*hh)
	q=x[4*hh+1]
	b=zeros(hh,2)
	b[1:hh,1]=x[2*hh+1:3*hh]
	b[1:hh,2]=x[3*hh+1:4*hh]
	hv=zeros(2,hh)
	hv[1,:]=lolk[1:hh]-x[1:hh]-q*b[1:hh,1]
	hv[2,:]=lolk[hh+1:2*hh]-x[hh+1:2*hh]-q*b[1:hh,2]
		  # first deri
	for h=1:hh
		J[h,h]=mut2(hv[1,h])
		J[hh+h,hh+h]=mut2(hv[2,h])
		J[2*hh+h,h]=hlam[h]*q*mut2(hv[1,h])
		J[3*hh+h,hh+h]=hlam[h]*q*mut2(hv[2,h])
		J[h,2*hh+h]=q*mut2(hv[1,h])
		J[hh+h,3*hh+h]=q*mut2(hv[2,h])
		J[2*hh+h,2*hh+h]=hlam[h]*q^2*mut2(hv[1,h])-(1-hlam[h])
		J[3*hh+h,3*hh+h]=hlam[h]*q^2*mut2(hv[2,h])-(1-hlam[h])
	end
		# wrtq and last equation
	for h=1:hh
		J[h,4*hh+1]=b[h,1]*mut2(hv[1,h])
		J[hh+h,4*hh+1]=b[h,2]*mut2(hv[2,h])
		J[2*hh+h,4*hh+1]=hlam[h]*(q*b[h,1]*mut2(hv[1,h])-mut(hv[1,h]))
		J[3*hh+h,4*hh+1]=hlam[h]*(q*b[h,2]*mut2(hv[2,h])-mut(hv[2,h]))
		J[4*hh+1,2*hh+h]=frac
		J[4*hh+1,3*hh+h]=(1-frac)
	end
#
	Jhx=zeros(hh+1,2,4*hh,nstates) # derivatives of cash-at hand wrt k1,k2,b1,b2, only h=1:hh
	help=zeros(hh+1,2,nstates)
	dhelp=zeros(2*hh,2*hh,nstates)
	hx=zeros(hh+1,2,nstates)
	for s=1:nstates
		r=mpk(nextk,s)
		wage=mpl(nextk,s)	
		hx[1,1,s]=en[1]*wage
		hx[1,2,s]=en[1]*wage
		for h=1:hh
			hx[h+1,1,s]=x[h]*r+en[h+1]*wage+b[h,1]*bond[s]
			hx[h+1,2,s]=x[hh+h]*r+en[h+1]*wage+b[h,2]*bond[s]
			Jhx[h+1,1,h,s]+=r
			Jhx[h+1,1,2*hh+h,s] += bond[s]
			Jhx[h+1,2,hh+h,s] += r
			Jhx[h+1,2,3*hh+h,s] += bond[s]
			for i=1:2*hh
				Jhx[h+1,1,i,s]+=x[h]*mpkk(nextk,s)*hfrac[i]+en[h+1]*mplk(nextk,s)*hfrac[i]
				Jhx[h+1,2,i,s]+=x[hh+h]*mpkk(nextk,s)*hfrac[i]+en[h+1]*mplk(nextk,s)*hfrac[i]
			end
		end
		for i=1:2*hh
			Jhx[1,1,i,s]+=en[1]*mplk(nextk,s)*hfrac[i]
			Jhx[1,2,i,s]+=en[1]*mplk(nextk,s)*hfrac[i]
		end
	end
	help[:,:,:]=evalp(hx)	
	dhelp[:,:,:]=devalp(hx)
	for s=1:nstates
		r=mpk(nextk,s)
		for h=1:2*hh # FOCS wrt k AND bond
			if h <= hh
				fbeta=beta[h]
				mu=mut(hx[h+1,1,s]-help[h+1,1,s])
				mu2=mut2(hx[h+1,1,s]-help[h+1,1,s])
				hl=hlam[h]
			else
				fbeta=beta[h-hh]
				mu=mut(hx[h+1-hh,2,s]-help[h+1-hh,2,s])
				mu2=mut2(hx[h+1-hh,2,s]-help[h+1-hh,2,s])
				hl=hlam[h-hh]
			end
			for i=1:2*hh # deri wrt k's
				if h <= hh
					hjac=Jhx[h+1,1,i,s]
				else
					hjac=Jhx[h-hh+1,2,i,s]
				end
				if h == hh || h==2*hh
					J[h,i] += fbeta*prob[s]*(hfrac[i]*mpkk(nextk,s)*mu+r*hjac*mu2) # price changes change cah
					J[2*hh+h,i] += hl*fbeta*prob[s]*bond[s]*hjac*mu2	
				else
					J[h,i] += fbeta*prob[s]*(hfrac[i]*mpkk(nextk,s)*mu+r*hjac*(1-dhelp[h,h+1,s])*mu2) # price changes change cah
					J[2*hh+h,i] += hl*fbeta*prob[s]*bond[s]*hjac*(1-dhelp[h,h+1,s])*mu2
				end
				if h+1<clevern
				#  mr i's capital effect mr j's cash at hand that effects mr h's optimal decision
				# need to loop over all j=1,2hh...		
					for j=1:2*hh 
						if j <= hh
							fhjac=Jhx[j,1,i,s]
						else
							fhjac=Jhx[j-hh,2,i,s]
						end	
						if j != (h+1)							
							J[h,i]+= -fbeta*prob[s]*r*fhjac*dhelp[h,j,s]*mu2  
							J[2*hh+h,i]+= -hl*fbeta*prob[s]*bond[s]*fhjac*dhelp[h,j,s]*mu2
						end
					end
				end
			end
				# wrt bond
			for i=2*hh+1:4*hh
				if i==2*hh+h
					if h<= hh
						hjac=Jhx[h+1,1,i,s]
					else
						hjac=Jhx[h-hh+1,2,i,s]
					end
					if h==hh || h==2*hh
						J[h,i] += fbeta*prob[s]*r*hjac*mu2
						J[2*hh+h,i]+= hl*fbeta*prob[s]*bond[s]*hjac*mu2
					else	
						J[h,i] += fbeta*prob[s]*r*hjac*(1-dhelp[h,h+1,s])*mu2
						J[2*hh+h,i]+= hl*fbeta*prob[s]*bond[s]*hjac*(1-dhelp[h,h+1,s])*mu2
					end
				elseif h+1 < clevern
					if i < 3*hh
						J[h,i] += -fbeta*prob[s]*r*Jhx[i-2*hh+1,1,i,s]*dhelp[h,i-2*hh+1,s]*mu2
						J[2*hh+h,i]+= -hl*fbeta*prob[s]*bond[s]*Jhx[i-2*hh+1,1,i,s]*dhelp[h,i-2*hh+1,s]*mu2
					elseif i > 3*hh && i < 4*hh
						J[h,i] += -fbeta*prob[s]*r*Jhx[i-3*hh+1,2,i,s]*dhelp[h,i-2*hh+1,s]*mu2
						J[2*hh+h,i]+= -hl*fbeta*prob[s]*bond[s]*Jhx[i-3*hh+1,2,i,s]*dhelp[h,i-2*hh+1,s]*mu2
					end
				end
			end
		end
	end
end

function dfdlolk(x,lolk)  # derivative  of first order conditions w.r.t. lolk
	nextk=sum(x[h] for h=1:2*hh)
	hlam=ones(hh)
	hlam[1:hh-1] .= lolk[2*hh+1]
	q=x[4*hh+1]
	b=zeros(hh,nn)
	b[1:hh,1]=x[2*hh+1:2*hh+hh]
	b[1:hh,2]=x[3*hh+1:4*hh]
	hv=ones(hh,nn)
  	hv[:,1]=lolk[1:hh]-x[1:hh]-q*b[1:hh,1]
	hv[:,2]=lolk[hh+1:2*hh]-x[hh+1:2*hh]-q*b[1:hh,2]
	F = zeros(4*hh+1,2*hh)
  for h=1:hh
	  F[h,h]=-mut2(hv[h,1])
	  F[hh+h,hh+h]=-mut2(hv[h,2])
	  F[2*hh+h,h]=-hlam[h]*q*mut2(hv[h,1])
	  F[3*hh+h,hh+h]=-hlam[h]*q*mut2(hv[h,2])
  end
  return F 
end

function dfdx(x,lolk) # wrapper for Jacobian, can be switched to finite differences 
  y=zeros(4*hh+1,4*hh+1)
  df!(y,x,lolk)
  # Still the finite difference version...
 # fdeps=0.000001
 # F=zeros(4*hh+1)
 # f!(F,x,lolk)
 #   h1=F
 # for j=1:4*hh+1
#	  x[j]=x[j]+fdeps
#	  F=zeros(4*hh+1)
#	  f!(F,x,lolk)
#	  h2=F
#	  for jj=1:4*hh+1
#		  y[jj,j]=(h2[jj]-h1[jj])/fdeps
#	  end
#	  x[j]=x[j]-fdeps
# end
  return y
end

function soleq(st,lolk,l,t) # wrapper for non-linear solver nlsolve
# wrapper for nonlinear solver
	stapo=st[1:4*hh+1,1]
	F=zeros(4*hh+1)
	J=zeros(4*hh+1,4*hh+1)
	heps=0.00005
#
	a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), stapo)
	a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), a.zero)
	if a.residual_norm > heps
		stapo=st[1:4*hh+1,2]
		a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), stapo)
		a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), a.zero)
	end	
#    a=nlsolve((F,x) ->f!(F,x,lolk), stapo)
	# println("Solved with error",a.residual_norm)
	i=1
	if a.residual_norm > heps
		println("First Error",a.residual_norm)
		println("I and T",l,t)
		for i=1:hh
			stapo[i]=lolk[i]*hpc[i]
			stapo[i+hh]=lolk[i+hh]*hpc[i]
		end
		stapo[2*hh+1:4*hh].=0.0
		stapo[4*hh+1]=0.95
#		a=nlsolve((F,x) ->f!(F,x,lolk), stapo)
#		a=nlsolve((F,x) ->f!(F,x,lolk), a.zero)
		a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), stapo)
		a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), a.zero)
	end
	if a.residual_norm > heps
		println("Ohje Error ",a.residual_norm)
 		while( a.residual_norm > heps)
		  stapo2=st[:,1] +(rand(4*hh+1).-0.5)/50
		  a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), stapo2)
		  a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), a.zero)
		  i+=1
		  if a.residual_norm > heps
			stapo2=st[:,2] +(rand(4*hh+1).-0.5)/50
			a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), stapo2)
			a=nlsolve((F,x) ->f!(F,x,lolk), (J,x) ->df!(J,x,lolk), a.zero)
		  end
		  if i>1
			heps=heps*1.5
			i=0
		  end
		end
	println("After some struggle, error: ",a.residual_norm,)
	end
if abs(a.residual_norm) > 0.001
	res=st[1:4*hh+1,1]
else
  res=a.zero
 end
 return res
end
#

function fiax(a,y,lflag)  # fit new coefficients, here simple ols
	beta=(a'*a) \ (a'*y)
	me=abs.(y-a*beta)
		println("Mr",lflag)
		println("maxerr",maximum(me))
		println("avgerr",norm(me)/sqrt(maxtime*maxsim))
	return beta,norm(me)/sqrt(maxtime*maxsim)
 end

 function newcoeff(lambda) # determine new coefficients, here via ols
	x1=zeros(hh+1,nn,maxsim*maxtime)
	y1=zeros(hh,nn,maxsim*maxtime)
	for sim=1:maxsim*maxtime
		for n=1:nn
			x1[1:hh,n,sim]=ah[:,n,sim]
			y1[:,n,sim]=yh[:,n,sim]
		end
	end
#	for h=1:hh
#		m=minimum(x1[h,:])
#		ll=maximum(x1[h,:])-m
#		global boundsadj[h,1]=m
#		global boundsadj[h,2]=ll
#		x1[h,:]=(x1[h,:].-m)./ll
#	end
	# Hopefully all variables correctly defined
		#Threads.@threads :dynamic
		err=zeros(hh,nn)
	for h=2:hh 
		for n=1:nn
			xhelp=zeros(2,maxtime*maxsim)
			for i=1:maxsim*maxtime
				xhelp[1,i]=1.0
				xhelp[2,i]=x1[h,n,i]
			end
			lincoeff[h,n,:].=0.0
			ikh,err[h,n]=fiax(xhelp[1:2,:]',y1[h,n,:],h)
			lincoeff[h,n,1]=ikh[1]
			lincoeff[h,n,2]=ikh[2]
			y1[h,n,:]=y1[h,n,:]-xhelp[1:2,:]'*ikh
			if n==1 && h<clevern
				cxhelp=zeros(dim[h],maxtime*maxsim)
				evec=eigvecs(magma[h,:,:])
				heig=eigvals(magma[h,:,:])
				println(heig[2*hh-3:2*hh])
		#		println(evec[:,2*hh-1:2*hh])
				global proj[h,:,1:dim[h]]=evec[:,(2*hh+1-dim[h]):2*hh]
				for i=1:maxsim*maxtime
					help=ones(2*hh)
					help[1:hh]=x1[1:hh,1,i]
					help[hh+1:2*hh]=x1[1:hh,2,i]
					cxhelp[:,i]=proj[h,:,1:dim[h]]'*help 
				end
				for j=1:2
					m=minimum(cxhelp[j,:])
					l=maximum(cxhelp[j,:])-m
					global bounds[1,h,j]=m
					global bounds[2,h,j]=l
					cxhelp[j,:]=(cxhelp[j,:].-m)./l
				end

	
				kstar = zeros(tdeg,maxsim*maxtime)
				for i=1:maxsim*maxtime
					kstar[:,i]=polyn(cxhelp[:,i])
				end
				co,err[h,1]=fiax(kstar', y1[h,1,:],h)
				global coeff[h,1:tdeg]=co*lambda
			end
		end
	end
	return err
end

function solvemodel() # Solve the model!
 global runflag=0 # 0 if no starting point
 cshocks=zeros(Int64,maxtime*maxsim)
 for i=1:maxtime*maxsim
	cshocks[i]=nstates
	h=0
	rn=rand(Float64)
	for s=1:nstates
		h+=prob[s]
		if rn<h
			cshocks[i]=s
			break
		end
	end
	global shock=cshocks[i]
 end

 global sol=zeros(4*hh+1,maxtime*maxsim)
 global oldk=zeros(hh+1,2,maxsim)
 global holdk=zeros(hh+1,2,maxsim)
 global holdb=zeros(hh+1,2,maxsim)
 global oldb=zeros(hh+1,2,maxsim)
 global oldk[2:hh+1,:,:] .=ones(hh)*ststk/(2*hh)
 global oldk[1,:,:] .=0.0
 global oldak=zeros(maxsim)
 global oldak[:]=frac*sum(oldk[h,1,:] for h=1:hh+1)+(1-frac)*sum(oldk[h,2,:] for h=1:hh+1)
 global magma=zeros(hh,2*hh,2*hh)
 global runflag=0
 for iter=1:maxiter
	lambda=0.0
	if iter >= 25
		lambda=(iter-25)/100
	end
	if lambda >= 1
		 lambda=1
	end
	if iter > 125 # after 125 iterations we fix initial conditions...
		global oldk[:,:,:]=holdk
		global oldb[:,:,:]=holdb
		global oldak[:]=frac*sum(oldk[h,1,:] for h=1:hh+1)+(1-frac)*sum(oldk[h,2,:] for h=1:hh+1)
	end
	global avdis=zeros(hh,2)
	global shock=[1,2,3,4,1,2,3,4,3,2] # This is hard-coded for maxsim=10, needs to be arrayo f size maxim
	println("Iteration number",iter)
	lmagma=zeros(hh,2*hh,2*hh,maxsim)
	# Start simulations...
	 	@time Threads.@threads :dynamic for i=1:maxsim
			for t=1:maxtime
				cah=ones(2*hh)	
				cah[1:hh]=oldk[1:hh,1,i]*mpk(oldak[i],shock[i])+
       				en[1:hh]*mpl(oldak[i],shock[i])+oldb[1:hh,1,i]*bond[shock[i]]
		  		cah[hh+1:2*hh]=oldk[1:hh,2,i]*mpk(oldak[i],shock[i])+
		    		en[1:hh]*mpl(oldak[i],shock[i])+oldb[1:hh,2,i]*bond[shock[i]]
				global ah[1:hh,1,t+(i-1)*maxtime]=cah[1:hh]
				global ah[1:hh,2,t+(i-1)*maxtime]=cah[hh+1:2*hh]
				st=zeros(4*hh+1,2)
				if iter <= 1
					for h=1:hh
						st[h,:] .=cah[h]*hpc[h]
						st[h+hh,:] .=cah[h+hh]*hpc[h]
					end
					st[4*hh+1]=1.0
				else
					st[:,1]=sol[:,(i-1)*maxtime+t]
					st[:,2]=sol[:,(i-1)*maxtime+t]
					if t>1
						st[:,2]=sol[:,(i-1)*maxtime+t-1]
					end
				end
				lolk=zeros(2*hh+1)
				lolk[1:2*hh]=cah
				lolk[2*hh+1]=lambda
				res=soleq(st,lolk,i,t)
				global sol[:,(i-1)*maxtime+t]=res[:]
				kap=res[1:2*hh]
				q=res[4*hh+1]
				b=zeros(hh,2)
				b[1:hh,1]=res[2*hh+1:3*hh]
				b[1:hh,2]=res[3*hh+1:4*hh]
				cons=zeros(hh,nn)
				cons[:,1]=cah[1:hh]-kap[1:hh]-q*b[1:hh,1]
				cons[:,2]=cah[hh+1:2*hh]-kap[hh+1:2*hh]-q*b[1:hh,2]
				for h=1:hh
					global avdis[h,1] += ut(cons[h,1])
					global avdis[h,2] += ut(cons[h,2])
				end
				hhelp=-(dfdx(res,lolk)) \ dfdlolk(res,lolk)  # is of dimension 4hh+1 X 2hh
				help=zeros(hh,2*hh) # 2hh-dim gradients for hh agents
				for h=1:hh
					help[h,:]=hhelp[h,:]+hhelp[2*hh+h,:]*q+hhelp[4*hh+1,:]*b[h,1] 
#					help[hh+h,:]=hhelp[hh+h,:]+hhelp[3*hh+h,:]*q+hhelp[4*hh+1,:]*b[h,2] 
					help[h,h]= help[h,h]-lincoeff[h,1,2]
#					help[hh+h,hh+h]= help[hh+h,hh+h]-lincoeff[h,2,2]
				end
				for h=1:hh
					lmagma[h,:,:,i]=lmagma[h,:,:,i]+help[h,:]*(help[h,:])'
				end
				if iter > maxiter -5
					println(yh[1:5,1,t+(i-1)*maxtime]-(kap[1:5]+b[1:5,1]*q))
				end
				for n=1:nn
					global yh[:,n,t+(i-1)*maxtime]=kap[1+hh*(n-1):hh+hh*(n-1)]+b[1:hh,n]*q
				end
				global oldak[i]=sum(kap[h]*hfrac[h] for h=1:2*hh)
				global shock[i]=cshocks[(i-1)*maxtime+t]
				global oldk[2:hh+1,1,i]=kap[1:hh]
				global oldk[2:hh+1,2,i]=kap[hh+1:2*hh]
				global oldk[1,:,i] .=0.0
				global oldb[2:hh+1,1,i]=b[1:hh,1]
				global oldb[2:hh+1,2,i]=b[1:hh,2]
				global oldb[1,:,i] .=0.0
			end
	   end
	global magma[:,:,:]=sum(lmagma[:,:,:,i] for i=1:maxsim)
	err=newcoeff(lambda)
	println("...................")
	global avdis=avdis/(maxtime*maxsim)
	gain=zeros(hh)
	for h=1:hh
		gain[h]=(utinv(avdis[h,1])-utinv(avdis[h,2]))/utinv(avdis[h,1])
	end
	println("Welfare difference", gain)
	println("Average cons.",utinv.(avdis[:,1]))
	println("Avg error1",err[:,1]./utinv.(avdis[:,1]))
	println("Avg error2",err[:,2]./utinv.(avdis[:,2]))
	global runflag=1
 #   male()
	if iter == 125
		global holdk[:,:,:]=oldk
		global holdb[:,:,:]=oldb
		println(holdk) 
	end
  end
end


# Start of main program...

solvemodel()

# Save coefficients and active subspaces if needed:

#cd("/Users/felixkuebler/Dropbox/JULIA/dataIM/2agents")
#rawrite(xstar,"xstarIM0.dat") 
#rawrite(coeff,"coeffIM0.dat")
#rawrite(lincoeff,"lincoeffIM0.dat")
#rawrite(proj,"proj0.dat")
#rawrite(sol,"sol0.dat")
#rawrite(anop,"anop0.dat")
#rawrite(oldk,"oldk0.dat")
#rawrite(oldb,"oldb0.dat")
