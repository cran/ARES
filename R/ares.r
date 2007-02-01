.packageName <- "ARES";

"aresCalc" <- function(rawdata, bootsize=NULL, maxsize=NULL)
{
	alri.est = 0;
	
	# bootsize = bootstrap size
	if(is.null(bootsize))
	{
		bootsize=200;
	}

	# how far to extrapolate
	if(is.null(maxsize))
	{
		maxsize=100;
	}

	# convert the data to binary matrix
	rawdata = rawdata>0;

	# convert the data to 'dens', which is used by the other functions
	dens = fmtbino(x=rowSums(rawdata),size=ncol(rawdata));
	size = length(dens);

	# Compute the 'allelic accumulation curve' from one to the nr of observed species,
	# with moment-based estimates.
	total = attr(dens, "total");
	est.mbe = aftbino.mbe(count = dens * total);

	if(size >= maxsize)
	{
		est.tot = est.mbe[1:maxsize,2:4];
	}
	else
	{
		# The calculations to extrapolate & get confidence bounds
		Q = vtbino(dens);

		# remove zeroprob from Q, to prevent an infinity returned by omtbino
		Qp = Q;
		Qp[,1] = pmax(1e-8,Qp[,1]);
		richness = round((omtbino(size, Qp) + 1) * total);

		# extra safety, but probably a richness of one million does not make sense
		# anyway - to be chaned in a future version.
		if (is.infinite(richness)|richness>1e6)
		{
			richness=1e6;
		}

		nvec = rbinom(bootsize, richness, total/richness);
		prob = dmtbino(1:size, size, Q);
		dmat = sapply(nvec, function(x) rmultinom(1, x, prob));
		dmat = dmat/matrix(nvec, size, bootsize, byrow = TRUE);
		Qlist = apply(dmat, 2, vtbino);
		while (length(Qlist)<bootsize) {
       Qlist = apply(dmat, 2, vtbino);
       }
		estimat = matrix(0, maxsize, bootsize);
		
		for (i in 1:bootsize) 
			estimat[, i] = aftbino(size, total = nvec[i], Q = Qlist[[i]], sizes = 1:maxsize)[,2];
		
		ci = t(apply(estimat, 1, quantile, c(0.025, 0.975)));
		esti = aftbino(size, total = total, Q = Q, sizes = 1:maxsize);
		est.mle = cbind(esti[, 2], ci);
		
		# glue the moment-based and maximum-likelyhood etsimation together
		alri.tot = rbind(est.mbe[,2:4], est.mle[(size+1):maxsize,1:3]);
	}

	return(alri.tot);

}

"aresPlot" <- function(output_aresCalc,T = NULL)
{

	if(is.null(T))
	{
		T = "";
	}
	
	nrind <- nrow(output_aresCalc);   # number of individuals
	xl <- c(1, nrind);
	yl <- c( min(output_aresCalc[,2]), max(output_aresCalc[,3]) );          # axis limits
	plot(1:nrind,output_aresCalc[,1], type = "l", col = "black", 
			ylab = "Allelic Richness", xlab = "Number of Individuals",
			xlim = xl, ylim = yl);
	title(main = T);
	lines(1:nrind,output_aresCalc[,2], type = "l", col = "red");
	lines(1:nrind,output_aresCalc[,3], type = "l", col = "red");
}


"fmtbino" <- function(count=NULL,x=NULL,size=NULL)
{
	if(is.null(count))
	{
		x=as.data.frame(table(x[x>0]));
		x[,1]=unclass(x[,1]);
		count=rep(0,size);
		count[x[,1]]=x[,2];
	}
	
	total=sum(count);
	attr(count,"total")=total;
	return(count/total);
}

"rpldtbino" <- function(Q,size,to.conditional=FALSE)
{
	x=1-(1-Q[,1])^size;
	
	if(to.conditional)
	{
		Q[,2]=Q[,2]*x;
	}
	else
	{
		Q[,2]=Q[,2]/x;
		Q[,2]=Q[,2]/sum(Q[,2]);
	}
	
	return(Q)
}

"dtbino" <- function(x,size,prob)
{
	z=suppressWarnings(dbinom(x,size,prob)/(1-dbinom(0,size, prob)));
	return(ifelse(is.finite(z),z,ifelse(x==0,Inf,ifelse(x==1,1,0))));
}

"dmtbino" <- function(x,size,Q) 
{
	return(drop(sapply(Q[,1],dtbino,x=x,size=size)%*%Q[,2]));
}

"omtbino" <- function(size,Q)
{
	return(sum(1/(1/(1-Q[,1])^size-1)*Q[,2]));
}

"ntbino" <- function(mu, size, neps=1e-8)
{
	out=rep(0,length(mu));
	out[mu>=size]=1;
	i=(mu>1)&(mu<size);
	
	if(any(i)) 
	{
		mu=mu[i];
		p=mu/size;
		d=rep(2, length(p));
		
		while(max(abs(d))>neps) 
		{
			d=mu*(1-p)^(size-1);
			d=(size*p-mu+(1-p)*d)/size/(1-d);
			p=p-d;
		}	
		
		out[i]=p;
	}
	return(out);
}

"gtbino" <- function(dens,Q,prob) 
{
	size=length(dens);
	i=(1:size)[dens>0];
	prob=sapply(i,dtbino,size=size,prob=prob);
	prob=prob%*%(dens[i]/dmtbino(x=i,size=size,Q=Q));
	return(drop(prob)-1);
}

"emtbino" <- function(dens,G=1,Q=NULL,control=list()) 
{

	ctrl=list(eeps=1e-06,emax=5000);
	ctrl[names(control)]=control;
	size=length(dens);
	i=which(dens>0);
	dens=dens[i];

	if(is.null(Q))
	{
		if(G>1)
		{
			Q=ntbino(mu=range(i),size=size);
			Q=cbind(sort(runif(G,Q[1],Q[2])),(G:1)/sum(1:G));
		}
		else
		{
			return(cbind(ntbino(mu=sum(dens*i),size=size),1));
		}
	}
		
	like=sum(dens*log(dmtbino(x=i,size=size,Q=Q)));
	
	for(j in 1:ctrl$emax)
	{
		mat=rbind(sapply(i,dtbino,prob=Q[,1],size=size));
		mat=t(mat*Q[,2]);
		mat=mat*(dens/rowSums(mat));
		Q[,2]=colSums(mat);
		Q[,1]=ntbino(drop(i%*%mat)/Q[,2],size);
		nlike=sum(dens*log(dmtbino(x=i,size=size,Q=Q)));
		if(nlike-like<ctrl$eeps)
		{
			break;
		}
		else
		{
			like = nlike;
		}
	}
	
	return(Q);
}

"sltbino" <- function(dens,Q,prob) 
{

	size=length(dens);
	i=which(dens>0);
	dens=dens[i];
	r=which.min(gtbino(dens=dens,Q=Q,prob=Q[,1]));
	f=dmtbino(x=i,size=size,Q=Q);
	fadd=dtbino(x=i,size=size,prob=prob);
	frem=dtbino(x=i,size=size,prob=Q[r,1]);
	a=(fadd-frem)*Q[r,2]/f;
	d1=c(sum(dens*a),sum(dens*a/(a+1)));
	d2=-c(sum(dens*a^2), sum(dens*(a/(a+1))^2));
	
	if(all(is.finite(c(d1,d2))))
	{
		if((d2[1]<=d2[2])||(d2[1]<=(d1[2]-d1[1])))
		{
			a=-d1[1]/d2[1];
		}
		else
		{
			a=-d1[1]/(d1[2]-d1[1]);
		}	
		
		a=ifelse(a>1,NA,ifelse(a<0,NA,a));
	}
	else
	{
		a=NA;
	}
	
	if(is.na(a))  
	{
		z=Q[r,2]*(fadd-frem);
		lik=function(p){ return(-sum(dens*log(f+p*z))); }
		a=suppressWarnings(optimize(lik,0:1)$minimum)
	}
	
	if(a>=1)
	{
		Q[r,1]=prob;
	}
	else
	{
		Q=cbind(c(prob,Q[r,1],Q[-r,1]),c(a*Q[r,2],(1-a)*Q[r,2],Q[-r,2]));
	}
	
	return(Q);
	
}

"vtbino" <- function(dens, Q=NULL, control=list()) 
{

	ctrl=list(geps=1e-3,gmax=100,length=1000,span=0.05,eeps=1e-06,emax=5000);
	ctrl[names(control)]=control;
	size=length(dens);
	if(is.null(Q)) 
	{
		Q=emtbino(dens=dens,G=floor((size-1)/2),control=ctrl);
		Q=mtbino(Q,span=ctrl$span);
		Q=emtbino(dens=dens,Q=Q,control=ctrl);
	}
	
	step=1;
	grid=range(which(dens>0));
	grid=seq(grid[1],grid[2],length=ctrl$length);
	grid=ntbino(mu=grid,size=size);
	grad=gtbino(dens=dens,Q=Q,prob=grid);
	i=which.max(grad);
	while((grad[i]>ctrl$geps)&&(step<=ctrl$gmax))
	{
		Q=sltbino(dens=dens,Q=Q,prob=grid[i]);
		Q=emtbino(dens=dens,Q=Q,control=ctrl);
		grad=gtbino(dens=dens,Q=Q,prob=grid);
		i=which.max(grad);
		step=step+1;
	}
	
	Q=mtbino(Q,span=0);
	attr(Q,"convergence")=(grad[i]<ctrl$geps);
	return(Q);
}

"batbino" <- function(dens,Q,total=NULL,control=list())
{
	ctrl=list(eeps=1e-06,emax=5000);
	ctrl[names(control)]=control;
	k=nrow(Q);
	Q=list(Q);
	if(k>1)
	{
		for (i in 2:k)
		{
			Q[[i]]=mpair(Q=Q[[(i-1)]]);
			Q[[i]]=emtbino(dens=dens,Q=Q[[i]],control=ctrl);
		}
		Q=rev(Q);
	}
	
	if(!is.null(total))
	{
		size=length(dens);
		attr(Q,"loglik")=total*sapply(Q,function(D) sum(dens*log(dmtbino(x=1:size,size=size,Q=D))))
	}
	
	return(Q);
}

"aftbino" <- function(size,total,Q,sizes=NULL) 
{
	Q=mtbino(Q,span=0);
	
	if(is.null(sizes))
	{
		sizes=1:(4*size);
	}
	
	if(Q[1,1]>0)
	{
		out=Q[,2]/(1-(1-Q[,1])^size);
		out=drop(out%*%(1-outer(1-Q[,1],sizes,"^")));
	}
	else
	{
		out=Q[-1,2]/(1-(1-Q[-1,1])^size);
		out=drop(out%*%(1-outer(1-Q[-1,1],sizes,"^")));
		out=out+Q[1,2]*sizes/size;
	}
	
	out=cbind(size=sizes,esti=total*out);
	
	return(invisible(out));
}

"aftbino.mbe" <- function(count,estimated.richness=NULL,conf.level=0.95) 
{
	size=length(count);
	if (is.null(estimated.richness))
	{
		estimated.richness=ifelse(all(count[1:2]>0), count[1]^2/count[2]/2/size*(size - 1)+sum(count), Inf);
	}
	posi=(1:size)[count>0];
	count=count[count>0];
	mat=1-outer(size-(1:size),posi, choose)/outer(rep(size,size),posi,choose);
	esti=drop(mat%*%count);
	se=drop(mat^2%*%count)-esti^2/estimated.richness;
	se=suppressWarnings(sqrt(se))*qnorm((1+conf.level)/2);
	out=cbind(size=1:size,esti=esti,lwr=esti-se, upr=esti+se);
	return(invisible(out));
}

"mtbino" <- function(Q,span)
{
	Q=Q[order(Q[,1]),,drop=FALSE];
	k=nrow(Q);

	if(k>1)
	{
		l=Q[-1,1]-Q[-k, 1];
		r=(1:k)[c(l,2)>span];
		l=(1:k)[c(2,l)>span];
		nQ=matrix(0,length(l),2);
		for(i in seq(along=l)) 
		{
			k=l[i]:r[i];
			nQ[i,2]=sum(Q[k,2]);
			nQ[i,1]=sum(Q[k,1]*Q[k,2])/nQ[i,2];
		}
		Q=nQ;
	}
	return(Q);
}

"mpair" <- function(Q) 
{
	Q=Q[order(Q[,1]),,drop=FALSE];
	i=nrow(Q);
	
	if(i>1) 
	{
		i=which.min(Q[-1,1]-Q[-i,1]);
		Q[i,1]=(Q[i,1]*Q[i,2]+Q[i+1,1]*Q[i+1,2])/(Q[i,2]+Q[i+1,2]);
		Q[i,2]=Q[i,2]+Q[i+1,2];
		Q=Q[-(i + 1),,drop=FALSE];
	}
	
	return(Q);
}

"removetsp" <- function(Q,epsilon=1e-6)
{
	Q=Q[Q[,2]>epsilon,,drop=FALSE];
	Q[,2]=Q[,2]/sum(Q[,2]);
	return(Q);
}
