########################################################################
############## Given distribution of interaction #######################
############## Interaction effect_Distribution Estimates ###############

###theta_c,Mp are not in percentages ######

effectPDF<- function( distribution=c("beta","normal","uniform","weibull"),
parameter1,parameter2=NULL,parameter3=NULL,p,N,theta_c,M,K=20000,nb_incr=50,range=NULL,verbose=FALSE)
{
   					
#############################################################
## Input the known parameters (variance explained ~ (0,1))
## Assume Gene x Environment 

### check inputs

  if (p > 0.5 | p <= 0) stop("minor allele frequency should be a number between 0 and 0.5")
	
	if (!(N > 0 | M > 0 | K > 0 | theta_c > 0)) stop( "negative input values are not allowed")
	
	if (! theta_c < 1 ) stop(" variance explained should be a number between 0 and 1 ")

  if (length(range)!=2) stop("range of interaction effect size should be a vector with length two")  
  
  N <- round(N)
  M <- round(M)
  K <- round(K)
  
    variance_exp_interaction <- seq(range[1],range[2],length.out=nb_incr)

########## Start of the loop
	output <- list()

#source("GEWIST.r")

		for ( j in 1:nb_incr) {
		output[[j]] <- gewistLevene(p, N, variance_exp_interaction[j], theta_c, M, K = 20000,verbose = TRUE )[,2]
			}

########## Densities

	if (distribution=="normal") { densities <- dnorm(100*variance_exp_interaction,mean=parameter1,sd=parameter2)/sum(dnorm(100*variance_exp_interaction,mean=parameter1,sd=parameter2))} 
	if (distribution=="beta") { densities <- dbeta(100*variance_exp_interaction,shape1=parameter1,shape2=parameter2,ncp=parameter3)/sum(dbeta(100*variance_exp_interaction,shape1=parameter1,shape2=parameter2,ncp=parameter3)) }
	if (distribution=="uniform") { densities <- dunif(100*variance_exp_interaction,min=parameter1,max=parameter2)/sum(dunif(100*variance_exp_interaction,min=parameter1,max=parameter2)) }
	if (distribution=="weibull"){ densities <- dweibull(100*variance_exp_interaction,shape=parameter1,scale=parameter2) /sum(dweibull(100*variance_exp_interaction,shape=parameter1,scale = parameter2))}

		p_matrix<-matrix(unlist(output),1000,(nb_incr))

		power_distri<-data.frame("p_cutoffs"=NA,"VP_power"=NA)

				for (k in 1:1000){		
					power_distri[k,]<-c( k/1000,sum(p_matrix[k,]*densities))
				}

############################################################################################################

	if (distribution=="beta"){ distri="Beta Distribution"}
	if (distribution=="normal"){ distri="Normal Distribution"}
	if (distribution=="uniform"){ distri="Uniform Distribution"}
	if (distribution=="weibull"){distri="Weibull Distribution"}

distribution <- c("\n########### Interaction Testing For GenexEnvironment ############\n\n",distri,"\n\n")

		optimal_VP_power <- max(power_distri[,2])
		optimal_VP_thres  <- which.max(power_distri[,2])/1000
		conventional_pow <- power_distri[1000,2]
		maximal <- list("Conventional_power"=conventional_pow,
					"Optimal_VP_power"=optimal_VP_power,
					"Optimal_pval_threshold"=optimal_VP_thres)
      
      
      
		cat(distribution)
			if (verbose) return(power_distri)
				else return(maximal)
					}




