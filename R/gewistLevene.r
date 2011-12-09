################################## GENE X ENVIRONMENT_FIXED EFFECT SIZES #################

# Outputs the optimal prioritization threshold for testing genetic interactions

	gewistLevene <- function(p, N, theta_gc, theta_c, M, K = 20000,verbose=FALSE){

## Input the known parameters (variance explained ~ (0,1))
## Assume Gene x Environment 

### check inputs

    
	if (p > 0.5 | p <= 0) stop("minor allele frequency should be a number between 0 and 0.5")
	
	if (!(N > 0 | M > 0 | K > 0 | theta_gc >0 | theta_c > 0)) stop( "negative input values are not allowed")
	
	if (!(theta_gc < 1 | theta_c < 1 )) stop(" variance explained should be a number between 0 and 1 ")

	N <- round(N)
  M <- round(M)
  K <- round(K)
  
### calculate beta coefficients
	
	b2 <- sqrt(( theta_c )/( 1 - theta_gc - theta_c ))

	b3 <- sqrt(theta_gc/(2*p*(1 - p)*( 1 - theta_gc - theta_c )))

### sample sizes 

	n1 <- round(N*(1 - p)^2)
	n2 <- round(N*(1 - p)*p*2)
	n3 <- N - n1 - n2

############ sample the variance of covariate and error per genotype

	var_C_G_1 <- rchisq(K,n1 - 1)/(n1 - 1)
	var_C_G_2 <- rchisq(K,n2 - 1)/(n2 - 1)
	var_C_G_3 <- rchisq(K,n3 - 1)/(n3 - 1)

	error_G_1 <- rchisq(K,n1 - 1)/(n1 - 1)
	error_G_2 <- rchisq(K,n2 - 1)/(n2 - 1)
	error_G_3 <- rchisq(K,n3 - 1)/(n3 - 1)

########### sample the b_x's 

	b_x1 <- rnorm(K,b2 + b3*(   - 2*p), sqrt(1/(var_C_G_1*(n1 - 1))))
	b_x2 <- rnorm(K,b2 + b3*(1 - 2*p),sqrt(1/(var_C_G_2*(n2 - 1))))
	b_x3 <- rnorm(K,b2 + b3*(2 - 2*p),sqrt(1/(var_C_G_3*(n3 - 1))))

########### the covariance of error term and covariate per genotype

	cov_error_C_1 <- (b_x1 - (b2 + b3*(   - 2*p)))*var_C_G_1
	cov_error_C_2 <- (b_x2 - (b2 + b3*(1 - 2*p)))*var_C_G_2
	cov_error_C_3 <- (b_x3 - (b2 + b3*(2 - 2*p)))*var_C_G_3	

	sum_of_C_square_G <- var_C_G_1*n1*( - 2*p) + var_C_G_2*n2*(1 - 2*p) + var_C_G_3*n3*(2 - 2*p)
	sum_of_C_G_square <- var_C_G_1*n1*( - 2*p)^2 + var_C_G_2*n2*(1 - 2*p)^2 + var_C_G_3*n3*(2 - 2*p)^2
	sum_C_2 <- var_C_G_1*n1 + var_C_G_2*n2 + var_C_G_3*n3

	total_cov_error <- ( - 2*p)*n1*cov_error_C_1 + (1 - 2*p)*n2*cov_error_C_2 + (2 - 2*p)*n3*cov_error_C_3
	total_cov	 <- n1*cov_error_C_1 + n2*cov_error_C_2 + n3*cov_error_C_3

############ observed b2s

	beta2 <- b2 + (sum_of_C_square_G*total_cov_error + total_cov*sum_of_C_G_square)/(sum_C_2*sum_of_C_G_square + sum_of_C_square_G^2)
	beta3 <- b3 + ((b2 - beta2)*sum_of_C_square_G + total_cov_error)/(sum_of_C_G_square)

########### Thus we can obtain the total variance per genotype

	var_group_1 <- abs(var_C_G_1*(b2 + b3*(   - 2*p))^2 + error_G_1 + 2*cov_error_C_1*(b2 + b3*(   - 2*p)))
	var_group_2 <- abs(var_C_G_2*(b2 + b3*(1 - 2*p))^2 + error_G_2 + 2*cov_error_C_2*(b2 + b3*(1 - 2*p)))
	var_group_3 <- abs(var_C_G_3*(b2 + b3*(2 - 2*p))^2 + error_G_3 + 2*cov_error_C_3*(b2 + b3*(2 - 2*p)))

	vaj_1 <- ((b2 + b3*(   - 2*p))^2 + 1)*((n1 - 1) - 2*(gamma(1/2)/beta(1/2,(n1 - 1)/2))^2)/(n1 - 1)
	vaj_2 <- ((b2 + b3*(1 - 2*p))^2 + 1)*((n2 - 1) - 2*(gamma(1/2)/beta(1/2,(n2 - 1)/2))^2)/(n2 - 1)
	vaj_3 <- ((b2 + b3*(2 - 2*p))^2 + 1)*((n3 - 1) - 2*(gamma(1/2)/beta(1/2,(n3 - 1)/2))^2)/(n3 - 1)

	sd_sam1 <- sqrt(abs((pi - 2)*(var_group_1)/(pi*(n1 - 1)) - 2/pi*vaj_1))
	sd_sam2 <- sqrt(abs((pi - 2)*(var_group_2)/(pi*(n2 - 1)) - 2/pi*vaj_2))
	sd_sam3 <- sqrt(abs((pi - 2)*(var_group_3)/(pi*(n3 - 1)) - 2/pi*vaj_3))

	z_value1 <- rnorm(K,sqrt(2*(var_group_1)/pi),sd_sam1)
	z_value2 <- rnorm(K,sqrt(2*(var_group_2)/pi),sd_sam2)
	z_value3 <- rnorm(K,sqrt(2*(var_group_3)/pi),sd_sam3)

	z_value <- (n1*z_value1 + n2*z_value2 + n3*z_value3)/N

      	denominator  <- 2*((n1 - 1)*var_group_1*(1 - 2/pi) + (n2 - 1)*var_group_2*(1 - 2/pi) + (n3 - 1)*var_group_3*(1 - 2/pi))
       	numerator   <- (N - 3)*((n1)*(z_value1 - z_value)^2 + (n2)*(z_value2 - z_value)^2 + (n3)*(z_value3 - z_value)^2)
    
    	levene_p <- 1 - pf(numerator/denominator,df1 = 2,df2 = N - 3)

############ Calculate interation p values

	RSS1 <- error_G_1 + (b2 + b3*(   - 2*p) - beta2 - beta3*(    -2*p))^2*var_C_G_1 + 2*(b2 + b3*(    -2*p) - beta2 - beta3*(   - 2*p))*cov_error_C_1
	RSS2 <- error_G_2 + (b2 + b3*(1 - 2*p) - beta2 - beta3*(1 - 2*p))^2*var_C_G_2 + 2*(b2 + b3*(1 - 2*p) - beta2 - beta3*(1 - 2*p))*cov_error_C_2
	RSS3 <- error_G_3 + (b2 + b3*(2 - 2*p) - beta2 - beta3*(2 - 2*p))^2*var_C_G_3 + 2*(b2 + b3*(2 - 2*p) - beta2 - beta3*(2 - 2*p))*cov_error_C_3

	f_stats <- 2*beta3^2*p*(1 - p)*(N - 4)^2/(RSS1*(n1 - 1) + RSS2*(n2 - 1) + RSS3*(n3 - 1))

	interaction_p <- 1 - pf(f_stats,df1 = 1,df2 = N - 4)

########## Calculate powers

		result <- data.frame("p-value_Cut-offs" = NA, "VP_power" = NA)

			for (i in 1:1000){
				power <- mean(levene_p < i/1000 & interaction_p<0.05/(M*i/1000))
				result[i,] <- c( i/1000,power)
				}


	conv_power <- power
	optimal_power <- max(result[,2])
	optimal_p_threshold <- which.max(result[,2])/1000


  				if (verbose){return(result)
				}else {
					
					return(list("Conventional_power"=conv_power,
					"Optimal_VP_power"=optimal_power,"Optimal_pval_threshold"=optimal_p_threshold))
					
					}
 				}


############## End of Script



