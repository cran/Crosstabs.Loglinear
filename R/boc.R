



# effect sizes for 2x2 tables
# https://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Basic-statistics/TwoByTwo.html

ES2x2table <- function(obsFreqs) {

	a <- a_orig <- obsFreqs[1,1]
	b <- b_orig <- obsFreqs[1,2]
	c <- c_orig <- obsFreqs[2,1]
	d <- d_orig <- obsFreqs[2,2]

	outext <- c()
	
	for (lupe1 in 1:4) {
	
		p1 <- a / (a+b) 
		p2 <- c / (c+d) 
		
		# Risk difference
		Risk_diff <- p1 - p2
		Risk_diff_SE <- sqrt(p1*(1-p1)/(a+b)+p2*(1-p2)/(c+d))
		Risk_diff_lb <- Risk_diff - 1.96 * Risk_diff_SE
		Risk_diff_ub <- Risk_diff + 1.96 * Risk_diff_SE
				
		# Risk ratio
		Risk_ratio = p1 / p2
		Risk_ratio_SE <- sqrt((1-p1)/a+(1-p2)/c)
		Risk_ratio_lb <- Risk_ratio * exp(- 1.96 * Risk_ratio_SE)
		Risk_ratio_ub <- Risk_ratio * exp(  1.96 * Risk_ratio_SE)
		
		# Odds ratio 1 v 2
		Odds_ratio <- (p1/(1-p1)) / (p2/(1-p2))
		Odds_ratio_SE <- sqrt((1/a+1/b+1/c+1/d))
		Odds_ratio_lb <- Odds_ratio * exp(- 1.96 * Odds_ratio_SE)
		Odds_ratio_ub <- Odds_ratio * exp(  1.96 * Odds_ratio_SE)
		
		# # Odds ratio 2 v 1
		# Odds_ratio21 <- (p2/(1-p2)) / (p1/(1-p1))
		# Odds_ratio_lb21 <- Odds_ratio21 * exp(- 1.96 * Odds_ratio_SE)
		# Odds_ratio_ub21 <- Odds_ratio21 * exp(  1.96 * Odds_ratio_SE)
		
		# Yule's Q     http://personality-project.org/r/html/Yule.html
		YulesQ <- (Odds_ratio-1) / (Odds_ratio+1)    #   same as  (a*d - b*c)/(a*d+b*c)
		YulesQ_lb <- (Odds_ratio_lb-1) / (Odds_ratio_lb+1)    
		YulesQ_ub <- (Odds_ratio_ub-1) / (Odds_ratio_ub+1)    
		
		ES_2x2_tab <- matrix(NA, 4, 3)
		ES_2x2_tab[1,] <- cbind(Risk_diff, Risk_diff_lb, Risk_diff_ub)
		ES_2x2_tab[2,] <- cbind(Risk_ratio, Risk_ratio_lb, Risk_ratio_ub)
		ES_2x2_tab[3,] <- cbind(Odds_ratio, Odds_ratio_lb, Odds_ratio_ub)
		# ES_2x2_tab[4,] <- cbind(Odds_ratio21, Odds_ratio_lb21, Odds_ratio_ub21),
		ES_2x2_tab[4,] <- cbind(YulesQ, YulesQ_lb, YulesQ_ub)
				
		colnames(ES_2x2_tab) <- c(' ','CI_lb','CI_ub')
		rownames(ES_2x2_tab) <- c('Risk difference','Risk ratio','Odds ratio','Yule\'s Q')
	
	
		if (lupe1 == 1) {ES_2x2_tab_1 <- ES_2x2_tab  
			             outext <- rbind(outext, paste('Effect Sizes for ', rownames(obsFreqs)[1], ' vs. ', 
			                             rownames(obsFreqs)[2], ' (',  names(dimnames(obsFreqs))[1],    
			                             ')  when ', names(dimnames(obsFreqs))[2], ' = ', colnames(obsFreqs)[1], sep='') )
			             a = c_orig; b = d_orig; c = a_orig; d = b_orig 
		}
		if (lupe1 == 2) {ES_2x2_tab_2 <- ES_2x2_tab
			             outext <- rbind(outext, paste('Effect Sizes for ', rownames(obsFreqs)[2], ' vs. ', 
			                             rownames(obsFreqs)[1], ' (',  names(dimnames(obsFreqs))[1],
			                             ')  when ', names(dimnames(obsFreqs))[2], ' = ', colnames(obsFreqs)[1], sep='') )
			             a = b_orig; b = a_orig; c = d_orig; d = c_orig
		}	
		if (lupe1 == 3) {ES_2x2_tab_3 <- ES_2x2_tab
			             outext <- rbind(outext, paste('Effect Sizes for ', rownames(obsFreqs)[1], ' vs. ', 
			                             rownames(obsFreqs)[2], ' (',  names(dimnames(obsFreqs))[1],
			                             ')  when ', names(dimnames(obsFreqs))[2], ' = ', colnames(obsFreqs)[2], sep='') )
			             a = d_orig; b = c_orig; c = b_orig; d = a_orig		}
		if (lupe1 == 4) {ES_2x2_tab_4 <- ES_2x2_tab
						 outext <- rbind(outext, paste('Effect Sizes for ', rownames(obsFreqs)[2], ' vs. ', 
			                             rownames(obsFreqs)[1], ' (',  names(dimnames(obsFreqs))[1],
			                             ')  when ', names(dimnames(obsFreqs))[2], ' = ', colnames(obsFreqs)[2], sep='') )		
		}
	}

	output <- list(ES_2x2_tab_1=ES_2x2_tab_1, ES_2x2_tab_2=ES_2x2_tab_2, ES_2x2_tab_3=ES_2x2_tab_3, ES_2x2_tab_4=ES_2x2_tab_4)
	
	names(output) <- unlist(outext)
	
	return(invisible(output))
}


# print(ES2x2table(obsFreqs))

# obsFreqsR <- cbind(obsFreqs[,2], obsFreqs[,1]); colnames(obsFreqsR) <- c(colnames(obsFreqs)[2], colnames(obsFreqs)[1])

# print(ES2x2table(obsFreqsR))









# rounds numeric columns in a matrix
# numeric columns named 'p' or 'plevel' or 'plevel.adj' are rounded to round_p places
# numeric columns not named 'p' are rounded to round_non_p places

round_boc <- function(donnes, round_non_p = 3, round_p = 5) {
	
	# identify the numeric columns
	#	numers <- apply(donnes, 2, is.numeric)  # does not work consistently 
	for (lupec in 1:ncol(donnes)) {

		if (is.numeric(donnes[,lupec]) == TRUE) 
		
			if (colnames(donnes)[lupec] == 'p' | colnames(donnes)[lupec] == 'plevel'  | 
			    colnames(donnes)[lupec] == 'p adj.')  {
				donnes[,lupec] <- round(donnes[,lupec],round_p)
			} else {
				donnes[,lupec] <- round(donnes[,lupec],round_non_p)				
			}		
		# if (is.numeric(donnes[,lupec]) == FALSE) numers[lupec] = 'FALSE'		
		# if (colnames(donnes)[lupec] == 'p') numers[lupec] = 'FALSE'		
	}
	
	# # set the p column to FALSE
	# numers_not_p <- !names(numers) %in% "p"
	
#	donnes[,numers_not_p] = round(donnes[,numers_not_p],round_non_p) 
	
#	if (any(colnames(donnes) == 'p'))  donnes[,'p'] = round(donnes[,'p'],round_p) 

	return(invisible(donnes))
}






squareTable <- function(var1, var2, grpnames, tabdimnames=c('variable 1','variable 2')) {
    var1fact <- factor(var1, labels = grpnames)
    var2fact <- factor(var2, labels = grpnames)
    table(var1fact, var2fact, dnn=tabdimnames)
}






kappa.cohen <- function (var1, var2, grpnames) {

	# the data for this function (kapdon) are the category values for each of 2 columns,
	# but the analyses are performed on the contingency table	

	kapdonCT <- squareTable(var1, var2, grpnames); kapdonCT

	# based on Valiquette 1994 BRM, Table 1
	n <- sum(kapdonCT)  # Sum of Matrix elements
	kapdonP <- kapdonCT / n  # table of proportions
	po <- sum(diag(kapdonP))
	c <- rowSums(kapdonP)
	r <- colSums(kapdonP)
	pe <- r %*% c
	num <- po - pe
	den <- 1 - pe
	kappa <- num / den

	# SE and variance of kappa from Cohen 1968
	sek <- sqrt((po*(1-po))/(n*(1-pe)^2))
		
	if (n < 100) var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	
	if (n >= 100) {
		s <- t(matrix(colSums(kapdonP))) # columns sum
		var <- (pe+pe^2-sum(diag(r%*%s) * (r+t(s)))) / (n*(1-pe)^2)
		# asymptotic kappa variance as reported by 
		# Fleiss, J. L., Lee, J. C. M., & Landis, J. R. (1979).
	    # The large sample variance of kappa in the case of different sets of raters. 
	    # Psychological Bulletin, 86, 974-977
	}

	zkappa <- kappa / sqrt(    abs(var)     )
	
	sig <- round(pnorm(abs(zkappa),lower.tail = FALSE),5) * 2 # 2-tailed test

#	print( c(kappa, sek, var, zkappa, sig))

	return(invisible(c(kappa, zkappa, sig)))

	# # based on kappa.m	
	# n <- sum(kapdon) # Sum of Matrix elements	
	# kapdonP <- kapdon/n # proportion		
	# r <- matrix(rowSums(kapdonP)) # rows sum
	# s <- t(matrix(colSums(kapdonP))) # columns sum	
	# Ex <- (r) %*% s # expected proportion for random agree
	# f <- diag(1,3)	
	# # pom <- apply(rbind(t(r),s),2,min)	
	# po <- sum(sum(kapdonP * f))  # sum(sum(x.*f))
	# pe <- sum(sum(Ex * f))
	# k <- (po-pe)/(1-pe)
	# # km <- (pom-pe)/(1-pe) # maximum possible kappa, given the observed marginal frequencies
	# # ratio <- k/km # observed as proportion of maximum possible	
	# # kappa standard error for confidence interval as reported by Cohen in his original work
	# sek <- sqrt((po*(1-po))/(n*(1-pe)^2))	
	# var <- pe/(n*(1-pe))  # kappa variance as reported by Cohen in his original work
	# # var <- (pe+pe^2-sum(diag(r*s).*(r+s')))/(n*(1-pe)^2)  # for N > 100	
	# zk  <-  k/sqrt(var)

}



# Fleiss's kappa

# source: https://en.wikipedia.org/wiki/Fleiss%27_kappa

# Fleiss's kappa is a generalisation of Scott's pi statistic, a
# statistical measure of inter-rater reliability. It is also related to
# Cohen's kappa statistic. Whereas Scott's pi and Cohen's kappa work for
# only two raters, Fleiss's kappa works for any number of raters giving
# categorical ratings (see nominal data), to a fixed number of items. It
# can be interpreted as expressing the extent to which the observed amount
# of agreement among raters exceeds what would be expected if all raters
# made their ratings completely randomly. Agreement can be thought of as
# follows, if a fixed number of people assign numerical ratings to a number
# of items then the kappa will give a measure for how consistent the
# ratings are. The scoring range is between 0 and 1. 

# Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328. 
# Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bul- letin, 76, 378-382. 
# Fleiss, J.L., Levin, B., & Paik, M.C. (2003). Statistical Methods for Rates and Proportions, 3rd Edition. New York: John Wiley & Sons. 

kappa.fleiss <- function(var1, var2) {

	kapdon <- cbind(var1, var2)
	
	# the data for this function (kapdon) are the category values for each column,
	# but the analyses are performed on a count matrix (not a contin table) = the fleissmat below	
	fleissmat <- matrix(0,nrow(kapdon),max(kapdon))
	for (luper in 1:nrow(kapdon)) {
		for (lupec in 1:ncol(kapdon)) {
			fleissmat[luper,kapdon[luper,lupec]] <- fleissmat[luper,kapdon[luper,lupec]] + 1				
		}
	}	
	n <- nrow(fleissmat) 
	m <- sum(fleissmat[1,]) 	
	a <- n * m		
	pj <- colSums(fleissmat) / a 
	b <- pj * (1-pj)
	c <- a*(m-1)
	d <- sum(b)
	kj <- 1-(colSums((fleissmat * (m-fleissmat))) / (c %*% b)) # the value of kappa for the j-th category
	# sekj <- sqrt(2/c) 
	# zkj <- kj / sekj
	# pkj <- round(pnorm(abs(zkj),lower.tail = FALSE),5) * 2  # 2-tailed test
	k <- sum(b*kj, na.rm=TRUE) / d  # Fleiss's (overall) kappa
	sek <- sqrt(2*(d^2-sum(b *(1-2 *pj))))/sum(b *sqrt(c)) 	
	#ci <- k+(c(-1,1) * (pnorm(zkj)) * sek) 	
	zk <- k / sek  # normalized kappa
	sig <- round(pnorm(abs(zk),lower.tail = FALSE),5) * 2  # 2-tailed test
	
#	print( c(k, sek, zk, sig))
	
	return(invisible(c(k, zk, sig)))
}



kappas <- function(var1, var2, grpnames) {
	
	kc <- kappa.cohen(var1, var2, grpnames)  
	kf <- kappa.fleiss(var1, var2)
	kappasOUT <- rbind(kc,kf)
	kappasOUT[,1:2] <- round(kappasOUT[,1:2],3)
	kappasOUT[,3] <- round(kappasOUT[,3],5)	
	rownames(kappasOUT) <- c( "Cohen's kappa", "Fleiss's kappa")
	colnames(kappasOUT) <- c( "    kappa", "        z", "         p" )
		
	# k2 <- irr::kappa2( grpdat )  # Unweighted Kappa for categorical data without a logical order
	# kf <- irr::kappam.fleiss( grpdat, exact = FALSE, detail = TRUE)
	# kl <- irr::kappam.light( grpdat )
	# kappasOUT <- matrix( -9999, 3, 4)
	# kappasOUT[1,] <- cbind( k2$subjects, k2$value, k2$statistic, k2$p.value)
	# kappasOUT[2,] <- cbind( kl$subjects, kl$value, kl$statistic, kl$p.value)
	# kappasOUT[3,] <- cbind( kf$subjects, kf$value, kf$statistic, kf$p.value)
	# rownames(kappasOUT) <- c( "Cohen's kappa", "Light's kappa", "Fleiss's kappa")
	# colnames(kappasOUT) <- c( "   N", "    kappa", "         z", "      p" )

	return (kappasOUT)
}






# Press' Q 

# When DFA is used to classify
# individuals in the second or holdout sample. The percentage of cases that are
# correctly classified reflects the degree to which the samples yield consistent
# information. The question, then is what proportion of cases should be correctly
# classified? This issue is more complex than many researchers acknowledge. To
# illustrate this complexity, suppose that 75% of individuals are Christian, 15%
# are Muslim, and 10% are Sikhs. Even without any information, you could thus
# correctly classify 75% of all individuals by simply designating them all as
# Christian. In other words, the percentage of correctly classified cases should
# exceed 75%.

# Nonetheless, a percentage of 76% does not indicate that classification is
# significantly better than chance. To establish this form of significance, you
# should invoke Press' Q statistic. 

# Compute the critical value, which equals the chi-square value at 1 degree of
# freedom. You should probably let alpha equal 0.05. When Q exceeds this critical
# value, classification can be regarded as significantly better than chance,
# thereby supporting cross-validation.

# The researcher can use Press's Q statistic to compare with the chi-square critical 
# value of 6.63 with 1 degree of freedom (p < .01). If Q exceeds this critical value, 
# the classification can be regarded as significantly better than chance. 

# Press'sQ = (N _ (n*K))^2 / (N * (K - 1))

# where
# N = total sample size
# n = number of observations correctly classified 
# K = number of groups_
# Given that Press's Q = 29.57 > 6.63, it can be concluded that the classification 
# results exceed the classification accuracy expected by chance at a statistically 
# significant level (p < .01). 


PressQ <- function(freqs) {
	
	N <- sum(colSums(freqs))
	n <- sum(diag(freqs))
	k <- ncol(freqs)
	PressQ <- (N - (n*k))^2 / (N * (k - 1))
	
	return(invisible(PressQ))
}







   
