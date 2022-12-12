



CROSSTABS <- function(data, data_type = 'raw', variables=NULL, Freq = NULL, verbose=TRUE) {


	if (is.null(variables) & (data_type == 'raw' | data_type == 'counts'))
		message('\n\nThe variables argument was not specified for this data_type, the analyses cannot proceed.')
	
	
	if (data_type == 'raw') {
		
		if (length(variables) != 2)  message('\n\nTwo, and only two, variables must be specified for the analyses to proceed.')
	
		tabel <- table(data[,variables[1]], data[,variables[2]])  # A will be rows, B will be columns 
	
	}
	
	if (data_type == 'cont.table')  tabel <- data 
	
	
	if (data_type == 'counts') {   # counts to contingency table
	
		if (length(variables) != 2)  message('\n\nTwo, and only two, variables must be specified for the analyses to proceed.')
	
		tabel <- xtabs(Freq ~ data[,variables[1]] + data[,variables[2]], data=data)
		
		names(dimnames(tabel)) <- variables
	}


# tabel <- cbind(tabel[,2], tabel[,1]);    colnames(tabel) <- c('Yes', 'No');   tabel

	# Pearson chi-square	
	Pearson <- chisq.test(tabel, correct=FALSE)
	X2_Pearson <- Pearson$statistic
	df <- Pearson$parameter
	pX2_Pearson  <- Pearson$p.value

	# Yates' corrected chi-square (Continuity Correction)
	Yates <- chisq.test(tabel, correct=TRUE)
	X2_Yates <- Yates$statistic
	pX2_Yates  <- Yates$p.value
	
	obsFreqs <- tabel
	expFreqs <- Pearson$expected                                                                                                                                                                                                      

	Nexplt5 <- sum(expFreqs < 5)  # cells with expected count < 5
	Pexplt5 <- (Nexplt5 / (nrow(expFreqs) * ncol(expFreqs))) * 100  # percentage

	# Likelihood Ratio   http://myweb.facstaff.wwu.edu/minerb2/biometrics/chi_conti_test.html
	G2 <- 2*sum(obsFreqs*log(obsFreqs/expFreqs))   
	pG2 <- 1 - pchisq(G2, 1)
	N <- sum(tabel)	

	# Fisher's exact test
	FisherExact <- fisher.test(tabel)$p.value
	
	# Linear-by-Linear Association test (Agresti, 2002, p 87)
	# also known as the Mantel-Haenszel Chi-Square Test    https://v8doc.sas.com/sashtml/stat/chap28/sect19.htm
	# phi -- SPSS    2012 Norusis 10 Crosstabulation p 173-174   phi can be > 1   SPSS adds a negative sign "if necessary"
	# phi <- sqrt( X2_Pearson / N) # not good for > 2 x 2 tables
	
	combins <- expand.grid(1:nrow(tabel), 1:ncol(tabel), stringsAsFactors = FALSE)  # adpated from DescTools - MHChisqTest
	# https://stackoverflow.com/questions/44120875/weighted-correlation-in-r
	weighted_corr <- cov.wt(combins, wt = as.vector(as.matrix(tabel)), cor = TRUE) 
	r_contin <- weighted_corr$cor[2,1]
	X2_MH <- (N - 1) * r_contin**2
	p_X2_MH <- 1 - pchisq(X2_MH, 1)
		
	
	X2tab <- rbind(cbind(X2_Pearson, df, pX2_Pearson), cbind(X2_Yates, df, pX2_Yates), cbind(G2, df, pG2),
	               cbind(NA,NA,FisherExact), cbind(X2_MH, df, p_X2_MH) )
	colnames(X2tab) <- c(' ','df','p')
	rownames(X2tab) <- c('Pearson Chi-Square','Yates Continuity Correction','Likelihood Ratio',
	                     'Fisher\'s Exact p','Linear-by-Linear Association')

	if (nrow(tabel) == 2 & ncol(tabel) == 2) {
		# McNemar's Chi-squared test with continuity correction  -- only for 2x2 tables
		McN <- mcnemar.test(as.matrix(tabel))
		X2_McN <- McN$statistic
		p_X2_McN <- McN$p.value

		X2tab <- rbind(X2tab, cbind(X2_McN, df, p_X2_McN))
		rownames(X2tab) <- c('Pearson Chi-Square','Yates Continuity Correction','Likelihood Ratio',
		                     'Fisher\'s Exact p','Linear-by-Linear Assocn.','McNemar Test')
	}




	# In Sheskin's Handbook of Parametric and Nonparametric Statistical Procedures, for each categorical 
	# association coefficient (e.g., Pearson's C, Pearson's phi, Cramer's V), the Author also takes into 
	# account Cohen's effect size measure. For example, for Cramer's V, Cohen's measure of 
	# effect size (namely, w) can be calculated as: w=Cramer's V * squareroot(k-1)

	# https://www.statisticshowto.com/contingency-coefficient/
	# The contingency coefficient (C) is a coefficient of association that tells whether two variables or data sets 
	# are independent or dependent of each other. It is also known as Pearsons Coefficient
	C <- sqrt((X2_Pearson / (N + X2_Pearson)))
	
	# phi -- SPSS    2012 Norusis 10 Crosstabulation p 173-174   phi can be > 1   SPSS adds a negative sign "if necessary"
	phi <- sqrt( X2_Pearson / N)

	# Cramer's V     https://www.statology.org/cramers-v-in-r/
	dmin_m1 <- min( (nrow(tabel)-1), (ncol(tabel)-1) )
	CramerV_X2 <- sqrt( (X2_Pearson / N) / dmin_m1)   # using X2
	CramerV_G2 <- sqrt( (G2 / N) / dmin_m1)   # using G2
	
	# Cohen's W from Cramer's V     https://www.spss-tutorials.com/effect-size/
	dmin <- min( (nrow(tabel)), (ncol(tabel)) )
	CohenW_X2 <- CramerV_X2 * sqrt(dmin - 1)
	CohenW_G2 <- CramerV_G2 * sqrt(dmin - 1)

	# Cohen's W from C   https://www.spss-tutorials.com/effect-size/
	# CohenW <- sqrt( C**2 / (1 - C**2))

	modEStab <- rbind(C, phi, CramerV_X2, CramerV_G2, CohenW_X2, CohenW_G2)
	colnames(modEStab) <- list(rep("", dim(modEStab)[2]))
	rownames(modEStab) <- c('Contingency coefficient C','Phi','Cramer\'s V (from X2)','Cramer\'s V (from G2)',
	                        'Cohen\'s W (from X2)','Cohen\'s W (from G2)')
	

	Output <- list(obsFreqs=obsFreqs, expFreqs=expFreqs, modEStab=modEStab, residuals=(obsFreqs - expFreqs),
	               stdresiduals=Pearson$residuals, adjresiduals=Pearson$stdres) 



	# # for the SPSS Directional Measures
	
	# v good = DescTools, and the chisquare packages
	
	# Nominal by Nominal  Lambda	
	# https://people.ohio.edu/ruhil/statsbookR/props.html
	# derived from the chisquare function in the chisquare package   Goodman-Kruskal's lambda		
	# E1 <- N - max(rowSums(tabel))
	# E2 <- sum(colSums(tabel) - apply(tabel, 2, max))
	# lambda.rowvar.dependent <- (E1 - E2) / E1	
	# E1 <- N - max(colSums(tabel))
	# E2 <- sum(rowSums(tabel) - apply(tabel, 1, max))
	# lambda.colvar.dependent <- (E1 - E2) / E1



	# effect sizes for 2x2 tables	
	if (nrow(tabel) == 2 & ncol(tabel) == 2) {
		EStab2x2 <- ES2x2table(obsFreqs)  # boc function
		Output <- c(Output, EStab2x2)
	}


	if (verbose == TRUE) {
			
		message('\n\nObserved Frequencies:\n'); print(tabel)
		
		message('\n\nRow Totals:\n');    print(rowSums(tabel))
	
		message('\n\nColumn Totals:\n'); print(colSums(tabel))
		
		message('\n\nN = ', N)
	
		message('\n\nExpected Frequencies:\n'); print(round_boc(expFreqs))
	
		message('\n\nNumber of cells with expected count < 5 = ', Nexplt5, ' (', Pexplt5,'%)')
	
		message('\n\nChi-Square Tests:\n'); print(round_boc(X2tab), print.gap=4)
	
		message('\n\nModel Effect Sizes:'); print(round_boc(modEStab), print.gap=4)
	
		message('\n\nResiduals:\n') #  (the Pearson residuals, (observed - expected) / sqrt(expected))
		print(round_boc(obsFreqs - expFreqs))
	
		message('\n\nStandardized Residuals:\n') #  (the Pearson residuals, (observed - expected) / sqrt(expected))
		print(round_boc(Pearson$residuals))
	
		message('\n\nAdjusted Residuals:\n') #  ((observed - expected) / sqrt(V), where V is the residual cell variance)
		print(round_boc(Pearson$stdres))

		if (nrow(tabel) == 2 & ncol(tabel) == 2) {
			message('\n\n2-by-2 Table Effect Sizes\n')
			print(lapply(EStab2x2,round,3))
		}
	}

return(invisible(Output))

}



# Interpreting Yules Q      https://www.statisticshowto.com/gamma-coefficient-goodman-kruskal/

# Yules Q is always a number between -1 and 1.

# Q is 0: no association between the variables.

# Q = 0 to  0.29: a negligible or very small association.

# Q = -0.30 to -0.49 or 0.30 to 0.49: a moderate association between the variables.

#		message('\n\nYule\'s Q (the number of pairs in agreement - the number in disagreement) / the total number of paired observations) = ', round(YulesQ,3) )




# # from SPSS:


# Phi is a chi-square-based measure of association that involves
# dividing the chi-square statistic by the sample size and taking the
# square root of the result

# Lambda. A measure of association that reflects the proportional
# reduction in error when values of the independent variable are used to
# predict values of the dependent variable. A value of 1 means that the
# independent variable perfectly predicts the dependent variable. A
# value of 0 means that the independent variable is no help in
# predicting the dependent variable.

# Uncertainty coefficient. A measure of association that indicates the
# proportional reduction in error when values of one variable are used
# to predict values of the other variable. For example, a value of 0.83
# indicates that knowledge of one variable reduces error in predicting
# values of the other variable by 83%. The program calculates both
# symmetric and asymmetric versions of the uncertainty coefficient.




# Correlations. For tables in which both rows and columns contain
# ordered values, Correlations yields Spearman's correlation
# coefficient, rho (numeric data only). Spearman's rho is a measure of
# association between rank orders. When both table variables (factors)
# are quantitative, Correlations yields the Pearson correlation
# coefficient, r, a measure of linear association between the variables.



# Nominal. 

# For nominal data (no intrinsic order, such as Catholic,
# Protestant, and Jewish), you can select Contingency coefficient, Phi
# (coefficient) and Cramr's V, Lambda (symmetric and asymmetric lambdas
# and Goodman and Kruskal's tau), and Uncertainty coefficient.

# Contingency coefficient. A measure of association based on chi-square.
# The value ranges between 0 and 1, with 0 indicating no association
# between the row and column variables and values close to 1 indicating
# a high degree of association between the variables. The maximum value
# possible depends on the number of rows and columns in a table. Phi and
# Cramer's V. Phi is a chi-square-based measure of association that
# involves dividing the chi-square statistic by the sample size and
# taking the square root of the result. Cramer's V is a measure of
# association based on chi-square. Lambda. A measure of association that
# reflects the proportional reduction in error when values of the
# independent variable are used to predict values of the dependent
# variable. A value of 1 means that the independent variable perfectly
# predicts the dependent variable. A value of 0 means that the
# independent variable is no help in predicting the dependent variable.
# Uncertainty coefficient. A measure of association that indicates the
# proportional reduction in error when values of one variable are used
# to predict values of the other variable. For example, a value of 0.83
# indicates that knowledge of one variable reduces error in predicting
# values of the other variable by 83%. The program calculates both
# symmetric and asymmetric versions of the uncertainty coefficient.




# Ordinal. 

# For tables in which both rows and columns contain ordered
# values, select Gamma (zero-order for 2-way tables and conditional for
# 3-way to 10-way tables), Kendall's tau-b, and Kendall's tau-c. For
# predicting column categories from row categories, select Somers' d.

# Gamma. A symmetric measure of association between two ordinal
# variables that ranges between -1 and 1. Values close to an absolute
# value of 1 indicate a strong relationship between the two variables.
# Values close to 0 indicate little or no relationship. For 2-way
# tables, zero-order gammas are displayed. For 3-way to n-way tables,
# conditional gammas are displayed. Somers' d. A measure of association
# between two ordinal variables that ranges from -1 to 1. Values close
# to an absolute value of 1 indicate a strong relationship between the
# two variables, and values close to 0 indicate little or no
# relationship between the variables. Somers' d is an asymmetric
# extension of gamma that differs only in the inclusion of the number of
# pairs not tied on the independent variable. A symmetric version of
# this statistic is also calculated. Kendall's tau-b. A nonparametric
# measure of correlation for ordinal or ranked variables that take ties
# into account. The sign of the coefficient indicates the direction of
# the relationship, and its absolute value indicates the strength, with
# larger absolute values indicating stronger relationships. Possible
# values range from -1 to 1, but a value of -1 or +1 can be obtained
# only from square tables. Kendall's tau-c. A nonparametric measure of
# association for ordinal variables that ignores ties. The sign of the
# coefficient indicates the direction of the relationship, and its
# absolute value indicates the strength, with larger absolute values
# indicating stronger relationships. Possible values range from -1 to 1,
# but a value of -1 or +1 can be obtained only from square tables.




# Nominal by Interval. 

# When one variable is categorical and the other is
# quantitative, select Eta. The categorical variable must be coded
# numerically.

# Eta. A measure of association that ranges from 0 to 1, with 0
# indicating no association between the row and column variables and
# values close to 1 indicating a high degree of association. Eta is
# appropriate for a dependent variable measured on an interval scale
# (for example, income) and an independent variable with a limited
# number of categories (for example, gender). Two eta values are
# computed: one treats the row variable as the interval variable, and
# the other treats the column variable as the interval variable. Kappa.
# Cohen's kappa measures the agreement between the evaluations of two
# raters when both are rating the same object. A value of 1 indicates
# perfect agreement. A value of 0 indicates that agreement is no better
# than chance. Kappa is based on a square table in which row and column
# values represent the same scale. Any cell that has observed values for
# one variable but not the other is assigned a count of 0. Kappa is not
# computed if the data storage type (string or numeric) is not the same
# for the two variables. For string variable, both variables must have
# the same defined length.

# Risk. For 2 x 2 tables, a measure of the strength of the association
# between the presence of a factor and the occurrence of an event. If
# the confidence interval for the statistic includes a value of 1, you
# cannot assume that the factor is associated with the event. The odds
# ratio can be used as an estimate or relative risk when the occurrence
# of the factor is rare.

# McNemar. A nonparametric test for two related dichotomous variables.
# Tests for changes in responses using the chi-square distribution.
# Useful for detecting changes in responses due to experimental
# intervention in "before-and-after" designs. For larger square tables,
# the McNemar-Bowker test of symmetry is reported.

# Cochran's and Mantel-Haenszel statistics. Cochran's and
# Mantel-Haenszel statistics can be used to test for independence
# between a dichotomous factor variable and a dichotomous response
# variable, conditional upon covariate patterns defined by one or more
# layer (control) variables. Note that while other statistics are
# computed layer by layer, the Cochran's and Mantel-Haenszel statistics
# are computed once for all layers.



# 2011 Howell - Fundamental Statistics for the Behavioral Sciences (7th ed) p 526
# Why do we have both odds and risk?
# Why do we have to complicate things by having both odds ratios and risk ratios? 
# That is a very good question, and it has some good answers. Risk is something 
# that I think most of us understand. When we say the risk of having a heart attack 
# in the No Aspirin condition is .0171, we are saying that 1.71% of the participants 
# in that condition had a heart attack, and that is pretty straightfor- ward. When 
# we say that the odds of a heart attack in that condition are .0174, we are saying 
# that the odds of having a heart attack are 1.74% of the odds of not having a heart 
# attack. That may be a popular way of setting bets on race horses, but it leaves me 
# dissatisfied. So why have an odds ratio in the first place?
# The odds ratio can be calculated in situations in which a true risk ratio cannot be. 
# In a retrospective study, where we find a group of people who had heart attacks and 
# another group of people who did not have heart attacks and we look back to see if 
# they took aspirin, we cant really calculate risk. Risk is future oriented. If we 
# give 1,000 people aspirin and withhold it from 1,000 oth- ers, we can look at these 
# people 10 years down the road and calculate the risk (and risk ratio) of heart attacks. 
# But if we take 1,000 people with (and without) heart attacks and look backward, we 
# cant really calculate risk because we have sampled heart attack patients at far 
# greater than their normal rate in the popu- lation (50% of our sample has had a heart 
# attack, but certainly 50% of the pop- ulation does not suffer from heart attacks). 
# But we can always calculate odds ratios. And, when we are talking about low probability 
# events, such as having a heart attack, the odds ratio is usually a very good estimate 
# of what the risk ratio would be if we could calculate it. The odds ratio is equally 
# valid for prospective and retrospective sampling designs. That is important.




	# message('\nPearson\'s Chi-squared = ', round(X2,3), '   df = ', df, '    p = ', round(pX2,5))

	# message('\nG2 = ', round(G2,3), '   df = ', df, '    p = ', round(pG2,5))

	# message('\nContingency coefficient C = ', round(C,3))

	# message('\nCramer\'s V (from X2) = ', round(CramerV_X2,5))

	# message('\nCramer\'s V (from G2) = ', round(CramerV_G2,5))

	# message('\nCohen\'s W (from X2) = ', round(CohenW_X22,5))

	# message('\nCohen\'s W (from G2) = ', round(CohenW_G2,5))
		


