


LOGLINEAR <- function(data, data_type = 'raw', variables=NULL, Freq='Freq', verbose=TRUE) {


if (is.null(variables))  message('\n\nThe variables argument was not specified, the analyses cannot proceed.')

if (data_type == 'raw') 	{ # produce a contingency table and then convert it to dataframe with the counts
	tabel <- xtabs( noquote(paste('~', paste(variables, collapse = ' + '), sep=' ')), data = data)		
}	
	
if (data_type == 'counts') { 	
	tabel <- xtabs( noquote(paste(Freq, '~', paste(variables, collapse = ' + '), sep=' ')), data = data)
}

if (data_type == 'cont.table') {	
	# making sure the order of the variable levels is the same as for xtabs (alphabetical)
	tabel <- as.data.frame(data)
	tabel <- xtabs( noquote(paste('Freq', '~', paste(variables, collapse = ' + '), sep=' ')), data = tabel)
}
 
tabel_noms <- list(names(dimnames(tabel)))

Nvars <- length(variables)


formule <- as.formula(paste('Freq', paste(paste(variables,sep=""),collapse=" * "),sep=" ~ "))
termsALL <- attr(terms.formula(formule), "term.labels")
modtermsList <- vector("list", Nvars)
modtermsList[[1]] <- variables
for (lupe1 in 2:Nvars) {
	for (lupe2 in 1:length(termsALL)) {
		if ( grepl(":", termsALL[lupe2], fixed=FALSE) == TRUE & 
		     lengths(gregexpr(':', termsALL[lupe2])) == (lupe1 -1)) 
			modtermsList[[lupe1]] <- c(modtermsList[[lupe1]], termsALL[lupe2])
	}	
}
#modtermsList



# K - Way and Higher-Order Effects

output_KwayHO <- matrix(NA, (Nvars+1), 7)

# null model
mod_null <- glm(Freq ~ 1, data = tabel, family = poisson); summary(mod_null)
df_null <- mod_null$df.residual
X2_null <- sum(residuals(mod_null, type = "pearson")^2)
pX2_null <- 1 - pchisq(X2_null, mod_null$df.residual)
G2_null <- mod_null$deviance
pG2_null <- 1 - pchisq(G2_null, mod_null$df.residual)
AIC_null <- mod_null$aic
output_KwayHO[1,] <- cbind(1, df_null, G2_null, pG2_null, X2_null, pX2_null, AIC_null)

# saturated model
formule_satd <- as.formula(paste('Freq', paste(modtermsList[[Nvars]]), sep=" ~ "))
mod_satd <- glm(formule_satd, data = tabel, family = poisson); summary(mod_satd)
AIC_satd <- mod_satd$aic
output_KwayHO[nrow(output_KwayHO),] <- cbind(0, 0, 0, 1, 0, 1, AIC_satd)


output_Kway <- matrix(NA, (Nvars), 7)
prevmod <- mod_satd
termsRED <- termsALL
# for (lupe in Nvars:2) {
for (lupe in Nvars:1) {

	if (lupe > 1)  termsRED <- termsRED[ - which(termsRED %in% modtermsList[[lupe]])];   # print(termsRED)

	if (lupe == 1) termsRED <- modtermsList[[lupe]];   # print(termsRED)
		
	formule_red <- as.formula(paste('Freq', paste(termsRED, collapse=" + "), sep=" ~ "))
	# print(formule_red)	
	mod_red <- glm(formule_red, data = tabel, family = poisson); summary(mod_red)

	if (lupe > 1) {
		
		df_KwayHO <- mod_red$df.residual
		
		X2_KwayHO <- sum(residuals(mod_red, type = "pearson")^2)
		pX2_KwayHO <- 1 - pchisq(X2_KwayHO, mod_red$df.residual)
		
		G2_KwayHO <- mod_red$deviance
		pG2_KwayHO <- 1 - pchisq(G2_KwayHO, mod_red$df.residual)
		
		AIC_KwayHO <- mod_red$aic
		
		output_KwayHO[lupe,] <- cbind(lupe, df_KwayHO, G2_KwayHO, pG2_KwayHO, 
		                              X2_KwayHO, pX2_KwayHO, AIC_KwayHO) 
	}
 
	if (lupe == 1)  mod_red <- mod_null

	# if (lupe == 1)  { mod_red <- mod_null; prevX2 <- X2_null }
	
	ano_Kway <- anova(prevmod, mod_red)

	G2_Kway <- abs(ano_Kway$Deviance[2])

	df_Kway <- abs(ano_Kway$Df[2])
	
	pG2_Kway <- 1 - pchisq(G2_Kway, df_Kway)
	
	AIC_Kway <- mod_red$aic - prevmod$aic
	
	X2_Kway <- sum(residuals(mod_red, type = "pearson")^2)

	if (lupe != Nvars & lupe != 1)  X2_Kway <- X2_Kway - prevX2

	if (lupe == 1)  X2_Kway <- output_KwayHO[(lupe),5] - output_KwayHO[(lupe+1),5] 

	pX2_Kway <- 1 - pchisq(X2_Kway, df_Kway)

	output_Kway[lupe,] <- cbind(lupe, df_Kway, G2_Kway, pG2_Kway, X2_Kway, pX2_Kway, AIC_Kway) 

	prevmod <- mod_red

	prevX2 <- X2_Kway
}

dimnames(output_KwayHO) <- list(rep("", dim(output_KwayHO)[1]))
colnames(output_KwayHO) <- c('K','df','LR Chi-Square','p','Pearson Chi-Square','p','AIC')

dimnames(output_Kway) <- list(rep("", dim(output_Kway)[1]))
colnames(output_Kway) <- c('K','df','LR Chi-Square','p','Pearson Chi-Square','p','AIC diff.')





# Partial Associations

outputPartAssocn <- StepEffects <- c()

termsRED <- termsALL[-length(termsALL)] # remove the highest order term, already have the info

for (lupe1 in (Nvars-1):1) {
 
	termsGENCLASS <-  modtermsList[[lupe1]]
	formuleGENCLASS <- as.formula(paste('Freq', paste(termsGENCLASS, collapse=" + "), sep=" ~ "))
	modGENCLASS <- glm(formuleGENCLASS, data = tabel, family = poisson); summary(modGENCLASS)

	StepEffects <- rbind(StepEffects, data.frame( Effect = NA) )

	outputPartAssocn <- rbind(outputPartAssocn, cbind(NA,NA,NA,NA))

	for (lupe2 in 1:length(modtermsList[[lupe1]])) {
		
		termsTEMP <- termsRED[ - which(termsRED %in% modtermsList[[lupe1]][lupe2])]

		formule_red <- as.formula(paste('Freq', paste(termsTEMP, collapse=" + "), sep=" ~ "))
		
		mod_red <- glm(formule_red, data = tabel, family = poisson); summary(mod_red)

		ano_part <- anova(modGENCLASS, mod_red)

		G2_part <- abs(ano_part$Deviance[2])

		df_part <- abs(ano_part$Df[2])
		
		pG2_part <- 1 - pchisq(G2_part, df_part)
		
		AIC_part <- mod_red$aic - modGENCLASS$aic
		
		outputPartAssocn <- rbind(outputPartAssocn, cbind(G2_part, df_part, pG2_part, AIC_part))

		StepEffects <- rbind(StepEffects, data.frame( Effect = modtermsList[[lupe1]][lupe2] ) )
	}
	
	termsRED <- termsRED[ - which(termsRED %in% modtermsList[[lupe1]])]
}

dimnames(outputPartAssocn) <- list(rep("", dim(outputPartAssocn)[1]))
colnames(outputPartAssocn) <- c('LR Chi-Square','df','p','AIC diff.')
outputPartAssocn <- round_boc(outputPartAssocn)

PartialAssociationsTab <- data.frame(StepEffects, data.frame(outputPartAssocn) , row.names=NULL)
PartialAssociationsTab[is.na(PartialAssociationsTab)] <- " "  # set NAs to blanks





# Parameter Estimates

# For saturated models, .500 has been added to all observed cells
# SPSS "Model Selection" (not "General") Parameter Estimates

# # using the stats::loglin function
# paramests <- loglin( (tabel + .5), margin=tabel_noms, start = rep(1, length(tabel)), fit = FALSE,  
                     # eps = 0.1, iter = 20, param = TRUE, print = FALSE)$param

# using glm   https://rstudio-pubs-static.s3.amazonaws.com/84177_4604ecc1bae246c9926865db53b6cc29.html
formule <- as.formula(paste('Freq', paste(paste(termsALL,sep=""),collapse=" + "),sep=" ~ "))

glmcommand <- paste(paste(modtermsList[[1]],sep=""), collapse="=contr.sum, ")
glmcommand <- noquote(paste('glm(formule,  family=poisson(link="log"), contrasts=list(', 
                             glmcommand, "= contr.sum), data=(tabel + .5))"))

glm_mod <- suppressWarnings(eval(parse(text = glmcommand)))
paramests <- coef(summary(glm_mod));  paramests
CI_ub <- paramests[,1] + 1.96 * paramests[,2]
CI_lb <- paramests[,1] - 1.96 * paramests[,2]
#CIs <- suppressMessages(suppressWarnings(confint(glm_mod)))  # defaults to MASS estimates & warnings
paramests <- cbind(paramests, CI_lb, CI_ub)






# Backward Elimination

keeper <- function(tokeep, todrop) {
	
	maineffs <- unique(unlist(strsplit(todrop[1],split = ":")))
	
	Xnlevelm1 <- length(maineffs) - 1
	
	lowers <- maineffs
	
	if (Xnlevelm1 > 1) {
			
		for (lupe in 2:length(Xnlevelm1))  
			lowers <- c(lowers, combn(maineffs, Xnlevelm1, FUN=paste, collapse=':', simplify=TRUE) )			
	}
	
	keepthese <- unique(c(tokeep, lowers))

	keepthese <- keepthese[ - which(keepthese %in% todrop)]

	return(keepthese)
}

# tokeep <- termsRED

# todrop <- termsRED[14]

# keepcela <- keeper(tokeep = tokeep, todrop = todrop)



keepallthese <- function(tokeep) {
	
	keepthese <- c()
	
	for (lupe in 1:length(tokeep)) {
		
		term <- tokeep[lupe]
	
		maineffs <- unique(unlist(strsplit(term[1],split = ":")))
		
		Xnlevelm1 <- length(maineffs) - 1
		
		lowers <- maineffs
		
		if (Xnlevelm1 > 1) {
				
			for (lupe in 2:length(Xnlevelm1))  
				lowers <- c(lowers, combn(maineffs, Xnlevelm1, FUN=paste, collapse=':', simplify=TRUE) )			
		}
		
		keepthese <- c(keepthese, unique(c(tokeep, lowers)) )

	}
	
	return(unique(keepthese))
}

# xx = c('gender:ethnic:cathelp', 'gender:income:cathelp')

# keepallthese(xx)



AddRowToDF <- function(dfname, newvalues, Nprior = 0, Npost = 0) {
	
	# add a row of NAs, into which the newvalues will be inserted, along with Nprior rows of NAs
	dfname2 <- dfname
	for (newrows in 1:(Nprior+1))	dfname2 <- rbind(dfname2, rep(NA, ncol(dfname)))
	
	# loop through newvalues, adding them to the correct columns in dfname
	for (lupe in 1:length(newvalues))  dfname2[names(newvalues[lupe])][nrow(dfname2),1] = newvalues[lupe]
	
	# add Npost rows with NAs
	if (Npost > 0) for (newrows in 1:Npost)	dfname2 <- rbind(dfname2, rep(NA, ncol(dfname)))
	
	return(dfname2)
}



Step_Summary <- data.frame(Step = NA, GenDel = NA, Effects = NA, LR_Chi_Square = NA, df = NA, p = NA, AIC = NA)

Step_num <- 0

Step_Summary <- AddRowToDF(dfname = Step_Summary, 
                           newvalues = list(Step = Step_num, GenDel = 'Generating Class', Effects = termsALL[length(termsALL)],  
                                            LR_Chi_Square = 0, df = 0, p = 1, AIC = output_KwayHO[nrow(output_KwayHO),7]), Nprior=0 )

Step_Summary <- AddRowToDF(dfname = Step_Summary, 
                           newvalues = list(GenDel = 'Deleted Effect', Effects = termsALL[length(termsALL)],
                                            LR_Chi_Square = output_KwayHO[(nrow(output_KwayHO)-1),3], 
                                            df = output_KwayHO[(nrow(output_KwayHO)-1),2], 
                                            p = output_KwayHO[(nrow(output_KwayHO)-1),4], 
                                            AIC = output_KwayHO[(nrow(output_KwayHO)-1),7] ), Nprior=1 )

Step_num <- Step_num + 1				
	
modtermsList_2 <- modtermsList

termsRED <- termsALL

# if removing the highest term does not reduce fit, continue with the Backward Elimination
if (output_KwayHO[(nrow(output_KwayHO)-1),4] > .05) {

	termsRED <- termsALL[-length(termsALL)] # remove the highest order term
	
	mustkeep <- c()

	for (lupe1 in (Nvars-1):1) {
	 
		previous_mod <- mod_satd  # for future comparison

		yesno <- 'yes'
		
		pvalues <- c()
	
		testset <- modtermsList_2[[lupe1]]

		testset <- setdiff(testset, mustkeep)

		while (yesno == 'yes' & length(testset) > 0) {

			terms_GENCLASS <- termsRED     # c(mustkeep, testset)    # modtermsTEMP     # modtermsList_2[[lupe1]]
			
			formule_GENCLASS <- as.formula(paste('Freq', paste(terms_GENCLASS, collapse=" + "), sep=" ~ "))
			
			mod_GENCLASS <- glm(formule_GENCLASS, data = tabel, family = poisson); summary(mod_GENCLASS)
			
			ano_GENCLASS <- anova(mod_satd, mod_GENCLASS)
	
			previous_mod <- mod_GENCLASS  # for future comparison
			
			G2_GENCLASS <- abs(ano_GENCLASS$Deviance[2])
	
			df_GENCLASS <- abs(ano_GENCLASS$Df[2])
	
			pG2_GENCLASS <- 1 - pchisq(G2_GENCLASS, df_GENCLASS)
			
			AIC_GENCLASS <- mod_GENCLASS$aic - AIC_satd
			
			X2_GENCLASS <- sum(residuals(mod_GENCLASS, type = "pearson")^2)
			
			pX2_GENCLASS <- 1 - pchisq(X2_GENCLASS, mod_GENCLASS$df.residual)
	
			output_GENCLASS <- cbind(G2_GENCLASS, df_GENCLASS, pG2_GENCLASS, AIC_GENCLASS)
			

			# Step_Summary <- AddRowToDF(dfname = Step_Summary, newvalues = c(GenDel = 'Generating Class', Effects = paste('All ',lupe1,'-ways',sep='')))
	
			Step_Summary <- AddRowToDF(dfname = Step_Summary, 
			                           newvalues = list(Step = Step_num, GenDel = 'Generating Class', Effects = "All of these terms:", LR_Chi_Square = G2_GENCLASS, 
			                                            df = df_GENCLASS, p = pG2_GENCLASS, AIC = AIC_GENCLASS), Nprior=2, Npost=0 )
			                                            			                                            
			if (length(terms_GENCLASS) > 1)	{
				
				for (lupe in 1:length(terms_GENCLASS))
				
				Step_Summary <- AddRowToDF(dfname = Step_Summary, newvalues = list(Effects = terms_GENCLASS[lupe]) )
				
				if ( lupe == length(terms_GENCLASS)) Step_Summary <- AddRowToDF(dfname = Step_Summary, newvalues = NA )
			}	
	

			testset <- setdiff(testset, mustkeep)
			
			for (lupe2 in 1:length(testset)) {
	
				termsTEMP <- termsRED[ - which(termsRED %in% testset[lupe2])]
		
				formule_red <- as.formula(paste('Freq', paste(termsTEMP, collapse=" + "), sep=" ~ "))
				
				mod_red <- glm(formule_red, data = tabel, family = poisson); summary(mod_red)
		
				ano_DEL <- anova(mod_GENCLASS, mod_red)
		
				G2_DEL <- abs(ano_DEL$Deviance[2])
		
				df_DEL <- abs(ano_DEL$Df[2])
	
				pG2_DEL <- 1 - pchisq(G2_DEL, df_DEL)
				
				AIC_DEL <- mod_red$aic - mod_GENCLASS$aic
				
				X2_DEL <- sum(residuals(mod_red, type = "pearson")^2)
				
				pX2_DEL <- 1 - pchisq(X2_DEL, mod_red$df.residual)
	
				Step_Summary <- AddRowToDF(dfname = Step_Summary, 
				                          newvalues = list(GenDel = 'Deleted Effect Test', Effects = testset[lupe2],
				                                           LR_Chi_Square = G2_DEL, df = df_DEL, p = pG2_DEL, AIC = AIC_DEL) )
				
				pvalues <- c(pvalues, pG2_DEL)  # save the p value
				
				# at the end of lupe 2, remove the term with the highest p value if it is  > .05  & reset lupe2
	
				if (lupe2 == length(testset)) {
					
					Step_num <- Step_num + 1				
								
					if (max(pvalues) > .05) {
											
						delthis <- which(pvalues == max(pvalues))
												
						todrop <- testset[delthis]    # modtermsTEMP[delthis]
						
						termsRED <- keeper(tokeep = termsRED, todrop = todrop)
																		
						testset <- testset[-delthis]
						
						if (length(testset) == 0)  yesno <- 'no'

						Step_Summary <- AddRowToDF(dfname = Step_Summary, 
						                           newvalues = list(GenDel = 'Deleted On This Step', Effects = todrop), Nprior=1, Npost=0 )
				
						pvalues <- c()
						
					}
					if (!is.null(pvalues) & (all(pvalues < .05) | length(testset)==1)) {
						
						mustkeep <- c(mustkeep, keepallthese(testset) )
						
						yesno <- 'no'
						
						Step_Summary <- AddRowToDF(dfname = Step_Summary, 
						                           newvalues = list(GenDel = 'Deleted On This Step', Effects = 'none deleted'), Nprior=1, Npost=0 )
					}
										
			}  # if (lupe2 == length(testset))
			
		} # lupe2
	
	} # while
		
}  # lupe1

} # if







# Final Model

formule_final <- as.formula(paste('Freq', paste(termsRED, collapse=" + "), sep=" ~ "))

mod_final <- glm(formule_final, data = tabel, family = poisson); summary(mod_final)

df_final <- mod_final$df.residual

X2_final <- sum(residuals(mod_final, type = "pearson")^2)
pX2_final <- 1 - pchisq(X2_final, df_final)

G2_final <- mod_final$deviance
pG2_final <- 1 - pchisq(G2_final, df_final)

AIC_final <- mod_final$aic

output_final <- cbind(df_final, G2_final, pG2_final, X2_final, pX2_final, AIC_final) 
dimnames(output_final) <- list(rep("", dim(output_final)[1]))
colnames(output_final) <- c('df','LR Chi-Square','p','Pearson Chi-Square','p','AIC')

# glm parameter estimates
mod_final_coefs <- summary(mod_final)$coefficients

# Cell Counts and Residuals
expfreqs <- fitted(mod_final)
resids <- mod_final$y - expfreqs
stdresids <- resids  / sqrt(expfreqs)  # Field p 825
adjresids <- resid(mod_final) # ??
tabel2 <- cbind(mod_final$y, expfreqs, resids, stdresids, adjresids)
colnames(tabel2) <- c('Obsd. Freq.','Exp. Freq.','Residuals','Std. Resid.','Adjusted Resid.')
tabel3 <- cbind(as.data.frame(tabel), tabel2)
# confirm that the Freq (from data2) & Obsd. Freq. (from tabel2) columns are identical, just in case
if (sum(tabel3[,'Freq'] - tabel3[,'Obsd. Freq.']) == 0)  { tabel3 <- subset(tabel3, select = -c(Freq) )  
	} else { message('\n\nProblem: the freqs in tabel2 are not = to those in tabel3') }
# place the rows in the same order as SPSS
tabel4 <- tabel3
for (lupe in Nvars:1)  tabel4 <- tabel4[ order(eval(parse(text = noquote(paste('tabel4$', variables[lupe], sep='')) ))), ]





if (verbose) {
	
	message('\n\nThe input data:\n')
	print(tabel, print.gap=4)
	
	message('\n\nK - Way and Higher-Order Effects\n')
	print(round_boc(output_KwayHO), print.gap=4)
	
	message('
	These are tests that K - Way and Higher-Order Effects are zero, i.e., tests
	of the hypothesis that Kth-order and higher interactions are zero
	If these effects and all higher order effects are removed from the model,
	then here are the consequences.
	
	The df values indicate the number of effects (model terms) that are removed.
	
	The first row, labeled as 1, shows the consequences of removing all of the main
	effects and all higher order effects (i.e., everything) from the model. This
	usually results in poor fit. A statistically significant chi-square indicates	
	that the prediction of the cell frequencies is significantly worse than the 
	prediction that is provided by the saturated model. It would suggest that at
	least one of the removed effects needs to be included in the model.
	
	The second row, labeled as 2, shows the consequences of removing all of the
	two-way and higher order effects from the model, while keeping the main effects.
	A statistically significant chi-square indicates a reduction in prediction success
	compared to the saturated model and that at least one of the removed effects needs
	to be included in the model.
	
	The same interpretation process applies if there is a K = 3 row, and so on.
	A K = 3 row in the table would show the consequences of removing all of the
	three-way and higher order effects from the model, while keeping the two-way
	interactions and main effects.
	
	A nonsignificant chi-square for a row would indicate that removing the
	model term(s) does not significantly worsen the prediction of the cell
	frequencies and the term(s) is nonessential and can be dropped from the model.
	
	The bottom row in the table, labeled as 0, is for the saturated mode. It
	includes all possible model terms and therefore provides perfect prediction
	of the cell frequencies. The AIC values for this model can be helpful in
	gaging the relative fit of models with fewer terms.')
	
	
	message('\n\nK-Way Effects\n')
	print(round_boc(output_Kway), print.gap=4)
	
	message('
	These are tests that the K - Way Effects are zero, i.e., tests whether
	interactions of a particular order are zero. The tests are for model
	comparisons/differences. For each K-way test, a model is fit with and then
	without the interactions and the change/difference in chi-square and
	likelihood ratio chi-square values are computed.
		
	For example, the K = 1 test is for the comparison of the model with
	all main effects and the intercept with the model with only the intercept.
	A statistically significant K = 1 test is (conventionally) considered to
	mean that the main effects are not zero and that they are needed in the model.
		
	The K = 2 test is for the comparison of the model with all two-way
	interactions, all main effects, and the intercept with the model with
	the main effects, and the intercept. A statistically significant K = 2 test
	is (conventionally) considered to mean that the two-way interactions are
	not zero and that they are needed in the model.
		
	The K = 3 test (if there is one) is for the comparison of the model
	with all three-way interactions, all two-way interactions, all main
	effects, and the intercept with the model with all two-way interactions,
	all main effects, and the intercept. A statistically significant K = 3 test
	is (conventionally) considered to mean that the three-way interactions
	are not zero and that they are needed in the model, and so on.
		
	The df values for the model comparisons are the df values associated
	with the K-way terms.
		
	The above "K - Way and Higher-Order Effects" and "K - Way" tests are for the
	ncollective importance of the effects at each value of K. There are not tests
	nof individual terms. For example, a significant K = 2 test means that the set
	nof two-way terms is important, but it does not mean that every two-way term is
	significant.')
	
	# \nFor example, the K = 1 test is for the comparison of the model with all main effects and
	# the intercept with the model with only the intercept. The K = 2 test is for the comparison
	# of the model with all two-way interactions, all main effects, and the intercept with the model
	# with the main effects, and the intercept. The K = 3 test (if there is one) is for the
	# comparison of the model with all three-way interactions, all two-way interactions, all main
	# effects, and the intercept with the model with all two-way interactions, all main effects,
	# and the intercept, and so on.')
	
	
	# # Interpreting the K-Way and Higher-Order Effects table in SPSS Statistics Model Selection (HILOGLINEAR) output
	
	# https://www.ibm.com/support/pages/interpreting-k-way-and-higher-order-effects-table-spss-statistics-model-selection-hiloglinear-output
	
	# The K-way and Higher Order Effects section of the K-Way and Higher-Order Effects table 
	# provides tests of all effects at level K and higher, while the K-way Effects section 
	# provides tests of just the K-way effects. For example, suppose you have four variables 
	# specified. The line for K=4 is the same for each section of the table because there 
	# is only one 4-way effect. In the bottom part of the table, the line for K=3 provides 
	# a test that all of the 3-way terms are null, while in the upper half of the table, 
	# the line for K=3 combines the test for 3-way terms with the one for the 4-way term, 
	# testing that the 3-way and all higher order effects are null. For K=2, the bottom 
	# half shows the test for the 2-way terms, while the top half shows the combined test 
	# for 2 through 4-way effects. Finally, the line for K=1 in the lower section shows 
	# the test for main effects, while the line in the top section shows the test for 
	# all 4 levels of effects combined (which is everything other than an intercept). 
	# Note that the df in the top half of the table for each line is the sum of the df in 
	# the lower half for the same K and higher values of K (which are lower in each section of the table).
	
	message('\n\nPartial Associations:\n')
	print(PartialAssociationsTab, print.gap=4)
	
	message('
	These are tests of individual terms in the model, with the restriction that
	higher-order terms at each step are excluded. The tests are for differences
	between models. For example, the tests of 2-way interactions are for the
	differences between the model with all 2-way interactions (and no higher-order
	interactions) and the model when each individual 2-way interaction is removed in turn.')

	
	message('\n\n\nParameter Estimates (SPSS "Model Selection", not "General", Parameter Estimates):

	For saturated models, .500 has been added to all observed cells:\n')
	colnames(paramests)[4] <- 'p'
	print(round_boc(paramests), print.gap=4)
	# print(lapply(paramests,round,3))   # for when param ests are from the loglin function


	message('\n\nBackward Elimination Statistics:\n')
	Step_Summary_show <- round_boc(Step_Summary)
	Step_Summary_show[is.na(Step_Summary_show)] <- " "  # set NAs to blanks
	print(round_boc(Step_Summary_show), print.gap=4, row.names = FALSE)
	
	message('
	The hierarchical backward elimination procedure begins with all possible
	terms in the model and then removes, one at a time, terms that do not
	satisfy the criteria for remaining in the model.
	A term is dropped only when it is determined that removing the term does
	not result in a reduction in model fit AND if the term is not involved in any
	higher order interaction. On each Step above, the focus is on the term that results
	in the least-significant change in the likelihood ratio chi-squre if removed.
	If the change is not significant, then the term is removed.')
	
	message('\n\n\nThe Final Model Formula:\n')
	message(noquote(paste('Freq', paste(termsRED, collapse=" + "), sep=" ~ ")))
	
	message('\n\nThe Final Model Goodness-of-Fit Tests:\n')
	print(round_boc(output_final), print.gap=4)
	
	message('\n\nGeneralized Linear Model Coefficients for the Final Model:\n')
	print(round_boc(mod_final_coefs), print.gap=4)
	
	message('\n\nCell Counts and Residuals:\n')
	print(round_boc(tabel4), print.gap=4)
}


Output <- list(KwayHO=output_KwayHO, Kway=output_Kway, PartialAssociations=PartialAssociationsTab,  
               paramests=paramests, Step_Summary=Step_Summary, FinalModeltests=output_final, 
               FinalModelcells=tabel4)

return(invisible(Output))

}

