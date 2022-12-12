\name{LOGLINEAR}
\alias{LOGLINEAR}
\title{Hierarchical loglinear analyses for two or more categorical variables}
\description{Provides 'SPSS'- and 'SAS'-like output for hierarchical loglinear analyses 
	of two or more categorical variables. The input can be raw data, a contingency
	table, or a dataframe with cell frequency counts. The output includes:
	(1) a table with the K-Way and higher-order effects;
	(2) a table with the K-Way effects;
	(3) a table with the the partial associations;
	(4) a table with the parameter estimates; 			
	(5) a table with the backward elimination statistics;
	(6) a table with the final model goodness of fit tests; and
	(7) a table with the final model observed and expected frequencies,
	standardized residuals, and adjusted residuals.
}
\usage{LOGLINEAR(data, data_type = 'raw', variables=NULL, Freq = 'Freq', verbose=TRUE) 
}
\arguments{
  \item{data}{
  \code{}The input data, which can be raw data or a dataframe 
  with cell frequency counts (see the Examples below).}

  \item{data_type}{
  \code{}The kind of input data. The options are 'raw' (for raw data), 'cont.table' 
  (for a two-dimensional contingency table), or 'counts' (for a dataframe with 
  the cell frequency counts).}

  \item{variables}{
  \code{}The variable names. Two or more variable names must be specified, 
  as in, variables=c('varA','varB', 'varC').}
   	    
  \item{Freq}{
  \code{}(optional) If data_type = 'counts', then Freq is the name of the column in data
  that has the frequency counts. If unspecified, it will be assumed that the column 
  is named 'Freq'.}
   	    
  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE.}
}
\details{The purpose of hierarchical loglinear procedures is to find a model that best 
	fits data given the model-fitting constraints, and to then provide the model 
	parameters. The analyses begin with the saturated model, which includes all 
	possible terms and for which there is a perfect fit to data. Terms are then 
	tested for possible exclusion, which occurs when removal of a term does not 
	result in a statistically significant reduction in fit \emph{and} when a term is 
	not involved in any higher order interactions. This function provides 
	statistics for the saturated model, for the hierarchal removal of the model 
	terms, for the backward elimination steps, and for the final model.	

	When data_type = 'cont.table', the data must be a two-dimensional contingency 
	table that has the names of the table dimensions/variables. See the Examples below.  
}
\value{A list with the following possible elements:
\item{KwayHO}{A table with the K-Way and higher-order effects.}
\item{Kway}{A table with the K-Way effects.}
\item{PartialAssociations}{A table with the partial associations.}
\item{paramests}{A table with the parameter estimates.}
\item{StepSummTab}{A table with the backward elimination statistics.}
\item{FinalModeltests}{A table with the final model goodness of fit tests.}
\item{FinalModelcells}{A table with the final model observed and expected frequencies and
	adjusted residuals.}
}
\references{
	{Agresti, A. (2013). \emph{Categorical data analysis (3rd ed).} Hobokon, NJ: John Wiley & Sons.}
\cr\cr {Ajzen, R., & Walker, C. M. (2021). \emph{Categorical data analysis for the behavioral and 
	    social sciences (2nd ed.).} New York, NY: Routledge.}
\cr\cr {Field, A. (2018). Chapter 18: Categorical data. 
		\emph{Discovering statistics using SPSS (5th ed.).} Los Angeles, CA: Sage.}
\cr\cr {Noursis, M. J. (2012). Chapter 1: Model selection loglinear analysis.
	\emph{IBM SPSS statistics 19: Advanced statistical 
	procedures Companion.} Upper Saddle River, NJ: Prentice Hall.}
\cr\cr {Nussbaum, E. M. (2015). \emph{Categorical and nonparametric data analysis 
		choosing the best statistical technique.} New York, NY: Routledge.}
\cr\cr {Stevens, J. P. (2009). Chapter 14: Categorical data analysis: The log linear model. 
	\emph{Applied multivariate statistics for the social sciences (5th ed.).} 
	New York, NY: Routledge.}
\cr\cr {Tabachnick, B. G., & Fidell, L. S. (2019). Chapter 16: Multiway 
	frequency analysis. \emph{Using multivariate statistics.} New York, NY: Pearson.}
\cr\cr {von Eye, A., & Mun, E. Y. (2013). \emph{Log-Linear modeling concepts, 
		interpretation, and application.} Hoboken, NJ: Wiley.}
}
\author{Brian P. O'Connor}
\examples{
# Field (2018). Chapter 19: Categorical data -- cats & dogs, entering raw data
LOGLINEAR(data = datasets$Field_2018, 
          data_type = 'counts', 
          variables=c('Animal', 'Training', 'Dance'), 
          Freq = 'Freq' )

# Field (2018). Chapter 19: Categorical data -- cats & dogs, entering raw counts	
LOGLINEAR(data = datasets$Field_2018_raw,  
          data_type = 'raw', 
          variables=c('Animal', 'Training', 'Dance'), 
          Freq = NULL )
          
# Field (2018). Chapter 19: Categorical data -- cats & dogs, entering a table
# example of creating and entering a two-dimensional contingency table for 'data'
food <- c(28, 10)
affection <- c(48, 114)
Field_2018_cats_conTable <- as.table(rbind(food, affection)) 
colnames(Field_2018_cats_conTable) <- c('danced', 'did not dance')
names(attributes(Field_2018_cats_conTable)$dimnames) <- c('Training','Dance') 
LOGLINEAR(data = Field_2018_cats_conTable, 
          data_type = 'cont.table', 
          variables=c('Training', 'Dance') )


# go to this web page to see many more examples of the LOGLINEAR function analyses:
# https://oconnor-psych.ok.ubc.ca/loglinear/LOGLINEAR_vignettes.html

}
