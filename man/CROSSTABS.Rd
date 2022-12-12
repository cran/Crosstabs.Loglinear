\name{CROSSTABS}
\alias{CROSSTABS}
\title{Cross tabulations of two categorical variables}
\description{Provides 'SPSS'- and 'SAS'-like output for cross tabulations of two
	categorical variables. The input can be raw data, a contingency
	table, or a dataframe with cell frequency counts. The output includes
	the contingency table, expected frequencies, Pearson's' Chi-Square, Yates's Chi-Square
	(continuity correction), the Likelihood Ratio, Fisher's Exact p, 
	the Linear-by-Linear Association.', the McNemar Test, the contingency coefficient C, phi,
	Cramer's V, Cohen's W, the residuals, standardized residuals, and adjusted 
	residuals. Additional output for 2-by-2 tables includes the
	risk difference, the risk ratio, the odds ratio, and Yule's Q.}
\usage{CROSSTABS(data, data_type = 'raw', variables=NULL, Freq = NULL, verbose=TRUE) 
}
\arguments{
  \item{data}{
  \code{}The input data, which can be raw data, a contingency table, or a dataframe 
  with cell frequency counts (see the Examples below).}

  \item{data_type}{
  \code{}The kind of input data. The options are 'raw' (for raw data), 'cont.table' 
  (for a two-dimensional contingency table), or 'counts' (for a dataframe with the cell frequency counts).}

  \item{variables}{
  \code{}(optional) The two variable names, which is required if data_type = 'raw' or 'counts', 
  e.g., variables=c('varA','varB'). Not required if data_type = 'cont.table'.}
   	    
  \item{Freq}{
  \code{}(optional) If data_type = 'counts', then Freq is the name of the column in data
  that has the frequency counts. If unspecified, it will be assumed that the column 
  is named 'Freq'.}
   	    
  \item{verbose}{
  \code{}(optional) Should detailed results be displayed in console? \cr The options are: 
  TRUE (default) or FALSE.}
}
% \details{This function }
\value{A list with the following possible elements:
\item{obsFreqs}{The observed frequencies.}
\item{expFreqs}{The expected frequencies.}
\item{modEStab}{Model test and effect size coefficients.}
\item{residuals}{The residuals.}
\item{stdresiduals}{The standardized residuals.}
\item{adjresiduals}{The adjusted residuals.}
\item{EStab2x2}{For a 2-by-2 contingency table, a list with the risk difference, 
	the risk ratio, the odds ratio, and Yule's Q values.}
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
\cr\cr {Stevens, J. P. (2009). Chapter 14: Categorical data analysis: The log linear model. 
	\emph{Applied multivariate statistics for the social sciences (5th ed.).} 
	New York, NY: Routledge.}
\cr\cr {Tabachnick, B. G., & Fidell, L. S. (2019). Chapter 16: Multiway 
	frequency analysis. \emph{Using multivariate statistics.} New York, NY: Pearson.}
}
\author{Brian P. O'Connor}
\examples{
# when 'data' is a raw data file (rather than counts/frequencies)
# Field (2018). Chapter 18: Categorical data -- cats only
CROSSTABS(data = subset(datasets$Field_2018_raw, Animal=='Cat'), 
          data_type = 'raw', 
          variables=c('Training','Dance') )


# when 'data' is a file with the counts/frequencies (rather than raw data points)
# Field (2018). Chapter 18: Categorical data -- cats only
CROSSTABS(data = subset(datasets$Field_2018, Animal=='Cat'), 
          data_type = 'counts', 
          variables=c('Training','Dance') )


# create and enter a two-dimensional contingency table for 'data'
# Field (2018). Chapter 18: Categorical data -- cats only
food <- c(28, 10)
affection <- c(48, 114)
Field_2018_cats_conTable <- rbind(food, affection) 
colnames(Field_2018_cats_conTable) <- c('danced', 'did not dance')
names(attributes(Field_2018_cats_conTable)$dimnames) <- c('Training','Dance') 
CROSSTABS(data = Field_2018_cats_conTable, data_type = 'cont.table')


# another way of creating the same two-dimensional contingency table for 'data'
# Field (2018). Chapter 18: Categorical data -- cats only
Field_2018_cats_conTable_2 <- matrix( c(28, 48, 10, 114), nrow = 2, ncol = 2)
colnames(Field_2018_cats_conTable_2) <- c('danced', 'did not dance')
rownames(Field_2018_cats_conTable_2) <- c('food', 'affection')
CROSSTABS(data = Field_2018_cats_conTable_2, data_type = 'cont.table')


# go to this web page to see many more examples of the CROSSTABS function analyses:
# https://oconnor-psych.ok.ubc.ca/loglinear/CROSSTABS_vignettes.html


}
