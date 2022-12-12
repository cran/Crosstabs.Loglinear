\name{datasets}
\alias{datasets}
\docType{data}
\title{datasets}
\description{A list with example data that were used in textbook
presentations of categorical data analyses
}
\usage{data(datasets)}
\details{
A list with example data that were used in the following textbook
presentations of categorical data analyses: 

\strong{datasets$Agresti_2019_Tab9.3} is tabled data from Agresti (2019, p. 346).

\strong{datasets$Agresti_2019_Tab9.8} is tabled data from Agresti (2019, p. 351).

\strong{datasets$Ajzen_2021_Tab7.11} is tabled data from Ajzen and Walker (2021, p. 178).

\strong{datasets$Ajzen_2021_Tab7.16} is tabled data from Ajzen (2021, p. 180).

\strong{datasets$Field_2018} is tabled data from Field (2018, Output 18.5 and Output 18.6).

\strong{datasets$Field_2018_raw} is raw data that simulates those from Field (2018, Output 18.5 and Output 18.6).

\strong{datasets$George_2019_26_Hierarchical} is tabled data from George (2019, pp. 346-347).

\strong{datasets$Gray_2012_2wqy} is tabled data from Gray and Kinnear (2012, p. 538).

\strong{datasets$Gray_2012_3wqy} is tabled data from Gray and Kinnear (2012, p. 551).

\strong{datasets$Green_2014} is tabled data from Green and Salkind (2014, p. 334).

\strong{datasets$Ho_2014} is tabled data from Ho (2014, p. 513).

\strong{datasets$Howell_2013} is tabled data from Howell (2013, p. 150).

\strong{datasets$Howell_2017} is tabled data from Howell (2019, p. 512).

\strong{datasets$Meyers_2013} is tabled data from Meyers (2013, p. 693).

\strong{datasets$Noursis_2012_marital} is tabled data from Noursis (2012a, p. 3).

\strong{datasets$Noursis_2012_voting_degree} is tabled data from Noursis (2012b, p. 513).

\strong{datasets$Noursis_2012_voting_degree_sex} is tabled data from Noursis (2012b, p. 527).

\strong{datasets$Stevens_2009_HeadStart_1} is tabled data from Stevens (2009, p. 472).

\strong{datasets$Stevens_2009_HeadStart_2} is tabled data from Stevens (2009, p. 474).

\strong{datasets$Stevens_2009_Inf_Survival} is tabled data from Stevens (2009, p. 481).

\strong{datasets$TabFid_2019_small} is tabled data from Tabachnick and Fidell (2019, p. 677).

\strong{datasets$Warner_2020_titanic} is tabled data from Warner (2020, p. 525).

\strong{datasets$Warner_2020_dog} is tabled data from Warner (2020, p. 530).
}

\references{
	{Agresti, A. (2013). \emph{Categorical data analysis (3rd ed).} Hobokon, NJ: John Wiley & Sons.}
\cr\cr {Ajzen, R., & Walker, C. M. (2021). \emph{Categorical data analysis for the behavioral and 
	    social sciences (2nd ed.).} New York, NY: Routledge.}
\cr\cr {Field, A. (2018). Chapter 18: Categorical data. 
		\emph{Discovering statistics using SPSS (5th ed.).} Los Angeles, CA: Sage.}
\cr\cr {George, D., & Mallery, P. (2019). Chapter 26: Hierarchical log-linear models. 
	\emph{IBM SPSS statistics for Windows, version 25.} IBM Corp., Armonk, N.Y., USA.}
\cr\cr {Gray, C. D., & Kinnear, P. R. (2012). Chapter 14: The analysis of multiway frequency tables.
	\emph{IBM SPSS statistics 19 made simple.} Psychology Press.}
\cr\cr {Green, S. B., Salkind, N. J. (2014). Chapter 41: Two-way contingency table analysis.
	\emph{Using SPSS for Windows and Macintosh: Analyzing and understanding data.} New York, NY: Pearson.}
\cr\cr {Ho, R. (2014). Chapter 19: Nonparametric tests.
	\emph{Handbook of univariate and multivariate data analysis with 
	IBM SPSS.} Boca Raton, FL: CRC Press.}
\cr\cr {Howell, D. C. (2013). Chapter 6: Categorical data and chi-square.
	\emph{Statistical methods for psychology (8th ed.).} Belmont, CA: Wadsworth Cengage Learning.}
\cr\cr {Howell, D. C. (2017). Chapter 19: Chi-square. 
	\emph{Fundamental statistics for the behavioral sciences} Belmont, CA: Wadsworth Cengage Learning.}
\cr\cr {Meyers, L. S., Gamst, G. C., & Guarino, A. J. (2013). Chapter 66: 
	Hierarchical loglinear analysis. \emph{Performing data analysis using IBM SPSS.}
	Hoboken, NJ: Wiley.}
\cr\cr {Noursis, M. J. (2012a). Chapter 22: General loglinear analysis.
	\emph{IBM SPSS statistics 19: Statistical 
	procedures companion.} Upper Saddle River, NJ: Prentice Hall.}
\cr\cr {Noursis, M. J. (2012b). Chapter 1: Model selection Loglinear analysis.
	\emph{IBM SPSS Statistics 19: Advanced statistical 
	procedures Companion.} Upper Saddle River, NJ: Prentice Hall.}
\cr\cr {Stevens, J. P. (2009). Chapter 14: Categorical data analysis: The log linear model. 
	\emph{Applied multivariate statistics for the social sciences (5th ed.).} 
	New York, NY: Routledge.}
\cr\cr {Tabachnick, B. G., & Fidell, L. S. (2019). Chapter 16: Multiway 
	frequency analysis. \emph{Using multivariate statistics.} New York, NY: Pearson.}
\cr\cr {Warner, R. M. (2021). Chapter 17: Chi-square analysis of contingency 
	tables. \emph{Applied statistics: Basic bivariate techniques (3rd ed.).} 
	Thousand Oaks, CA: SAGE Publications.}
}

\examples{
names(datasets)

datasets$Agresti_2019_Tab9.3

datasets$Agresti_2019_Tab9.8

datasets$Ajzen_2021_Tab7.11

datasets$Ajzen_2021_Tab7.16

datasets$Field_2018

head(datasets$Field_2018_raw)

datasets$George_2019_26_Hierarchical

datasets$George_2019_27_Nonhierarchical

datasets$Gray_2012_2way

datasets$Gray_2012_3way

datasets$Green_2014

datasets$Ho_2014

datasets$Howell_2013

datasets$Howell_2017

datasets$Meyers_2013

datasets$Noursis_2012_marital

datasets$Noursis_2012_voting_degree

datasets$Noursis_2012_voting_degree_sex

datasets$Stevens_2009_HeadStart_1

datasets$Stevens_2009_HeadStart_2

datasets$Stevens_2009_Inf_Survival

datasets$TabFid_2019_small

datasets$Warner_2020_titanic

datasets$Warner_2020_dog
}
