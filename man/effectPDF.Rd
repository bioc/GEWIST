\name{effectPDF}
\alias{effectPDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genetic interaction testing given effect size distribution
}
\description{
Compute the optimal Variance Prioritization power and corresponding Levene's test p-value threshold for prioritization given the interaction effect size distribution using GEWIST.
}
\usage{
effectPDF(distribution = c("beta", "normal", "uniform", "weibull"),
parameter1, parameter2 = NULL, parameter3 = NULL, p, N, theta_c, M, 
K = 20000, nb_incr = 50,  range = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distribution}{
distribution of interaction effect size. Possible distributions are:

		"beta" for beta distribution
  
		"normal" for normal distribution
  
		"uniform" for uniform distribution
  
		"weibull" for weibull distribution
}
  \item{parameter1}{
the first parameter used in the corresponding distribution
}
  \item{parameter2}{
the second parameter used in the corresponding distribution, could set to be null
}
  \item{parameter3}{
the third parameter used in the corresponding distribution, could set to be null 
}
  \item{p}{
minor allele frequency of the SNP, a number between 0 and 0.5
}
  \item{N}{
sample size
}
  \item{theta_c}{
proportion of quantitative trait variance explained by the covariate, should be a number between 0 and 1
}
  \item{M}{
total number of SNPs to be tested
}
  \item{K}{
number of GEWIST procedures, by default, set to be 20,000
}
  \item{nb_incr}{
number of effect size points in the range to be prioritized using GEWIST; by default set to be 50.  
}
  \item{range}{
range of variance explained by interaction effect sizes, a vector of length 2
}
  \item{verbose}{
logical; if TRUE, for each interaction effect size, function returns a data.frame class object listing the VP power at each p-value, from 0.001 to 1 with 0.001 incremental increase.
}
}
%\details{
%%
%}

\value{
A list with three components:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{Optimal_VP_power}{VP power to detect interactions at the optimal VP p-value threshold}
\item{Conventional_power}{power to detect interactions without prioritization, i.e, VP power at Levene' test p-value of 1}
\item{Optimal_pval_threshold}{levene'e test p-value at which optimal VP power is achieved}
%% ...
}
\references{
 Deng W.Q, Pare G. (2011) 
 A fast algorithm to optimize SNP prioritization for gene-gene and gene-environment interactions. 
 \emph{Genetic Epidemiology}. 35: 729-738. doi: 10.1002/gepi.20624

 Pare G, Cook NR, Ridker PM, Chasman DI (2010) 
 On the Use of Variance per Genotype as a Tool to Identify Quantitative Trait Interaction Effects: A Report 
 from the Women's Genome Health Study. 
 \emph{PLoS Genet} 6(6): e1000981. doi:10.1371/journal.pgen.1000981

 Levene H. (1960) 
 Robust tests for equality of variances. In \emph{Contributions to Probability and Statistics: Essays in Honor of Harold 
 Hotelling} 
 eds:I. Olkin, S.G. Ghurye, W. Hoeffding, W.G. Madow & H.B.Mann, pp.278-292. Stanford: Stanford University Press.
}
\author{
Wei Q. Deng <dengwq@mcmaster.ca>
Guillaume Pare <pareg@mcmaster.ca>
}
%\note{
%%  ~~further notes~~
%}
\section{Warning}{Computational time is directly proportional to \emph{nb_incr}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{


# Given a SNP with minor allele frequency of 10\% and a sample 
# of 10,000 individuals, we are interested in testing interactions 
# between this SNP and a covariate of effect size 10\%. The 
# total number of SNP is 500,000. Assume the unknown interaction 
# effect size has a Weibull distribution in the range of 0.05\% 
# and 0.3\% variance explained with 50 increments.Repeat GEWIST 
# for each of the 50 interaction effect sizes.  


library(GEWIST)
effectPDF(distribution = "weibull", parameter1 = 0.8, parameter2 = 0.3,
parameter3 = NULL, p = 0.1 ,N = 10000, theta_c = 0.1, M = 350000,
K = 20000, nb_incr = 50, range = c(0.05/100,0.3/100), verbose = FALSE)

## End of script

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetic interactions }
\keyword{ quantitative trait }
\keyword{ variance prioritization }% __ONLY ONE__ keyword per line
\keyword{ Levene's Test}
\keyword{ Bonferroni correction}