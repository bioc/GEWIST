\name{gewistLevene}
\alias{gewistLevene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genome Environment Wide Interaction Search Threshold with Levene's test}
\description{
Compute the optimal Variance Prioritization power and the corresponding Levene's test prioritization p-value threshold for a given interaction effect size
}
\usage{
gewistLevene(p, N, theta_gc, theta_c, M, K = 20000, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
minor allele frequency of the SNP, a number between 0 and 0.5
}
  \item{N}{
sample size
}
  \item{theta_gc}{
proportion of quantitative trait variance explained by the interaction, should be a  number between 0 and 1
}
  \item{theta_c}{
proportion of quantitative trait variance explained by the covariate, should be a  number between 0 and 1
}
\item{M}{
total number of SNPs to be tested
}
  \item{K}{
number of procedures, by default, set to be 20,000
}
  \item{verbose}{
logical; if TRUE, function returns a data.frame class object listing the VP power at each p-value, from 0.001 to 1 with 0.001 incremental increase.
}
}					

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list with three components:
\item{Optimal_VP_power}{VP power to detect interactions at the optimal Levene's test p-value threshold}
\item{Conventional_power}{power to detect interactions without prioritization, i.e,
VP power at Levene's test p-value of 1}
\item{Optimal_pval_threshold}{levene'e test p-value at which optimal VP power is achieved}
%% ...
}
\references{
 Deng W.Q, Pare G. (2011) 
 A fast algorithm to optimize SNP prioritization for gene-gene and gene-environment interactions. 
 \emph{Genetic Epidemiology} 35: 729-738. doi: 10.1002/gepi.20624

 
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

%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

# Given a SNP with minor allele frequency of 10\% and a sample 
# of 15,000 individuals, we are interested in testing interactions
# between this SNP and a covariate of effect size #10\%. The 
# total number of SNP is 500,000. Assume the interaction 
# explains 0.1\% of the quantitative trait variance. 


library(GEWIST)
gewistLevene(p = 0.1, N = 15000, theta_gc =0.1/100, theta_c = 0.1 , M = 500000, 
K = 20000, verbose=FALSE)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ genetic interactions }
\keyword{ quantitative trait }
\keyword{ variance prioritization }
\keyword{ Levene's Test}
\keyword{ Bonferroni correction}