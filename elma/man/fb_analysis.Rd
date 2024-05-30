\name{fb_analysis}
\alias{fb_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Makes a forward and backward analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function makes forward and backward analysis, where it checks in interval that is defined in the call. One can also define if the alarm should go off when there is two consecutive HPI values are greater than or equal to 85, or if it should be when there are two consecutive MAP values that are less than or equal to 72.
}
\usage{
fb_analysis(case, table, minutes, hpi=HPI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case}{ Case number.
}
  \item{table}{ A dataframe.
}
  \item{minutes}{ Number of minutes in the interval.
}
  \item{hpi}{ If the alarm should be defined by HPI.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For the input the case number should be the number of the case that is being analysed. The dataframe should consist of two coloumns, one with the MAP data and the other with the HPI data. The number of minutes defines how long we look ahead and back for the forward and backward analysis. By default the alarm will be defined by HPI, if 'HPI = False' is specified the alarm will be defined by MAP.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The function returnes the case number and the count of how many events the analysis qualifies as the different defined events.

  \item{case }{ Case number for the case being analysed}
  \item{LA }{Number of late alarms}
  \item{OA }{Number of ongoing alarms}
  \item{MI }{Number of medical interventions}
  \item{TP_f }{Number of true positives in the forward analysis}
  \item{TP_b }{Number of true positives in the backward analysis}
  \item{FP }{Number of false positives}
  \item{FN }{Number of false negatives}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Nour Mamoun Abdul Rahman Al Said,
Anne Sofie Gam Husted,
Sif Michelsen Jonsdottir
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
##function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
