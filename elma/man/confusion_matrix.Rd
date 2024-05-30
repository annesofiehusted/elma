\name{confusion_matrix}
\alias{confusion_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Loops through files in directory and runs them through the function 'fb_analysis'.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function loops through the case files in the given directory, and passes them through the forward backward analysis function 'fb_analysis', and makes a dataframe with all returns.
}
\usage{
confusion_matrix(directory, minutes, HPI)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{directory}{
%%     ~~Describe \code{x} here~~
The chosen directory.
}
 \item{minutes}{ Number of minutes in the interval.}
 \item{HPI}{ A boolean value.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
The directory for the cases to be analysised. The number of minutes defines how long we look ahead and back for the forward and backward analysis and the boolean value determines if the alarm should be defined by HPI (TRUE) or MAP (FALSE).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A dataframe is returned with all values from the analysis for the different cases.
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
function (x)
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
