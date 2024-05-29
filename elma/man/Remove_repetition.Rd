\name{remove_repetition}
\alias{remove_repetition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Removes repetition
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Removes the repetition that maybe on the file from the monitor, it does so by identifying how long the file should be, and then cuts it to be this length.

When the data is played from the computer to the moniter via the soundcard, the computer replays the data when it is done. This means that the data repeats itself and the datafile from the moniter also has this repetition. We use this function to get the data without repetition
}
\usage{
remove_repetition(case, table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{case}{ Case number.
%%     ~~Describe \code{x} here~~
}
 \item{table}{ A dataframe
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The dataframe should contain two coloumns. One with the MAP data and the other with the HPI data.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a dataframe with the datapoints, where each datapoint only occurs once.
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
function (x, table)
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
