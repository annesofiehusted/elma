\name{remove_noise}
\alias{remove_noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Removes noise
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function removes the static noise that are at the beginning and end of the VitalDB data.
}
\usage{
remove_noise(table)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{ A dataframe.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The dataframe should consist of 2 coloumns, one being time and the other shoud be the artirial pressure at that time.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The function returns a dataframe where tthe noise has been removed at the start and end, so it will be smaller than the original dataframe.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Nour Mamoun Abdul Rahman Al Said,
Anne Sofie Gam Husted,
Sif Michelen Jonsdottir
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
