\name{biggest_dif}
\alias{biggest_dif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Finds the biggest positive difference in a two minute window.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

Finds the biggest positive difference between two MAP values in a two minute interval, and if this is greater than 8 TRUE is returned, if not then FALSE will be returned. If there a value below 65 it also checks if this is a part of a hypotensive event.

}
\usage{
biggest_dif(table, i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{ A dataframe.
%%     ~~Describe \code{x} here~~
}
  \item{i}{ An index.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Input should be a dataframe with two coloumns, one with MAP data and the other with HPI data, and the index input should be the index of the point in the table we have reached.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A boolean value is returned, depending on the biggest difference found in the interval.
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
function (table, i)
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
