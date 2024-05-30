\name{new_alarm_forward}
\alias{new_alarm_forward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Finds the end of an event for a forward analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds the end of an event and returns the index of when the event is over, which is defined as when HPI is less than 85.
}
\usage{
new_alarm_forward(table, i)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{ A dataframe.
%%     ~~Describe \code{x} here~~
}
 \item{i}{ An index.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Input should be a dataframe with two coloumns, one with MAP data and the other with HPI data. The index input should be the index for when the event starts.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the index for when the event is over.
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
