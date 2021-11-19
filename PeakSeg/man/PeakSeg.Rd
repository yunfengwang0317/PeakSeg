\name{PeakSeg}
\alias{PeakSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Segment for numeric vectors
%%  ~~function to do ... ~~
}
\description{
	This package is designed for finding breakpoints from any numeric vectors. The breakpoints can be located based on "Peaks". A "Peak" refers to the position with the largest deviation to either left or right flanking region.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PeakSeg(vec, winnum = 3, step = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
	  The vec represents the input numeric vector (1-demision).
%%     ~~Describe \code{vec} here~~
}
  \item{winnum}{
	  The winnum is the value of the "lag" parameter in diff function. This function is used to compute the differences of each window to the previous windows.
%%     ~~Describe \code{winnum} here~~
}
  \item{step}{
	  The step is the number of elements in the input vector to create a window. When step is set to 20, every 20 elements would be combined to generate one window.
%%     ~~Describe \code{step} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Segmented}{The complete segmented signal}
  \item{SegmentsTable}{Segments result table: (segment start index, segmentsize, segment value)}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
	PhD YF.WANG
%%  ~~who you are~~
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
function (vec, winnum = 3, step = 20) 
{
    quan = function(oo, alpha) {
        a1 = as.numeric(quantile(sort(oo), probs = c(alpha, 1 - 
            alpha))[1])
        b1 = as.numeric(quantile(sort(oo), probs = c(alpha, 1 - 
            alpha))[2])
        result = oo[oo >= a1 & oo <= b1]
        return(result)
    }
    moveave <- function(v) {
        for (i in 3:(length(v) - 2)) {
            tmp = v[(i - 2):(i + 2)]
            v[i] = mean(tmp, na.rm = T)
        }
        v[1] = mean(v[1:4], na.rm = T)
        v[2] = mean(v[2:5], na.rm = T)
        v[length(v) - 1] = mean(v[(length(v) - 4):(length(v) - 
            1)])
        v[length(v)] = mean(v[(length(v) - 3):length(v)])
        return(v)
    }
    supermedian <- function(startend, zscore) {
        if (all(is.na(zscore[startend[1]:startend[2]]))) {
            return(median(zscore, na.rm = T))
        }
        else if (sum(is.na(zscore[startend[1]:startend[2]]))/(startend[2] - 
            startend[1]) < 0.67) {
            return(median(zscore[startend[1]:startend[2]], na.rm = T))
        }
        else {
            return(median(zscore[startend[1]:startend[2]], na.rm = T) * 
                0.5)
        }
    }
    compress <- function(vec, step = 20) {
        init_pos = c(seq(1, length(vec), by = step))
        if (!length(vec) \%in\% init_pos) 
            init_pos = c(init_pos, length(vec))
        mat = vector()
        for (i in 1:(length(init_pos) - 1)) {
            mat = rbind(mat, c(init_pos[i], init_pos[i + 1]))
        }
        medvec = apply(mat, 1, supermedian, zscore = vec)
        return(medvec)
    }
    cumsumwin <- function(medvec, winnum) {
        diff_init = diff(c(rep(median(medvec), winnum), medvec, 
            rep(median(medvec), winnum)), lag = winnum)
        diff_wins = c()
        for (i in 1:(length(diff_init) - winnum)) diff_wins = c(diff_wins, 
            median(diff_init[i:(i + winnum)]))
        return(diff_wins)
    }
    findpeak <- function(v, step) {
        peak = c()
        valley = c()
        v[is.na(v)] = median(v, na.rm = T)
        trimv = quan(v, 0.05)
        CIup = median(trimv) + 1 * sd(trimv)
        CIdown = median(trimv) - 1 * sd(trimv)
        for (i in 2:(length(v))) {
            if (i >= 5 & i <= (length(v) - 5)) {
                if (v[i] > max(v[(i - 4):(i - 1)]) & max(v[(i + 
                  1):(i + 4)]) < v[i] & v[i] > CIup) {
                  peak = c(peak, i)
                }
                else if (v[i] < min(v[(i - 4):(i - 1)]) & v[i] < 
                  min(v[(i + 1):(i + 4)]) & v[i] < CIdown) {
                  valley = c(valley, i)
                }
            }
            else if (i < 5) {
                if (v[i] > max(v[-i][1:(i + 4)]) & v[i] > CIup & 
                  !any(1:(i + 4) \%in\% peak)) {
                  peak = c(peak, i)
                }
                else if (v[i] < min(v[-i][1:(i + 4)]) & v[i] < 
                  CIdown & !any(1:(i + 4) \%in\% valley)) {
                  valley = c(valley, i)
                }
            }
            else if (i > length(v) - 5) {
                if (v[i] > max(v[-i][(i - 4):(length(v) - 1)]) & 
                  v[i] > CIup & !any((i - 4):length(v) \%in\% peak)) {
                  peak = c(peak, i)
                }
                else if (v[i] < min(v[-i][(i - 4):(length(v) - 
                  1)]) & v[i] < CIdown & !any((i - 4):length(v) \%in\% 
                  valley)) {
                  valley = c(valley, i)
                }
            }
        }
        bp = sort(c(peak, valley))
        if (length(bp) > 2) {
            rmbp = c()
            for (i in 2:length(bp)) {
                if (!all(c(bp[i - 1], bp[i + 1]) \%in\% peak) & 
                  !all(c(bp[i - 1], bp[i + 1]) \%in\% valley)) {
                  rmbp = c(rmbp, i)
                }
            }
            bp = bp[-rmbp]
        }
        return(sort(unique(c(c(peak - 1, valley - 1) * step, 
            1, length(v) * step))))
    }
    medvec = compress(vec, step)
    medvec = moveave(medvec)
    forward = cumsumwin(medvec, winnum)
    backward = -rev(cumsumwin(rev(medvec), winnum))
    backward = c(rep(NA, 5), backward)
    residule = 999
    slide = 1
    for (i in 1:5) {
        res_ = sum(abs(forward - backward[i:(i + length(forward) - 
            1)]), na.rm = T)
        if (res_ < residule) {
            residule = res_
            slide = i
        }
    }
    backward = backward[slide:(slide + length(forward) - 1)]
    waves = forward + backward
    waves[is.na(waves)] = forward[is.na(waves)] * 2
    breakpoints = findpeak(waves, step)
    SegmentsTable = c()
    Segmented = c()
    for (i in 1:(length(breakpoints) - 1)) {
        from = breakpoints[i]
        to = min(length(vec), breakpoints[i + 1])
        SegmentsTable = rbind(SegmentsTable, c(from, to - from + 
            1, median(vec[from:to])))
        Segmented = c(Segmented, rep(median(vec[from:to]), to - 
            from + 1))
    }
    seg.data = list(SegmentsTable = SegmentsTable, Segmented = Segmented)
    return(seg.data)
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
