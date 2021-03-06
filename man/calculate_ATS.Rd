% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_ATS.R
\name{calculate_ATS}
\alias{calculate_ATS}
\title{Calculate ATS}
\usage{
calculate_ATS(
  chart_matrix,
  time_matrix,
  nobs,
  starttime,
  endtime,
  design_interval,
  n_time_units,
  time_unit,
  CL,
  no_signal_action = "omit"
)
}
\arguments{
\item{chart_matrix}{charting statistic values arranged as a numeric matrix. \cr
\code{chart_matrix[i,j]} is the jth charting statistic of the ith subject.}

\item{time_matrix}{observation times arranged as a numeric matrix. \cr
\code{time_matrix[i,j]} is the jth observation time of the ith subject,
corresponding to the time the charting statistic \code{chart_matrix[i,j]} is computed.}

\item{nobs}{number of observations arranged as an integer vector. \cr
\code{nobs[i]} is the number of observations for the ith subject.}

\item{starttime}{a numeric vector that gives the start times. \cr
\code{starttime[i]} is the time that the ith subject starts to be monitored.}

\item{endtime}{a numeric vector that gives the end times. \cr
\code{endtime[i]} is the time that the ith subject is lost to be monitored.}

\item{design_interval}{a numeric vector of length two that 
gives the left- and right- limits of the design interval. 
By default, \code{design_interval=range(time_matrix,na.rm=TRUE)}.}

\item{n_time_units}{an integer value that gives the number of basic time units
in the design time interval. \cr
The design interval will be discretized to \cr
\code{seq(design_interval[1],design_interval[2],length.out=n_time_units)}}

\item{time_unit}{an optional numeric value of basic time unit. Only used when \code{n_time_units} is missing. \cr
The design interval will be discretized to \cr
\code{seq(design_interval[1],design_interval[2],by=time_unit)}}

\item{CL}{a numeric value specifying the control limit. \cr
\code{CL} is the control limit, signals will be given if charting statistics are greater than the control limit.}

\item{no_signal_action}{a character specifying the method to use when a signal is not given to a process.
If \code{no_signal_action="omit"} take averages by omitting the processes with no signals, namely, average only the processes with signals. \cr
If \code{no_signal_action="maxtime"} impute the signal times by the maximum time, which is the right limit of design time interval. \cr
If \code{no_signal_action="endtime"} impute the signal times by the end times.}
}
\value{
a numeric value, the ATS given the charting statistics and the control limit.
}
\description{
The function \code{calculate_ATS} calculates the average time to signals (ATS) given 
a control chart matrix and a specified control limit (CL). ATS is defined as the average time from
the start of process monitoring to signal times.
}
\details{
Calculate ATS
}
\examples{
data("data_example_long_1d")

result_pattern<-estimate_pattern_long_1d(
  data_matrix=data_example_long_1d$data_matrix_IC,
  time_matrix=data_example_long_1d$time_matrix_IC,
  nobs=data_example_long_1d$nobs_IC,
  design_interval=data_example_long_1d$design_interval,
  n_time_units=data_example_long_1d$n_time_units,
  estimation_method="meanvar",
  smoothing_method="local linear",
  bw_mean=0.1,
  bw_var=0.1)

result_monitoring<-monitor_long_1d(
  data_matrix_new=data_example_long_1d$data_matrix_OC,
  time_matrix_new=data_example_long_1d$time_matrix_OC,
  nobs_new=data_example_long_1d$nobs_OC,
  pattern=result_pattern,
  side="upward",
  chart="CUSUM",
  method="standard",
  parameter=0.5)

result_ATS<-calculate_ATS(
  chart_matrix=result_monitoring$chart,
  time_matrix=data_example_long_1d$time_matrix_OC,
  nobs=data_example_long_1d$nobs_OC,
  starttime=rep(0,nrow(data_example_long_1d$time_matrix_OC)),
  endtime=rep(1,nrow(data_example_long_1d$time_matrix_OC)),
  design_interval=data_example_long_1d$design_interval,
  n_time_units=data_example_long_1d$n_time_units,
  CL=2.0)
  
}
\references{
Qiu, P. and Xiang, D. (2014). Univariate dynamic screening system: an approach for identifying individuals with irregular longitudinal behavior. Technometrics, 56:248-260. \cr
Qiu, P., Xia, Z., and You, L. (2020). Process monitoring roc curve for evaluating dynamic screening methods. Technometrics, 62(2).
}
