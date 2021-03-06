% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_dacs.R
\docType{data}
\name{dacs}
\alias{dacs}
\alias{sys1}
\alias{sys2}
\alias{sys3}
\alias{sys4}
\alias{sys5}
\alias{sys6}
\alias{sys14c}
\alias{sys17}
\alias{sys27}
\alias{sys40}
\alias{ss1a}
\alias{ss1b}
\alias{ss1c}
\alias{ss2}
\alias{ss3}
\alias{ss4}
\alias{sys1g}
\alias{sys2g}
\alias{sys3g}
\alias{sys4g}
\alias{sys5g}
\alias{sys6g}
\alias{sys14cg}
\alias{sys17g}
\alias{sys27g}
\alias{sys40g}
\alias{ss1ag}
\alias{ss1bg}
\alias{ss1cg}
\alias{ss2g}
\alias{ss3g}
\alias{ss4g}
\alias{tohma}
\title{Software failure data}
\description{
Software failure data (time data, grouped data) collected by John Musa and Yoshihiro Tohma
}
\details{
The follwoing data are software failure data in real application, which was
collected by John Musa. The name without \code{g} (e.g. \code{sys1}, \code{sys2}) indicates
the interval time of failures in CPU time, and the last negative value
means the time in which there is no failure from the last failure. The name with
the suffix \code{g} (e.g. \code{sys1g}, \code{sys2g}) means the grouped fauilure data, i.e., the sequence
of the number of failures per woring day.
\describe{
 \item{sys1, sys1g}{136 failures have occured during system test and operations in the program with 21,700 instructions.}
 \item{sys2, sys2g}{54 failures have occured during system test and operations in the program with 27,700 instructions.}
 \item{sys3, sys3g}{38 failures have occured during system test and operations in the program with 23,400 instructions.}
 \item{sys4, sys4g}{53 failures have occured during system test and operations in the program with 33,500 instructions.}
 \item{sys5, sys5g}{831 failures have occured during system test in the program with 2,445,000 instructions.}
 \item{sys6, sys6g}{73 failures have occured during subsystem test in the program with 5,700 instructions.}
 \item{sys14c, sys14cg}{36 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{sys17, sys17g}{38 failures have occured during system test in the program with 61,900 instructions.}
 \item{sys27, sys27g}{41 failures have occured during system test in the program with 126,100 instructions.}
 \item{sys40, sys40g}{101 failures have occured during system test in the program with 180,000 instructions.}
 \item{ss1a, ss1ag}{112 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{ss1b, ss1bg}{375 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{ss1c, ss1cg}{277 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{ss2, ss2g}{192 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{ss3, ss3g}{278 failures have occured during operations in the program with hundreds of thousands instructions.}
 \item{ss4, ss4g}{196 failures have occured during operations in the program with hundreds of thousands instructions.}
}
The data \code{tohma} is the faliure data collected by Yoshihiro Tohma. This is the grouped failure data consisting of
the number of failures per a test.
}
