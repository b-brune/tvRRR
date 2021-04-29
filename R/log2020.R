#' Logarithmized Stock Index Time Series
#'
#' @details
#' `log2020` contains **2019 - 2021 Logarithmized Stock Index Time Series**
#'
#' A dataset containing logarithmized closing prices for different
#' stock indices with data ranging from June 04, 2019 to March 01, 2021.
#' Missing values have been imputed using \code{\link[imputeTS]{na_kalman}()}
#' from the \code{imputeTS} R-package.
#'
#'
#' @format
#' A data frame with 454 rows and 8 variables:
#'         \describe{
#'         \item{Date}{Date}
#'         \item{DJI}{Dow Jones Industrial Average}
#'         \item{GDAXI}{DAX Performance Index}
#'         \item{HSI}{Hang Seng Index}
#'         \item{N225}{Nikkei 225}
#'         \item{NASDAQ}{NASDAQ Composite}
#'         \item{S+P500}{S+P 500}
#'         \item{TWII}{TSEC Weighted Index}
#'         }
#'
#' @source The data was retrieved March 01, 2021 from
#' \describe{
#' \item{DAX}{\url{https://finance.yahoo.com/quote/\%5EGDAXI/history?p=\%5EGDAXI}}
#' \item{DJI}{\url{https://finance.yahoo.com/quote/\%5EDJI/history?p=\%5EDJI}}
#' \item{HSI}{\url{https://sg.finance.yahoo.com/quote/\%5EHSI/history?p=\%5EHSI}}
#' \item{N225}{\url{https://finance.yahoo.com/quote/\%5EN225/history?p=\%5EN225}}
#' \item{NASDAQ}{\url{https://finance.yahoo.com/quote/\%5EIXIC/history?p=\%5EIXIC}}
#' \item{S+P 500}{\url{https://finance.yahoo.com/quote/\%5EGSPC/history?p=\%5EGSPC}}
#' \item{TWII}{\url{https://finance.yahoo.com/quote/\%5ETWII/history?p=\%5ETWII}}
#' }
#'
#' @rdname stockdata

"log2020"
