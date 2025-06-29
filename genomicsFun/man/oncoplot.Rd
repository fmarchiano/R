% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomicsFun.r
\name{oncoplot}
\alias{oncoplot}
\title{Title: Draw Oncoplot}
\usage{
oncoplot(mat_filtered, anno_df, title, rm_empty_rows = FALSE)
}
\arguments{
\item{mat_filtered}{matrix containing relevant variants}

\item{anno_df}{annotation data frame containing sample annotations}

\item{title}{Title of the oncoplot}

\item{rm_empty_rows}{Logical, if TRUE remove empty rows (default FALSE)}
}
\description{
Create annotated complex heatmap (oncoplot)
}
