% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stanmodcode_function preR3.R,
%   R/stanmodcode_function.R
\name{get_marinesurvival_stancode}
\alias{get_marinesurvival_stancode}
\title{Choose model and generate approiate code}
\usage{
get_marinesurvival_stancode(input)

get_marinesurvival_stancode(input)
}
\description{
Becasue the code is an indexing nightmare need three versions
1. One species, one datasets
2. One species, many datasets
3. Many species, one or many datasets

Becasue the code is an indexing nightmare need three versions
1. One species, one datasets
2. One species, many datasets
3. Many species, one or many datasets
}
\details{
All are very similar but with slightly differetn indexing structure (and no
latent variable if only one species).

This code chooses the right code for a given dataset

All are very similar but with slightly differetn indexing structure (and no
latent variable if only one species).

This code chooses the right code for a given dataset
}
