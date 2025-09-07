# Copyright 2022 UNAIDS, World Health Organization, and individual contributors.
# Author: Eduard Grebe
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(tidyverse)
library(arrow)

filter_mdri_data_assay <- function(df, assay_name, threshold, viral_load_threshold, time_cutoff) {
  #browser()
    filtered_df <- df %>%
        filter(
            assay == assay_name,
            ia_threshold == threshold,
            vl_threshold == viral_load_threshold,
            bigT == time_cutoff
        )
    return(filtered_df)
}

filter_frr_data_assay <- function(df, assay_name, threshold, viral_load_threshold, time_cutoff) {
  #browser()
  filtered_df <- df %>%
    filter(
      assay == assay_name,
      ia_threshold == threshold,
      vl_threshold == viral_load_threshold,
      bigT == time_cutoff,
      subtype == "All"
    )
  return(filtered_df)
}

filter_mdri_data_subtype <- function(df, subtypes) {
    filtered_df <- df %>%
        filter(
            subtype %in% subtypes
        )
    return(filtered_df)
}

read_bs_params = function(runs, n = 10000) {
    bsdat <- arrow::read_parquet("data/mdri_bs_params.parquet") %>%
        filter(
            run %in% runs
        ) %>%
      group_by(run) %>%
      slice(1:n) %>%
      ungroup()
    return(bsdat)
}

read_bs_params_dataset = function(runs, n = 10000) {
  bsdat <- arrow::open_dataset("data/mdri_fits_bs_dataset/") %>%
    filter(
      run %in% runs
    ) %>%
    group_by(run) %>%
    collect() %>%
    slice(1:n) %>%
    ungroup()
  return(bsdat)
}

read_bs_params_frr = function(runs, n = 10000) {
  bsdat <- arrow::open_dataset("data/prt_fits_bs_params.parquet") %>%
    filter(
      run %in% runs
    ) %>%
    group_by(run) %>%
    collect() %>%
    slice(1:n) %>%
    ungroup()
  return(bsdat)
}

read_bs_params_frr_dataset = function(runs, n = 10000) {
  bsdat <- arrow::open_dataset("data/prt_fits_bs_dataset/") %>%
    filter(
      run %in% runs
    ) %>%
    group_by(run) %>%
    collect() %>%
    slice(1:n) %>%
    ungroup()
  return(bsdat)
}


weighted_mean = function(mu, sigma, weights) {
    w = weights / sum(weights)
    mean = sum(mu * w)
    sd = sqrt(sum(sigma^2 * w^2))
    return(c(mean, sd))
}
