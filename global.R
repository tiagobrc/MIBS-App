# global.R is sourced by Shiny before ui.R and server.R.
#
# The app needs MIBS >= 0.0.0.9004: earlier versions fail in group.dist() on
# R >= 4.3 (the "All" plot group) and do not import stringdist.
min_mibs <- "0.0.0.9004"

if (!requireNamespace("MIBS", quietly = TRUE) ||
    utils::packageVersion("MIBS") < min_mibs) {
  stop("MIBS-app needs the MIBS package, version >= ", min_mibs, ".\n",
       "Install or update it with:\n",
       "  remotes::install_github(\"victoraLab/MIBS@package-conversion\")",
       call. = FALSE)
}
