userlib <- Sys.getenv("R_LIBS_USER")
if(!dir.exists(userlib)) dir.create(userlib, recursive=TRUE)
options(repos = c(CRAN = "https://cloud.r-project.org/"))
packages <- c("dplyr", "readr", "stringr", "furrr", "future", "betareg", "VGAM", "mixsqp")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if(length(to_install)) install.packages(to_install, lib=userlib)
