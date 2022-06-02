# Install rCMEM

# Uninstall and reinstall developer branch from GitHub

# 1. If rCTM is loaded and in the memory, forget rCTM
if ("rCMEM" %in% (.packages())){
  detach("package:rCMEM", unload=TRUE) 
}

# 2. If remotes is not already installed, install it
if (! ("remotes" %in% installed.packages())) {
  install.packages("remotes")
}

# 3. Install package from developer branch of GitHub
devtools::install_github("https://github.com/tilbud/rCMEM")

# 4. Load version into memory
library(rCMEM)
