#Function borrowed from reticulate (accessed 09.12.2024, partly updated with 02.09.2026):
#https://github.com/rstudio/reticulate/blob/main/R/miniconda.R
#This is used to specify a version of miniconda instead of installing the latest to ensure compatibility with the VertexWiseR requirements

miniconda_installer_py310url <- function(version = "3") {
    
  #miniconda does not yet have arm64 binaries for macOS,
  #use miniforge instead
  info <- as.list(Sys.info())
  if (info$sysname == "Darwin" && info$machine == "arm64") {
      
    url <- getOption("reticulate.miniconda.url")
    if (!is.null(url))
      return(url)
    
    base <- "https://github.com/conda-forge/miniforge/releases/latest/download"
    
    info <- as.list(Sys.info())
    arch <- miniconda_installer_arch(info)
    version <- as.character(version)
    name <- sprintf("Miniforge%s-MacOSX-%s.sh", version, arch)
    
    return(file.path(base, name))
  } 
  
  #Windows and Linux
  base <- "https://repo.anaconda.com/miniconda"
  
  info <- as.list(Sys.info())
  arch <- miniconda_installer_arch(info)
  version <- as.character(version)
  name <- if (is_windows())
    sprintf("Miniconda%s-py310_26.7.1-1-Windows-%s.exe", version, arch)
  else if (is_osx())
    sprintf("Miniconda%s-py310_26.7.1-1-MacOSX-%s.sh", version, arch)
  else if (is_linux())
    sprintf("Miniconda%s-py310_26.7.1-1-Linux-%s.sh", version, arch)
  else
    stopf("unsupported platform %s", shQuote(Sys.info()[["sysname"]]))
  
  return(file.path(base, name))
  
}


miniconda_installer_arch <- function(info) {
  
  # allow user override
  arch <- getOption("reticulate.miniconda.arch")
  if (!is.null(arch))
    return(arch)
  
  # miniconda url use x86_64 not x86-64 for Windows
  if (info$machine == "x86-64")
    return("x86_64")
  
  # otherwise, use arch as-is
  info$machine
  
}



is_windows <- function() {
  identical(.Platform$OS.type, "windows")
}

is_osx <- function() {
  Sys.info()["sysname"] == "Darwin"
}

is_linux <- function() {
  identical(tolower(Sys.info()[["sysname"]]), "linux")
}

stopf <- function(fmt = "", ...) {
  stop(sprintf(fmt, ...))
}
