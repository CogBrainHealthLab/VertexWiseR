#This script is a modified version of reticulate's install_python(), which installs pyenv in a user-specified directory instead of the "r-reticulate" root folder.

#As opposed to install_python() it does NOT check and handle a list of the different pyenv versions. It only installs the latest pyenv version at the user-specified path. 

#It also, as in the default install_python,does not force reinstallation, and assumes the 'optimized' option to be enabled

#The script contains a modified version of codes from pyenv.R within the reticulate package (v.1.38.0) https://github.com/rstudio/reticulate/blob/main/R/pyenv.R

python_custominstall=function(custompath) {

#This installation requires git
if (Sys.which("git") == "") { stop("Git is needed for your system to install Python instead of Miniconda. Please install Git first and ensure it is on your PATH. You may set git in R by typing:\n Sys.setenv(PATH = paste(\"[path\\to\\git]\\bin\", Sys.getenv(\"PATH\"), sep=\";\"))") }
  
  
if (missing(custompath))  #default installation
{
  
root=reticulate::install_python(version = "3.10:latest")

} else #a custom path was specified
{

#example of custom path
#root='C:/Users/john.doe/Desktop/Python'
root=custompath

##########################Installs Pyenv in the custom path##
#(modified from pyenv_bootstrap_windows() and pyenv_bootstrap_unix()


if (.Platform$OS.type=='windows') 
{

    #clone if necessary
    if (!file.exists(root)) {
      url <- "https://github.com/pyenv-win/pyenv-win"
      system(paste("git","clone", shQuote(url), shQuote(root)),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
  
    oldwd <- getwd()
    on.exit(setwd(oldwd)) #will restore user's working directory path on function break
    
    # ensure up-to-date
    setwd(root)
    system("git pull", ignore.stdout = TRUE)
    
    # path to pyenv binary
    pyenv <- file.path(root, "pyenv-win/bin/pyenv")
    
    # running 'update' after install on windows is basically required
    # https://github.com/pyenv-win/pyenv-win/issues/280#issuecomment-1045027625
    system(paste(pyenv, "update"), ignore.stdout = TRUE)
    
    # return path to pyenv binary
    Sys.setenv(PYENV_ROOT=pyenv)

} else { #########linux and mac installation##########################

    
    # move to tempdir
    owd <- setwd(tempdir())
    on.exit(setwd(owd), add = TRUE)
  
    # reticulate's pyenv installation pipeline states: 'pyenv python builds are substantially faster on macOS if we pre-install some dependencies (especially openssl) as pre-built but "untapped kegs" (i.e., unlinked to somewhere on the PATH but tucked away under $BREW_ROOT/Cellar).'
    if (.Platform$OS.type == "mac") {
      brew <- Sys.which("brew")
      
      if(brew == "" && file.exists(brew <- "/opt/homebrew/bin/brew"))
      { #adds homebrew to path, replacing withr::local_path()
        Sys.setenv(PATH=paste(normalizePath("/opt/homebrew/bin", winslash = "/", mustWork = TRUE), Sys.getenv("PATH"), sep=";"))}
  
      if(file.exists(brew)) {
        system2(paste(brew, c("install -q openssl readline sqlite3 xz zlib tcl-tk")))
        system2(paste(brew, c("install --only-dependencies pyenv python")))
      }}
    
    # download the installer
    url <- "https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer"
    name <- basename(url)
    download.file(url, destfile = name, mode = "wb")
    
    # ensure it's runnable
    Sys.chmod(name, mode = "0755")
    
    #Set custompath instead of default root
    Sys.setenv(PYENV_ROOT=root)
    
    # run the script -- for some reason, the pyenv installer will return
    # a non-zero status code even on success?
    message("Installing pyenv ...")
    
    suppressWarnings(system("./pyenv-installer", intern = TRUE))
    path <- file.path(root, "bin/pyenv")
    if (!file.exists(path)) {stop("installation of pyenv failed")}
    
    pyenv=path
  
}

###########################Now that pyenv is installed, builds Python environment with it (modified from pyenv_install)

#As if `optimized = TRUE` (the default), Python is build with:
  Sys.setenv(PYTHON_CONFIGURE_OPTS = "--enable-shared --enable-optimizations --with-lto");
  Sys.setenv(PYTHON_CFLAGS = "-march=native -mtune=native")


  # build install arguments
  if (.Platform$OS.type=='windows')
  {  
    #latest does not work on windows so have to get latest version manually
    available = system(paste0(pyenv, " install --list"), 
                       intern = TRUE)
    available <- available[startsWith(available, '3.10')]
    out <- as.character(max(numeric_version(available, strict = FALSE),
                            na.rm = TRUE))
    system(paste(pyenv, "install", out))
    
    #May get a '.vbs' file permission error on windows. Try uninstalling
    #all python versions in the pyenv environment and attempt again.
    #or see https://github.com/pyenv-win/pyenv-win/issues/236
    
  } else 
  {  args <- c("install", "--skip-existing", "3.10:latest");
     system2(pyenv, args)
  }
  
  
#################Saves Python path in a local R environment##############

} #the path is saved for both default and custom path installations
  
#will store path in .Renviron in tools::R_user_dir() 
#location specified by CRAN, create it if not existing:
envpath=tools::R_user_dir(package='VertexWiseR')
if (!dir.exists(envpath)) {dir.create(envpath) }
  
#make .Renviron file there and set conda/python paths in it:
renviron_path <- file.path(envpath, ".Renviron")


if (!missing(custompath)) 
{
   #list python executables and save path to the one 
   #in the main installation directory
    pythonfile=list.files(path = root, recursive = TRUE, full.names = TRUE)
    pattern <- "\\d/python.exe"; 
    pythonpath=pythonfile[stringr::str_detect(pythonfile, pattern)]
    
    if (!length(pythonpath) > 0) { stop('Python installation failed.') 
    } else { message('Installation successful!')}

} else #the path to .exe is directly returned by default reticulate install_python() function
{ 
    pythonpath=root
    if (!file.exists(pythonpath)) { stop('Python installation failed.') 
    } else { message('Installation successful!')}
}

# Write to the .Renviron file
env_vars <- paste0('RETICULATE_PYTHON="',pythonpath,'"')
cat(paste(env_vars, collapse = "\n"), file = renviron_path, 
    sep = "\n", append = TRUE);
fallback <- paste0('RETICULATE_PYTHON_FALLBACK="',pythonpath,'"')
cat(paste(fallback, collapse = "\n"), file = renviron_path, 
    sep = "\n", append = TRUE)


message(paste0("Your Python path is set in ", 
               renviron_path, ' \n'))

}
