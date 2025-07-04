

############################################################################################################################
############################################################################################################################
#' @title VertexWiseR system requirements installation
#'
#' @description Helps the user verify if VertexWisrR's system requirements are present and install them (a suitable 'Python' or 'Miniconda' environment, 'BrainStat' toolbox and libraries). If they are installed already, nothing will be overwritten. 
#'
#' @details VertexWiseR imports and makes use of the R package 'reticulate.' 'reticulate' is a package that allows R to borrow or translate Python functions into R. Using 'reticulate', the package calls functions from the 'BrainStat' Python module. For 'reticulate' to work properly with VertexWiseR, a Python environment needs to be installed with it — the default choice offered by VWRfirstrun is to let reticulate (version 1.41.0) create an ephemeral Python virtual environment using UV and py_require(). 
#' If for a reason this is not desirable, VWRfirstrun() also gives the choice to install a 'Miniconda' (lightweight version of Python) or Python environment in a reticulate default path or a specified path.
#' Vertex-wise statistical analyses of cortical surface require fsaverage and parcellation templates as imported by default in 'BrainStat'. 
#' The decode_surf_data() function also requires the 'Neurosynth' database to be downloaded.
#' @param requirement String that specifies a requirement to enquire about: 
#' - For only Python/Conda installation: 'python/conda only'
#' - For Python/Conda and Brainstat installation: 'conda/brainstat'
#' - For specific 'BrainStat' libraries: 'fsaverage5', 'fsaverage6', 'fslr32k', 'yeo_parcels'
#' - For the neurosynth database: 'neurosynth'. 
#' Default is 'any' and checks everything.
#' @param n_vert Numeric vector indicating the number of vertices of a given surface data so that only the required templates are asked for. It will modify the requirement argument accordingly.
#' @param promptless A boolean object specifying whether to prompt the user for action when system requirements are missing. If TRUE, VWRfirstrun() will simply inform of what is missing and will not prompt for action. Default is FALSE.
#' @return No returned value in interactive session. In non-interactive sessions, a string object informing that system requirements are missing.
#' @examples
#' VWRfirstrun()
#' @importFrom reticulate conda_binary py_module_available miniconda_path py_require
#' @importFrom fs path_home
#' @importFrom methods is 
#' @importFrom utils menu
#' @importFrom rappdirs user_data_dir
#' @importFrom tools R_user_dir
#' @export

VWRfirstrun=function(requirement="any", n_vert=0, promptless=FALSE) 
{
  #First checks the n_vert argument. This ensures only the necessary fsaverage data is demanded:
  #are fsaverage5 templates in brainstat_data?
  if (n_vert==20484)
  {requirement='fsaverage5'}  
  #are fsaverage6 templates in brainstat_data?
  if  (n_vert==81924)
  {requirement='fsaverage6'}
  #are fslr32k templates in brainstat_data?
  if  (n_vert==64984)
  {requirement='fslr32k'}
  #is yeo parcellation data in brainstat_data?
  #if (n_vert>0 & n_vert!=20484 & n_vert!=81924 & n_vert!=64984)
  #{requirement='yeo_parcels'} 
  
  # If custom installation paths have been defined by the user, source
  # them from the package directory:
  Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
  if (file.exists(Renvironpath)) {
    readRenviron(Renvironpath)
  
  #if cache was cleared and UV environment no longer exists, remove it
  if (!file.exists(Sys.getenv('VIRTUAL_ENV')))
  {
    Sys.unsetenv('VIRTUAL_ENV')
    lines = readLines(Renvironpath) 
    lines = lines[!grepl("^VIRTUAL_ENV", lines)]
    writeLines(lines, Renvironpath)
  }
  }
  
  #default time limit to download is 60s which can be too short:
  options(timeout=500); #set to 500s instead
  on.exit(options(timeout=60))
  
  if (interactive()==TRUE & promptless==FALSE) { 
    #can only run interactively as it requires user's action
    
    ##################################################################
    ##################################################################
    #check if miniconda (or suitable python environment) is installed
    message('Checking for Miniconda or Python environment...')
    
    if(reticulate::py_available()==FALSE) #fast, but unreliable check if Python is in a custom location. It will work if VWRfirstrun() was run once already
    {
      
      if (is.null(tryCatch(reticulate::py_discover_config(), error = function(e) NULL))) #slow but reliable check if first time 
      {
        missingobj=1
        
        prompt = utils::menu(title="A suitable Python environment for reticulate or Miniconda could not be found in the environment. \n Do you want Miniconda or Python to be installed now?",
                  choices = c("Yes, install a virtual Python environment (Recommended)",
                              "Yes, install Miniconda",
                              "Yes, install Python", 
                              "No"), 
                             )
        if (prompt==1) #Install ephemeral virtual environment via UV 
        { 
          #reset variable if RETICULATE_PYTHON was set to NA
          Sys.unsetenv("RETICULATE_PYTHON")
          
          #check if virtual environment available with the right settings and initialize it.
          message('Installing Python ephemeral environment via reticulate\'s py_require() and UV...\n')
          reticulate::py_require(packages=c("numpy<=1.26.4","matplotlib","brainstat==0.4.2","vtk==9.3.1"), python_version = "<3.11")
          reticulate::py_config()
          
          #will store cache path in .Renviron in tools::R_user_dir() 
          #location specified by CRAN, create it if not existing:
          envpath=tools::R_user_dir(package='VertexWiseR')
          if (!dir.exists(envpath)) {dir.create(envpath, 
                                                recursive = TRUE)} 
          #make .Renviron file there and set conda/python paths:
          renviron_path <- file.path(envpath, ".Renviron")
          env_vars <- c(paste0('VIRTUAL_ENV="',
                               Sys.getenv('VIRTUAL_ENV'),'"'))
          cat(paste0(env_vars, "\n", collapse = "\n"), 
              file = renviron_path, sep = "\n", append = TRUE)
          
          message(paste0('\nThis virtual environment is cached in the following path:\n', tools::R_user_dir("reticulate", "cache"), 
                         '\n\nIf it gets cleared, VWRfirstrun() will prompt you again for a new installation. Unless the cache is removed, it will NOT use a new virtual environment created via py_require(), as the virtual path is fixated there:\n',
                         paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
                         ,'\n'))
        } 
        else if (prompt==2) #Install Miniconda
        {
          
          #define default miniconda directory 
          defaultpath=reticulate::miniconda_path()
          
          #give the choices to specify path
          choice = utils::menu(c("Default", "Custom"), 
                               title=paste0("Miniconda's default installation path is ", defaultpath,". Type \"1\" or \"Default\" if you want to install Miniconda in the default Path. \nYou can, alternatively, type your own path (note that the Miniconda installer does not support paths containing spaces)."))
          
          if (choice==1) #Install Miniconda within default path
          { 
            message('Installing Miniconda (v24.9.2)...')
            Sys.setenv(RETICULATE_PYTHON=defaultpath)
            
            #custom url to get version 24.9.2
            on.exit(options(reticulate.miniconda.url=NULL))
            options(reticulate.miniconda.url=miniconda_installer_py39url())
            reticulate::install_miniconda(update = FALSE)
            message("Installing dependency packages with appropriate versions...")
            reticulate::py_install("numpy==1.26.4", pip=TRUE, 
                                   envname = defaultpath)
            reticulate::py_install("vtk==9.3.1",pip = TRUE, 
                                   envname = defaultpath) # latest vtk==9.4.0 causes problems
            
            #will store path in .Renviron in tools::R_user_dir() 
            #location specified by CRAN, create it if not existing:
            envpath=tools::R_user_dir(package='VertexWiseR')
            if (!dir.exists(envpath)) {dir.create(envpath, 
                                                  recursive = TRUE) }
            
            #get path to python executable (will differ across OS)
            #silenced with sink and con
            defaultpathexe <- tryCatch({
              tmpfile <- tempfile()
              con <- file(tmpfile, open = "wt")
              sink(con, type = "message")
              result <- reticulate::py_discover_config()$python  
              # Close sink and connection 
              sink(type = "message")
              close(con)
              file.remove(tmpfile)
              result
            }, error = function(e) {
              # Ensure cleanup in case of error
              sink(type = "message")
              close(con)
              file.remove(tmpfile)
              NULL
            })
            
            
            #make .Renviron file there and set conda/python paths:
            renviron_path <- file.path(envpath, ".Renviron")
            env_vars <- c(
              paste0('RETICULATE_PYTHON="',
                     defaultpathexe,'"'),
              paste0('RETICULATE_MINICONDA_PATH="',
                     defaultpath,'"'),
              paste0('RETICULATE_PYTHON_FALLBACK="',
                     defaultpath,'"'))
            # Write to the .Renviron file
            cat(paste(env_vars, "\n", collapse = "\n"), 
                file = renviron_path, sep = "\n", append = TRUE)
            #now everytime VWRfirstrun() is called, the .Renviron file
            #is read and the custom path accessed:
            readRenviron(renviron_path)
 
          }
          else {        #Install Miniconda within custom path
            
            userpath <- readline("Enter the full path to the directory:")
            Sys.setenv(RETICULATE_PYTHON=userpath)
            #will store path in .Renviron in tools::R_user_dir() 
            #location specified by CRAN, create it if not existing:
            envpath=tools::R_user_dir(package='VertexWiseR')
            if (!dir.exists(envpath)) {dir.create(envpath, recursive = TRUE) }
            #make .Renviron file there and set conda/python paths:
            renviron_path <- file.path(envpath, ".Renviron")
            env_vars <- c(
              paste0('RETICULATE_MINICONDA_PATH="',
                     userpath,'"'),
              paste0('RETICULATE_PYTHON_FALLBACK="',
                     userpath,'"'))
            # Write to the .Renviron file
            cat(paste(env_vars, "\n", collapse = "\n"), 
                file = renviron_path, sep = "\n", append = TRUE)
            #now everytime VWRfirstrun() is called, the .Renviron file
            #is read and the custom path accessed:
            readRenviron(renviron_path)
            message(paste0("Your custom Miniconda path is set in ", 
                           renviron_path, ' \n'))
            
            #Install miniconda in the new path
            message('Installing Miniconda (v24.9.2)...')
            #custom url to get version 24.9.2
            on.exit(options(reticulate.miniconda.url=NULL))
            options(reticulate.miniconda.url=miniconda_installer_py39url())
            #install_miniconda will use miniconda_path() which relies on RETICULATE_MINICONDA_PATH defined above
            reticulate::install_miniconda(update = FALSE)
            message("Installing dependency packages with appropriate versions...")
            #set environment variable to make sure packages 
            #arrive at the same place, not in 'r-miniconda'
            reticulate::py_install("numpy==1.26.4", pip=TRUE, 
                                   envname=userpath)
            reticulate::py_install("vtk==9.3.1",pip = TRUE, 
                                   envname=userpath) # latest vtk==9.4.0 causes problems
            
            #get path to python executable (will differ across OS)
            #silenced with sink and con
            userpathexe <- tryCatch({
              tmpfile <- tempfile()
              con <- file(tmpfile, open = "wt")
              sink(con, type = "message")
              result <- reticulate::py_discover_config()$python  
              # Close sink and connection 
              sink(type = "message")
              close(con)
              file.remove(tmpfile)
              result
            }, error = function(e) {
              # Ensure cleanup in case of error
              sink(type = "message")
              close(con)
              file.remove(tmpfile)
              NULL
            })
            
            #add user path to the .Renviron:
            env_vars <- paste0('RETICULATE_PYTHON="',
                               userpathexe,'"')
            # Write to the .Renviron file and read again
            cat(paste(env_vars, "\n", collapse = "\n"), 
                file = renviron_path, sep = "\n", append = TRUE)
            readRenviron(renviron_path)
            
          }
          
        message('Please restart R after Miniconda installation for its environment to be properly detected by reticulate.')
          
        }
        else if (prompt==3) #Install Python instead of Miniconda
        { 
          
          #give the choices to specify path
          choice = utils::menu(title=paste0("Reticulate's Python default installation is done within 'r-reticulate' in your OS libraries, via install_python(). Type \"1\" or \"Default\" if you want the default installation. \nYou can, alternatively, specify the path where Python will be installed."),
                               choices=c("Default", "Custom"))
          
          if (choice==1) #Install Python within default path
          {
            python_custominstall()
          }
          else 
          {        #Install Python within custom path
            userpath <- readline("Enter the full path:")
            message(paste('Installing Python in',userpath,'...\n'))
            python_custominstall(userpath)
          }
          
          
          #Read new python enviroment
          Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
          if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
          
          #vtk installation
          message('A specific version of the vtk package (v9.3.1) will be installed.\n')
          #latest vtk==9.4.0 causes problems
          status <- system(paste0(Sys.getenv('RETICULATE_PYTHON')," -m pip install vtk==9.3.1"));
          #if failed then try pip3 
          if (status != 0) { 
            message('Could not install package with pip, trying pip3...\n')
            status2 = system(paste0(Sys.getenv('RETICULATE_PYTHON')," -m pip3 install vtk==9.3.1")); 
          if (status2 != 0) {message("vtk 9.3.1 was not installed successfully. This may cause further issues.")}
          }
          
          message('Please restart R after Python installation for its environment to be properly detected by reticulate.')
          
      }
      else { 
        stop('VertexWiseR will not work properly without Miniconda or a suitable version of Python for reticulate.\n\n')
      }
      
    } else 
    {
    reticulate::py_config() #if py_discover_config() successful, initialise for py_available to return TRUE next time
    }} 
    
    #################################################################
    ##################################################################
    ###check if BrainStat is installed
    if (requirement!="python/conda only")
    {message('Checking for BrainStat package...')}
    
    if( !reticulate::py_module_available("brainstat") 
        & requirement!="python/conda only")
    {
      missingobj=1
      
      prompt = utils::menu(c("Yes", "No"), title="\nThe Brainstat package could not be found in your Python/Conda environment. It is needed for vertex-wise linear models and the surface plotter to work. \n Do you want Brainstat (v0.4.2) to be installed now (~1.65 MB)? The NiMARE (~20.4 MB) and Brainspace (~84.2 MB) libraries and other BrainStat dependencies will automatically be installed within your Python library.")
      if (prompt==1)
      {	
        
        #if miniconda exists, use reticulate to install
        if (tryCatch(file.exists(reticulate::conda_binary()), error = function(e) FALSE)) {
          reticulate::py_install("brainstat==0.4.2",pip=TRUE) 
        }
        else { #if only Python, install via pip
        
          status <- system("pip install brainstat==0.4.2");
          #if failed then try pip3 
          if (status != 0) {
            message('Could not install package with pip, trying pip3...\n')
            system("pip3 install brainstat==0.4.2");}
          
          #reticulate might not search again for the list of modules
          #so R needs to be restarted
          message('Restarting R may be needed for the installation to take effect.')
        }
        
      } 
      else {
        stop('VertexWiseR will not work properly without BrainStat.\n\n')}
    }
    
    ###############################################################################################################################
    #check if BrainStat fsaverage/parcellation templates are installed (stops only if function needs it)
    
#can be avoided with if the requirement is only to check python conda, and/or brainstat installations
if (requirement!="python/conda only" & requirement!='conda/brainstat')
{   
    message('Checking BrainStat\'s analysis data...')
  
    #if no path to BrainStat's data path is defined yet
    #and no default $HOME/brainstat_data folder exists,
    #prompt for user to define a path or set default
    if (Sys.getenv('BRAINSTAT_DATA')=="" & 
        !dir.exists(paste0(fs::path_home(),'/brainstat_data/'))) 
    {
      missingobj=1
      
      choice = utils::menu(c("Default", "Custom"), title=paste0('By Default, BrainStat stores data required for analyses in ', fs::path_home(), '/brainstat_data/. Alternatively, you may want to specify your own custom path for the BrainStat data. Where do you want BrainStat data to be saved?\n'))
      if (choice==1)  #set path to $HOME_DIR by default
      {
        message(paste('The brainstat_data directory will be located in',
                      fs::path_home()))
        brainstat_data_path=fs::path_home()
      }
      else if (choice==2) #set custom path
      {
        brainstat_data_usrpath <- readline("Enter the full path to store Brainstat data:")
        message(paste('The path to BrainStat data is now',
                      brainstat_data_usrpath,'.'))
        dir.create(brainstat_data_usrpath, showWarnings = FALSE,
                   recursive = TRUE)
        
        if (!dir.exists(brainstat_data_usrpath))
        {stop('No directory could be found or created via dir.create()  within the given path. Please try typing it again.')}
        
        brainstat_data_path=brainstat_data_usrpath
      }
      
      #Write it in a local .Renviron file 
      envpath=tools::R_user_dir(package='VertexWiseR')
      if (!dir.exists(envpath)) {dir.create(envpath, recursive = TRUE) }
      #make .Renviron file there and set conda/python paths:
      renviron_path <- file.path(envpath, ".Renviron")
      env_vars <- paste0('BRAINSTAT_DATA="',brainstat_data_path,'"')
      # Write to the .Renviron file
      cat(paste0(env_vars, "\n", collapse = "\n"), 
          file = renviron_path, sep = "\n", append = TRUE)
      #now everytime VWRfirstrun() is called, the .Renviron file 
      #is read and the custom path accessed:
      readRenviron(renviron_path)
      
      message(paste0("The path where BrainStat data will be stored is set in ", renviron_path))
      
      
      #if no path set but brainstat_data folder is already here, it stays as default:
    } else if (Sys.getenv('BRAINSTAT_DATA')=="" & 
               dir.exists(paste0(fs::path_home(),'/brainstat_data/')))
    {
      brainstat_data_path=fs::path_home()
    } else if (!Sys.getenv('BRAINSTAT_DATA')=="") 
      #if a path has been set, it is read in the .Renviron local file as BRAINSTAT_DATA variable
    {
      brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')
    }
    

    #create brainstat_data folder substructure if not already in 
    #the default or custom path 
    #(this allows later downloads to be specifically located)
    dir.create(paste0(brainstat_data_path,'/brainstat_data/'),
               showWarnings = FALSE)
    dir.create(paste0(brainstat_data_path, '/brainstat_data/surface_data/'), showWarnings = FALSE)
    dir.create(paste0(brainstat_data_path,'/brainstat_data/parcellation_data/'), showWarnings = FALSE)
    
    
    #for each required data, it will now check the set path (default or custom) and prompt for download if they are missing
  
    #fsaverage5 data
    if ((requirement=="any" | requirement=='fsaverage5')==TRUE 
        & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
    { 
      missingobj=1
      
      prompt = utils::menu(c("Yes", "No"), title=paste("VertexWiseR could not find BrainStat's fsaverage5 templates in", brainstat_data_path, ". They are needed if you want to analyse cortical surface in fsaverage5 space. \n  Do you want the fsaverage5 templates (~7.81 MB) to be downloaded now?"))
      
      if (prompt==1){    
        brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
        brainstat.datasets.base$fetch_template_surface(template = "fsaverage5", data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))
        
      } else if (requirement=='fsaverage5') {
        stop('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')
      } else if (requirement=='any') {
        warning('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')}
    } 
    
    #Fsaverage6 data
    if ((requirement=="any" | requirement=='fsaverage6')==TRUE 
        & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6'))) 
    { 
      missingobj=1
      
      prompt = utils::menu(c("Yes", "No"), title=paste("VertexWiseR could not find BrainStat's fsaverage6 templates in", brainstat_data_path, ".  They are needed if you want to analyse cortical surface in fsaverage6 space. \n Do you want the fsaverage6 templates (~31.2 MB) to be downloaded now?"))
      
      if (prompt==1)
      { brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
      brainstat.datasets.base$fetch_template_surface(template = "fsaverage6", data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))
      
      } else if (requirement=='fsaverage6') { 
        stop('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')
      } else if (requirement=="any") {
        warning('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')}
    } 
    
    #fsLR32k data
    if ((requirement=="any" | requirement=='fslr32k')==TRUE 
        & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-conte69'))) 
    { 
      missingobj=1
      
      prompt = utils::menu(c("Yes", "No"), title=paste("VertexWiseR could not find BrainStat's fslr32k templates in", brainstat_data_path, ".  They are needed if you want to analyse cortical surface in fslr32k space. \n Do you want the fslr32k templates (~4.12 MB) to be downloaded now?"))
      
      if (prompt==1)
      { brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
      brainstat.datasets.base$fetch_template_surface(template = "fslr32k", data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))
      
      } else if (requirement=='fslr32k') { 
        stop('VertexWiseR will not be able to analyse fslr32k data without the brainstat templates.\n\n')
      } else if (requirement=="any") {
        warning('VertexWiseR will not be able to analyse fslr32k data without the brainstat templates.\n\n')}
    } 
    
    #Yeo parcellatio data
#    if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='fslr32k' | requirement=='yeo_parcels')==TRUE 
#        & !file.exists(paste0(brainstat_data_path,'/brainstat_data/parcellation_data/__MACOSX/')))
#    {
#      missingobj=1
#      
#      prompt = utils::menu(c("Yes", "No"), title=paste("VertexWiseR could not find BrainStat's yeo parcellation data in", brainstat_data_path, ". They are fetched by default by BrainStat for vertex-wise linear models to run and cannot be ignored. \n Do you want the yeo parcellation data (~1.01 MB) to be downloaded now?"))
#      
#      if (prompt==1){    
#        brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
#        try(brainstat.datasets.base$fetch_parcellation(template="fsaverage",atlas="yeo", n_regions=7, data_dir = paste0(brainstat_data_path,'/brainstat_data/parcellation_data/')), 
#            silent=TRUE)}  
#      
#      else if  (requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='fslr32k' | requirement=='yeo_parcels') 
#      {
#        stop('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')}
#      else if (requirement=="any") 
#      {
#        warning('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')
#      }
#    }
} 
    
    #####################################################################
    #####################################################################
    #Check if neurosynth database is present and download
    if ((requirement=="any" | requirement=='neurosynth')==TRUE 
        & !file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR')))
    {
      missingobj=1
      
      prompt = utils::menu(c("Yes", "No"), title=paste0(
        "\nneurosynth_dataset.pkl is not detected inside VertexWiseR's installed package directory (", 
        system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'), 
        "). It is needed to be able to run decode_surf_data(). It can be downloaded from the github VertexWiseR directory.\n\nDo you want the neurosynth database (7.5 MB) to be downloaded now?"))
      if (prompt==1) {
        
        #function to check if url exists
        #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
        #https://stackoverflow.com/a/60627969
        valid_url <- function(url_in,t=2){
          con <- url(url_in)
          check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=TRUE)[1])
          suppressWarnings(try(close.connection(con),silent=TRUE))
          ifelse(is.null(check),TRUE,FALSE)}
        
        #Check if URL works and avoid returning error but only print message as requested by CRAN:
        url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz"
        if(valid_url(url)) {
          download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz",destfile = paste0(system.file(package='VertexWiseR'),'/extdata/neurosynth_dataset.pkl.gz'))
        } else { 
          warning("The neurosynth database (neurosynth_dataset.pkl.gz) failed to be downloaded from the github VertexWiseR directory. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata and download the object manually.") #ends function
        } 
        
        #if user refuses, stops if required, just returns a message if optionnal at this stage
      } else if (requirement=="neurosynth") {
        stop("\ndecode_surf_data() can only work with the neurosynth database.\n") }       
      else if (requirement=="any") {
        warning("\ndecode_surf_data() can only work with the neurosynth database.\n")}
    }
    
    if(!exists('missingobj'))
    { message('No system requirements are missing. \u2713 \n') }
  }
  
  #############################################################
  #############################################################
  #############################################################
  #If the session is non-interactive and any required file is missing, the script will stop and require the user to run VWRfirstrun() interactively.
  #non-interactive checks proceed when VWRfirstrun() is called by a function with the argument VWR_check defined.
  
  #If the session is interactive but has the promptless option set as TRUE, the same checks will be run to inform of what is missing while not prompting for any action from the user.
  
  else if (exists("VWR_check") | 
           (interactive()==TRUE & promptless==TRUE)) 
  { 
    #creates the 'non_interactive' message to inform upper functions that the non-interactive session cannot proceed when files are missing
    non_interactive="Run VWRfirstrun() in an interactive R session to check for the missing system requirements and to install them."
    
    #miniconda or python missing?
    if(reticulate::py_available()==FALSE) #fast, but unreliable check if Python is in a custom location. It will work if VWRfirstrun() was run once already
    {
      if (is.null(reticulate::py_discover_config()) |
          is(tryCatch(reticulate::py_discover_config(), error=function(e) e))[1] == 'simpleError') #slow but reliable check if first time 
      {
        missingobj="Miniconda or a suitable version of Python for reticulate could not be found in the environment.\n"
        
        if (interactive()==FALSE)
        { non_interactive=paste0(missingobj,non_interactive)
        return(non_interactive)
        } else {message(missingobj)}
      } else {
        reticulate::py_config #needed for py_available to return TRUE at the next check if py_discover_config was successful
      }
    } 
    
    #brainstat missing
    if (!reticulate::py_module_available("brainstat") 
        & (requirement!="python/conda only" | requirement=="conda/brainstat")) 
    {
      missingobj="The Brainstat package could not be found in your Python/Conda environment. It is needed for vertex-wise linear models and the surface plotter to work.\n";
      
      if (interactive()==FALSE)
      { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
      } else {message(missingobj)}
    } 
    
    #For brainstat data, it will look either in default $HOME path or 
    #custom if it's been set
    if (Sys.getenv('BRAINSTAT_DATA')=="")
    { 
      brainstat_data_path=fs::path_home()
    } 
    else if (!Sys.getenv('BRAINSTAT_DATA')=="") 
    {
      brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')
    }
    
    
    
    #fsaverage5 missing
    if ((requirement=="any" | requirement=='fsaverage5')==TRUE & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
    {
      missingobj=paste0("VertexWiseR could not find brainstat fsaverage5 templates in the ",brainstat_data_path,"/brainstat_data/ directory. They are needed if you want to analyse cortical surface in fsaverage5 space.\n");
      
      if (interactive()==FALSE)
      { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
      } else {message(missingobj)}
    } 
    
    #fsaverage6 missing
    if ((requirement=="any" | requirement=='fsaverage6')==TRUE & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6')))  
    {
      missingobj=paste0("VertexWiseR could not find brainstat fsaverage6 templates in the ",brainstat_data_path,"/brainstat_data/ directory. They are needed if you want to analyse cortical surface in fsaverage6 space.\n");
      
      if (interactive()==FALSE)
      { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
      } else {message(missingobj)}
    } 
    
    #fslr32k missing
    if ((requirement=="any" | requirement=='fslr32k')==TRUE 
        & !file.exists(paste0(brainstat_data_path,'/brainstat_data/surface_data/tpl-conte69'))) 
    {
      missingobj=paste0("VertexWiseR could not find brainstat fslr32k templates in the ",brainstat_data_path,"/brainstat_data/ directory. They are needed if you want to analyse cortical surface in fslr32k space.\n");
      
      if (interactive()==FALSE)
      { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
      } else {message(missingobj)}
    } 
    
    #yeo parcels missing
    #if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels')==TRUE 
    #    & !file.exists(paste0(brainstat_data_path,
    #                          '/brainstat_data/parcellation_data/__MACOSX/'))) 
    #{
    #  missingobj=paste0("VertexWiseR could not find brainstat yeo parcellation data in the ",brainstat_data_path,"/brainstat_data/ directory. They are fetched by default by brainstat for vertex-wise linear models to run and cannot be ignored.\n");
    #  
    #  if (interactive()==FALSE)
    #  { non_interactive=paste0(missingobj,non_interactive)
    #  return(non_interactive)
    #  } else {message(missingobj)}
    #} 
    
    #neurosynth data missing
    if ((requirement=="any" | requirement=='neurosynth')==TRUE & !file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))) 
    {
      missingobj=paste0("neurosynth_dataset.pkl is not detected inside VertexWiseR's installed package directory (", system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'), "). It is needed to be able to run decode_surf_data().\n");
      
      if (interactive()==FALSE)
      { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
      } else {message(missingobj)}
    } 
    
    #If nothing is missing, missingobj will not have been created
    #thus allowing the function to inform that no (specified) requirements are missing
    if(!exists('missingobj'))
    { message('No system requirements are missing. \u2713 \n') }
    
  }
  
}

