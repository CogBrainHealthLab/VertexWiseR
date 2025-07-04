#' @title Getting started with tutorials
#'
#' @description Guides the user through the various tutorials and assist with downloading the demo data
#'
#' @param demo_data A boolean object specifying whether to prompt user about downloading the demo data (default is TRUE)
#' @examples
#' getting_started()
#' @importFrom utils menu untar unzip vignette
#' @export
#' 
getting_started=function(demo_data=TRUE)
{
if (interactive()==TRUE) #only works in interactive session
{
  cat("Welcome to VertexWiseR's demo!\n")
  
  if (demo_data==TRUE)
  {
   prompt = utils::menu(title="Before starting, you may want to download the demo data folder used by the various tutorials. Would you like to download it now?",
     choices = c(paste0("Yes, in a temporary R directory (default path is tempdir())"),
                 paste0("Yes, in the current working directory (",getwd(),")"),
                 "Yes, let me type a path where it will go", 
                 "No"), 
  )
  if (prompt!=4) 
  {
    if (prompt==1) {out_dir=tempdir()}
    else if (prompt==2){out_dir=getwd()}
    else if (prompt==3)
    {out_dir <- readline("Enter the full path to the directory:")}
    
    # Folder name
    out_dir <- paste0(out_dir,"/demo_data")
    
    # Create demo_data directory if it doesn't exist
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    
    # GitHub repository URL
    base_url <- "https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/refs/heads/main/inst/demo_data/"
    
    files_to_download <- c(
      "FINK_Tv_ses13.rds",
      "FINK_behdata_ses13.rds",
      "SPRENG_CTv_site1.rds",
      "SPRENG_behdata_site1.rds",
      "fink_surf_data_hippunfold.zip",
      "model1_TFCE.rds",
      "spreng_surf_data_cat12.zip",
      "spreng_surf_data_fmriprep.zip",
      "spreng_surf_data_freesurfer.zip"
    )
    
    # Download each file
    for (file in files_to_download) {
      # Construct raw GitHub URL (different from the tree URL)
      raw_url <- paste0(base_url, file)
      
      # Destination path
      dest_path <- file.path(out_dir, file)
      
      tryCatch({
        # Download the file
        download.file(url = raw_url, destfile = dest_path, mode = "wb")
        cat("Downloaded:", file, "\n")
        
        # Check if file is a zip file and extract it
        if (grepl("\\.zip$", file, ignore.case = TRUE)) {
        #unzip or untar as unzip struggles with long paths
          zip_result <- suppressWarnings(try(unzip(zipfile = dest_path, exdir = out_dir),silent=T))
          if (is.null(zip_result) | inherits(zip_result, "try-error")==T) {untar(tarfile = dest_path, exdir = out_dir)}
          cat("Unzipped:", file, "\n")
        }
      }, error = function(e) {
        cat("Failed to download:", file, "\n")
        cat("Error:", e$message, "\n")
      })
    }
    
    cat("\nDownload and extraction process completed.\n")
    cat("\nFiles were saved in:", normalizePath(out_dir), "\n")
  }
  }

  prompt_2 = utils::menu(title="\nAvailable tutorials (will open in Rstudio's Help panel):",
            choices = c("Extracting surface data in VertexWiseR", "Example analyses with VertexWiseR - Example 1", "Example analyses with VertexWiseR - Example 2"))
  if (prompt_2==1) { vignette(package='VertexWiseR', topic='VertexWiseR_surface_extraction')}
  if (prompt_2==2) { vignette(package='VertexWiseR', topic='VertexWiseR_Example_1')}
  if (prompt_2==3) { vignette(package='VertexWiseR', topic='VertexWiseR_Example_3')}

} 
}