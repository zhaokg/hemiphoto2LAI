.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  requireNamespace("utils",quitely=TRUE)
  if (interactive())
  {
    v = utils::packageVersion("hemiphoto2LAI")
    packageStartupMessage("hemiphoto2LAI", v, ". For help or bug reporting, please contact Kaiguang Zhao at lidar.rs@gmail.com.")
  }
   
   
}

.onLoad <- function(libname, pkgname) {
   #library.dynam("beast", pkgname, libname )
   utils::data(LAD_list, package=pkgname,         envir=parent.env(environment())) 
   utils::data(sampleGapData, package=pkgname,         envir=parent.env(environment())) 
   #utils::data(modis_ohio, package=pkgname,      envir=parent.env(environment())) 
   #utils::data(simAnnualData01, package=pkgname,    envir=parent.env(environment())) 
}

.onUnload <- function(libpath) {
  library.dynam.unload("hemiphoto2LAI", libpath)
}

 