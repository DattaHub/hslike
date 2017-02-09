

.onLoad <- function(libname, pkgname) {
  if (exists(".hslikeStarted", .GlobalEnv)) 
    return()

  invisible(assign('.hslikeStarted', TRUE, pos=.GlobalEnv))
}

.onUnload <- function(libname, pkgname) {
  if (!exists(".hslikeStarted", .GlobalEnv)) 
    return()

  invisible(remove('.hslikeStarted', pos=.GlobalEnv))
}


