.onUnload <- function (libpath) {
  library.dynam.unload("hdInference", libpath)
}
