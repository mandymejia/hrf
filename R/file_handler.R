#' Save a complex R object to a `.qs` file (temporary or persistent).
#'
#' Uses the `qs` package for high-performance serialization.
#'
#' @param object R object to save.
#' @param label Optional deterministic label for reproducible filename.
#' @param prefix Optional filename prefix (default: "qs_obj_").
#' @param preset Compression level: "fast", "balanced", "high", etc.
#' @param tmp Logical; TRUE = save to `tempdir()`, FALSE = save to `./work/`
#' @inheritParams work_dir_Param
#'
#' @return Full file path to the saved object.
#' @export
save_object <- function(object, label = NULL, prefix = "qs_obj_", compress_level = 1, tmp = TRUE, work_dir = "work") {
  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop("The 'qs2' package is required. Please install it using install.packages('qs2').")
  }

  # Decide save location
  save_dir <- if (tmp) tempdir() else normalizePath(work_dir, mustWork = FALSE)
  if (!tmp && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  if (!tmp) message("Work folder is located at: ", save_dir)
  # Generate path from label or tempfile
  file_path <- if (is.null(label)) {
    tempfile(pattern = prefix, tmpdir = save_dir, fileext = ".qs")
  } else {
    file.path(save_dir, paste0(prefix, label, ".qs"))
  }

  # Save object using qs
  qs2::qs_save(object, file_path, compress_level = compress_level)
  return(file_path)
}

#' Load a `.qs` object from disk (temporary or persistent).
#'
#' @param label Optional label used when saving (must match exactly).
#' @param file_path Optional full path to a `.qs` file. Overrides label logic.
#' @param prefix Prefix used when saving (default: "qs_obj_").
#' @param tmp Logical; only used if loading by label. TRUE = `tempdir()`, FALSE = `./work/`
#' @param delete_after_load Whether to delete the file after reading. Defaults to TRUE.
#'
#' @return Deserialized R object.
#' @export
load_object <- function(label = NULL, file_path = NULL, prefix = "qs_obj_", tmp = TRUE, delete_after_load = TRUE) {
  # Case 1: full path provided
  if (!is.null(file_path)) {
    if (!file.exists(file_path)) stop("File does not exist: ", file_path)
    obj <- qs2::qs_read(file_path)
    if (delete_after_load) unlink(file_path)
    return(obj)
  }

  # Case 2: using label + tmp
  if (is.null(label)) stop("Must provide either 'label' or 'file_path'")
  load_dir <- if (tmp) tempdir() else file.path(getwd(), "work")
  file_path <- file.path(load_dir, paste0(prefix, label, ".qs"))

  if (!file.exists(file_path)) stop("File does not exist: ", file_path)
  obj <- qs2::qs_read(file_path)
  if (delete_after_load) unlink(file_path)
  return(obj)
}

#' Try saving an object and print any serialization errors.
#'
#' @param object R object to try saving.
#' @param prefix Optional prefix for debug filename.
#'
#' @return TRUE if successful, FALSE if it fails.
#' @keywords internal
try_save_debug <- function(object, prefix = "debug_qs_") {
  tryCatch({
    path <- save_object(object, label = NULL, prefix = prefix, tmp = TRUE)
    message("Saved object successfully at: ", path)
    TRUE
  }, error = function(e) {
    message("Failed to save object: ", e$message)
    FALSE
  })
}

#' Delete the persistent 'work' directory and all saved `.qs` files.
#'
#' @return TRUE if deleted successfully, FALSE otherwise.
#' @export
delete_work_folder <- function() {
  work_dir <- file.path(getwd(), "work")
  if (dir.exists(work_dir)) {
    unlink(work_dir, recursive = TRUE, force = TRUE)
    return(TRUE)
  }
  return(FALSE)
}
