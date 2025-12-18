# Kitchen Utensils -
# This file contains functions for the pipeline.
# The functions are sorted alphabetically

#' process_in_chunks (parallel)
#'
#' Run an action over chunks of a large file in parallel.
#'
#' @param large_file_path Path to a large csv file.
#' @param chunk_size Number of rows per chunk.
#' @param action Function to apply to each chunk.
#' @param ... Additional parameters to pass to the action.
#' @param workers Number of parallel workers (default: 6)
process_in_chunks <- function(large_file_path, chunk_size = 1e6, action, ...,
                              workers = 6, tmp_dir = "/tmp") {
  library(data.table)
  library(furrr)
  library(fs)
  library(future)

  ## -------- Decompress once if .gz -------- ##
  local_path <- large_file_path
  staged <- FALSE

  if (grepl("\\.gz$", large_file_path, ignore.case = TRUE)) {
    staged <- TRUE
    tmp_file <- file.path(tmp_dir, paste0("decompressed_", basename(tools::file_path_sans_ext(large_file_path))))
    message("Decompressing ", large_file_path, " to ", tmp_file, " ...")
    system2("zcat", args = shQuote(large_file_path), stdout = tmp_file)
    message("Decompression complete.")
    local_path <- tmp_file

    on.exit({
      if (fs::file_exists(tmp_file)) {
        fs::file_delete(tmp_file)
        message("Temporary file ", tmp_file, " deleted.")
      }
    }, add = TRUE)
  }

  ## -------- Parallel plan setup -------- ##
  current_plan <- tryCatch(future::plan("list"), error = function(e) NULL)
  if (is.null(current_plan) || length(current_plan) == 0) {
    plan_type <- if (.Platform$OS.type == "unix") multicore else multisession
    plan(plan_type, workers = workers)
  }

  ## -------- Get header and total lines -------- ##
  header_col_names <- fread(local_path, header = TRUE, sep = ",", nrows = 0) |> colnames()
  wc_out <- system(paste("wc -l", shQuote(local_path)), intern = TRUE)
  total_lines <- as.numeric(strsplit(wc_out, " ")[[1]][1])
  message("Detected ", total_lines, " total lines (including header).")

  message("Starting parallel chunk processing with ", workers, " workers...")
  chunk_idx <- 1
  chunk_collector <- list()

  repeat {
    chunk_indices <- chunk_idx:(chunk_idx + workers - 1)
    chunk_results <- future_map(
      chunk_indices,
      function(i) {
        skip <- chunk_size * (i - 1) + 1
        if (skip > total_lines) return(NULL)  # stop overshoot
        dt <- fread(local_path, header = FALSE, sep = ",",
                    col.names = header_col_names,
                    nrows = chunk_size, skip = skip)
        if (nrow(dt) == 0) return(NULL)
        message("Working on chunk ", i, " with ", nrow(dt), " rows.")
        do.call(action, c(list(uncollapsed_raw_counts = dt), list(...)))
      },
      .progress = TRUE
    )

    chunk_results <- Filter(Negate(is.null), chunk_results)
    if (length(chunk_results) == 0) break
    chunk_collector <- c(chunk_collector, chunk_results)
    chunk_idx <- chunk_idx + workers
  }

  plan(sequential)
  chunk_collector
}



#' Read a CSV file with enforced data types from a master schema using data.table
#'
#' This function reads a CSV file into a data.table, applying column types
#' specified in a central schema.yaml file. It is designed for high performance
#' by leveraging data.table::fread. It uses the 'here' package to reliably
#' locate the schema file relative to the project root.
#'
#' @param csv_path The path to the CSV file to be read.
#' @param schema_path The path to the schema YAML file. Defaults to
#'   a path constructed from the project root.
#'
#' @return A data.table with data types enforced according to the schema.
#'
#' @import data.table
#' @import yaml
#' @import here
#'
#' @examples
#' \dontrun{
#' df <- read_data_table("path/to/your/data.csv")
#' }
read_data_table <- function(csv_path, schema_path = NULL, nrows = Inf) {

  # Ensure required packages are available
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required but not installed.")
  }
  if (!requireNamespace("here", quietly = TRUE)) {
    stop("Package 'here' is required but not installed.")
  }

  # --- Define the mapping from abstract types to data.table class names ---
  # Note: For performance, datetimes are read as character and converted after.
  type_mapping <- c(
    "integer" = "integer",
    "float" = "numeric",
    "string" = "character",
    "boolean" = "logical",
    "datetime" = "character" # Read as character first for speed
  )

  # --- Determine the schema path ---
  if (is.null(schema_path)) {
    # Automatically find the schema file from the project root
    schema_path <- here::here("sushilib", "sushi_io", "schema.yml")
  }

  if (!file.exists(schema_path)) {
    stop(paste("Schema file not found at the expected path:", schema_path))
  }

  # --- Read the master schema file ---
  schema <- yaml::read_yaml(schema_path)

  # --- Create the colClasses vector for fread ---
  # This creates a named character vector, e.g., c(col_a = "integer", col_b = "character")
  col_classes <- vapply(schema, function(dtype) type_mapping[[dtype]], FUN.VALUE = character(1))

  # --- Filter col_classes to only include columns present in the CSV ---
  # First, read only the header of the target file.
  file_headers <- names(data.table::fread(csv_path, nrows = 0))

  # Find the intersection of columns in the schema and columns in the file.
  valid_cols <- intersect(names(col_classes), file_headers)

  # Create a new col_classes vector containing only the columns that actually exist.
  final_col_classes <- col_classes[valid_cols]


  cat("Applying data.table Schema to existing columns:\n")
  print(final_col_classes)

  # --- Read the CSV using fread with the filtered schema ---
  dt <- data.table::fread(csv_path, colClasses = final_col_classes, sep=",", header = TRUE, nrows = nrows)

  # --- Post-processing for datetime columns ---
  # Find which columns were originally specified as 'datetime' that also exist in the data
  datetime_cols <- names(schema)[which(unlist(schema) == "datetime")]
  datetime_cols_in_data <- intersect(datetime_cols, names(dt)) # Ensure we only try to convert existing columns

  if (length(datetime_cols_in_data) > 0) {
    cat("\nConverting datetime columns:", paste(datetime_cols_in_data, collapse = ", "), "\n")
    # Use data.table's efficient `:=` to convert columns by reference
    for (col in datetime_cols_in_data) {
      # The get() is no longer needed in modern data.table when using (col)
      dt[, (col) := as.POSIXct(dt[[col]])]
    }
  }

  return(dt)
}

#' validate_columns_entries
#'
#' This function checks that for a list of columns, all entries are filled in.
#' It checks all column entries against a list of potential empty values.
#'
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @param empty_values Optional vector of values that equate to empty. Defaults to NA, "NA", "", and " ".
#' @return Boolean
validate_columns_entries= function(selected_columns, df, empty_values= c(NA, 'NA', '', ' ')) {
  # Check for rows in selected_columns that equate to predefined empty values.
  missing_rows= df %>% dplyr::filter(if_any(all_of(selected_columns), ~ . %in% empty_values))
  if(nrow(missing_rows) > 0) {
    print('The following rows in the sample meta are not filled out for the sequencing index columns.')
    print(missing_rows) # show the empty rows
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' validate_columns_exist
#'
#' This function checks that a list of columns are present in a dataframe.
#' Columns that were not found in the dataframe are printed out.
#'
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_columns_exist= function(selected_cols, df) {
  # Check that all of selected_columns are in df - base::setdiff(A, B) = A[!A %in% B].
  unmatched_cols= base::setdiff(selected_cols, colnames(df))

  if(length(unmatched_cols) > 0) {
    print('The following columns are missing: ')
    print(unmatched_cols)
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' validate_unique_samples
#'
#' This function checks that a list of columns uniquely identifies every row of a dataframe.
#'
#' @param selected_columns A vector of strings each representing a column name
#' @param df A dataframe to check against
#' @return Boolean
validate_unique_samples= function(selected_columns, df) {
  unique_column_values= df %>% dplyr::distinct(pick(all_of(selected_columns)))
  if(nrow(unique_column_values) != nrow(df)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' check_file_exists
#'
#' This function checks that a given file exists, and if not, stops execution.
#' This is meant to be used at the end of a module to validate that outputs have been generated.
#'
#' @param file_path Path to a file to check
check_file_exists <- function(file_path) {
  if (file.exists(file_path)) {
    message(paste("File ", file_path, " was successfully generated."))
  } else {
    stop(paste("Error: File ", file_path, " does not exist, please check the output of this module."))
  }
}

#' Filter out control barcodes
#'
#' This function filters out control barcodes from a dataframe.
#' tryCatch is used in the event that the `cb_name` column does not exist, it then returns the original dataframe.
#'
#' @param df A dataframe to filter.
#' @return A filtered dataframe.
filter_control_barcodes <- function(df) {
  df %>% dplyr::filter(tryCatch(is.na(cb_name), error = function(e) TRUE))
}

#' Append a print statement to a file to track critical console output
#'
#' This function, when applied to a given print statement, will append the statement to a file in order to
#' track critical console outputs for human review. It will also write the statement to the console as usual.
#'
#' @param print_statement A string to print and append to a file
append_critical_output <- function(statement, out) {
  # If the out/logs directory does not exist, create it
    if (!dir.exists(paste0(out, "/logs"))) {
        dir.create(paste0(out, "/logs"), recursive = TRUE)
    }
  # Convert data frames or lists into a string if necessary
  if (!is.character(statement)) {
    statement <- capture.output(print(statement))
  }
  # Print to console
  cat(statement, sep = "\n")
  # Append to file
  cat(statement, file = paste0(out, "/logs/critical_output.txt"), append = TRUE, sep = "\n")
}

#' Delete existing files when re-running a module
#'
#' This function deletes existing files in a given output directory that match a given pattern.
#' This is useful when re-running a module to ensure that old files do not interfere with
#' the new output.
#' @param out_dir The output directory to check for existing files.
#' @param pattern A pattern to match files against. This can be a regular expression.
delete_existing_files <- function(out_dir, pattern) {
  files_to_delete <- list.files(path = out_dir, pattern = pattern, full.names = TRUE)
  if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
    message(paste("Deleted", length(files_to_delete), "existing files in", out_dir))
    message("Deleted files:")
    print(files_to_delete)
  } else {
    message(paste("No existing files to delete in", out_dir))
  }
}

#' Writes out a table to a csv
#'
#' This function uses readr's write_csv to write out a table to a given path.
#'
#' @param table Table of data
#' @param path Path to the output csv file.
write_out_table = function(table, path) {
  # Defaults should be no rownames and quote only when necessary
  readr::write_csv(x = table, file = path)
}