#' @title Runs a STRUCTURE analysis using a genlight object
#'
#' @export
run_structure <- function(
    input_path,
    k.range = 1:5,
    numrep = 3,
    burnin = 1000,
    numreps = 1000,
    noadmix = FALSE,
    structure_path = "/usr/local/bin/structure",
    output_base_dir = tempdir(),
    plot_dir = file.path(tempdir(), "evanno_plots"),              
    plot.out = TRUE,
    delete.files = TRUE
) {
  dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  genind_obj <- to_genind(input_path)
  str_file <- genind_to_structure_v2(genind_obj, file = "structure_input.str", dir = output_base_dir)
  
  check_structure_format(str_file)
  
  result <- running_structure(
    input_file = str_file,
    k.range = k.range,
    numrep = numrep,
    burnin = burnin,
    numreps = numreps,
    noadmix = noadmix
  )
  
  return(result)
}

##############

to_genind <- function(input_path) {
  ext <- tools::file_ext(input_path)
  raw <- if (ext == "csv") {
    readr::read_csv(input_path)
  } else if (ext == "xlsx") {
    readxl::read_excel(input_path)
  } else {
    stop("Unsupported file type. Please provide a CSV or XLSX.")
  }
  
  raw <- as.data.frame(lapply(raw, function(x) gsub("|", "/", x, fixed = TRUE)))
  raw[is.na(raw)] <- "N"
  raw <- dplyr::mutate(raw, across(everything(), ~ dplyr::case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))) |> 
    dplyr::rename(Ind = 1, Pop = 2)
  
  ind <- as.character(raw$Ind)
  pop <- as.character(raw$Pop)
  geno <- raw[, 3:ncol(raw)]
  
  genind_obj <- adegenet::df2genind(geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
  genind_obj@pop <- as.factor(pop)
  return(genind_obj)
}

############

genind_to_structure_v2 <- function(genind_obj, file = "structure_input.str", include_pop = TRUE, dir = tempdir()) {
  out_path <- file.path(dir, file)
  # Get basic info
  ind <- adegenet::indNames(genind_obj)
  pop <- if (include_pop) as.character(genind_obj@pop) else rep(1, length(ind))
  ploidy <- max(genind_obj@ploidy)
  loci <- adegenet::locNames(genind_obj)
  n_loc <- length(loci)
  
  # Convert allele matrix to integer format
  allele_matrix <- matrix(-9, nrow = length(ind), ncol = n_loc * ploidy)
  colnames(allele_matrix) <- paste(rep(loci, each = ploidy), paste0(".a", 1:ploidy), sep = "")
  
  for (i in seq_along(ind)) {
    for (j in seq_along(loci)) {
      allele_values <- genind_obj@tab[i, grep(paste0("^", loci[j], "\\."), colnames(genind_obj@tab))]
      allele_values[is.na(allele_values)] <- -9
      allele_matrix[i, ((j - 1) * ploidy + 1):(j * ploidy)] <- allele_values
    }
  }
  
  final_data <- data.frame(ID = ind, POP = pop, allele_matrix, stringsAsFactors = FALSE)
  
  write.table(final_data, file = out_path, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  return(out_path)
}

###########

check_structure_format <- function(file_path, ploidy = 2, verbose = TRUE) {
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  df <- tryCatch(
    read.table(file_path, sep = " ", header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE),
    error = function(e) stop("Unable to read file: ", e$message)
  )
  
  # Minimum expected columns: ID, POP, alleles
  expected_cols <- 2 + ploidy * floor((ncol(df) - 2) / ploidy)
  
  if (ncol(df) < expected_cols) {
    stop("The file has fewer columns than expected. Check ploidy or locus formatting.")
  }
  
  issues <- list()
  
  # Check for missing IDs or POP assignments
  if (any(df[[1]] == "" | is.na(df[[1]]))) issues$id <- "Missing individual IDs"
  if (any(df[[2]] == "" | is.na(df[[2]]))) issues$pop <- "Missing population assignments"
  
  # Check allele columns: should be integers or -9
  allele_matrix <- df[, -(1:2)]
  invalid <- sapply(allele_matrix, function(col) {
    all(grepl("^(-9|\\d+)$", col)) == FALSE
  })
  
  if (any(invalid)) {
    issues$alleles <- paste(sum(invalid), "allele columns have invalid values (must be integers or -9).")
  }
  
  # Summary report
  if (verbose) {
    message("✅ STRUCTURE Format Check — Report")
    message("File: ", file_path)
    message("Individuals: ", nrow(df))
    message("Populations: ", length(unique(df[[2]])))
    message("Loci estimated: ", (ncol(df) - 2) / ploidy)
    if (length(issues) == 0) {
      message("All format checks passed.")
    } else {
      for (item in issues) message("❌ ", item)
    }
  }
  
  return(invisible(length(issues) == 0))
}

###########

running_structure <- function(
    input_file,
    k.range,
    numrep,
    burnin,
    numreps,
    noadmix = FALSE,
    structure_path = "/usr/local/bin/structure",
    output_dir = tempdir(),
    plot_dir = file.path(output_dir, "evanno_plots"),
    plot.out = TRUE,
    save_plots_dir = tempdir()
) {

  # Validate K range
  if (length(k.range) < 2) stop("Provide at least two K values for Evanno analysis.")
  if (!dir.exists(save_plots_dir)) dir.create(save_plots_dir, recursive = TRUE)
  
  numinds <- length(readLines(input_file))
  numloci <- ncol(read.table(input_file, header = FALSE, sep = " ")) - 2
  base_label <- tools::file_path_sans_ext(basename(input_file))
  
  rep.df <- expand.grid(rep = 1:numrep, k = k.range)
  rownames(rep.df) <- paste0(base_label, ".k", rep.df$k, ".r", rep.df$rep)
  
  out_files <- lapply(rownames(rep.df), function(run_label) {
    k <- rep.df[run_label, "k"]
    out_path <- file.path(output_dir, run_label)
    
    mainparams <- paste0(out_path, "_mainparams")
    extraparams <- paste0(out_path, "_extraparams")
    
    writeLines(paste("#define", c(
      paste("MAXPOPS", k),
      paste("BURNIN", burnin),
      paste("NUMREPS", numreps),
      paste("INFILE", input_file),
      paste("OUTFILE", out_path),
      paste("NUMINDS", numinds),
      paste("NUMLOCI", numloci),
      "MISSING -9", "LABEL 1", "POPDATA 1",
      "POPFLAG 0", "LOCDATA 0", "PHENOTYPE 0",
      "EXTRACOLS 0", "MARKERNAMES 0"
    )), con = mainparams)
    
    # change alpha value divide 1 with max K value #######################################################
    
    writeLines(paste("#define", c(
      paste("NOADMIX", ifelse(noadmix, 1, 0)), 
      "FREQSCORR 1", "INFERALPHA 1",
      "ALPHA 1.0", "COMPUTEPROB 1",
      paste("SEED", sample(1e6, 1)),
      "METROFREQ 10", "REPORTHITRATE 0",
      "STARTATPOPINFO 0"
    )), con = extraparams)
    
    cmd <- paste(structure_path,
                 "-m", mainparams,
                 "-e", extraparams,
                 "-i", input_file,
                 "-o", out_path)
    
    message("Running STRUCTURE: ", cmd)
    log <- tryCatch(system(cmd, intern = TRUE), error = function(e) e$message)
    writeLines(log, paste0(out_path, "_log.txt"))
    
    final_out <- paste0(out_path, "_f")
    if (!file.exists(final_out) && file.exists(out_path)) final_out <- out_path
    if (!file.exists(final_out)) {
      warning("Missing output file for run: ", run_label)
      return(NULL)
    }
    
    result <- tryCatch(structureRead(final_out), error = function(e) {
      warning("Failed to parse STRUCTURE output: ", run_label)
      return(NULL)
    })
    result$label <- run_label
    result
  })
  
  run.result <- Filter(Negate(is.null), out_files)
  message("Successful STRUCTURE runs: ", length(run.result))
  names(run.result) <- sapply(run.result, `[[`, "label")
  class(run.result) <- c("structure.result", class(run.result))
  
  if (length(run.result) < 3) stop("Evanno analysis needs at least three successful runs.")
  
  ev <- utils.structure.evanno(run.result, plot = FALSE)
  
  # 
  lapply(names(ev$plots), function(pname) {
    plot_obj <- ev$plots[[pname]]
    png_path <- file.path(save_plots_dir, paste0("Evanno_", pname, ".png"))
    grDevices::png(filename = png_path, width = 1600, height = 1200, res = 200)
    print(plot_obj)
    grDevices::dev.off()
    message("Saved PNG plot: ", png_path)
  })
  
  if (plot.out && "delta.k" %in% names(ev$plots)) {
    suppressMessages(print(ev$plots$delta.k))
  }
  
  if (delete.files) unlink(output_dir, recursive = TRUE)
  
  invisible(list(
    results = run.result,
    evanno = ev,
    plot.paths = list.files(save_plots_dir, full.names = TRUE)
  ))
}

########