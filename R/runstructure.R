#' @title Runs a STRUCTURE analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a STRUCTURE analysis based on
#' functions from \code{strataG}
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param ... Parameters to specify the STRUCTURE run (check \code{structureRun}
#'  within strataG.
#' for more details). Parameters are passed to the \code{structureRun} function.
#' For example you need to set the k.range and the type of model you would like
#' to run (noadmix, locprior) etc. If those parameter names do not tell you
#' anything, please make sure you familiarize with the STRUCTURE program
#' (Pritchard 2000).
#' @param exec Full path and name+extension where the structure executable is
#' located. E.g. \code{'c:/structure/structure.exe'} under Windows. For Mac and
#' Linux it might be something like \code{'./structure/structure'} if the
#' executable is in a subfolder 'structure' in your home directory
#' [default working directory "."].
#' @param plot.out Create an Evanno plot once finished. Be aware k.range needs
#' to be at least three different k steps [default TRUE].
#' @param plot_theme Theme for the plot. See details for options
#' [default theme_dartR()].
#' @param save2tmp If TRUE, saves any ggplots and listings to the session
#' temporary directory (tempdir) [default FALSE].
#' @param verbose Set verbosity for this function (though structure output
#' cannot be switched off currently) [default NULL]
#' @details The function is basically a convenient wrapper around the beautiful
#' strataG function \code{structureRun} (Archer et al. 2016). For a detailed
#' description please refer to this package (see references below).
#' @author Bernd Gruber (Post to \url{https://groups.google.com/d/forum/dartr})
#' @export
structure_analysis <- function(
    input_path,
    output_dir = tempdir(),
    save_plots_dir = file.path(output_dir, "evanno_plots"),
    structure_exec = "/usr/local/bin/structure",
    k.range = 1:5,
    numrep = 3,
    burnin = 1000,
    numreps = 1000,
    plot.out = TRUE,
    delete.files = TRUE
) {
  # Load input file
  file_ext <- tools::file_ext(input_path)
  if (file_ext == "csv") {
    raw_data <- readr::read_csv(input_path)
  } else if (file_ext == "xlsx") {
    raw_data <- readxl::read_excel(input_path)
  } else {
    stop("Unsupported file type. Use .csv or .xlsx.")
  }
  
  # Clean genotype data
  raw_data <- lapply(raw_data, function(x) gsub("|", "/", x, fixed = TRUE)) |> as.data.frame()
  raw_data[is.na(raw_data)] <- "N"
  
  raw_data <- dplyr::mutate(
    raw_data,
    across(everything(), ~ dplyr::case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))
  ) |> dplyr::rename(Ind = 1, Pop = 2)
  
  ind <- as.character(raw_data$Ind)
  pop <- as.character(raw_data$Pop)
  fsnps_geno <- raw_data[, 3:ncol(raw_data)]
  
  fsnps_gen <- adegenet::df2genind(fsnps_geno, ind.names = ind, pop = pop,
                                   sep = "/", NA.char = "N", ploidy = 2, type = "codom")
  fsnps_gen@pop <- as.factor(pop)
  fsnps_gen_sub <- fsnps_gen  # If subsetting is needed later
  
  # Generate STRUCTURE input file
  structure_input <- file.path(output_dir, "structure_input.str")
  genind_to_structure_clean(fsnps_gen_sub, file = structure_input)
  
  # Run STRUCTURE across user-defined K range and analyze Evanno
  result <- run_structure_with_evanno(
    input_file = structure_input,
    k.range = k.range,
    numrep = numrep,
    burnin = burnin,
    numreps = numreps,
    structure_path = structure_exec,
    output_dir = file.path(output_dir, "structure_runs"),
    save_plots_dir = save_plots_dir,
    plot.out = plot.out,
    delete.files = delete.files
  )
  
  return(result)
}


run_structure <- function(
    input_file,
    k.range = 1:5,                
    numrep = 3,                   
    burnin = 1000,                
    numreps = 1000,               
    structure_path = "/usr/local/bin/structure",
    output_dir = tempdir(),
    save_plots_dir = tempdir(),
    delete.files = TRUE,
    plot.out = TRUE
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
    
    writeLines(paste("#define", c(
      "NOADMIX 0", "FREQSCORR 1", "INFERALPHA 1",
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
  names(run.result) <- sapply(run.result, `[[`, "label")
  class(run.result) <- c("structure.result", class(run.result))
  
  if (length(run.result) < 3) stop("Evanno analysis needs at least three successful runs.")
  
  ev <- utils.structure.evanno(run.result, plot = FALSE)
  
  # 📸 Save Evanno plots as PNG
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