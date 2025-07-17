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
gl.run.structure <- function(x,
                             ...,
                             plot.out = TRUE,
                             plot_theme = theme_dartR(),
                             exec = "/usr/local/bin/structure",
                             save2tmp = FALSE,
                             verbose = NULL) {
  
  pkg <- "purrr"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }
  
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "Jody",
                   verbose = verbose)
  
  # CHECK DATATYPE
  #datatype <- utils.check.datatype(x, verbose = verbose)
  
  #if (datatype != "SNP") {
  # stop(error(
  #   "You need to provide a SNP genlight object (ploidy=2)!"
  # ))
  #}
  
  if (tools::file_ext(x) == "csv") {
    file <- readr::read_csv(x)
  } else if (tools::file_ext(x) == "xlsx") {
    file <- readxl::read_excel(x)
  } else {
    stop("Input file should be in csv or xlsx format.")
  }
  
  file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
  file <- as.data.frame(file)
  
  library(dplyr)
  
  file[is.na(file)] <- "N"
  file <- file %>%
    mutate(across(everything(), ~ case_when(
      . == "N/A" ~ "N", 
      . == "NA" ~ "N",
      TRUE ~ .x))) %>%
    rename(Ind = 1, Pop = 2)
  
  ind <- as.character(file$Ind)
  pop <- as.character(file$Pop)
  fsnps_geno <- file[, 3:ncol(file)]
  
  fsnps_gen <- adegenet::df2genind(fsnps_geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
  fsnps_gen@pop <- as.factor(file$Pop)
  
  ### ================
  
  gg <- dartR::gi2gl(fsnps_gen, verbose = 0)
  #gg <- fsnps_gen
  
  # check that Structure is installed
  #structure <- Sys.which(exec)
  
  
  # DO THE JOB
  gg <- dartR::utils.structure.genind2gtypes(dartR::gl2gi(gg, verbose = 0)) 
  
  sr <- utils.structure.run(gg, exec = exec, ...)
  
  ev <- dartR::utils.structure.evanno(sr)
  
  pa <- ((ev$plots$mean.ln.k + ev$plots$mean.ln.k) / 
           (ev$plots$ln.ppk + ev$plots$delta.k)) + plot_theme
  
  # PRINTING OUTPUTS
  if (plot.out) {
    suppressMessages(print(pa))
  }
  
  # SAVE INTERMEDIATES TO TEMPDIR
  if (save2tmp & plot.out) {
    # check for '/' in match.call
    mc <- gsub("/", ":", as.character(funname))
    mc <- gsub(":", "-", mc)
    nmc <- gsub("/", "_over_", names(funname))
    nmc <- gsub(":", "-", nmc)
    
    # creating temp file names
    temp_plot <-
      tempfile(pattern = paste0("Plot", paste0(nmc, "_", mc,
                                               collapse = "_")))
    
    # saving to tempdir
    saveRDS(pa, file = temp_plot)
    if (verbose >= 2) {
      message(
        "  Saving the plot in ggplot format to the tempfile as",
        temp_plot
      )
      message(
        "  NOTE: Retrieve output files from tempdir using 
                        gl.list.reports() and gl.print.reports()"
      )
      
    }
  }
  
  # FLAG SCRIPT END
  if (verbose >= 1) {
    message("Completed:", funname, "\n\n")
  }
  
  # RETURN
  return(sr)
  
}

#' @export
get_structure_exec <- function(exec = Sys.which("structure")) {
  return(exec)
}

###=============== VERBOSITY ================###
gl.check.verbosity <- function(x = NULL) {
  # SET VERBOSITY or GET it from global
  if (is.null(x)) {
    if (is.null(options()$dartR_verbose)) {
      verbose <- 2
    } else {
      verbose <- options()$dartR_verbose
    }
  } else {
    if (is.numeric(x) & x >= 0 & x <= 5) {
      verbose <- x
    } else {
      cat(
        warn(
          "Warning: Parameter verbose must be an integer in the range 
                    0 to 5, set to 2\n"
        )
      )
      verbose <- 2
    }
  }
  
  return(verbose)
}

###=============FLAG START=============###
utils.flag.start <- function(func = NULL,
                             build = NULL,
                             verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  if (is.null(func)) {
    stop(paste("Fatal Error: The calling function must be specified.\n"))
  }
  if (verbose >= 1) {
    if (verbose == 5) {
      if (!is.null(build)) {
        message("[dartR] Starting ", func)
      } else {
        message("Starting ", func)
      }
    } else {
      message("Starting ", func)
    }
  }
  invisible(func)
}

### =================== UTILS.STRUCTURE.RUN =============###
utils.structure.run <- function (g,
                                 k.range = NULL,
                                 num.k.rep = 1,
                                 label = NULL,
                                 delete.files = TRUE,
                                 exec = "/usr/local/bin/structure",
                                 ...) {

################################################################
  
    .structureParseQmat <- function (q.mat.txt, 
                                     pops) {
      
      q.mat.txt <- sub("[*]+", "", q.mat.txt)
      q.mat.txt <- sub("[(]", "", q.mat.txt)
      q.mat.txt <- sub("[)]", "", q.mat.txt)
      q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
      cols1to4 <- c("row", "id", "pct.miss", "orig.pop")
      
      strsplit(q.mat.txt, " ") %>% 
        purrr::map(function(q) {
        q <- q[!q %in% c("", " ", ":")] %>% 
          as.character() %>% 
          rbind() %>% 
          as.data.frame(stringsAsFactors = FALSE)
        
        stats::setNames(q, 
                        c(cols1to4,
                          paste("Group", 1:(ncol(q) - 4), sep = ".")))
        
      }) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate_at(dplyr::vars("row",
                                     "pct.miss",
                                     "orig.pop",
                                     dplyr::starts_with("Group.")), as.numeric) %>% 
        dplyr::mutate(orig.pop = if (!is.null(pops)) {
                                                       pops[.data$orig.pop]
                                                     } else {
                                                       .data$orig.pop
                                                     })
    }
    
################################################################
    
    .alleles2integer <- function (g,
                                  min.val = 0) {
      g$data %>% 
        dplyr::group_by(.data$locus) %>% 
        dplyr::mutate(allele = min.val - 1 + as.integer(factor(.data$allele))) %>%
        dplyr::ungroup()
    }    
    
################################################################

    .stackedAlleles <- function (g, 
                                 alleles2integer = FALSE, 
                                 na.val = NULL, 
                                 ...) {
      
      x <- if (alleles2integer){ 
        .alleles2integer(g, ...)
      } else {
        g$data
        }
      
      if (!is.null(na.val)){ 
        x$allele[is.na(x$allele)] <- na.val
      }
      
      x %>%
        dplyr::arrange(.data$id, .data$locus) %>% 
        dplyr::mutate(a = rep(1:g$ploidy, dplyr::n()/g$ploidy)) %>% 
        tidyr::spread(.data$locus, .data$allele) %>%
        dplyr::rename(allele = "a") %>% 
        dplyr::select(.data$id, .data$stratum, .data$allele, dplyr::everything())
    }    
        
####################################################    
    .getFileLabel <- function (g, 
                               label = NULL){
      desc <- g$description
      label <- if (!is.null(label)) {
        label
      }else if (!is.null(desc)) {
        desc
      }else{
        "strataG.gtypes"
        }
      gsub("[[:punct:]]", ".", label)
    }
    
####################################################
    
    structureWrite <- function (g,  
                                label = NULL,  
                                maxpops = 1:(dplyr::n_distinct(g$data$stratum)),
                                burnin = 1000, 
                                numreps = 1000, 
                                noadmix = TRUE, 
                                freqscorr = FALSE, 
                                randomize = TRUE, 
                                seed = 0, 
                                pop.prior = NULL, 
                                locpriorinit = 1, 
                                maxlocprior = 20, 
                                gensback = 2,
                                migrprior = 0.05, 
                                pfrompopflagonly = TRUE, 
                                popflag = NULL, 
                                inferalpha = FALSE,
                                alpha = 1, 
                                unifprioralpha = TRUE, 
                                alphamax = 20,
                                alphapriora = 0.05, 
                                alphapriorb = 0.001)  {
     
      if (!is.null(pop.prior)) {
        if (!pop.prior %in% c("locprior", "usepopinfo")) {
          stop(paste("'pop.prior' must be 'locprior' or 'usepopinfo'."))
        }
      }
      
      if (is.null(popflag)){ 
        popflag <- rep(1, length(unique(g$data$id)))
      }
      
      if (length(popflag) != length(unique(g$data$id))) {
        stop(sprintf("  'popflag' should be the same length as the number of individuals in 'g'."))
      }
      if (!all(popflag %in% c(0, 1))) {
        stop(paste("  All values in 'popflag' must be 0 or 1."))
      }
      
      if (is.null(names(popflag))){ 
        names(popflag) <- unique(g$data$id)
      }
      
      in.file <- ifelse(is.null(label), "data", paste(label, "data", sep = "_"))
      out.file <- ifelse(is.null(label), "out", paste(label, "out", sep = "_"))
      main.file <- ifelse(is.null(label), "mainparams", paste(label, "mainparams", sep = "_"))
      extra.file <- ifelse(is.null(label), "extraparams", paste(label, "extraparams", sep = "_"))
      mat <- .stackedAlleles(g, alleles2integer = TRUE, na.val = -9) %>% 
        dplyr::select(-.data$allele) %>%
        dplyr::mutate(id = gsub(" ", "_", .data$id), 
                      stratum = as.numeric(factor(.data$stratum)), 
                      popflag = popflag[.data$id]) %>% 
        dplyr::select(.data$id, .data$stratum, .data$popflag, dplyr::everything()) %>% 
        as.matrix()
      
      write(paste(sort(unique(g$data[["locus"]])), collapse = " "), file = in.file)
      
      for (i in 1:nrow(mat)) {
        write(paste(mat[i, ], collapse = " "), file = in.file, 
              append = TRUE)
      }
      
      main.params <- c(paste("MAXPOPS", as.integer(maxpops)), 
                       paste("BURNIN", as.integer(burnin)), 
                       paste("NUMREPS", as.integer(numreps)), 
                       paste("INFILE", in.file), 
                       paste("OUTFILE", out.file), 
                       paste("NUMINDS", length(unique(g$data$id))), 
                       paste("NUMLOCI", length(unique(g$data$locus))), 
                       "MISSING -9", "LABEL 1", "POPDATA 1", 
                       "POPFLAG 1", "LOCDATA 0", "PHENOTYPE 0", 
                       "EXTRACOLS 0", "MARKERNAMES 1")
      
      main.params <- paste("#define", main.params)
      
      write(main.params, file = main.file)
      
      extra.params <- c(paste("NOADMIX", as.integer(noadmix)), 
                        paste("FREQSCORR", as.integer(freqscorr)),
                        paste("INFERALPHA", as.integer(inferalpha)), 
                        paste("ALPHA", as.numeric(alpha)), 
                        "FPRIORMEAN 0.01", "FPRIORSD 0.05", "LAMBDA 1.0", 
                        paste("UNIFPRIORALPHA", as.integer(unifprioralpha)), 
                        paste("ALPHAMAX", as.numeric(alphamax)), 
                        paste("ALPHAPRIORA", as.numeric(alphapriora)), 
                        paste("ALPHAPRIORB", as.numeric(alphapriorb)), 
                        "COMPUTEPROB 1", 
                        paste("ADMBURNIN", max(0, as.integer(burnin/2))), 
                        "ALPHAPROPSD 0.025", "STARTATPOPINFO 0", 
                        paste("RANDOMIZE", as.integer(randomize)), 
                        paste("SEED", as.integer(seed)),
                        "METROFREQ 10", "REPORTHITRATE 0")
      
      if (!is.null(pop.prior)) {
        pop.prior <- tolower(pop.prior)
        prior.params <- if (pop.prior == "locprior") {
          c("LOCPRIOR 1", "LOCISPOP 1", paste("LOCPRIORINIT", 
                                              locpriorinit), paste("MAXLOCPRIOR", maxlocprior))
        } else if (pop.prior == "usepopinfo") {
          c("USEPOPINFO 1", paste("GENSBACK", trunc(gensback)), 
            paste("MIGRPRIOR", migrprior), paste("PFROMPOPFLAGONLY", 
                                                 as.integer(pfrompopflagonly)))
        }
        extra.params <- c(extra.params, prior.params)
      }
      
      extra.params <- extra.params[!is.na(extra.params)]
      extra.params <- paste("#define", extra.params)
      
      write(extra.params, file = extra.file)
      
      invisible(list(files = c(data = in.file, 
                               mainparams = main.file, 
                               extraparams = extra.file, 
                               out = out.file),
                     pops = sort(unique(g$data$stratum))))
    }
    
###########################################################
    structureRead <- function(file,
                              pops = NULL) {
      
      if (!file.exists(file)) {
        stop(sprintf(paste("the file '", file, "' can't be found.", 
                   sep = "")))
      }
      
      result <- scan(file, "character", quiet = TRUE)
      loc <- grep("Estimated", result, ignore.case = FALSE, 
                  value = FALSE)
      est.ln.prob <- as.numeric(result[loc[1] + 6])
      loc <- grep("likelihood", result, ignore.case = FALSE, 
                  value = FALSE)
      mean.lnL <- as.numeric(result[loc[1] + 2])
      var.lnL <- as.numeric(result[loc[2] + 2])
      loc <- grep("MAXPOPS", result, value = F)
      maxpops <- result[loc]
      maxpops <- sub("MAXPOPS=", "", maxpops)
      maxpops <- as.integer(sub(",", "", maxpops))
      loc <- grep("GENSBACK", result, value = F)
      gensback <- result[loc]
      gensback <- sub("GENSBACK=", "", gensback)
      gensback <- as.integer(sub(",", "", gensback))
      smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, 
                var.lnL = var.lnL)
      result <- scan(file, "character", sep = "\n", 
                     quiet = TRUE)
      first <- grep("(%Miss)", result, value = FALSE) + 1
      last <- grep("Estimated Allele", result, value = FALSE) - 
        1
      tbl.txt <- result[first:last]
      tbl.txt <- sub("[*]+", "", tbl.txt)
      tbl.txt <- sub("[(]", "", tbl.txt)
      tbl.txt <- sub("[)]", "", tbl.txt)
      tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
      prior.lines <- grep("[|]", tbl.txt)
      
      no.prior <- if (length(prior.lines) < length(tbl.txt)) {
        no.prior.q.txt <- if (length(prior.lines) == 0){ 
          tbl.txt
        }else{
          tbl.txt[-prior.lines]
        }
        .structureParseQmat(no.prior.q.txt, pops)
      } else { 
        NULL
      }
      
      if (maxpops == 1) {
        no.prior$row <- NULL
        return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
      }
      
      has.prior <- if (length(prior.lines) > 0) {
        prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
        prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
        df <- .structureParseQmat(prior.q.txt, pops)
        
        prior.anc <- purrr::map(prior.txt, function(x) {
          anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
          rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat),  sep = ".")
          colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
          
          x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
            y <- y[y != ""]
            y[-1]
          })
          
          for (i in 1:ncol(x)) {
            pop <- as.numeric(x[1, i])
            anc.mat[pop, ] <- as.numeric(x[-1, i])
          }
          anc.mat
        }) %>% stats::setNames(df$id)
        prob.mat <- t(sapply(1:nrow(df), function(i) {
          pop.probs <- rowSums(prior.anc[[i]])
          pop.probs[is.na(pop.probs)] <- df$Group.1[i]
          pop.probs
        }))
        colnames(prob.mat) <- paste("Group", 1:ncol(prob.mat), 
                                    sep = ".")
        df$Group.1 <- NULL
        df <- cbind(df, prob.mat)
        list(df = df, prior.anc = prior.anc)
      }else{
         NULL
      }
      
      has.prior.df <- if(is.null(has.prior)){ 
        NULL
      }else{ 
        has.prior$df
        }
      
      q.mat <- rbind(no.prior, has.prior.df)
      q.mat <- q.mat[order(q.mat$row), ]
      q.mat$row <- NULL
      rownames(q.mat) <- NULL
      q.mat[, -(1:3)] <- t(apply(q.mat[, -(1:3)], 1, function(i) i/sum(i)))
      prior.anc <- if (is.null(has.prior)){ 
        NULL
      }else {
        has.prior$prior.anc
        }
      
      list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
    }

###########################################################
    
    # Create output base directory in a safe sandbox
    base_label <- file.path(tempdir(), gsub("[[:space:]:]", ".", g$description), "structureRun")
    run_label <- file.path(base_label, "run")
    dir.create(run_label, recursive = TRUE, showWarnings = FALSE)
    
    # Verify folder creation
    if (!utils::file_test("-d", run_label)) {
      stop(sprintf("'%s' is not a valid folder.", run_label))
    }
    
    # Write permission check
    test_file <- file.path(run_label, "test_write_check.txt")
    tryCatch({
      writeLines("Write permission confirmed.", test_file)
      message("Write test succeeded: ", test_file)
    }, error = function(e) {
      stop(sprintf("Write test failed: %s\nDirectory '%s' is not writeable.", e$message, run_label))
    })
    
    # 🧬 Set up K range and replicate structure
    if (is.null(k.range)) {
      k.range <- 1:(dplyr::n_distinct(g$data$stratum))
    }
    rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
    rownames(rep.df) <- paste(run_label, ".k", rep.df$k, ".r", rep.df$rep, sep = "")
    
    # Run STRUCTURE for each replicate
    out.files <- lapply(rownames(rep.df), function(x) {
      sw.out <- structureWrite(g, label = x, maxpops = rep.df[x, "k"])
      files <- sw.out$files
      
      cmd <- paste(exec,
                   "-m", files["mainparams"],
                   "-e", files["extraparams"],
                   "-i", files["data"],
                   "-o", files["out"])
      
      message("Running STRUCTURE with command:\n", cmd)
      err.code <- system(cmd)
      message("STRUCTURE exited with code: ", err.code)
      
      # Inspect output directory contents
      output_dir <- dirname(files["out"])
      message("Checking output directory: ", output_dir)
      print(list.files(output_dir, full.names = TRUE))
      
      # Check for output file and fallback
      base_out <- files["out"]
      alt_out <- paste0(base_out, "_f")
      
      if (file.exists(alt_out)) {
        message("Using STRUCTURE output with '_f' suffix.")
        files["out"] <- alt_out
      } else if (file.exists(base_out)) {
        message("Using fallback STRUCTURE output without '_f' suffix.")
        files["out"] <- base_out
      } else {
        stop(sprintf("Neither expected STRUCTURE output file '%s' nor fallback '%s' found.\nCheck if STRUCTURE executed correctly or verify input formatting.", alt_out, base_out))
      }
      
      # Read STRUCTURE results
      result <- structureRead(files["out"], sw.out$pops)
      if (file.exists("seed.txt")) file.remove("seed.txt")
      files <- if (delete.files) NULL else files
      
      result <- c(result, list(files = files, label = basename(x)))
      fname <- paste(x, ".ws.rdata", sep = "")
      save(result, file = fname)
      fname
    })
    
    # Compile final results
    run.result <- lapply(out.files, function(f) {
      result <- NULL
      load(f)
      result
    })
    names(run.result) <- sapply(run.result, function(x) x$label)
    class(run.result) <- c("structure.result", class(run.result))
    
    # Cleanup (optional)
    if (delete.files) unlink(base_label, recursive = TRUE, force = TRUE)
    
    run.result }