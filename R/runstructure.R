#' @title Run STRUCTURE analysis
#' @description Adapted from the dartR package which was based on strataG. Build for docker deployment.
#' @export
run_structure <- function(
    file,
    k.range = 1:5,
    num.k.rep = 3,
    burnin = 1000,
    numreps = 1000,
    noadmix = FALSE,
    phased = FALSE,
    ploidy = 2,
    linkage = FALSE,
    #structure_path = "/usr/local/bin/structure",
    output_base_dir = tempdir(), # changed to structure_files
    plot_dir = file.path(output_base_dir, "evanno_plots"),
    clumpp_plots = file.path(output_base_dir, "str_plots")
    #plot.out = TRUE,
    #delete.files = TRUE
){
  dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(clumpp_plots, recursive = TRUE, showWarnings = FALSE)
  
  # Get STRUCTURE path
  structure_path <- Sys.which("structure")
  if(structure_path == ""){
    stop("STRUCTURE not available in the PATH")
  }
  
  ext <- tools::file_ext(file)
  
  raw <- if (ext == "csv") {
    readr::read_csv(file)
  } else if (ext == "xlsx") {
    readxl::read_excel(file)
  } else {
    stop("Unsupported file type. Please provide a CSV or XLSX.")
  }
  
  library(dplyr) # to make sure it works
  
  # Homogenize file format
  # check allele format
  raw <- as.data.frame(lapply(raw, function(x) gsub("|", "/", x, fixed = TRUE)))
  raw[is.na(raw)] <- "N"
  raw <- dplyr::mutate(raw, across(everything(), ~ dplyr::case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))) |> 
    dplyr::rename(Ind = 1, Pop = 2)
  
  # For Plotting
  populations_df <- raw[,1:2]
  colnames(populations_df) <- c("Label", "Population")
  populations_df$Label <- rownames(populations_df)
  
  ### Change pop to numeric - for STRUCTURE
  pop_df <- as.data.frame(raw$Pop) %>%
    dplyr::rename(pops = 1)
  
  # Get total no. of pops
  pop_df_unique <- as.data.frame(pop_df[!duplicated(pop_df), ]) %>%
    dplyr::rename(pops = 1)
  # add row names as numbers
  pop_df_unique$num <- rownames(pop_df_unique)
  write.csv(pop_df_unique, file = "population_order.csv")
  
  # replace the pops in the original df (pop_df) with the numbers
  pop_df_corr <- left_join(pop_df, pop_df_unique, by = "pops") 
  pops <- pop_df_corr$num
  
  ### Change Ind to numeric
  ind_only <- as.data.frame(raw$Ind)
  ind_only$num <- rownames(ind_only)
  
  ind <- as.character(ind_only$num)
  pop <- as.character(pops)
  geno <- raw[, 3:ncol(raw)]
  
  genind_obj <- adegenet::df2genind(geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
  genind_obj@pop <- as.factor(pop)
  
  str.file <- to_structure(genind_obj, file = "structure_input.str", dir = output_base_dir)
  
  result <- running_structure(
    input_file = str.file,
    k.range = k.range,
    num.k.rep = num.k.rep,
    burnin = burnin,
    numreps = numreps,
    noadmix = noadmix,
    structure_path = structure_path
  )
  
  #devtools::install_github("sa-lee/starmie")
  str.dir <- output_base_dir
  str.files <- list.files(str.dir, pattern = "\\_f$", full.names = TRUE)
  str.data <- lapply(str.files, starmie::loadStructure)
  
  for (i in str.data){
    #path_plot <- paste(clumpp_plots, i$K)
    file_name <- paste0(clumpp_plots, "/", i$K, "plot.png")
    plotQ(i, populations_df, outfile = file_name)
  }

}


## Adapted from the dartR package
to_structure <- function(genind_obj, file = "structure_input.str", include_pop = TRUE, dir = output_base_dir) {
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

## Directly from the dartR package, revised sections with comments
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


## Directly from the dartR package, revised sections with comments
structureRead <- function(file,
                          pops = NULL) {
  
  if (!file.exists(file)) {
    stop(error(paste("the file '", file, "' can't be found.", 
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
  
  # Initial section's error: no function to return from, jumping to top level
  # Resolved by wrapping in a function 'maxpops_check'
  maxpops_check <- function(maxpops, smry, no.prior){   
    if (maxpops == 1) {
      no.prior$row <- NULL
      return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
    }}
  
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
  q.mat <- q.mat[order(q.mat$id), ] # removed q.mat$row since row is nonexistent
  #q.mat$row <- NULL
  rownames(q.mat) <- NULL
  
  # Initial single-line code did not work, revised as ff::
  subset <- q.mat[,-(1:3)]
  subset <- as.data.frame(lapply(subset, as.numeric))
  normalized <- subset/rowSums(subset)
  q.mat <- cbind(q.mat[,1:3], normalized)
  
  #q.mat[, -c(1:3)] <- t(apply(q.mat[, -c(1:3)], 1, function(i) i/sum(i)))
  
  prior.anc <- if (is.null(has.prior)){ 
    NULL
  }else {
    has.prior$prior.anc
  }
  
  list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
}


## Directly from the dartR package, revised sections with comments
utils.structure.evanno <- function (sr, plot = TRUE) 
{
  if (!"structure.result" %in% class(sr)) {
    stop(error("'sr' is not a result from 'structure.run'."))
  }
  k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
  
  
  #f (length(k.tbl) < 3) 
  # stop(error("must have at least two values of k."))
  sr.smry <- t(sapply(sr, function(x) x$summary))
  ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[,"k"], mean)
  sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[,"k"], sd)
  
  ln.pk <- diff(ln.k)
  ln.ppk <- abs(diff(ln.pk))
  delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
    abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1])/sd.ln.k[i]
  })
  df <- data.frame(k = as.numeric(names(ln.k)), 
                   reps = as.numeric(table(sr.smry[, "k"])),
                   mean.ln.k = as.numeric(ln.k), 
                   sd.ln.k = as.numeric(sd.ln.k), 
                   ln.pk = c(NA, ln.pk), 
                   ln.ppk = c(NA, ln.ppk, NA), 
                   delta.k = c(NA, delta.k, NA))
  
  rownames(df) <- NULL
  df$sd.min <- df$mean.ln.k - df$sd.ln.k
  df$sd.max <- df$mean.ln.k + df$sd.ln.k
  plot.list <- list(mean.ln.k = 
                      ggplot2::ggplot(df, 
                                      ggplot2::aes_string(x = "k", 
                                                          y = "mean.ln.k")) + 
                      ggplot2::ylab("mean LnP(K)") + 
                      ggplot2::geom_segment(ggplot2::aes_string(x = "k", 
                                                                xend = "k", 
                                                                y = "sd.min", 
                                                                yend = "sd.max")), 
                    ln.pk = ggplot2::ggplot(df[!is.na(df$ln.pk), ], 
                                            ggplot2::aes_string(x = "k", 
                                                                y = "ln.pk")) +
                      ggplot2::ylab("LnP'(K)"), 
                    ln.ppk = ggplot2::ggplot(df[!is.na(df$ln.ppk), ],
                                             ggplot2::aes_string(x = "k", 
                                                                 y = "ln.ppk")) +
                      ggplot2::ylab("LnP''(K)"))
  
  if (!all(is.na(df$delta.k))) {
    plot.list$delta.k <- ggplot2::ggplot(df[!is.na(df$delta.k), 
    ], ggplot2::aes_string(x = "k", y = "delta.k")) + 
      ggplot2::ylab(expression(Delta(K)))
  }
  for (i in 1:length(plot.list)) {
    plot.list[[i]] <- plot.list[[i]] + ggplot2::geom_line() + 
      ggplot2::geom_point(fill = "white", shape = 21, 
                          size = 3) + ggplot2::xlim(c(1, max(df$k))) + 
      ggplot2::theme(axis.title.x = ggplot2::element_blank())
  }
  if (plot) {
    p <- plot.list %>% purrr::map(function(x) {
      ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
    })
    maxWidth <- do.call(grid::unit.pmax, purrr::map(p, function(x) x$widths[2:3]))
    for (i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
    p$bottom <- "K"
    p$ncol <- 2
    do.call(gridExtra::grid.arrange, p)
  }
  df$sd.min <- df$sd.max <- NULL
  invisible(list(df = df, plots = plot.list))
}


### Adapted from the dartR package

running_structure <- function(input_file,
                              k.range,
                              num.k.rep,
                              burnin,
                              numreps,
                              noadmix = FALSE,
                              phased = FALSE,
                              ploidy = 2,
                              linkage = FALSE,
                              structure_path = structure_path,
                              output_dir = tempdir(), # replaced "structure_files" with tempdir()
                              plot_dir = file.path(output_dir, "evanno_plots")){
  
  
  # Validate K range
  if (length(k.range) < 2) stop("Provide at least two K values for Evanno analysis.")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  numinds <- length(readLines(input_file))
  numloci <- (ncol(read.table(input_file, header = FALSE, sep = " ")) - 2)/2
  
  base_label <- tools::file_path_sans_ext(basename(input_file))
  
  rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
  rep.df$run <- paste0(base_label, ".k", rep.df$k, ".r", rep.df$rep)
  
  
  out_files <- lapply(1:nrow(rep.df), function(run_label) {
    #k <- rep.df[i, "k"]
    #####################
    out_path <- file.path(output_dir, rep.df[run_label, "run"])
    
    mainparams <- paste0(output_dir, "/mainparams")
    extraparams <- paste0(output_dir, "/extraparams")
    
    #mainparams <- paste0(dir, "mainparams")
    #extraparams <- paste0(dir, "extraparams")
    
    if(is.null(ploidy)){
      ploidy = "2"
    } else {
      ploidy = ploidy
    }
    
    writeLines(paste("#define", c(
      paste("MAXPOPS", rep.df[run_label, "k"]),
      paste("BURNIN", burnin),
      paste("NUMREPS", numreps),
      paste("INFILE", input_file),
      #paste("OUTFILE", out_path),
      paste("NUMINDS", numinds),
      paste("NUMLOCI", numloci),
      paste("PLOIDY", ploidy),
      "MISSING -9", 
      "ONEROWPERIND 1",
      "LABEL 1", 
      "POPDATA 1",
      "POPFLAG 0", 
      "LOCDATA 0", 
      "PHENOTYPE 0",
      "EXTRACOLS 0", 
      "MARKERNAMES 0",
      "RECESSIVEALLELES 0",
      "MAPDISTANCES 0",
      paste("PHASED", ifelse(phased, 1, 0)),
      "MARKOVPHASE 0",
      "NOTAMBIGUOUS -999"
    )), con = mainparams)
    
    # change alpha value divide 1 with max K value #######################################################
    alpha = 1/max(k.range)
    
    writeLines(paste("#define", c(
      paste("NOADMIX", ifelse(noadmix, 1, 0)), 
      paste("LINKAGE", ifelse(linkage, 1, 0)),
      "USEPOPINFO 0",
      "FREQSCORR 1",
      "ONEFST 0",
      "INFERALPHA 1",
      "POPALPHAS 0",
      paste("ALPHA", alpha),
      "INFERLAMBDA 0",
      "POPSPECIFICLAMBDA 0",
      "LAMBDA 1.0",
      "FPRIORMEAN 0.01",
      "FPRIORSD 0.05",
      "UNIFPRIORALPHA 1",
      "ALPHAMAX 10.0",
      "ALPHAPRIORA 1.0",
      "ALPHAPRIORB 2.0",
      "LOG10RMIN -4.0",
      "LOG10RMAX 1.0",
      "LOG10RSTAT -2.0",
      "GENSBACK 2",
      "MIGRPRIOR 0.01",
      "PFROMPOPFLAGONLY 0",
      "LOCISPOP 1",
      "LOCPRIORINIT 1.0",
      "MAXLOCPRIOR 20.0",
      "PRINTNET 1",
      "PRINTLAMBDA 1",
      "PRINTQSUM 1",
      "SITEBYSITE 0",
      "PRINTQHAT 0",
      "UPDATEFREQ 100",
      "PRINTLIKES 0",
      "INTERMEDSAVE 0",
      "ECHODATA 1",
      "ANCESTDIST 0",
      "NUMBOXES 1000",
      "ANCESTPINT 0.90",
      "COMPUTEPROB 1",
      "ADMBURNIN 500",
      "ALPHAPROPSD 0.025",
      "STARTATPOPINFO 0",
      "RANDOMIZE 0",
      "SEED 2245",
      "METROFREQ 10",
      "REPORTHITRATE 0"
    )), con = extraparams)
    
    #cmd <- paste(structure_path,
    #             "-i", input_file,
    #             "-m", mainparams,
    #             "-e", extraparams,
    #             "-o", out_path)
    
    # for(i in 1:nrow(rep.df)){
    cmd <- paste(structure_path,
                 "-i", input_file,
                 "-K", rep.df[run_label, "k"],
                 "-m", mainparams,
                 "-e", extraparams,
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
    png_path <- file.path(plot_dir, paste0("Evanno_", pname, ".png"))
    grDevices::png(filename = png_path, width = 1600, height = 1200, res = 200)
    print(plot_obj)
    grDevices::dev.off()
    message("Saved PNG plot: ", png_path)
  })
  
  #if (plot.out && "delta.k" %in% names(ev$plots)) {
  suppressMessages(print(ev$plots$delta.k))
  #}
  
  #if (delete.files) unlink(output_dir, recursive = TRUE)
  
  invisible(list(
    results = run.result,
    evanno = ev,
    plot.paths = list.files(plot_dir, full.names = TRUE)
  ))
}

# Plot STRUCTURE
plotQ <- function(qmat, populations_df, outfile = outfile) {
  
  # Revised to be compatible with large list of matrices
  facet = FALSE
  K <- qmat$K
  Label <- seq_len(nrow(qmat$ancest_df))
  
  clusters <- qmat$ancest_df %>%
    as.data.frame() %>%
    select(contains("Cluster"))
  colnames(clusters) <- gsub(" ", ".", colnames(clusters))
  
  df <- data.frame(ID = as.data.frame(Label), clusters)
  
  if (is.null(populations_df)){
    #Generate a plot without any family information
    #Q_melt <- reshape2::melt(Q, variable.name="Cluster")
    Q_melt <- stats::reshape(
      df, 
      varying = list(names(clusters)),
      v.names = "Value",
      timevar = "Cluster",
      times = paste0("K", seq_len(K)),
      direction = "long"
    )
    
    colnames(Q_melt) <- c("Label", "Cluster", "Value")
  } else{
    if (length(populations_df$Label) != length(rownames(df)))
      stop("Unequal sample size.")
    colnames(populations_df) <- c("Label", "Population")
    Q_merge <- merge(populations_df, df, by="Label")
    #Q_melt <- melt(Q_merge, id.vars=c("Label", "Population"), variable.name="Cluster")
    
    Q_melt <- stats::reshape(
      Q_merge, 
      varying = list(names(clusters)),
      v.names = "Value",
      timevar = "Cluster",
      times = paste0("K", seq_len(K)),
      direction = "long"
    )
  }
  
  #Generate plot
  # Updated to proper ggplot syntax
  Q_melt <- Q_melt[order(Q_melt$Cluster),]
  Q_melt$Label <- factor(Q_melt$Label)
  
  gg <- ggplot(Q_melt, aes(x=Label, y=Value, fill=Cluster))
  if (!is.null(populations_df)){
    if (facet){
      gg <- gg + facet_grid( Cluster ~ Population, scales = "free_x", space = "free_x")
    } else{
      gg <- gg + facet_grid( . ~ Population, scales = "free_x", space = "free_x")
    }
  } else{
    if (facet){
      gg <- gg + facet_grid( Cluster ~ ., scales = "free_x", space = "free_x")
    }
  }
  gg <- gg + geom_bar(stat = "identity", width=1)
  gg <- gg + scale_y_continuous(expand=c(0,0), breaks=c(0.25,0.75))
  gg <- gg + coord_cartesian(ylim=c(0,1))
  gg <- gg + xlab("Sample ID") + ylab("Proportion of cluster") 
  gg <- gg + theme_bw()
  gg <- gg + guides(fill=guide_legend(title="Cluster"))
  gg <- gg + theme(axis.text.x = element_text(angle = 90))
  
  ggplot2::ggsave(outfile, plot = gg, width = 12, height = 10, dpi = 600)
}
