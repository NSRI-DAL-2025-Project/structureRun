#' @title Runs a faststructure analysis using a genlight object
#'
#' @description
#' This function takes a genlight object and runs a faststructure analysis.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param k.range Range of the number of populations [required].
#' @param num.k.rep Number of replicates [default 1].
#' @param exec Full path and name+extension where the fastStructure executable
#' is located [default working directory "./fastStructure"].
#' @param exec.plink path to plink executable [default working directory].
#' @param output Path to output file [default getwd()].
#' @param tol Convergence criterion [default 10e-6].
#' @param prior Choice of prior: simple or logistic [default "simple"].
#' @param cv Number of test sets for cross-validation, 0 implies no CV step
#'  [default 0].
#' @param seed Seed for random number generator [default NULL].
#' @details
#' Download faststructure binary for your system from here (only runs on Mac or
#' Linux):

#' @author Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
#'
#' @examples
#' \dontrun{
#' # Please note: faststructure needs to be installed
#' # Please note: faststructure is not available for windows
#' t1 <- gl.filter.callrate(platypus.gl, threshold = 1)
#' res <- gl.run.faststructure(t1,
#'   exec = "./fastStructure", k.range = 2:3,
#'   num.k.rep = 2, output = paste0(getwd(), "/res_str")
#' )
#' qmat <- gl.plot.faststructure(res, k.range = 2:3)
#' gl.map.structure(qmat, K = 2, t1, scalex = 1, scaley = 0.5)
#' }
#' @export
gl.run.faststructure <- function(x,
                                 k.range,
                                 num.k.rep = 1,
                               #  exec = "./fastStructure",
                                # exec.plink = getwd(),
                                 output = getwd(),
                                 tol = 10e-6,
                                 prior = "simple",
                                 cv = 0,
                                 seed = NULL) {
  pkg <- "gsubfn"
  if (!(requireNamespace(pkg, quietly = TRUE))) {
    cat(error(
      "Package",
      pkg,
      " needed for this function to work. Please install it.\n"
    ))
    return(-1)
  }

  exec <- Sys.which("fastStructure-master")
  exec.plink <- Sys.which("plink")
  
  ############################
  
  if (tools::file_ext(x) == "csv") {
    file <- readr::read_csv(x)
  } else if (tools::file_ext(x) == "xlsx") {
    file <- readxl::read_excel(x)
  } else {
    stop("Input file should be in csv or xlsx format.")
  }
  
  
  file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
  file <- as.data.frame(file)
  
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
  
  x <- dartR::gi2gl(fsnps_gen, verbose = 0)

  ###########################
  
  
  dartR.base::gl2plink(x,
    bed.files = TRUE,
    outpath = output,
    verbose = 0,
    plink.bin.path = exec.plink
  )

  for (k_n in k.range) {
    
    for (rep_n in 1:num.k.rep) {
      print(paste("Running K =",k_n,"Replicate =",rep_n))
      if (is.null(seed)) {
        system(
          paste0(
            exec,
            " -K ",
            k_n,
            " --input=",
            output,
            "/gl_plink",
            " --output=",
            output,
            "/genotypes_output",
            "_",
            rep_n,
            " --tol=",
            tol,
            " --prior=",
            prior,
            " --cv=",
            cv
          )
        )
      } else {
        system(
          paste0(
            exec,
            " -K ",
            k_n,
            " --input=",
            output,
            "/gl_plink",
            " --output=",
            output,
            "/genotypes_output",
            "_",
            rep_n,
            " --tol=",
            tol,
            " --prior=",
            prior,
            " --cv=",
            cv,
            " --seed=",
            seed
          )
        )
      }
    }
  }

  files_structure <- list.files(
    path = output,
    pattern = "^genotypes_output.+log"
  )
  files_structure_2 <- paste0(output, "/", files_structure)
  n_first_line <- 23
  df_likelihood <-
    as.data.frame(matrix(nrow = length(k.range), ncol = num.k.rep))
  for (i in 1:length(files_structure_2)) {
    file_name <- files_structure[i]
    file_name <-
      strsplit(sub("(^[^_]+_[^_]+)_(.*)$", "\\2", file_name), " ")
    file_name <- strsplit(file_name[[1]], "\\.")
    k_replicate <- as.numeric(file_name[[1]][1])
    k_run <- as.numeric(file_name[[1]][2])
    likelihood <- as.character(unname(unlist(
      gsubfn::read.pattern(
        file = files_structure_2[i],
        pattern = "^Marginal Likelihood = .*"
      )
    )))
    likelihood_2 <-
      as.numeric(substr(likelihood, n_first_line, nchar(likelihood)))
    df_likelihood[k_run, k_replicate] <- likelihood_2
  }
  df_likelihood <-
    as.data.frame(df_likelihood[stats::complete.cases(df_likelihood), ])
  df_likelihood_res <- rowMeans(df_likelihood)

  p3 <- ggplot() +
    geom_line(aes(x = k.range, y = df_likelihood_res), size = 1) +
    geom_point(aes(x = k.range, y = df_likelihood_res),
      size = 2,
      color = "blue"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "bottom") +
    xlab("K") +
    ylab("Marginal Likelihood") +
    scale_x_continuous(breaks = round(seq(1, 10, by = 1), 1))

  print(p3)

  files_structure <- list.files(
    path = output,
    pattern = "^genotypes_output.+log"
  )
  files_structure_2 <- paste0(output, "/", files_structure)

  files_q <- list.files(path = output, pattern = "*meanQ")
  files_q_2 <- paste0(output, "/", files_q)
  q_list <-
    rep(list(as.list(rep(NA, num.k.rep))), length(k.range) + 1)
  for (i in 1:length(files_q)) {
    file_name <- files_q[i]
    file_name <-
      strsplit(sub("(^[^_]+_[^_]+)_(.*)$", "\\2", file_name), " ")
    file_name <- strsplit(file_name[[1]], "\\.")
    k_replicate <- as.numeric(file_name[[1]][1])
    k_run <- as.numeric(file_name[[1]][2])
    q_df <- read.table(files_q_2[i])
    q_df <- cbind(id = x$ind.names, orig.pop = pop(x), q_df)
    q_list[[k_run]][[k_replicate]] <- q_df
  }

  names(q_list) <- 1:(length(k.range) + 1)

  q_list <- lapply(q_list, function(y) {
    names(y) <- 1:num.k.rep
    return(y)
  })

  q_list <- q_list[-(1)]

  return(list(q_list=q_list,plot=p3))
}
