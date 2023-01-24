#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(knitr, warn.conflicts = F, quietly = T)))
suppressWarnings(suppressMessages(library(argparse, warn.conflicts = F, quietly = T)))
suppressWarnings(suppressMessages(library(pandoc, warn.conflicts = F, quietly = T)))
suppressWarnings(suppressMessages(library(R.utils, warn.conflicts = F, quietly = T)))


#' current script file (in full path)
#' @description current script file (in full path)
#' @examples
#' works with Rscript, source() or in RStudio Run selection, RStudio Console
#' @export
ez.csf <- function() {
  # http://stackoverflow.com/a/32016824/2292993
  cmdArgs = commandArgs(trailingOnly = FALSE)
  needle = "--file="
  match = grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript via command line
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    ls_vars = ls(sys.frames()[[1]])
    if ("fileName" %in% ls_vars) {
      # Source'd via RStudio
      return(normalizePath(sys.frames()[[1]]$fileName))
    } else {
      if (!is.null(sys.frames()[[1]]$ofile)) {
        # Source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
      } else {
        # RStudio Run Selection
        # http://stackoverflow.com/a/35842176/2292993
        pth = rstudioapi::getActiveDocumentContext()$path
        if (pth!='') {
          return(normalizePath(pth))
        } else {
          # RStudio Console
          tryCatch({
            pth = rstudioapi::getSourceEditorContext()$path
            pth = normalizePath(pth)
          }, error = function(e) {
            # normalizePath('') issues warning/error
            pth = ''
          }
          )
          return(pth)
        }
      }
    }
  }
}


cargs <- commandArgs(TRUE)


if (!argparse:::detects_python()) {
  invokeRestart("abort")
}

parser <- ArgumentParser(description='benchplot')
parser$add_argument('--input', dest='input', type = "character",
                    help='Input csv file')
parser$add_argument('--output', dest='output', type = "character",
                    help='Output html file')
parser$add_argument('--tools', dest='tools', type ="character", default='',
                    help='Only for listed tools (comma separated)')
parser$print_help()
  
args <- parser$parse_args(cargs)

print(paste("output: ", args$output))
pandoc_version()
getwd()

input <- args$input
output <- args$output
print(paste("output: ", output))
print(paste("input:  ", input))
if (!isAbsolutePath(input)) {
  input <- file.path(getwd(), sub('./', '', args$output))
}
if (!isAbsolutePath(output)) {
  output <- file.path(getwd(), sub('./', '', output))
}

print(paste("output: ", output))
print(paste("input:  ", input))

args_list = list(
  input = args$input,
  tools = args$tools
)
rmd_path <- "/usr/users/QIB_fr017/fritsche/Projects/benchpro/benchpro.Rmd"
out_path <- "/usr/users/QIB_fr017/fritsche/Projects/benchpro/benchplot/output.html"


sourceDir <- dirname(ez.csf())
print(sourceDir)
print(paste("Sourcedir: " , sourceDir))
path_to_rmd <- file.path(sourceDir, "benchpro.Rmd")

print(paste("args_list: ", args_list))
print(paste("Rmd: ", path_to_rmd))

rmarkdown::render(path_to_rmd, output_file=output, params=args_list)
