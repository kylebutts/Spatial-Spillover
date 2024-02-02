# REPRO_ORDER <- list()
# repro_init <- function() { 
#   REPRO_ORDER <<- list()
# }
# repro_section <- function(str) {
#   if (!exists(str, where = REPRO_ORDER)) {
#     REPRO_ORDER[str] = c()
#   }
# }

#' Code layout:
#' 1. Functions `preprocess_R/py/jl_file` that will do small changes before
#' rendering. Mainly converting # Header ---- to #' ## Header.
#' This will save the file as a `index.R/py/jl`
#'
#' 2. Function that takes `index.R/py/jl`, runs `quarto render` to produce
#' `index.html.md`, using `keep-md` true. This is the intermediate step
#' before producing a html file that has run the code, but not converted.
#' This will be moved to the logbook as `index.md`. Additionally, a
#' `readme.md` will be made from rendering to `gfm`.
#' This function is also in charge of cleanup.
#

#' Take file and create a default title
#' @param file A file name
#' @examples 
#'  file_to_title("render_file.R")
#' @return String of title
file_to_title <- function(file) {
  default_title <- fs::path_ext_remove(fs::path_file(file))
  default_title <- gsub("_", " ", default_title)
  default_title <- gsub("-", " ", default_title)
  # title case
  default_title <- gsub(
    "(^|\\s)([a-z])",
    "\\1\\U\\2",
    default_title,
    perl = TRUE
  )
  return(default_title)
}
preprocess_R_file <- function(file, index_file_name) {
  txt <- xfun::read_utf8(here::here(file))

  # Add default yaml frontmatter, if needed
  # Look for `#' ---` at the start of the file
  # matches https://github.com/quarto-dev/quarto-cli/blob/770c47fbfeefe404fa43799116922bdb7c47a3c1/src/execute/rmd.ts#L395
  has_yaml <- grepl("^\\s*#'\\s*---", txt[1])
  if (!has_yaml) {
    txt <- c(
      "#' ---",
      paste0("#' title: '", file_to_title(file), "'"),
      "#' ---",
      txt
    )
  }

  # Convert code headers
  # `# Header ----` -> `#' ## Header`
  regex_header <- "^(#{1,}) (.*?)-{4,}"
  txt <- gsub(
    regex_header,
    "#' \\1# \\2",
    txt
  )

  # Save as index file
  xfun::write_utf8(txt, index_file_name)
  return(invisible(NULL))
}
preprocess_jl_file <- function(file, index_file_name) {
  txt <- xfun::read_utf8(here::here(file))

  # Add default yaml frontmatter, if needed
  # Look for `# %% [markdown]` at the start of the file
  # Then either `# ---` or `"""\n---`
  starts_with_md <- grepl("^#\\s*%%\\s*\\[markdown\\]", txt[1])
  has_yaml <- grepl("^#\\s*---", txt[2]) |
    grepl("^\\s*\"\"\"", txt[2]) & grepl("^\\s*---", txt[3])

  if (!(starts_with_md & has_yaml)) {
    txt <- c(
      "# %% [markdown]",
      "# ---",
      paste0("# title: '", file_to_title(file), "'"),
      "# ---",
      txt
    )
  }

  # Save as index file
  xfun::write_utf8(txt, index_file_name)
  return(invisible(NULL))
}
preprocess_py_file <- function(file, index_file_name) {
  txt <- xfun::read_utf8(here::here(file))

  # Add default yaml frontmatter, if needed
  # Look for `# %% [markdown]` at the start of the file
  # Then either `# ---` or `"""\n---`
  starts_with_md <- grepl("^#\\s*%%\\s*\\[markdown\\]", txt[1])
  has_yaml <- grepl("^#\\s*---", txt[2]) |
    grepl("^\\s*\"\"\"", txt[2]) & grepl("^\\s*---", txt[3])

  if (!(starts_with_md & has_yaml)) {
    txt <- c(
      "# %% [markdown]",
      "# ---",
      paste0("# title: '", file_to_title(file), "'"),
      "# ---",
      txt
    )
  }

  # Save as index file
  xfun::write_utf8(txt, index_file_name)
  return(invisible(NULL))
}
preprocess_qmd_file <- function(file, index_file_name) {
  txt <- xfun::read_utf8(here::here(file))

  # Save as index file
  xfun::write_utf8(txt, index_file_name)
  return(invisible(NULL))
}

delete_if_exists <- function(file) {
  if (fs::file_exists(file)) {
    fs::file_delete(file)
  }
}

#' Create a log of script files (in `.R`/`.py`/`.jl`)
#'
#' @description
#' This takes a plain `.R`/`.py`/`.jl` script and records all output to a log. 
#' The function uses `quarto render` under the hood to create markdown logs. 
#' 
#' @details
#' This function produces two markdown files. First is `readme.md` which will
#' display nicely in github and `index.md` which will be rendered to html. 
#' This function renders the file only one time and points to the same files 
#' (e.g. any plots created), so the memory and compute costs of two versions
#' is minimal. 
#' 
#' The function takes a path to the `.R` file and creates an entry to `out_dir` using the same folder strucutre as the `.R file`. 
#' For example, if you have `file = "code/analysis/estimate_regressions.R"`,
#' `out_base_dir = "logbook"` and `code_base_dir = "code"` it will create a logbook entry at `logbook/analysis/create_land_price_index.Rmd`.
#'
#' File can specify their own yaml frontmatter following 
#' <https://quarto.org/docs/prerelease/1.4/script.html>, 
#' or a default title block can be used if none are detected. 
#' 
#' In R, you can specify frontmatter as
#' ```
#' #' ---
#' #' title: "Good Title"
#' #' ---
#' ```
#' In python/julia, you can specify frontmatter as
#' ```
#' # %% [markdown]
#' # ---
#' # title: "Good Title"
#' # ---
#' ```
#'
#' @param file A `.R`/`.py`/`.jl` script that you want to log results from.
#'   Can either be relative to `here::here()` or an absolute path.
#' @param out_base_dir The base directory for the logbook. 
#'   Default is `logbook`. 
#' @param code_base_dir The base directory for the code files. 
#'   Paths in `out_base_dir` are based on the location of the script. 
#'   This lets you specify the base directory for the scripts.
#'   Default is `code`.
#'   E.g. `{code_base_dir}/analysis/estimate_regressions.R` would become 
#'   `{out_base_dir}/analysis/estimate_regressions/` in the logbook.
#' @param force_rerun The default behavior is to rerender each script.
#'   Set this option to `FALSE` to only rerun the file if the script is 
#'   modified more recently than the logbook entry.
#'
#' @return NULL, but creates a log of the file in the logbook.
#' 
render_file <- function(
  file, out_base_dir = "logbook", code_base_dir = "code", 
  force_rerun = TRUE
) {
  ## Process paths
  # In case absolute path is given, take relative to `here::here()`
  file <- fs::path_rel(file, here::here())
  file_name <- fs::path_file(file)
  file_dir <- dirname(file)
  file_ext <- fs::path_ext(file)

  # Check if `index.R/py/jl` exists to avoid overwriting a file
  index_file_name <- here::here(
    file_dir, paste0("index.", file_ext)
  )
  if (fs::file_exists(index_file_name)) {
    stop(paste0(
      "Stopping to avoid overwriting file.\n",
      "File already exists: ", index_file_name
    ))
  }

  # Rendering files
  knit_file <- here::here(file_dir, "index.html.md")
  html_file <- here::here(file_dir, "index.html")
  readme_file <- here::here(file_dir, "index.html-gfm.md")
  md_file_folder <- here::here(file_dir, "index_files/figure-html")
  index_files_dir <- here::here(file_dir, "index_files")
  
  # This entry's logbook folder
  out_dir <- here::here(
    out_base_dir,
    fs::path_ext_remove(fs::path_rel(file, code_base_dir))
  )

  ## Check if needed to rerun
  if (fs::file_exists(here::here(out_dir, "index.md")) & force_rerun == FALSE) {
    file_info <- fs::file_info(here::here(file))
    out_info <- fs::file_info(here::here(out_dir, "index.md"))

    if (file_info$modification_time < out_info$modification_time) {
      message(paste0(
        "The logbook entry is more recent than ", file, ".\n",
        "Skipping for efficiency reasons. To override and rerender, set `force_rerun = TRUE`."
      ))
      return(invisible(NULL))
    }
  }

  ## Clean up 
  ## Doing it this way in case quarto fails
  on.exit({
    if (fs::dir_exists(index_files_dir)) {
      fs::dir_delete(index_files_dir)
    }
    delete_if_exists(knit_file)
    delete_if_exists(readme_file)
    delete_if_exists(html_file)
    delete_if_exists(index_file_name)
    # For julia/python
    delete_if_exists(here::here(file_dir, "index.ipynb"))
  })

  ## Preprocess file, creating `index_file_name`
  switch(file_ext,
    "R" = preprocess_R_file(file, index_file_name),
    "py" = preprocess_py_file(file, index_file_name),
    "jl" = preprocess_jl_file(file, index_file_name),
    "qmd" = preprocess_qmd_file(file, index_file_name),
    stop(paste0(
      "File extension not supported: ", file_ext
    ))
  )

  ## Render
  # Render to html, keeping .html.md
  quarto_path <- quarto::quarto_path()
  system(paste0(
    quarto_path,
    " render '", index_file_name, "'",
    " --to html",
    " -M keep-md:true",
    " --execute-dir '", here::here(), "'",
    collapse = ""
  ))
  
  # Render to gfm
  quarto::quarto_render(
    knit_file,
    output_format = "gfm",
    execute_dir = here::here()
  )

  ## Copy file over to logbook
  fs::dir_create(out_dir)

  fs::file_copy(
    knit_file,
    here::here(out_dir, "index.md"),
    overwrite = TRUE
  )
  fs::file_copy(
    readme_file,
    here::here(out_dir, "readme.md"),
    overwrite = TRUE
  )

  if (fs::dir_exists(md_file_folder)) {
    out_dir_files <- here::here(out_dir, "index_files/figure-html")
    if (fs::dir_exists(out_dir_files)) fs::dir_delete(out_dir_files)
    fs::dir_create(out_dir_files)

    fs::file_move(
      md_file_folder,
      here::here(out_dir, "index_files")
    )
  }

  return(invisible(TRUE))
}
