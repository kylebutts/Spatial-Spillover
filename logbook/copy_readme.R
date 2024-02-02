readme_file <- here::here("README.md")
out_file <- here::here("logbook", "index.md")

if (fs::file_exists(readme_file)) {
  fs::file_copy(readme_file, out_file, overwrite = TRUE)
} else {
  cat(
    "# TODO: Create a README.md file in the project directory",
    file = here::here("logbook", "index.md")
  )
}
