
# Allow only letters, numbers, underscore, period for IDs
is_valid_identifier <- function(x) {
  grepl("^[A-Za-z0-9_.]+$", x)
}
