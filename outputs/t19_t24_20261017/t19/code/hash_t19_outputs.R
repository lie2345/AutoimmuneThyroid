#!/usr/bin/env Rscript
# Record and immediately verify SHA256 for every deliverable except this manifest.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[1], mustWork = TRUE)
stopifnot(endsWith(out, "/outputs/t19_t24_20261017/t19"))
if (!requireNamespace("digest", quietly = TRUE)) stop("Installed digest package is required")
manifest <- file.path(out, "validation/output_sha256.tsv")
files <- list.files(out, full.names = TRUE, recursive = TRUE)
files <- files[files != manifest]
sha <- function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE)
hashes <- vapply(files, sha, "")
records <- data.frame(path = substring(files, nchar(out) + 2L), bytes = file.info(files)$size,
                      sha256 = unname(hashes))
write.table(records, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(identical(records$sha256, unname(vapply(file.path(out, records$path), sha, ""))))
cat("Output SHA256 roundtrip: ", nrow(records), "/", nrow(records), " PASS\n", sep = "")
