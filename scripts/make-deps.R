#scans your scripts/ and utils/ folders
#finds all used packages
#resolves their versions
#prints a nice Markdown-formatted dependency list to the console (or saves it to a file if you want)

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# 1. Detect packages in your project
deps <- renv::dependencies(paths = c("scripts", "utils"), progress = FALSE)
direct_pkgs <- sort(unique(deps$Package))

# 2. Remove base/recommended packages
base_pkgs <- rownames(installed.packages(priority = "base"))
direct_pkgs <- setdiff(direct_pkgs, base_pkgs)

# 3. Get installed versions
ip <- as.data.frame(installed.packages()[, c("Package", "Version")])
pkg_info <- ip[ip$Package %in% direct_pkgs, , drop = FALSE]

# 4. Sort alphabetically
pkg_info <- pkg_info[order(pkg_info$Package), ]

# 5. Print Markdown list
cat("## Dependencies\n\n")
cat("Detected by scanning scripts and utils:\n\n")
for (i in seq_len(nrow(pkg_info))) {
  cat(sprintf("- `%s` (>= %s)\n", pkg_info$Package[i], pkg_info$Version[i]))
}

# 6. Optional: also write to a file (uncomment next line)
# write.table(pkg_info, "dependencies.txt", row.names = FALSE, quote = FALSE, sep = "\t")
