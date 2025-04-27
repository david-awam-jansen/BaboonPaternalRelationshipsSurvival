# build_readme.R

# Read in README.md (up to an insertion point)
readme <- readLines("../README.md")

# Read Rpackages.md
rpackages <- readLines("Rpackages.md")

# Find a placeholder in README.md like <!-- RPACKAGES_START -->
insert_point <- grep("<!-- RPACKAGES_START -->", readme)

# Build new README
new_readme <- c(
    readme[1:(insert_point)],
    "",
    rpackages,
    "",
    readme[(insert_point + 1):length(readme)]
)

# Write back to README.md
writeLines(new_readme, "../README.md")
