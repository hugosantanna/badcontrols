# =============================================================================
# Title: Run revdeplite
# Description: Run lightweight reverse dependency checks for badcontrols.
# Author: Brant Callaway
# Last update: 2026-11-24
# Date created: 2026-07-27
# =============================================================================

# No reverse dependencies on CRAN or GitHub as of 2026-11-24; add
# github_deps = c("owner/repo", ...) here if downstream packages appear.
revdeplite::revdeplite(
  check_dir = ".revdeplite"
)
