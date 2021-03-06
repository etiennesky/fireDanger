#' @importFrom utils packageDescription

.onAttach <- function(...) {
    pkgname <- "fireDanger"
    ver <- packageDescription(pkgname)$Version
    builddate <- packageDescription(pkgname)$Date
    mess <- paste(pkgname, " version ", ver, " (", builddate,") is loaded", sep = "")
    packageStartupMessage(mess)
    url <- paste0("https://raw.githubusercontent.com/SantanderMetGroup/", pkgname, "/master/DESCRIPTION")
    b <- tryCatch(suppressWarnings(readLines(url)), error = function(er) {
        er <- NULL
        return(er)
    })
    if (!is.null(b)) {
        latest.ver <- package_version(gsub("Version: ", "", b[grep("^Version", b)]))
        if (ver < latest.ver) {
            ver.mess1 <- paste0("WARNING: Your current version of ", pkgname, " (v", ver, ") is not up-to-date")
            ver.mess <- paste0("Get the latest stable version (", latest.ver,
                               ") using <devtools::install_github('SantanderMetGroup/", pkgname, "')>")
            packageStartupMessage(ver.mess1)
            packageStartupMessage(ver.mess)
        }
    }
    packageStartupMessage("Type <vignette(\"Climate_Services_2017\", package = \"fireDanger\")> for package overview and examples\nPlease use citation(\"fireDanger\") to cite this package.")
}
# End

