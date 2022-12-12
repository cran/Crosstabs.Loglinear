
.onAttach<- function(libname, pkgname){
  packageStartupMessage("**************************************************************************************************\n",
                        pkgname," ",packageDescription("Crosstabs.Loglinear")$Version,
                        "\n\nPlease contact Brian O'Connor at brian.oconnor@ubc.ca if you have questions or suggestions.\n",
                        "**************************************************************************************************", 
                        appendLF = TRUE)
}

NULL


# cat("Package: Crosstabs.Loglinear    Version: ", packageDescription("Crosstabs.Loglinear")$Version,"\n")


# # .onAttach<- function(libname, pkgname){
  # packageStartupMessage("**************************************************************************************************\nWelcome to ",
                        # pkgname," ",packageDescription("Crosstabs.Loglinear")$Version,
                        # "\n\nPlease contact Brian O'Connor at brian.oconnor@ubc.ca if you have questions or suggestions.\n",
                        # "**************************************************************************************************", 
                        # appendLF = TRUE)
# }

