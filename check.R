# automatic check for R package
library(devtools)
# R CMD check
check()
# check package on windows
check_win_devel()
check_win_release()
check_win_oldrelease()
# check package on mac
check_mac_release() 
# check package on Rhub
check_rhub()

# submit to CRAN
submit_cran()
