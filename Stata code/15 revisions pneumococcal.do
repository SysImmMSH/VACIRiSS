* 15 revisions pneumococcal.do
* Author: David Harrison
* Date: 12/12/2025
* Response to reviewers - Check for pneumococcal infections in HES
version 17
set more off
clear
* base working directory stored in global
cd "$base"

log using "revisions pneumococcal.log", replace
use "1 Raw data\hes_stacked", clear

preserve

* list of all 20 diag suffixes
local diag = "01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20"

* create infection variables for each episode
gen byte pneumococcal=0
foreach i of local diag {
	replace diag_`i'=substr(diag_`i',1,4)
	replace diag_`i'=subinstr(diag_`i', "-","",.) 
	replace diag_`i'=subinstr(diag_`i', "X","",.) 
	replace pneumococcal=1 if inlist(diag_`i',"A403","B953")
}
tab pneumococcal
log close
