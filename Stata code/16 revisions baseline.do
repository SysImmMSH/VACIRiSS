* 16 revisions baseline.do
* Author: David Harrison
* Date: 12/12/2025
* Response to reviewers - Check for pneumococcal infections in HES
version 17
set more off
clear
* base working directory stored in global
cd "$base"

log using "revisions baseline.log", replace
use working, clear
bys alloc: su age, d
gen yulosd=yulos/24
label var yulosd "ICU length of stay in days"
bys alloc: su yulosd, d
log close
