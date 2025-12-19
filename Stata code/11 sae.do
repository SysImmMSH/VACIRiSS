* 11 sae.do
* Author: David Harrison
* Date: 02/03/2023
* 
version 18
set more off
clear
cap log close
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"
use working, clear
count if alloc==0
local N0=r(N)
count if alloc==1
local N1=r(N)
use "1 Raw data\forms\AE", clear
drop if AE_01==0
drop if AE_01>=777
keep if inlist(AE_10,1,2,3,4)
merge m:1 Label using working, keep(3) keepusing(alloc RAN_02)
replace AE_03=8 if AE_01==6 & mi(AE_03)
replace AE_03=4 if AE_01==11 & mi(AE_03)
replace AE_03=16 if AE_03==999|AE_03==777
gen byte events=1
preserve
egen patients=tag(Label)
collapse (sum) events patients, by(alloc)
local n0=patients[1]
local n1=patients[2]
gen AE_03=99
reshape wide events patients, i(AE_03) j(alloc)
tempfile sae
save `sae'
restore
egen patients=tag(Label AE_03)
collapse (sum) events patients, by(alloc AE_03)
reshape wide events patients, i(AE_03) j(alloc)
mvencode events0 patients0 events1 patients1, mv(0)
append using `sae'
sort AE_03
label define AE_03_ 99 "Any", add
forvalues i=0/1 {
	gen Events`i'=string(events`i')
	gen Patients`i'=string(patients`i')+" ("+trim(string(100*patients`i'/`N`i'',"%4.1f"))+")"
}
label var AE_03 "Body system"
label var Events0 "Placebo events"
label var Events1 "PCV13 events"
label var Patients0 "Placebo patients (%)"
label var Patients1 "PCV13 patients (%)"
export excel AE_03 Events0 Patients0 Events1 Patients1 using VACIRiSS_Final_`today'.xlsx, first(varl) sheet(SAE, replace)
log using sae_`today'.log, replace
tabi `n1' `=`N1'-`n1'' \\  `n0' `=`N0'-`n0'', row exact
log close
