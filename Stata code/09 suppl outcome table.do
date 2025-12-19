* 09 suppl outcome table.do
* Author: David Harrison
* Date: 22/07/2024
* 
version 18
set more off
clear
cap log close
set scheme s2color
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"
use primary, clear
append using subgroups
* row order
gen byte y=.
replace y=1 if outcome=="covid"
replace y=3 if outcome=="SSIP0"
replace y=4 if outcome=="SSIP5"
replace y=5 if outcome=="SSIP7"
replace y=6 if outcome=="SSIP11"

* column 1 - outcome names
gen col1=""
replace col1="Sensitivity analysis - pre-Covid" if outcome=="covid"
replace col1="  Risk stratum 1 (low, 0-4 points)" if outcome=="SSIP0"
replace col1="  Risk stratum 2 (5-6 points)" if outcome=="SSIP5"
replace col1="  Risk stratum 3 (7-10 points)" if outcome=="SSIP7"
replace col1="  Risk stratum 4 (highest, ≥11 points)" if outcome=="SSIP11"
set obs `=_N+1'
replace y=2 in l
replace col1="Subgroup analysis by SSiP category" in l
sort y
* column 2 - percent or rate - PCV-13
gen col2=""
replace col2=string(events1)+"/"+string(denom1/365.25,"%4.1f")+" ("+string(365.25*events1/denom1,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col2=string(events1)+"/"+string(round(denom1))+" ("+string(100*events1/denom1,"%4.1f")+")" if effect=="RR"
* column 3 - percent or rate - placebo
gen col3=""
replace col3=string(events0)+"/"+string(denom0/365.25,"%4.1f")+" ("+string(365.25*events0/denom0,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col3=string(events0)+"/"+string(round(denom0))+" ("+string(100*events0/denom0,"%4.1f")+")" if effect=="RR"
* column 5 - relative effect
gen col5=string(estimate,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if !mi(effect)
* column 6 - p value
gen col6="<0.001" if pval<.001
replace col6=string(pval,"%5.3f") if inrange(pval,.001,.01)
replace col6=string(pval,"%4.2f") if pval>=.01 & !mi(pval)
replace col6=">0.99" if pval>.99 & !mi(pval)
replace col6=col6+"a" if inlist(outcome,"primary - 65 and under","SSIP0")
drop if mi(y)
tostring N1 N0, replace
replace N1="" if N1=="."
replace N0="" if N0=="."

label var col1 "Analysis"
label var N1 "N PCV13"
label var col2 "PCV13"
label var N0 "N Placebo"
label var col3 "Placebo"
label var col5 "Hazard ratio (95% CI)"
label var col6 "P value"
export excel col1 N1 col2 N0 col3 col5 col6 using VACIRiSS_Final_`today'.xlsx if !mi(y), first(varl) sheet(Suppl_outcome, replace)
