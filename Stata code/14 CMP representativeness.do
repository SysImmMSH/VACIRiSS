* 14 CMP representativeness.do
* Author: David Harrison
* Date: 22/07/2024
* 
version 18
set more off
clear
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"

use CMP_sepsis_survivors.dta, clear
tempname pf
postfile `pf' row str200 desc str30(allunits trialunits) using representativeness, replace
local row 1
count
local N1=r(N)
count if vaciriss_site==1
local N2=r(N)
gen byte grp1=1
gen byte grp2=vaciriss_site
post `pf' (`row++') ("Characteristic") ("All CMP (N = `N1')") ("VACIRiSS sites (N = `N2')")
* Age
forvalues i=1/2 {
	su calage if grp`i'
	local val`i'=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+")"
}
post `pf' (`row++') ("Age (years), mean (SD)") ("`val1'") ("`val2'")
* Sex
forvalues i=1/2 {
	count if !mi(sex) & grp`i'
	local N`i'=r(N)
}
post `pf' (`row++') ("Sex, n (%)") ("(N = `N1')") ("(N = `N2')")
forvalues i=1/2 {
	count if !mi(sex) & grp`i'
	local N=r(N)
	count if sex=="F" & grp`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
}
post `pf' (`row++') ("  Female") ("`val1'") ("`val2'")
forvalues i=1/2 {
	count if !mi(sex) & grp`i'
	local N=r(N)
	count if sex=="M" & grp`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
}
post `pf' (`row++') ("  Male") ("`val1'") ("`val2'")
* Ethnicity
recode ethng (1 = 1) (2 = 4) (3 = 2) (4 = 3) (5 = 5) (9 = .), gen(eth)
local eth1 "White"
local eth2 "Asian"
local eth3 "Black"
local eth4 "Mixed"
local eth5 "Other"
forvalues i=1/2 {
	count if !mi(eth) & grp`i'
	local N`i'=r(N)
}
post `pf' (`row++') ("Ethnicity, n (%)") ("(N = `N1')") ("(N = `N2')")
forvalues e=1/5 {
	forvalue i=1/2 {
		count if !mi(eth) & grp`i'
		local N=r(N)
		count if eth==`e' & grp`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
	}
	post `pf' (`row++') ("  `eth`e''") ("`val1'") ("`val2'")
}
* Surgical status
recode admtype (1 = 3) (2 = 1) (3 = 2), gen(ss)
local surg1 "Elective surgery"
local surg2 "Emergency surgery"
local surg3 "Medical"
forvalues i=1/2 {
	count if !mi(ss) & grp`i'
	local N`i'=r(N)
}
post `pf' (`row++') ("Surgical status, n (%)") ("(N = `N1')") ("(N = `N2')")
forvalues s=1/3 {
	forvalue i=1/2 {
		count if ss==`s' & grp`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' (`row++') ("  `surg`s''") ("`val1'") ("`val2'")
}
* Dependency
forvalues i=1/2 {
	count if !mi(dep_cat) & grp`i'
	local N`i'=r(N)
}
post `pf' (`row++') ("Pre-admission dependence, n (%)") ("(N = `N1')") ("(N = `N2')")
forvalues i=1/2 {
	count if dep_cat==0 & grp`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' (`row++') ("  None") ("`val1'") ("`val2'")
forvalues i=1/2 {
	count if dep_cat==1 & grp`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' (`row++') ("  Moderate (some assistance with ADLs)") ("`val1'") ("`val2'")
forvalues i=1/2 {
	count if dep_cat==2 & grp`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' (`row++') ("  Severe (total assistance with ADLs)") ("`val1'") ("`val2'")
* Residence
gen byte res=1 if inlist(resa,"M","R")|inlist(resa_v4,"M","A")
replace res=3 if inlist("H",resa,resa_v4)
replace res=4 if resa=="O"
replace res=6 if inlist("N",resa,resa_v4)
forvalues i=1/2 {
	count if !mi(res) & grp`i'
	local N`i'=r(N)
}
post `pf' (`row++') ("Pre-admission residence, n (%)") ("(N = `N1')") ("(N = `N2')")
local res1 "Home"
local res2 "Nursing home"
local res3 "Health-related institution"
local res4 "Non-health-related institution"
local res5 "Hospice or equivalent"
local res6 "No fixed address/abode or temporary abode"
foreach r in 1 3 4 6 {
	forvalues i=1/2 {
		count if res==`r' & grp`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' (`row++') ("  `res`r''") ("`val1'") ("`val2'")
}
* APACHE II score
forvalues i=1/2 {
	su AP2score if grp`i'
	local val`i'=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+") [N = "+string(r(N))+"]"
}
post `pf' (`row++') ("APACHE II score, mean (SD)") ("`val1'") ("`val2'")
postclose `pf'

*** Column for VACIRiSS participants
tempfile tf
postfile `pf' row str30 trialpatients using `tf', replace
local row 1
use working, clear
count
local N=r(N)
post `pf' (`row++') ("VACIRiSS participants (N = `N')")
* Age
su age
local val=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+")"
post `pf' (`row++') ("`val'")
* Sex
count if !mi(REG_05)
local N=r(N)
post `pf' (`row++') ("(N = `N')")
count if REG_05==2
local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
post `pf' (`row++') ("`val'")
count if REG_05==1
local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
post `pf' (`row++') ("`val'")
* Ethnicity
count if !mi(DEM_01)
local N=r(N)
post `pf' (`row++') ("(N = `N')")
forvalues e=1/5 {
	count if DEM_01==`e'
	local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
	post `pf' (`row++') ("`val'")
}
* Surgical status
count if inlist(ISC_04,1,2,3)
local N=r(N)
post `pf' (`row++') ("(N = `N')")
forvalues s=1/3 {
	count if ISC_04==`s'
	local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
	post `pf' (`row++') ("`val'")
}
* Dependency
count if inlist(ISC_05,1,2,3,4)
local N=r(N)
post `pf' (`row++') ("(N = `N')")
count if ISC_05==1
local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
post `pf' (`row++') ("`val'")
count if inlist(ISC_05,2,3)
local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
post `pf' (`row++') ("`val'")
count if ISC_05==4
local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
post `pf' (`row++') ("`val'")
* Residence
count if inlist(ISC_06,1,2,3,4,5,6)
local N=r(N)
post `pf' (`row++') ("(N = `N')")
foreach r in 1 3 4 6 {
	count if ISC_06==`r'
	local val=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
	post `pf' (`row++') ("`val'")
}
* APACHE II score
su ISC_11 if ISC_11<72
local val=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+") [N = "+string(r(N))+"]"
post `pf' (`row++') ("`val'")
postclose `pf'

use representativeness, clear
merge 1:1 row using `tf'
drop row _merge
compress
save, replace
export excel using VACIRiSS_Final_`today'.xlsx, sheet(Representativeness, replace)

