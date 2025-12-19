* 03 table 1.do
* Author: David Harrison
* Date: 01/03/2023
* 
version 18
set more off
clear
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"
use working, clear
tempname pf
postfile `pf' str200 desc str30(intervention control) using table1, replace
count if alloc==0
local N0=r(N)
count if alloc==1
local N1=r(N)
post `pf' ("Characteristic") ("PCV13 (n = `N1')") ("Placebo (n = `N0')")
forvalues i=0/1 {
	su age if alloc==`i'
	local val`i'=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+")"
}
post `pf' ("Age (years), mean (SD)") ("`val1'") ("`val0'")
post `pf' ("Sex, No. (%)") ("") ("")
forvalues i=0/1 {
	count if !mi(REG_05) & alloc==`i'
	local N=r(N)
	count if REG_05==2 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
}
post `pf' ("  Female") ("`val1'") ("`val0'")
forvalues i=0/1 {
	count if !mi(REG_05) & alloc==`i'
	local N=r(N)
	count if REG_05==1 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
}
post `pf' ("  Male") ("`val1'") ("`val0'")
local eth1 "White"
local eth2 "Asian"
local eth3 "Black"
local eth4 "Mixed"
local eth5 "Other"
post `pf' ("Self-reported ethnicity, No. (%)") ("") ("")
forvalues e=1/5 {
	forvalue i=0/1 {
		count if !mi(DEM_01) & alloc==`i'
		local N=r(N)
		count if DEM_01==`e' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N',"%4.1f")+")"
	}
	post `pf' ("  `eth`e''") ("`val1'") ("`val0'")
}
* Comorbidity
forvalues i=0/1 {
	count if !mi(charlson) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("Charlson comorbidity index, No. (%)") ("(n = `N1')") ("(n = `N0')")
forvalues j=0/3 {
	forvalue i=0/1 {
		count if charlson==`j' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' ("  `j'") ("`val1'") ("`val0'")	
}
forvalue i=0/1 {
	count if charlson>=4 & !mi(charlson) & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  4+") ("`val1'") ("`val0'")
* Vaccination history
post `pf' ("Pneumococcal vaccination any time prior, No. (%)") ("") ("")
forvalue i=0/1 {
	count if VH_01==1 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Yes") ("`val1'") ("`val0'")
forvalue i=0/1 {
	count if VH_01==0 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  No") ("`val1'") ("`val0'")
forvalue i=0/1 {
	count if !inlist(VH_01,0,1) & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Unknown") ("`val1'") ("`val0'")
post `pf' ("Influenza vaccination last 12 months, No. (%)") ("") ("")
forvalue i=0/1 {
	count if VH_02==1 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Yes") ("`val1'") ("`val0'")
forvalue i=0/1 {
	count if VH_02==0 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  No") ("`val1'") ("`val0'")
forvalue i=0/1 {
	count if !inlist(VH_02,0,1) & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Unknown") ("`val1'") ("`val0'")

* Surgical status
forvalues i=0/1 {
	count if inlist(ISC_04,1,2,3) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("Surgical status, No. (%)") ("") ("")
local surg1 "Elective surgery"
local surg2 "Emergency surgery"
local surg3 "Medical"
forvalues s=1/3 {
	forvalue i=0/1 {
		count if ISC_04==`s' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' ("  `surg`s''") ("`val1'") ("`val0'")
}
* Dependency
forvalues i=0/1 {
	count if inlist(ISC_05,1,2,3,4) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("Pre-admission dependence, No. (%)") ("") ("")
forvalues i=0/1 {
	count if ISC_05==1 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  None") ("`val1'") ("`val0'")
forvalues i=0/1 {
	count if inlist(ISC_05,2,3) & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Moderate (some assistance with ADLs)") ("`val1'") ("`val0'")
forvalues i=0/1 {
	count if ISC_05==4 & alloc==`i'
	local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
}
post `pf' ("  Severe (total assistance with ADLs)") ("`val1'") ("`val0'")
* Residence
forvalues i=0/1 {
	count if inlist(ISC_06,1,2,3,4,5,6) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("Pre-admission residence, No. (%)") ("") ("")
local res1 "Home"
local res3 "Health-related institution"
local res4 "Non-health related institution"
local res6 "No fixed address/abode or temporary abode"
foreach r in 1 3 4 6 {
	forvalue i=0/1 {
		count if ISC_06==`r' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' ("  `res`r''") ("`val1'") ("`val0'")
}
* Primary site of infection
forvalues i=0/1 {
	count if !mi(primary_site) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("Primary site of infection, No. (%)") ("(n = `N1')") ("(n = `N0')")
local site1 "Respiratory"
local site2 "Gastrointestinal"
local site3 "Genitourinary"
local site4 "Haematological"
local site5 "Musculoskeletal"
local site6 "Dermatological"
local site7 "Bloodstream infection (including endocarditis and catheter-related infection)"
local site8 "Neurological"
local site9 "Surgical wound infection"
local site10 "Other*"
local site11 "Unknown"
foreach s in 1 2 3 4 5 6 7 8 9 11 {
	forvalue i=0/1 {
		count if primary_site==`s' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' ("  `site`s''") ("`val1'") ("`val0'")
}
* SSIP score
egen SSIPcat=cut(SSIP), at(0 5 7 11 23)
forvalues i=0/1 {
	count if !mi(SSIPcat) & alloc==`i'
	local N`i'=r(N)
}
post `pf' ("SSIP score, No. (%)") ("") ("")
local ssip0 "0-4"
local ssip5 "5-6"
local ssip7 "7-10"
local ssip11 "11+"
foreach s in 0 5 7 11 {
	forvalue i=0/1 {
		count if SSIPcat==`s' & alloc==`i'
		local val`i'=string(r(N),"%3.0f")+" ("+string(100*r(N)/`N`i'',"%4.1f")+")"
	}
	post `pf' ("  `ssip`s''") ("`val1'") ("`val0'")
}
* APACHE II score
forvalues i=0/1 {
	su ISC_11 if ISC_11<72 & alloc==`i'
	local val`i'=string(r(mean),"%4.1f")+" ("+string(r(sd),"%3.1f")+") [n = "+string(r(N))+"]"
}
post `pf' ("APACHE II score, mean (SD)") ("`val1'") ("`val0'")
* SOFA score
forvalues i=0/1 {
	su SOFAscore if alloc==`i', d
	local val`i'=string(r(p50))+" ("+string(r(p25))+", " +string(r(p75))+") [n = "+string(r(N))+"]"
}
post `pf' ("SOFA score, median (IQR)") ("`val1'") ("`val0'")

* White blood cell count
forvalues i=0/1 {
	su BL_01 if BL_01<1000 & alloc==`i', d
	local val`i'=string(r(p50),"%4.1f")+" ("+string(r(p25),"%4.1f")+", " +string(r(p75),"%4.1f")+") [n = "+string(r(N))+"]"
}
post `pf' ("White blood cell count, median (IQR)") ("`val1'") ("`val0'")

* CRP
forvalues i=0/1 {
	su BL_10 if BL_10<1000 & alloc==`i', d
	local val`i'=string(r(p50))+" ("+string(r(p25))+", " +string(r(p75))+") [n = "+string(r(N))+"]"
}
post `pf' ("CRP, median (IQR)") ("`val1'") ("`val0'")

postclose `pf'
use table1, clear
compress
save, replace
export excel using VACIRiSS_Final_`today'.xlsx, sheet(Table1, replace)
br
