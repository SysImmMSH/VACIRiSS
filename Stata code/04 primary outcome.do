* 04 primary outcome.do
* Author: David Harrison
* Date: 02/03/2023
* 
version 18
set more off
clear
cap log close
set scheme s2color
* base working directory stored in global
cd "$base"
use working, clear
tempname pf
postfile `pf' str25 outcome str3 effect double(N1 N0 events1 events0 denom1 denom0 estimate lb ub abs abs_lb abs_ub pval) using primary, replace
log using "primary outcome.log", replace
foreach event in primary death infrel {
	gen Tyears=T_`event'/365.25
	bys alloc: ci means F_`event' if Tyears>0, poisson exp(Tyears)
	forvalues y=0/1 {
		qui count if F_`event'==1 & T_`event'>0 & alloc==`y'
		local e`y'=r(N)
		su T_`event' if T_`event'>0 & alloc==`y', meanonly
		local d`y'=r(sum)
		local n`y'=r(N)
	}
	ir F_`event' alloc Tyears, by(age65) istandard ird
	local a=r(ird)
	local al=r(lb_ird)
	local au=r(ub_ird)
	ir F_`event' alloc Tyears, by(age65) istandard
	local r=r(irr)
	local l=r(lb_irr)
	local u=r(ub_irr)
	post `pf' ("`event'") ("IRR") (`n1') (`n0') (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (`a') (`al') (`au') (.)
	*M-H for test of homogeneity only
	ir F_`event' alloc Tyears, by(age65)
	local p=chi2tail(1,r(chi2_mh))
	forvalues y=0/1 {
		qui count if F_`event'==1 & T_`event'>0 & alloc==`y' & age65==0
		local e`y'=r(N)
		su T_`event' if T_`event'>0 & alloc==`y' & age65==0, meanonly
		local d`y'=r(sum)
		local n`y'=r(N)
	}
	ir F_`event' alloc Tyears if age65==0
	local a=r(ird)
	local al=r(lb_ird)
	local au=r(ub_ird)
	local r=r(irr)
	local l=r(lb_irr)
	local u=r(ub_irr)
	post `pf' ("`event' - 65 and under") ("IRR") (`n1') (`n0') (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (`a') (`al') (`au') (`p')
	forvalues y=0/1 {
		qui count if F_`event'==1 & T_`event'>0 & alloc==`y' & age65==1
		local e`y'=r(N)
		su T_`event' if T_`event'>0 & alloc==`y' & age65==1, meanonly
		local d`y'=r(sum)
		local n`y'=r(N)
	}
	ir F_`event' alloc Tyears if age65==1
	local a=r(ird)
	local al=r(lb_ird)
	local au=r(ub_ird)
	local r=r(irr)
	local l=r(lb_irr)
	local u=r(ub_irr)
	post `pf' ("`event' - over 65") ("IRR") (`n1') (`n0') (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (`a') (`al') (`au') (.)
	drop Tyears
}
stset T_primary, f(F_primary)
encode Site, gen(site)
forvalues y=0/1 {
	qui count if F_primary==1 & T_primary>0 & alloc==`y'
	local e`y'=r(N)
	su T_primary if T_primary>0 & alloc==`y', meanonly
	local d`y'=r(sum)
	local n`y'=r(N)
}
stcox alloc age SSIP, shared(site) technique(dfp)
local r=exp(_b[alloc])
local l=exp(_b[alloc]-invnormal(.975)*_se[alloc])
local u=exp(_b[alloc]+invnormal(.975)*_se[alloc])
local p=2*(1-normal(abs(_b[alloc]/_se[alloc])))
post `pf' ("primary") ("HR") (`n1') (`n0') (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (.) (.) (.) (`p')
* check proportional hazards
stphplot, by(alloc) adjustfor(age SSIP)
graph save stphplot, replace
graph export stphplot.pdf, replace
* sensitivity analysis pre-covid
gen T_covid=T_primary
gen F_covid=F_primary
replace F_covid=0 if RAN_02+T_primary>d(29feb2020)
replace T_covid=d(29feb2020)-RAN_02 if RAN_02+T_primary>d(29feb2020)
replace F_covid=. if RAN_02>d(1feb2020)
replace T_covid=. if RAN_02>d(1feb2020)
stset T_covid, f(F_covid)
forvalues y=0/1 {
	qui count if F_covid==1 & T_covid>0 & alloc==`y'
	local e`y'=r(N)
	su T_covid if T_covid>0 & alloc==`y', meanonly
	local d`y'=r(sum)
	local n`y'=r(N)
}
stcox alloc age SSIP, shared(site) technique(dfp)
local r=exp(_b[alloc])
local l=exp(_b[alloc]-invnormal(.975)*_se[alloc])
local u=exp(_b[alloc]+invnormal(.975)*_se[alloc])
local p=2*(1-normal(abs(_b[alloc]/_se[alloc])))
post `pf' ("covid") ("HR") (`n1') (`n0') (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (.) (.) (.) (`p')
log close
postclose `pf'
