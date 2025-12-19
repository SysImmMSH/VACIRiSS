* 05 subgroups.do
* Author: David Harrison
* Date: 02/03/2023
* 
version 18
set more off
clear
cap log close
* base working directory stored in global
cd "$base"
use working, clear
tempname pf
postfile `pf' str7 outcome str3 effect double(N1 N0 events1 events0 denom1 denom0 estimate lb ub pval) using subgroups, replace
log using "subgroups.log", replace
stset T_primary, f(F_primary)
encode Site, gen(site)
egen SSIPcat=cut(SSIP), at(0 5 7 11 100)
forvalues y=0/1 {
	foreach s in 0 5 7 11 {
		qui count if F_primary==1 & T_primary>0 & alloc==`y' & SSIPcat==`s'
		local e`y'`s'=r(N)
		su T_primary if T_primary>0 & alloc==`y' & SSIPcat==`s', meanonly
		local d`y'`s'=r(sum)
		local n`y'`s'=r(N)
	}
}
stcox i.alloc##i.SSIPcat age, shared(site) technique(dfp)
test 1.alloc#5.SSIPcat 1.alloc#7.SSIPcat 1.alloc#11.SSIPcat
local p=r(p)
foreach s in 0 5 7 11 {
	lincom 1.alloc+1.alloc#`s'.SSIPcat, eform
	local r=r(estimate)
	local l=r(lb)
	local u=r(ub)
	post `pf' ("SSIP`s'") ("HR") (`n1`s'') (`n0`s'') (`e1`s'') (`e0`s'') (`d1`s'') (`d0`s'') (`r') (`l') (`u') (`p')
	local p=.
}
log close
postclose `pf'
