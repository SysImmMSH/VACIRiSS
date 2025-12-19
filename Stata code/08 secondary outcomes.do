* 08 secondary outcomes.do
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
postfile `pf' str14 outcome str3 effect double(events1 events0 denom1 denom0 estimate lb ub abs abs_lb abs_ub pval) using secondary, replace
encode Site, gen(site)
log using "secondary outcomes.log", replace
foreach endpoint in rehosp rehospinf reinfection {
	foreach t in 30 90 180 365 {
		tab alloc `endpoint'`t', row nokey
		cs `endpoint'`t' alloc
		local a=r(rd)
		local al=r(lb_rd)
		local au=r(ub_rd)
		forvalues y=0/1 {
			qui count if `endpoint'`t'==1 & alloc==`y'
			local e`y'=r(N)
			qui count if !mi(`endpoint'`t') & alloc==`y'
			local d`y'=r(N)
		}
		glm `endpoint'`t' alloc age SSIP, family(bin) link(log) cluster(site) eform
		local r=exp(_b[alloc])
		local l=exp(_b[alloc]-invnormal(.975)*_se[alloc])
		local u=exp(_b[alloc]+invnormal(.975)*_se[alloc])
		local p=2*(1-normal(abs(_b[alloc]/_se[alloc])))
		post `pf' ("`endpoint'`t'") ("RR") (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (`a') (`al') (`au') (`p')
	}
}
foreach event in rehosp ab {
	stset T_`event', failure(F_`event')
	forvalues y=0/1 {
		qui count if F_`event'==1 & T_`event'>0 & alloc==`y'
		local e`y'=r(N)
		su T_`event' if T_`event'>0 & alloc==`y', meanonly
		local d`y'=r(sum)
	}
	stcox alloc age SSIP, shared(site) technique(dfp)
	local r=exp(_b[alloc])
	local l=exp(_b[alloc]-invnormal(.975)*_se[alloc])
	local u=exp(_b[alloc]+invnormal(.975)*_se[alloc])
	local p=2*(1-normal(abs(_b[alloc]/_se[alloc])))
	post `pf' ("`event'") ("HR") (`e1') (`e0') (`d1') (`d0') (`r') (`l') (`u') (.) (.) (.) (`p')
}
log close
postclose `pf'
