* 06 km primary.do
* Author: David Harrison
* Date: 06/12/2024
* 
version 18
set more off
clear
cap log close
set scheme s2color
* base working directory stored in global
cd "$base"
use working, clear
label define alloc 1 "PCV13", modify
stset T_primary, f(F_primary)
sts graph, by(alloc) failure ///
	risktable(,order(2 "PCV13" 1 "Placebo")) ///
	legend(order(2 "PCV13" 1 "Placebo") col(1) pos(5) ring(0) region(lc(none))) ///
	yla(0 "0" .1 "10" .2 "20" .3 "30" .4 "40" .5 "50", gmax angle(0)) ///
	yti("Cumulative incidence" "of rehospitalisation or death (%)") ///
	xla(0 30 90 180 365) ///
	xti("Days from randomisation") ///
	title("") graphr(c(white))
gr save km_noci, replace
gr export km_noci.png, replace width(2800)
