* 07 km secondary.do
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
stset T_ab, f(F_ab)
sts graph, by(alloc) failure ///
	risktable(,order(2 "PCV13" 1 "Placebo")) ///
	legend(order(2 "PCV13" 1 "Placebo") col(1) pos(5) ring(0) region(lc(none))) ///
	yla(0 "0" .2 "20" .4 "40" .6 "60" .8 "80", gmax angle(0)) ///
	yti("Cumulative incidence" "of antibiotics in general practice (%)") ///
	xla(0 30 90 180 365) ///
	xti("Days from randomisation") ///
	title("") graphr(c(white))
gr save km_ab_noci, replace
gr export km_ab_noci.png, replace width(2800)
