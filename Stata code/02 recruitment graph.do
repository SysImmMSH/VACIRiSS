* 02 recruitment graph.do
* Author: David Harrison
* Date: 22/07/2024
* 
version 17
set more off
clear
cap log close
set scheme s2color
* base working directory stored in global
cd "$base"
use "1 Raw data\recruitment", clear
twoway (pci 220 716 0 716, lc(gs8) lp(dash)) ///
	(pci 220 724 0 724, lc(gs14) lw(1.1cm)) ///
	(line expected observed month, pstyle(p1 p2)), ///
	legend(order(3 "Expected" 4 "Actual") col(1) pos(5) ring(0) region(lc(none))) ///
	xla(703(6)747, angle(45) format(%tmm_CY)) ///
	yla(0(50)200 214, angle(0) gmax) ///
	graphr(c(white)) plotr(margin(b 0 t 0)) ///
	yti(Cumulative number of participants) xti("") ///
	text(207 716 "Internal pilot ", place(w) size(small)) ///
	text(90 726 "Sponsor paused recruitment" "Mar 20, 2020 to Jul 6, 2020", place(e) size(small))
gr export recruitment.png, replace width(2800)

