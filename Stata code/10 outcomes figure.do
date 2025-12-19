* 10 outcomes figure.do
* Author: David Harrison
* Date: 06/12/2024
* 
version 18
set more off
clear
cap log close
set scheme s2color
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"
* First figure - primary outcome, components and age subgroups
use primary, clear
append using secondary
* row order
gen byte y=.
replace y=2 if outcome=="primary - over 65"
replace y=4 if outcome=="primary - 65 and under"
replace y=7 if outcome=="death"
replace y=9 if outcome=="infrel"
replace y=11 if outcome=="primary" & effect=="HR"
replace y=13 if outcome=="primary" & effect=="IRR"
drop if mi(y)
* add rows for absolute effects
set obs `=_N+1'
replace y=1 in l
set obs `=_N+1'
replace y=3 in l
set obs `=_N+1'
replace y=6 in l
set obs `=_N+1'
replace y=8 in l
set obs `=_N+1'
replace y=12 in l
replace effect="IRD" if mi(effect)
gsort -y
replace estimate=abs[_n-1] if effect=="IRD"
replace lb=abs_lb[_n-1] if effect=="IRD"
replace ub=abs_ub[_n-1] if effect=="IRD"
* column 1 - outcome names
gen col1=""
replace col1="  Time to death" if outcome=="death"
replace col1="  Time to first infection-related rehospitalisation" if outcome=="infrel"
replace col1="  >65 years" if outcome=="primary - over 65"
replace col1="  ≤65 years" if outcome=="primary - 65 and under"
replace col1="death (primary outcome)" if y==12
replace col1="Time to first infection-related rehospitalisation or" if outcome=="primary" & effect=="IRR"
set obs `=_N+1'
replace y=5 in l
replace col1="Primary outcome by age strata" in l
set obs `=_N+1'
replace y=10 in l
replace col1="Components of primary outcome" in l
set obs `=_N+1'
replace y=14 in l
replace col1="{bf:Outcome}" in l
set obs `=_N+1'
replace y=15 in l
gsort -y
* column 2 - percent or rate - PCV-13
gen col2=""
replace col2=string(events1)+"/"+string(denom1/365.25,"%4.1f")+" ("+string(365.25*events1/denom1,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col2=string(events1)+"/"+string(round(denom1))+" ("+string(100*events1/denom1,"%4.1f")+")" if effect=="RR"
replace col2="{bf:PCV13}" if y==14
* column 3 - percent or rate - placebo
gen col3=""
replace col3=string(events0)+"/"+string(denom0/365.25,"%4.1f")+" ("+string(365.25*events0/denom0,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col3=string(events0)+"/"+string(round(denom0))+" ("+string(100*events0/denom0,"%4.1f")+")" if effect=="RR"
replace col3="{bf:Placebo}" if y==14
* column 4 - treatment effect
gen col4=effect+" "+string(estimate,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if inlist(effect,"IRR","HR","RR","IRD")
replace col4=effect+" "+string(100*estimate,"%4.1f")+" ("+string(100*lb,"%4.1f")+", "+string(100*ub,"%4.1f")+")" if effect=="RD"
replace col4="{bf:(95% CI)}" if y==14
replace col4="{bf:Treatment effect}" if y==15
* x positions for text columns
gen x1=0
replace x1=.25 if substr(col1,1,1)=="d"
replace x1=.5 if substr(col1,1,1)==" "
gen x2=8.36
gen x3=11.08
gen x4=14.04
*rescale points on log scale so 1/5 = 15.7 and 5 = 24.74
local m=(24.74-15.7)/(2*ln(5))
local c=(24.74+15.7)/2
gen x=`m'*ln(estimate)+`c' if inlist(effect,"IRR","HR","RR")
gen xl=`m'*ln(lb)+`c' if inlist(effect,"IRR","HR","RR")
gen xu=`m'*ln(ub)+`c' if inlist(effect,"IRR","HR","RR")
drop if mi(y)
twoway (pci 15.5 0 15.5 24.8 -2.3 0 -2.3 24.8 13.5 0 13.5 15.64, lc(black) lw(medium)) ///
	(pci 14.5 7 14.5 12.44, lc(black) lw(thin)) ///
	(pci 10.5 0 10.5 15.64, lc(gs8) lw(thin)) ///
	(pci 5.5 0 5.5 15.64, lc(gs8) lw(thin)) ///
	(pci 0 15.7 0 24.74, lc(black) lw(thin)) ///
	(pci 0 15.7 -.3 15.7 0 18.27 -.3 18.27 0 20.22 -.3 20.22 0 22.17 -.3 22.17 0 24.74 -.3 24.74, lc(black) lw(thin)) ///
	(pci 0 20.22 14.5 20.22, lc(black) lw(thin) lp("-")) ///
	(scatter y x1, ms(i) mlabel(col1) mlabp(3) mlabgap(0) mlabc(black)) ///
	(scatter y x2, ms(i) mlabel(col2) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x3, ms(i) mlabel(col3) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x4, ms(i) mlabel(col4) mlabp(0) mlabgap(0) mlabc(black)) ///
	(rspike xl xu y, hor lc(black) lw(thin)) ///
	(scatter y x, pstyle(p1) ms(S)) ///
	, text(15 9.72 "{bf:No. / person-yrs at risk (rate)}" -1.8 20.22 "Treatment effect (95% CI)", size(small)) ///
	text(-.8 15.7 "0.2" -.8 18.27 "0.5" -.8 20.22 "1" -.8 22.17 "2" -.8 24.74 "5", size(small)) ///
	text(14 20.18 "Favours PCV13", size(small) place(w)) ///
	text(14 20.26 "Favours placebo", size(small) place(e)) ///
	scale(*1.25) legend(off) yscale(off) yla(none) xscale(off) xla(none) ///
	graphr(col(white) margin(zero)) ///
	ysize(10) xsize(28) ///
	title("{bf:A}", pos(11) col(black)) name(A, replace)

* Second figure - All-cause rehospitalisation and time to 
use primary, clear
append using secondary
* row order
gen y=.
replace y=1 if outcome=="rehosp"
replace y=3 if outcome=="rehosp365"
replace y=5 if outcome=="rehosp180"
replace y=7 if outcome=="rehosp90"
replace y=9 if outcome=="rehosp30"
drop if mi(y)
* add rows for absolute effects
set obs `=_N+1'
replace y=2 in l
set obs `=_N+1'
replace y=4 in l
set obs `=_N+1'
replace y=6 in l
set obs `=_N+1'
replace y=8 in l
replace effect="RD" if mi(effect)
gsort -y
replace estimate=abs[_n-1] if effect=="RD"
replace lb=abs_lb[_n-1] if effect=="RD"
replace ub=abs_ub[_n-1] if effect=="RD"
* column 1 - outcome names
gen col1=""
replace col1="Time to first all-cause rehospitalisation" if outcome=="rehosp"
replace col1="  By 365 days" if index(outcome,"365")
replace col1="  By 180 days" if index(outcome,"180")
replace col1="  By 90 days" if index(outcome,"90")
replace col1="  By 30 days" if index(outcome,"30")
set obs `=_N+1'
replace y=10 in l
replace col1="All-cause rehospitalisation" in l
set obs `=_N+1'
replace y=11 in l
replace col1="{bf:Outcome}" in l
set obs `=_N+1'
replace y=12 in l
gsort -y
* column 2 - percent or rate - PCV-13
gen col2=""
replace col2=string(events1)+"/"+string(denom1/365.25,"%4.1f")+" ("+string(365.25*events1/denom1,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col2=string(events1)+"/"+string(round(denom1))+" ("+string(100*events1/denom1,"%4.1f")+")" if effect=="RR"
replace col2="{bf:PCV13}" if y==11
* column 3 - percent or rate - placebo
gen col3=""
replace col3=string(events0)+"/"+string(denom0/365.25,"%4.1f")+" ("+string(365.25*events0/denom0,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col3=string(events0)+"/"+string(round(denom0))+" ("+string(100*events0/denom0,"%4.1f")+")" if effect=="RR"
replace col3="{bf:Placebo}" if y==11
* column 4 - treatment effect
gen col4=effect+" "+string(estimate,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if inlist(effect,"IRR","HR","RR","IRD")
replace col4=effect+" "+string(100*estimate,"%4.1f")+" ("+string(100*lb,"%4.1f")+", "+string(100*ub,"%4.1f")+")" if effect=="RD"
replace col4="{bf:(95% CI)}" if y==11
replace col4="{bf:Treatment effect}" if y==12
* x positions for text columns
gen x1=0
replace x1=.25 if substr(col1,1,1)=="d"
replace x1=.5 if substr(col1,1,1)==" "
gen x2=8.36
gen x3=11.08
gen x4=14.04
*rescale points on log scale so 1/5 = 15.7 and 5 = 24.74
local m=(24.74-15.7)/(2*ln(5))
local c=(24.74+15.7)/2
gen x=`m'*ln(estimate)+`c' if inlist(effect,"IRR","HR","RR")
gen xl=`m'*ln(lb)+`c' if inlist(effect,"IRR","HR","RR")
gen xu=`m'*ln(ub)+`c' if inlist(effect,"IRR","HR","RR")
drop if mi(y)
twoway (pci 13.5 0 13.5 24.8 -2.3 0 -2.3 24.8 10.5 0 10.5 15.64, lc(black) lw(medium)) ///
	(pci 11.5 7 11.5 12.44, lc(black) lw(thin)) ///
	(pci 1.5 0 1.5 15.64, lc(gs8) lw(thin)) ///
	(pci 0 15.7 0 24.74, lc(black) lw(thin)) ///
	(pci 0 15.7 -.3 15.7 0 18.27 -.3 18.27 0 20.22 -.3 20.22 0 22.17 -.3 22.17 0 24.74 -.3 24.74, lc(black) lw(thin)) ///
	(pci 0 20.22 11.5 20.22, lc(black) lw(thin) lp("-")) ///
	(scatter y x1, ms(i) mlabel(col1) mlabp(3) mlabgap(0) mlabc(black)) ///
	(scatter y x2, ms(i) mlabel(col2) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x3, ms(i) mlabel(col3) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x4, ms(i) mlabel(col4) mlabp(0) mlabgap(0) mlabc(black)) ///
	(rspike xl xu y, hor lc(black) lw(thin)) ///
	(scatter y x, pstyle(p1) ms(S)) ///
	(scatteri 14.5 0 -3.3 0, ms(i)) ///
	, text(13 9.72 "{bf:No. / total No. or person-yrs at risk}" 12 9.72 "{bf:(% or rate)}" -1.8 20.22 "Treatment effect (95% CI)", size(small)) ///
	text(-.8 15.7 "0.2" -.8 18.27 "0.5" -.8 20.22 "1" -.8 22.17 "2" -.8 24.74 "5", size(small)) ///
	text(11 20.18 "Favours PCV13", size(small) place(w)) ///
	text(11 20.26 "Favours placebo", size(small) place(e)) ///
	scale(*1.25) legend(off) yscale(off) yla(none) xscale(off) xla(none) ///
	graphr(col(white) margin(zero)) ///
	ysize(10) xsize(28) ///
	title("{bf:B}", pos(11) col(black)) name(B, replace)

* Third figure - Infection-related rehospitalisation
use primary, clear
append using secondary
* row order
gen y=.
replace y=2 if outcome=="rehospinf365"
replace y=4 if outcome=="rehospinf180"
replace y=6 if outcome=="rehospinf90"
replace y=8 if outcome=="rehospinf30"
drop if mi(y)
* add rows for absolute effects
set obs `=_N+1'
replace y=1 in l
set obs `=_N+1'
replace y=3 in l
set obs `=_N+1'
replace y=5 in l
set obs `=_N+1'
replace y=7 in l
replace effect="RD" if mi(effect)
gsort -y
replace estimate=abs[_n-1] if effect=="RD"
replace lb=abs_lb[_n-1] if effect=="RD"
replace ub=abs_ub[_n-1] if effect=="RD"
* column 1 - outcome names
gen col1=""
replace col1="  By 365 days" if index(outcome,"365")
replace col1="  By 180 days" if index(outcome,"180")
replace col1="  By 90 days" if index(outcome,"90")
replace col1="  By 30 days" if index(outcome,"30")
set obs `=_N+1'
replace y=9 in l
replace col1="Infection-related rehospitalisation" in l
set obs `=_N+1'
replace y=10 in l
replace col1="{bf:Outcome}" in l
set obs `=_N+1'
replace y=11 in l
gsort -y
* column 2 - percent or rate - PCV-13
gen col2=""
replace col2=string(events1)+"/"+string(denom1/365.25,"%4.1f")+" ("+string(365.25*events1/denom1,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col2=string(events1)+"/"+string(round(denom1))+" ("+string(100*events1/denom1,"%4.1f")+")" if effect=="RR"
replace col2="{bf:PCV13}" if y==10
* column 3 - percent or rate - placebo
gen col3=""
replace col3=string(events0)+"/"+string(denom0/365.25,"%4.1f")+" ("+string(365.25*events0/denom0,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col3=string(events0)+"/"+string(round(denom0))+" ("+string(100*events0/denom0,"%4.1f")+")" if effect=="RR"
replace col3="{bf:Placebo}" if y==10
* column 4 - treatment effect
gen col4=effect+" "+string(estimate,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if inlist(effect,"IRR","HR","RR","IRD")
replace col4=effect+" "+string(100*estimate,"%4.1f")+" ("+string(100*lb,"%4.1f")+", "+string(100*ub,"%4.1f")+")" if effect=="RD"
replace col4="{bf:(95% CI)}" if y==10
replace col4="{bf:Treatment effect}" if y==11
* x positions for text columns
gen x1=0
replace x1=.25 if substr(col1,1,1)=="d"
replace x1=.5 if substr(col1,1,1)==" "
gen x2=8.36
gen x3=11.08
gen x4=14.04
*rescale points on log scale so 1/5 = 15.7 and 5 = 24.74
local m=(24.74-15.7)/(2*ln(5))
local c=(24.74+15.7)/2
gen x=`m'*ln(estimate)+`c' if inlist(effect,"IRR","HR","RR")
gen xl=`m'*ln(lb)+`c' if inlist(effect,"IRR","HR","RR")
gen xu=`m'*ln(ub)+`c' if inlist(effect,"IRR","HR","RR")
drop if mi(y)
twoway (pci 11.5 0 11.5 24.8 -2.3 0 -2.3 24.8 9.5 0 9.5 15.64, lc(black) lw(medium)) ///
	(pci 10.5 7 10.5 12.44, lc(black) lw(thin)) ///
	(pci 0 15.7 0 24.74, lc(black) lw(thin)) ///
	(pci 0 15.7 -.3 15.7 0 18.27 -.3 18.27 0 20.22 -.3 20.22 0 22.17 -.3 22.17 0 24.74 -.3 24.74, lc(black) lw(thin)) ///
	(pci 0 20.22 10.5 20.22, lc(black) lw(thin) lp("-")) ///
	(scatter y x1, ms(i) mlabel(col1) mlabp(3) mlabgap(0) mlabc(black)) ///
	(scatter y x2, ms(i) mlabel(col2) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x3, ms(i) mlabel(col3) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x4, ms(i) mlabel(col4) mlabp(0) mlabgap(0) mlabc(black)) ///
	(rspike xl xu y, hor lc(black) lw(thin)) ///
	(scatter y x, pstyle(p1) ms(S)) ///
	(scatteri 13.5 0 -4.3 0, ms(i)) ///
	, text(11 9.72 "{bf:No. / total No. (%)}" -1.8 20.22 "Treatment effect (95% CI)", size(small)) ///
	text(-.8 15.7 "0.2" -.8 18.27 "0.5" -.8 20.22 "1" -.8 22.17 "2" -.8 24.74 "5", size(small)) ///
	text(10 20.18 "Favours PCV13", size(small) place(w)) ///
	text(10 20.26 "Favours placebo", size(small) place(e)) ///
	scale(*1.25) legend(off) yscale(off) yla(none) xscale(off) xla(none) ///
	graphr(col(white) margin(zero)) ///
	ysize(10) xsize(28) ///
	title("{bf:C}", pos(11) col(black)) name(C, replace)

* Fourth figure - Reinfection and time to abx
use primary, clear
append using secondary
* row order
gen y=.
replace y=1 if outcome=="ab"
replace y=3 if outcome=="reinfection365"
replace y=5 if outcome=="reinfection180"
replace y=7 if outcome=="reinfection90"
replace y=9 if outcome=="reinfection30"
drop if mi(y)
* add rows for absolute effects
set obs `=_N+1'
replace y=2 in l
set obs `=_N+1'
replace y=4 in l
set obs `=_N+1'
replace y=6 in l
set obs `=_N+1'
replace y=8 in l
replace effect="RD" if mi(effect)
gsort -y
replace estimate=abs[_n-1] if effect=="RD"
replace lb=abs_lb[_n-1] if effect=="RD"
replace ub=abs_ub[_n-1] if effect=="RD"
* column 1 - outcome names
gen col1=""
replace col1="Time to first antibiotic therapy in general practice" if outcome=="ab"
replace col1="  By 365 days" if index(outcome,"365")
replace col1="  By 180 days" if index(outcome,"180")
replace col1="  By 90 days" if index(outcome,"90")
replace col1="  By 30 days" if index(outcome,"30")
set obs `=_N+1'
replace y=10 in l
replace col1="Reinfection" in l
set obs `=_N+1'
replace y=11 in l
replace col1="{bf:Outcome}" in l
set obs `=_N+1'
replace y=12 in l
gsort -y
* column 2 - percent or rate - PCV-13
gen col2=""
replace col2=string(events1)+"/"+string(denom1/365.25,"%4.1f")+" ("+string(365.25*events1/denom1,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col2=string(events1)+"/"+string(round(denom1))+" ("+string(100*events1/denom1,"%4.1f")+")" if effect=="RR"
replace col2="{bf:PCV13}" if y==11
* column 3 - percent or rate - placebo
gen col3=""
replace col3=string(events0)+"/"+string(denom0/365.25,"%4.1f")+" ("+string(365.25*events0/denom0,"%4.2f")+")" if effect=="IRR"|(effect=="HR" & outcome!="primary")
replace col3=string(events0)+"/"+string(round(denom0))+" ("+string(100*events0/denom0,"%4.1f")+")" if effect=="RR"
replace col3="{bf:Placebo}" if y==11
* column 4 - treatment effect
gen col4=effect+" "+string(estimate,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if inlist(effect,"IRR","HR","RR","IRD")
replace col4=effect+" "+string(100*estimate,"%4.1f")+" ("+string(100*lb,"%4.1f")+", "+string(100*ub,"%4.1f")+")" if effect=="RD"
replace col4="{bf:(95% CI)}" if y==11
replace col4="{bf:Treatment effect}" if y==12
* x positions for text columns
gen x1=0
replace x1=.25 if substr(col1,1,1)=="d"
replace x1=.5 if substr(col1,1,1)==" "
gen x2=8.36
gen x3=11.08
gen x4=14.04
*rescale points on log scale so 1/5 = 15.7 and 5 = 24.74
local m=(24.74-15.7)/(2*ln(5))
local c=(24.74+15.7)/2
gen x=`m'*ln(estimate)+`c' if inlist(effect,"IRR","HR","RR")
gen xl=`m'*ln(lb)+`c' if inlist(effect,"IRR","HR","RR")
gen xu=`m'*ln(ub)+`c' if inlist(effect,"IRR","HR","RR")
drop if mi(y)
twoway (pci 13.5 0 13.5 24.8 -2.3 0 -2.3 24.8 10.5 0 10.5 15.64, lc(black) lw(medium)) ///
	(pci 11.5 7 11.5 12.44, lc(black) lw(thin)) ///
	(pci 1.5 0 1.5 15.64, lc(gs8) lw(thin)) ///
	(pci 0 15.7 0 24.74, lc(black) lw(thin)) ///
	(pci 0 15.7 -.3 15.7 0 18.27 -.3 18.27 0 20.22 -.3 20.22 0 22.17 -.3 22.17 0 24.74 -.3 24.74, lc(black) lw(thin)) ///
	(pci 0 20.22 11.5 20.22, lc(black) lw(thin) lp("-")) ///
	(scatter y x1, ms(i) mlabel(col1) mlabp(3) mlabgap(0) mlabc(black)) ///
	(scatter y x2, ms(i) mlabel(col2) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x3, ms(i) mlabel(col3) mlabp(0) mlabgap(0) mlabc(black)) ///
	(scatter y x4, ms(i) mlabel(col4) mlabp(0) mlabgap(0) mlabc(black)) ///
	(rspike xl xu y, hor lc(black) lw(thin)) ///
	(scatter y x, pstyle(p1) ms(S)) ///
	(scatteri 14.5 0 -3.3 0, ms(i)) ///
	, text(13 9.72 "{bf:No. / total No. or person-yrs at risk}" 12 9.72 "{bf:(% or rate)}" -1.8 20.22 "Treatment effect (95% CI)", size(small)) ///
	text(-.8 15.7 "0.2" -.8 18.27 "0.5" -.8 20.22 "1" -.8 22.17 "2" -.8 24.74 "5", size(small)) ///
	text(11 20.18 "Favours PCV13", size(small) place(w)) ///
	text(11 20.26 "Favours placebo", size(small) place(e)) ///
	scale(*1.25) legend(off) yscale(off) yla(none) xscale(off) xla(none) ///
	graphr(col(white) margin(zero)) ///
	ysize(10) xsize(28) ///
	title("{bf:D}", pos(11) col(black)) name(D, replace)
graph combine A B C D, imargin(tiny) cols(1) iscale(*.45) graphr(col(white) margin(tiny)) ysize(40) xsize(28)
gr save outcome_figure_4panelv3, replace
gr export outcome_figure_4panelv3.png, replace width(1200)
gr export outcome_figure_4panelv3.pdf, replace
