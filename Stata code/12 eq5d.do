* 12 eq5d.do
* Author: David Harrison
* Date: 02/03/2023
* 
version 18
set more off
clear
cap log close
set scheme s2color
local today=string(date(c(current_date),"DMY"),"%tdCCYY-NN-DD")
* base working directory stored in global
cd "$base"
tempname pf
postfile `pf' str5 outcome str2 timepoint double(mean1 mean0 sd1 sd0 n1 n0 dif lb ub pval) using eq5d, replace
use "1 Raw data\forms\EQ5", clear
merge 1:1 Label using working, keep(2 3) keepusing(RAN_02 date_of_withdrawal date_of_death date_of_follow_up alloc Site age REG_05 SSIP)
drop _merge
gen byte male=REG_05==1 if inlist(REG_05,1,2)
encode Site, gen(site)
forvalues i=1/5 {
	gen byte item`i'=EQ5_0`i' if EQ5_0`i'<6
}
gen EQ5_vas=EQ5_06 if EQ5_06<777 
eq5dmap EQ5_sindex, covariates(age male) items(item1 item2 item3 item4 item5)
gen EQ5_index=EQ5_sindex
local t1 30
local t2 90
local t3 180
local t4 365
forvalues j=1/4 {
	forvalues i=1/5 {
		gen byte item`i'`j'=EQ5_0`i'`j' if EQ5_0`i'`j'<6
	}
	gen EQ5_vas`j'=EQ5_06`j' if EQ5_06`j'<777 
	eq5dmap EQ5_sindex`j', covariates(age male) items(item1`j' item2`j' item3`j' item4`j' item5`j')
	gen EQ5_index`j'=EQ5_sindex`j'
	replace EQ5_index`j'=0 if date_of_death-RAN_02<=`t`j''
}
log using "07 eq5d.log", replace
tabstat EQ5_sindex*, by(alloc) s(n mean sd min q max)
tabstat EQ5_index*, by(alloc) s(n mean sd min q max)
tabstat EQ5_vas*, by(alloc) s(n mean sd min q max)
forvalues i=0/1 {
	su EQ5_vas if alloc==`i'
	local m`i'=r(mean)
	local s`i'=r(sd)
	local n`i'=r(N)
}
post `pf' ("VAS") ("t0") (`m1') (`m0') (`s1') (`s0') (`n1') (`n0') (.) (.) (.) (.)
forvalues i=0/1 {
	su EQ5_index if alloc==`i'
	local m`i'=r(mean)
	local s`i'=r(sd)
	local n`i'=r(N)
}
post `pf' ("Index") ("t0") (`m1') (`m0') (`s1') (`s0') (`n1') (`n0') (.) (.) (.) (.)
forvalues j=1/4 {
	forvalues i=0/1 {
		su EQ5_vas`j' if alloc==`i'
		local m`i'=r(mean)
		local s`i'=r(sd)
		local n`i'=r(N)
	}
	post `pf' ("VAS") ("t`j'") (`m1') (`m0') (`s1') (`s0') (`n1') (`n0') (.) (.) (.) (.)
	forvalues i=0/1 {
		su EQ5_index`j' if alloc==`i'
		local m`i'=r(mean)
		local s`i'=r(sd)
		local n`i'=r(N)
	}
	xtreg EQ5_index`j' alloc EQ5_index age SSIP, i(site) 
	local d=_b[alloc]
	local l=_b[alloc]-invnormal(.975)*_se[alloc]
	local u=_b[alloc]+invnormal(.975)*_se[alloc]
	local p=2*(1-normal(abs(_b[alloc]/_se[alloc])))
	post `pf' ("Index") ("t`j'") (`m1') (`m0') (`s1') (`s0') (`n1') (`n0') (`d') (`l') (`u') (`p')
}
log close
postclose `pf'

forvalues i=1/5 {
	ren item`i' item`i'0
}
gen byte n=1
preserve
tempfile f
forvalues i=1/5 {
	forvalues j=0/4 {
		drop if mi(item`i'`j')
		collapse (sum) n, by(alloc item`i'`j')
		su n if alloc==0, meanonly
		gen pc=100*n/r(sum) if alloc==0
		su n if alloc==1, meanonly
		replace pc=100*n/r(sum) if alloc==1
		ren item`i'`j' response
		ren n n`i'`j'
		ren pc pc`i'`j'
		cap nois merge 1:1 alloc response using `f', nogen
		save `f', replace
		restore, preserve
	}
}
restore
use `f', clear
forvalues i=1/5 {
	forvalues j=0/4 {
		replace n`i'`j'=0 if mi(n`i'`j')
		replace pc`i'`j'=0 if mi(pc`i'`j')
	}
}
reshape long n1 pc1 n2 pc2 n3 pc3 n4 pc4 n5 pc5, i(alloc response) j(timepoint)
reshape wide n1 pc1 n2 pc2 n3 pc3 n4 pc4 n5 pc5, i(alloc timepoint) j(response)
save eq5dresponses, replace
use eq5dresponses, clear
* Plot
local title1 "Mobility"
local item11 "No problems in walking about"
local item12 "Slight problems in walking about"
local item13 "Moderate problems in walking about"
local item14 "Severe problems in walking about"
local item15 "Unable to walk about"
local title2 "Self-care"
local item21 "No problems washing or dressing"
local item22 "Slight problems washing or dressing"
local item23 "Moderate problems washing or dressing"
local item24 "Severe problems washing or dressing"
local item25 "Unable to wash or dress"
local title3 "Usual activities"
local item31 "No problems doing usual activities"
local item32 "Slight problems doing usual activities"
local item33 "Moderate problems doing usual activities"
local item34 "Severe problems doing usual activities"
local item35 "Unable to do usual activities"
local title4 "Pain / Discomfort"
local item41 "No pain or discomfort"
local item42 "Slight pain or discomfort"
local item43 "Moderate pain or discomfort"
local item44 "Severe pain or discomfort"
local item45 "Extreme pain or discomfort"
local title5 "Anxiety / Depression"
local item51 "Not anxious or depressed"
local item52 "Slightly anxious or depressed"
local item53 "Moderately anxious or depressed"
local item54 "Severely anxious or depressed"
local item55 "Extremely anxious or depressed"
gen zero=0
gen tot=100
gen l=timepoint-.2
gen r=timepoint+.2
forvalues i=1/5 {
	gen cut1=pc`i'1+pc`i'2
	gen cut2=cut1+pc`i'3
	gen cut3=cut2+pc`i'4
	twoway (rbar pc`i'1 zero l if alloc==1, barw(.4) fintens(10)) ///
		(rbar pc`i'1 zero r if alloc==0, barw(.4) fintens(10)) ///
		(rbar cut1 pc`i'1 l if alloc==1, barw(.4) pstyle(p1) fintens(30)) ///
		(rbar cut1 pc`i'1 r if alloc==0, barw(.4) pstyle(p2) fintens(30)) ///
		(rbar cut2 cut1 l if alloc==1, barw(.4) pstyle(p1) fintens(55)) ///
		(rbar cut2 cut1 r if alloc==0, barw(.4) pstyle(p2) fintens(55)) ///
		(rbar cut3 cut2 l if alloc==1, barw(.4) pstyle(p1) fintens(80)) ///
		(rbar cut3 cut2 r if alloc==0, barw(.4) pstyle(p2) fintens(80)) ///
		(rbar tot cut3 l if alloc==1, barw(.4) pstyle(p1) fintens(100)) ///
		(rbar tot cut3 r if alloc==0, barw(.4) pstyle(p2) fintens(100)), ///
		xla(0 "Baseline" 1 "30 days" 2 "90 days" 3 "180 days" 4 "1 year") xti("") ///
		yla(, angle(0)) yti(Percentage of survivors) ///
		legend(order(- "`item`i'5'" 9 "" 10 "" - "`item`i'4'" 7 "" 8 "" - "`item`i'3'" 5 "" 6 "" - "`item`i'2'" 3 "" 4 "" - "`item`i'1'" 1 "" 2 "" - " ") ///
		col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) rowgap(0) colgap(2) size(small)) ///
		plotr(margin(b 0)) graphr(color(white)) title("`=word(c(ALPHA),`i')'", pos(10) size(huge)) subtitle("`title`i''", size(medium) color(black)) name(`=word(c(ALPHA),`i')', replace)
	drop cut1 cut2 cut3
}
graph combine A B C D E, col(2) imargin(tiny) iscale(*.75) graphr(color(white)) ysize(12) xsize(9)
gr save eq5d, replace
gr export eq5d.svg, replace
gr export eq5d.eps, replace
gr export eq5d.pdf, replace
gr export eq5d.png, replace width(2800)

* Create EQ5D table
use eq5d, clear
*row order
gen byte y=.
local y 2
forvalues j=0/4 {
	replace y=`y++' if outcome=="Index" & timepoint=="t`j'"
}
local y=`y'+1
forvalues j=0/4 {
	replace y=`y++' if outcome=="VAS" & timepoint=="t`j'"
}
* column 1 - Labels
gen col1=""
replace col1="  Baseline" if timepoint=="t0"
replace col1="  30 days" if timepoint=="t1"
replace col1="  90 days" if timepoint=="t2"
replace col1="  180 days" if timepoint=="t3"
replace col1="  365 days" if timepoint=="t4"
set obs `=_N+1'
replace y=1 in l
replace col1="EQ-5D utility" in l
set obs `=_N+1'
replace y=7 in l
replace col1="EQ-5D visual analog scale" in l
sort y
* column 2 - mean (SD) [N] - PCV-13
gen col2=string(mean1,"%4.2f")+" ("+string(sd1,"%4.2f")+") ["+string(n1)+"]" if !mi(mean1) & outcome=="Index"
replace col2=string(mean1,"%4.1f")+" ("+string(sd1,"%4.1f")+") ["+string(n1)+"]" if !mi(mean1) & outcome=="VAS"
* column 3 - mean (SD) [N] - Placebo
gen col3=string(mean0,"%4.2f")+" ("+string(sd0,"%4.2f")+") ["+string(n0)+"]" if !mi(mean0) & outcome=="Index"
replace col3=string(mean0,"%4.1f")+" ("+string(sd0,"%4.1f")+") ["+string(n0)+"]" if !mi(mean1) & outcome=="VAS"
* column 4 - mean difference (95% CI)
gen col4=string(dif,"%4.2f")+" ("+string(lb,"%4.2f")+", "+string(ub,"%4.2f")+")" if !mi(dif)
* column 5 - p value
gen col5="<0.001" if pval<.001
replace col5=string(pval,"%5.3f") if inrange(pval,.001,.01)
replace col5=string(pval,"%4.2f") if pval>=.01 & !mi(pval)
replace col5=">0.99" if pval>.99 & !mi(pval)

label var col1 "Outcome"
label var col2 "PCV13"
label var col3 "Placebo"
label var col4 "Mean difference (95% CI)"
label var col5 "P value"
export excel col1 col2 col3 col4 col5 using VACIRiSS_Final_`today'.xlsx if !mi(y), first(varl) sheet(EQ5D, replace)
