* 13 vaccine diary.do
* Author: David Harrison
* Date: 01/03/2023
version 18
set more off
clear
set scheme s2color
* base working directory stored in global
cd "$base"
foreach sheet in "GSTT" "UCLH" "Portsmouth " "KCH" "Edinburgh" "Oxford" "Newport" " Surrey" "Cambridge" "Belfast" {
		import excel using "1 Raw data\Vaccine diary collated V1.xlsx", sheet("`sheet'") first clear
		ren VACIRiSSID PIN
		replace PIN=PIN[_n-1] if mi(PIN)
		gen site=trim("`sheet'")
		if "`sheet'"!="GSTT" {
			append using vaccine_diary
		}
		save vaccine_diary, replace
}
egen nmissing=rowmiss(aPainseverity bSwellingandorrednesssever CLimitationofarmmovementse DFeverseverity EMuscleandjointpainseverit fHeadacheseverity)
replace aPainseverity="1 Mild pain " if mi(aPainseverity) & nmissing==1
replace bSwellingandorrednesssever="None  " if mi(bSwellingandorrednesssever) & nmissing==1
replace fHeadacheseverity="Does not interfere with physical activity" if mi(fHeadacheseverity) & nmissing==1
foreach v of varlist aPainseverity bSwellingandorrednesssever CLimitationofarmmovementse DFeverseverity EMuscleandjointpainseverit fHeadacheseverity {
	replace `v'=trim(`v')
}
drop if nmissing>=5
drop nmissing
merge m:1 PIN using "1 Raw data\alloc", keep(1 3) assert(2 3)
drop _merge
gen a_pain=1 if aPainseverity=="1 Mild pain"
replace a_pain=2 if aPainseverity=="2 Moderate pain"
gen b_swelling=0 if bSwellingandorrednesssever=="None"
replace b_swelling=1 if bSwellingandorrednesssever=="Less than 20mm"
replace b_swelling=2 if bSwellingandorrednesssever=="More than 20mm"
gen c_arm=0 if CLimitationofarmmovementse=="No limitation"
replace c_arm=1 if CLimitationofarmmovementse=="Minimum limitation"
replace c_arm=2 if CLimitationofarmmovementse=="Moderate limitation"
replace c_arm=3 if CLimitationofarmmovementse=="Severe limitation"
gen d_fever=0 if DFeverseverity=="Does not interfere with physical activity"
replace d_fever=1 if DFeverseverity=="Impairs physical activity"
replace d_fever=2 if DFeverseverity=="Impairs physical activity and requires medications"
gen e_muscle=0 if EMuscleandjointpainseverit=="Does not interfere with physical activity"
replace e_muscle=1 if EMuscleandjointpainseverit=="Impairs physical activity"
replace e_muscle=2 if EMuscleandjointpainseverit=="Impairs physical activity and requires medications"
gen f_headache=0 if fHeadacheseverity=="Does not interfere with physical activity"
replace f_headache=1 if fHeadacheseverity=="Impairs physical activity"
replace f_headache=2 if fHeadacheseverity=="Impairs physical activity and requires medications"
save vaccine_diary, replace

use vaccine_diary, clear
gen byte n=1
preserve
foreach v of varlist a_pain-f_headache {
	levelsof `v', local(lev)
	collapse (sum) n, by(`v' alloc Day)
	bys alloc Day: egen N=total(n)
	gen percent=100*n/N
	reshape wide n percent, i(alloc Day) j(`v')
	save daily_`v', replace
	restore, preserve
	collapse (max) `v', by(PIN alloc)
	gen byte n=1
	bys alloc: egen N=total(n)
	collapse (sum) n (max) N, by(`v' alloc)
	gen Day=8
	gen percent=100*n/N
	reshape wide n percent, i(alloc Day) j(`v')	
	append using daily_`v'
	foreach l of local lev {
		replace percent`l'=0 if mi(percent`l')
	}
	save daily_`v', replace
	restore, preserve
}

use daily_a_pain, clear
gen zero=0
gen cut=100 /* no patients reported severe pain */
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent1 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent1 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut percent1 l if alloc==1, barw(.4) pstyle(p1) fintens(inten50)) ///
	(rbar cut percent1 r if alloc==0, barw(.4) pstyle(p2) fintens(inten50)) ///
	(rbar tot cut l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "Severe pain" 5 "" 6 "" - "Moderate pain" 3 "" 4 "" - "No pain or mild pain" 1 "" 2 "" - " ") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2) bmargin(b+.7)) ///
	plotr(margin(b 0)) graphr(color(white)) title("A", pos(10) size(huge)) subtitle("Pain at injection site", size(medium) color(black)) name(A, replace)

use daily_b_swelling, clear
gen zero=0
gen cut=percent0+percent1
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent0 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent0 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten50)) ///
	(rbar cut percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten50)) ///
	(rbar tot cut l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "More than 20 mm" 5 "" 6 "" - "Less than 20 mm" 3 "" 4 "" - "None" 1 "" 2 "" - " ") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2) bmargin(b+.7)) ///
	plotr(margin(b 0)) graphr(color(white)) title("B", pos(10) size(huge)) subtitle("Swelling/redness at injection site", size(medium) color(black)) name(B, replace)

use daily_c_arm, clear
gen zero=0
gen cut1=percent0+percent1
gen cut2=cut1+percent2
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent0 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent0 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut1 percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten30)) ///
	(rbar cut1 percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten30)) ///
	(rbar cut2 percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten70)) ///
	(rbar cut2 percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten70)) ///
	(rbar tot cut2 l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut2 r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "Severe limitation" 7 "" 8 "" - "Moderate limitation" 5 "" 6 "" - "Minimum limitation" 3 "" 4 "" - "No limitation" 1 "" 2 "") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2) bmargin(b+.7)) ///
	plotr(margin(b 0)) graphr(color(white)) title("C", pos(10) size(huge)) subtitle("Limitation of arm movement", size(medium) color(black)) name(C, replace)


use daily_d_fever, clear
gen zero=0
gen cut=percent0+percent1
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent0 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent0 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten50)) ///
	(rbar cut percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten50)) ///
	(rbar tot cut l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "Impairs physical activity" "and requires medications" 5 "" 6 "" - "Impairs physical activity" 3 "" 4 "" - "Does not interfere with" "physical activity" 1 "" 2 "") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2)) ///
	plotr(margin(b 0)) graphr(color(white)) title("D", pos(10) size(huge)) subtitle("Fever", size(medium) color(black)) name(D, replace)

use daily_e_muscle, clear
gen zero=0
gen cut=percent0+percent1
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent0 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent0 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten50)) ///
	(rbar cut percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten50)) ///
	(rbar tot cut l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "Impairs physical activity" "and requires medications" 5 "" 6 "" - "Impairs physical activity" 3 "" 4 "" - "Does not interfere with" "physical activity" 1 "" 2 "") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2)) ///
	plotr(margin(b 0)) graphr(color(white)) title("E", pos(10) size(huge)) subtitle("Muscle and joint pain", size(medium) color(black)) name(E, replace)

use daily_f_headache, clear
gen zero=0
gen cut=percent0+percent1
gen tot=100
gen l=Day-.2
gen r=Day+.2
twoway (rbar percent0 zero l if alloc==1, barw(.4) fintens(inten10)) ///
	(rbar percent0 zero r if alloc==0, barw(.4) fintens(inten10)) ///
	(rbar cut percent0 l if alloc==1, barw(.4) pstyle(p1) fintens(inten50)) ///
	(rbar cut percent0 r if alloc==0, barw(.4) pstyle(p2) fintens(inten50)) ///
	(rbar tot cut l if alloc==1, barw(.4) pstyle(p1) fintens(inten100)) ///
	(rbar tot cut r if alloc==0, barw(.4) pstyle(p2) fintens(inten100)), ///
	xla(1(1)7 8 "Worst") xti("Day          ") ///
	yla(, angle(0)) yti(Percentage of patients) ///
	legend(order(- "Impairs physical activity" "and requires medications" 5 "" 6 "" - "Impairs physical activity" 3 "" 4 "" - "Does not interfere with" "physical activity" 1 "" 2 "") col(3) subti("PCV13      Placebo  ", pos(1) size(medsmall)) region(lc(none)) span symx(*.7) colgap(2)) ///
	plotr(margin(b 0)) graphr(color(white)) title("F", pos(10) size(huge)) subtitle("Headache", size(medium) color(black)) name(F, replace)

graph combine A B C D E F, col(2) imargin(tiny) iscale(*.75) graphr(color(white)) ysize(12) xsize(9)
gr save vaccine_diary, replace
gr export vaccine_diary.svg, replace
gr export vaccine_diary.eps, replace
gr export vaccine_diary.png, replace width(2800)
