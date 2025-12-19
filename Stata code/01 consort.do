* 01 consort.do
* Author: David Harrison
* Date: 01/03/2023
* 
version 18
set more off
clear
cap log close
* base working directory stored in global
cd "$base"
log using "consort.log", replace
use working, clear
* allocation
tab alloc
* received intervention
tab IMP_01 alloc, col
* reasons for not receiving intervention
sort alloc Label
li Site Label RAN_02 alloc IMP_04 IMP_05 WD_02 WD_03 if IMP_01==0, noo sepby(alloc)
* loss to follow-up (withdrawals - excluding due to death)
tab WD_02 alloc if WD_02!=3
* only excluded from analysis if T=0 (withdrawn on day of randomisation)
li Site Label RAN_02 alloc WD_02 WD_03 T_primary if inlist(WD_02,1,5), noo sepby(alloc)
log close
