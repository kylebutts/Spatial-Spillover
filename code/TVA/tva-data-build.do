********************************************************************************
* This file merges all the datasets used 
********************************************************************************

clear all
set maxvar 8000
global dir "/Users/kylebutts/Documents/Projects/Spatial Spillover/data/tva"


********************************************************************************
* TVA DUMMIES
********************************************************************************

infile fips using "${dir}/tvacounties.txt", clear
g tva = 1
sort fips
save tmp, replace


********************************************************************************
* READ STATE-LEVEL DATA
* AND TAKE 1930 INCOME VARIABLE
********************************************************************************

use "${dir}/state_level_data", clear
keep state pcp_income_30
sort state
save tmp7, replace


********************************************************************************
* WAGE DATA
* THESE ARE THE
* CORRECT MANUF AND TRADE (RETAIL+WHOLESALE) WAGES.
********************************************************************************

use "${dir}/tva1", clear
keep fips mwage* pcmwage* pctwage* var88_county72 var89_county72
duplicates report fips
duplicates drop fips, force
sort fips
save tmp76, replace

* mwage is the total county level manufacturing wages in thousands of dollars.
* For 1900, 1920, 1930, it is from the Census of Population
* and it is given in dollars, so divided by 1,000. For 1940 defined as 1939 wages, from the Census of Population. For 1950 defined as
* 1954 wages from Economic Census. For 1960 defined as 1963 wages from Economic Census. 
* For 1970 is defined as 1972 wages from Economic Census, given in millions
* of dollars, so multiplied by 1,000. For 1980 defined as 1982 wages from Economic Census, given in millions of dollars, so
* multiplied by 1,000. For 1990 defined as 1987 wages from Economic Census, given in millions of dollars, so multiplied
* by 1,000. For 2000, defined as 1997 wages from Economic Census.

* pcmwage is equal to mwage divided by the relevant number of workers in thousand
* For 1900, 1920, 1930, and 1940 
* mage is divided by number of manufacturing workers in the county,  from the Census of Population. 
* For 1950 mwage is divided by 
* the 1954 number of manufacturing workers, from the Economic Census. For 1960 mwage is divided by the 1963 number of manufacturing workers 
* from Economic Census. 
* For 1970 mwage is divided by the 1972 number of manufacturing workers, from Economic Census.
* For 1980 mwage is divided by 1982 number of manufacturing workers, from Economic Census.
* For 1990 mwage is divided by the 1987 number of manufacturing workers, from Economic Census
* For 2000, mwage is divided by the 1997 number of manufacturing workers, from Economic Census.


********************************************************************************
* COUNTY-LEVEL DATA
********************************************************************************
//this file is built by tva_update.do (certified by Valentina 7/14/2011)

use "${dir}/county_level_data", clear
duplicates report fips
drop mwage* twage*
drop _merge
sort fips

merge 1:1 fips using tmp76
tab _merge
keep if _merge ==3
drop _merge

merge 1:1 fips using tmp
tab _merge
drop if _merge ==2
replace tva = 0 if tva ==.


********************************************************************************
* Merge on Fishback Data 
********************************************************************************

drop _merge
drop if county==.|state==.

merge 1:1 county state using "${dir}/fishback"
tab _merge
drop if _merge==2
drop N10


********************************************************************************
* Merge on Ag Land Values
********************************************************************************

drop _merge
merge 1:1 fips using "${dir}/ag_data"
tab _merge
drop if _merge==2
duplicates report fips


********************************************************************************
* Merge on TVA 'Donut'   
********************************************************************************

drop _merge
joinby fips using "${dir}/tva_donut", unmatched(both)
tab _merge
drop if _merge==2


********************************************************************************
* Merge on Correct Employment Variables
********************************************************************************

drop _merge
merge 1:1 fips using "${dir}/enrico_jobs"

tab _merge
drop if _merge==2
duplicates report fips


********************************************************************************
* Merge on Housing Val/Rents 
********************************************************************************

drop _merge
merge 1:1 fips using "${dir}/housingvals"
tab _merge
drop if _merge==2
duplicates report fips

drop _merge
merge 1:1 fips using "${dir}/county62_1", keepusing(var61_county62 var63_county62)

ren var61_county62 medhsval60
ren var63_county62 medrnt60
tab _merge
drop if _merge==2
duplicates report fips



********************************************************************************
* Merge on weather variables 
********************************************************************************

drop _merge

merge 1:1 fips using "${dir}/july_weather_1968_2002.dta"
foreach var of varlist tmin tmean tmax{
        ren `var' `var'_jul
	}
drop if _merge==2
drop _merge

merge 1:1 fips using "${dir}/jan_weather_1968_2002.dta"
foreach var of varlist tmin tmean tmax{
		ren `var' `var'_jan
	}
drop if _merge==2


********************************************************************************
* Merge extra manuf vars  
********************************************************************************

drop _merge
merge 1:1 fips using "${dir}/maunf_data_1.dta"
tab _merge
drop if _merge==2

drop _merge
merge 1:1 fips using "${dir}/manuf_data_2.dta"
tab _merge
drop if _merge==2

drop _merge
merge 1:1 fips using "${dir}/manuf_data_3.dta"
tab _merge
drop if _merge==2
summ *merge*
drop _merge*


********************************************************************************
* PROPOSED AUTHORITIES
********************************************************************************

sort fips
merge fips using "${dir}/proposed_authorities.dta"
tab _merge
drop if _merge==2
drop _merge
tab aut_euclidean1
tab aut_euclidean2
tab aut_euclidean3
tab aut_cent
tab aut_euclidean5
tab aut_euclidean6
g       aut4 =0
replace aut4 =1  if max(aut_euclidean1,aut_euclidean2,aut_euclidean3,aut_cent,aut_euclidean5,aut_euclidean6) ==1
summ fips tva aut_*
tab tva aut4
drop aut_eucl* aut_cent*



********************************************************************************
* Cleanup and Definition 
********************************************************************************

* State and region
drop state
g state = int(fips/1000)
drop if state ==51 | state ==52 | state ==2 | state ==3 | state ==15
g northeast =0
g midwest=0
g south=0
g west=0

replace northeast =1 if state == 9 | state == 23 | state == 25 | state == 33 | state == 44 | state  == 50 | state == 34 | state == 36 | state == 42
replace midwest=1    if state == 17 | state == 18 | state == 26 | state == 39 | state == 55 | state ==19 | state == 20 | state == 27 | state == 29 | state == 31 | state == 38 | state == 46
replace south=1      if state ==10 | state == 11 | state == 12 | state == 13 | state == 24 | state  == 37 | state == 45 | state == 51   | state == 54 | state == 1 | state == 21 | state == 28  | state == 47 | state == 5 | state == 22 | state == 40 | state == 48
replace west=1       if state == 4 | state == 8 | state == 16 | state == 30 | state == 32 | state == 35 | state == 49 | state == 56 | state == 2 | state == 6 | state == 15 | state == 41 | state ==53
g       region =1 if northeast ==1
replace region =2 if midwest ==1
replace region =3 if south ==1
replace region =4 if west ==1


* CPI
g cpi890= 5 
g cpi0  = 7
g cpi10 = 9
g cpi20 = 20
g cpi30 = 16.7
g cpi40 = 14
g cpi50 = 24.1
g cpi60 = 29.6
g cpi70 = 38.8
g cpi80 = 82.4
g cpi90 = 130.7
g cpi2000 = 172.2


* Merge the state-level variable
sort state
merge state using tmp7
tab _merge
drop if tva ==.
rm tmp.dta
rm tmp7.dta
rm tmp76.dta


********************************************************************************
* Drop Donut and counties with missing lat/long
********************************************************************************

* drop if border_county==1|latitude==.|longitud==.


********************************************************************************
* Drop counties with very low populations in any year
********************************************************************************

drop if pop0<1000|pop10<1000|pop20<1000|pop30<1000|pop40<1000|pop50<1000|pop60<1000|pop70<1000|pop80<1000|pop90<1000|pop2000<1000



********************************************************************************
* from evaluation.do, clean data
********************************************************************************

* Per capita wage in manufacturing, in real dollars
* Manufacturing payroll of production workers over average number of production workers 
g wage0    = (pcmwage00)/(cpi0/100)
g wage20   = (pcmwage20)/(cpi20/100)
g wage30   = (pcmwage30)/(cpi30/100)
g wage40   = (pcmwage39)/(cpi40/100)
g wage50   = (pcmwage47)/(cpi50/100)
g wage60   = (pcmwage58)/(cpi60/100)
g wage70   = (pcmwage67)/(cpi70/100)
g wage80   = (pcmwage82)/(cpi80/100)
g wage90   = (pcmwage87)/(cpi90/100)
g wage2000 = (pcmwage97)/(cpi2000/100)


* Per capita wage in trade, in real dollars
* total payroll in wholesale establishments + retail establishments / (total workers in retail+ in wholsale)
g twage30  = pctwage30/(cpi30/100) 
g twage40  = pctwage40/(cpi40/100)
g twage50  = pctwage54/(cpi50/100)
g twage60  = pctwage63/(cpi60/100)
g twage70  = pctwage72/(cpi70/100)
g twage80  = pctwage82/(cpi80/100)
g twage90  = pctwage87/(cpi90/100)
g twage2000  = pctwage97/(cpi2000/100)

* Agricultural values

g lnfaval0 = ln(faval900/(cpi0/100))
g lnfaval10 = ln(faval910/(cpi10/100))
g lnfaval20 = ln(faval920/(cpi20/100))
g lnfaval30 = ln(faval930/(cpi30/100))
g lnfaval40 = ln(faval940/(cpi40/100))
g lnfaval50 = ln(faval950/(cpi50/100))
g lnfaval60 = ln(faval959/(cpi60/100))
g lnfaval70 = ln(faval1969/(cpi70/100))
g lnfaval80 = ln(faval1982/(cpi80/100))
g lnfaval90 = ln(faval1992/(cpi90/100))
g lnfaval2000=ln(faval2002/(cpi2000/100))


* Median family income
gen lnmedfaminc50  = ln(medfaminc50/(cpi50/100)) 
gen lnmedfaminc60  = ln(medfaminc60/(cpi60/100)) 
gen lnmedfaminc70  = ln(medfaminc70/(cpi70/100)) 
gen lnmedfaminc80  = ln(medfaminc80/(cpi80/100)) 
gen lnmedfaminc90  = ln(medfaminc80/(cpi90/100)) 
gen lnmedfaminc2000= ln(medfaminc2000/(cpi2000/100))

* Farm production
gen lnvfprod30   = log(vfprod30/(cpi30/100))
gen lnvfprod40   = log(vfprod40/(cpi40/100))
gen lnvfprod50   = log(vfprod50/(cpi50/100))
gen lnvfprod60   = log(vfprod60/(cpi60/100))
gen lnvfprod70   = log(vfprod70/(cpi70/100))
gen lnvfprod80   = log(vfprod80/(cpi80/100))
gen lnvfprod90   = log(vfprod90/(cpi90/100))
gen lnvfprod2000 = log(vfprod2000/(cpi2000/100))

*Foreign born
gen fb0=fbwmtot00 + fbwftot00 
gen fb10=fbwtot10
gen fb20=fbwmtot20 + fbwftot20 
gen fb30=fbwmtot + fbwftot

gen fbshr20=fb20/(wmtot20 + wftot20)
gen fbshr30=fb30/(wmtot + wftot)

*Housing Values/Rents

ren medrnt30_NHGIS medrnt30
ren medhsval30_NHGIS medhsval30
ren var88_county72 medhsval70
ren var89_county72 medrnt70

foreach var in medhsval medrnt{
	foreach yr in 30 40 50 60 70 80 90 2000{
		cap replace `var'`yr'=`var'`yr'/(cpi`yr'/100)
	}
}

*various
drop other60

*drop counties experiencing big changes in area
gen d=(b1_lnd01_county00 - area)/(b1_lnd01_county00 + area)/2
drop if abs(d)>.03
replace area=(b1_lnd01_county00 + area)/2



********************************************************************************
* standardize variable names
********************************************************************************


ren emp00 emp0
ren manuf_jobs_00 manuf_jobs_0

foreach yr in 0 10 20 30 40 50 60 70 80 90 2000{
	cap drop manuf`yr'
	cap drop agr`yr'
	ren manuf_jobs_`yr' manuf`yr'
	cap ren ag_jobs`yr' agr`yr'
}

foreach yr in 0 10 20 30 60 80 90 2000{
	gen other`yr'=emp`yr'-agr`yr'-manuf`yr'
}


********************************************************************************
*   Make Share    
********************************************************************************

foreach var in manuf agr{
	foreach yr in 0 10 20 30 40 50 60 70 80 90 2000{
		cap gen `var'shr`yr'=`var'`yr'/emp`yr'
	}
}

********************************************************************************
* prepare outcomes
********************************************************************************


foreach var in pop emp house wage twage agr manuf other medhsval medrnt fb{
	cap gen ln`var'0=ln(`var'0)
	cap gen ln`var'10=ln(`var'10)
	cap gen ln`var'20=ln(`var'20)
	cap gen ln`var'30=ln(`var'30)
	cap gen ln`var'40=ln(`var'40)
	cap gen ln`var'50=ln(`var'50)
	cap gen ln`var'60=ln(`var'60)
	cap gen ln`var'70=ln(`var'70)
	cap gen ln`var'80=ln(`var'80)
	cap gen ln`var'90=ln(`var'90)
	cap gen ln`var'2000=ln(`var'2000)
}

drop lntwage20  lnmedhsval20 lnmedrnt20 //stata confuses 20 and 2000


********************************************************************************
* Generate Controls
********************************************************************************
g urbshare0=popurb0/pop0
g urbshare10=popurb10/pop10
g urbshare20=popurb20/pop20
g urbshare30=popurb30/pop20
gen popdens0=pop0/b1_lnd01_county00


*fix zeros in covariate quantities

foreach var in agr manuf other{
	foreach yr in 10 20 30{
		cap replace ln`var'`yr'=0 if `var'`yr'<=0
		gen no`var'`yr'dum=`var'`yr'<=0
	}
}

ren mfgcap00 mfgcap0
gen lnmfgcap0=ln(mfgcap0)
replace lnmfgcap0=0 if mfgcap0==0



*fix wages

foreach var in wage twage{
	foreach yr in 20 30{
		cap replace ln`var'`yr'=-1 if `var'`yr'==.
		cap gen no`var'`yr'dum=`var'`yr'==.
	}
}



*transformations
foreach var of varlist elevmax elevrang popdens0 tillit10 tillit1020 tillit1010 retsales radiorep totunemp area{
	gen ln`var'=ln(`var')
}


center tmean* lnelevmax

foreach var of varlist c_lnelevmax white0 white20 white30 c_tmean* lnmanuf0 lnmanuf10 lnmanuf20 lnradiorep lnemp20 lnemp30 lnwage0 lnwage20 lnwage30 lntwage30 lnpop0 lnpop20 lnpop30 lnfaval0 lnfaval20 lnfaval30 tmean*{
	gen `var'sq=`var'^2
	gen `var'cub=`var'^3
}

gen popdifsq=(lnpop30-lnpop20)^2
gen empdifsq=(lnemp30-lnemp20)^2
gen manufdifsq=(lnmanuf20-lnmanuf0)^2
gen urbsharedifsq=(urbshare20-urbshare0)^2
gen whitedifsq=(white20-white0)^2
gen agrdifsq=(lnagr20-lnagr0)^2
gen wagedifsq=(lnwage20-lnwage0)^2
gen favaldifsq=(lnfaval20-lnfaval0)^2

gen agrshr30sq=agrshr30^2
gen agrshr20sq=agrshr20^2

gen lnagr20sq=lnagr20^2
gen lnagr0sq=lnagr0^2
gen fbdifsq=(lnfb20-lnfb0)^2

gen wagedif2sq=(lnwage30-lnwage20)^2
gen tmean_jan_jul=tmean_jan*tmean_jul

gen pctil20=tillit1020/t10tot20
gen pctil30=tillit10/t10tot
replace PRADIO=PRADIO/100
gen urate30=totunemp/(totunemp+emp30)



save "${dir}/build.dta", replace





