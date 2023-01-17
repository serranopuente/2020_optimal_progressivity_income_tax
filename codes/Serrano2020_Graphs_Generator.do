/*
Optimal Progressivity of Personal Income Tax
A General Equilibrium Evaluation for Spain
Output file for graphs
Darío Serrano Puente (2020)

Last update: September 20th, 2020
*/


********************************************************************************
// 0. Preamble
********************************************************************************
clear matrix
clear all
set more off

// Define input/output directories
global basePath "C:\Darío\Research\Optimal Progressivity Spain"
global figurePath "$basePath\Codes\Output\Figures"
global tablePath "$basePath\Codes\Output\Tables"
global dataPath "$basePath\Data\AEAT\Output\Households 2015"

// ssc install grstyle // If not installed
// ssc install vlookup // If not installed
grstyle init
grstyle set graphsize 14cm 24cm
grstyle set legend 2, inside nobox


********************************************************************************
// 1. Data with every reform evaluation
********************************************************************************
import excel "$tablePath\OutputRef_actual.xlsx", sheet("Hoja1") clear

// 1.1 Transpose the imported data to have it in Stata structure
xpose, clear

// 1.2 Rename variables
rename v1-v39 (tau r lambda omega w Y K L C H_elle K_L L_H Y_H C_H K_Y C_Y I_Y G_Y Tr_Y T_Y ///
ratio_old_young corr_fat_sons gini_income p40_income p40p60_income p60p80_income p80p100_income ///
 p90p95_income p95p99_income p99p100_income gini_wealth p40_wealth p40p60_wealth p60p80_wealth ///
 p80p100_wealth p90p95_wealth p95p99_wealth p99p100_wealth Agg_welf_change)

********************************************************************************
// 1. Aggregate welfare changes graph
********************************************************************************
// Find optimal progressivity
sort Agg_welf_change
scalar opt_tau = tau[_N]
sort tau

twoway line Agg_welf_change tau, lcolor("218 108 122") msize(small) lwidth(medthick) ///
legend(on label(1 "Aggregate Welfare Change") region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("% change in consumption",margin(small) size(medlarge)) xline(0.1146, lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lpattern("---") lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
xlabel(0 .05 .10 .15 .20 .25 .30 .35 .40 .45 .50, format(%3.2f) labsize(medlarge)) ylabel(,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black) lwidth(thin)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(-6 `=scalar(0.1146+0.032)' "{&tau}=0.1146", size(medlarge)) text(-6 `=scalar(opt_tau+0.024)' "{&tau}=0.23", size(medlarge)) ///
text(-8 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(-8 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Average_Welfare_Gain.pdf", replace

********************************************************************************
// 2. Macro Aggregates change
********************************************************************************
// Create Indexes
sort Agg_welf_change
qui gen Y_index = Y/Y[16]*100
qui gen K_index = K/K[16]*100
qui gen L_index = L/L[16]*100
qui gen C_index = C/C[16]*100
sort tau

twoway (line Y_index tau, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern("l")) ///
(line K_index tau, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(shortdash_dot)) ///
(line L_index tau, lcolor("255 190 93") msize(small) lwidth(medthick) lpattern(dash)) ///
(line C_index tau, lcolor("26 140 255") msize(small) lwidth(medthick) lpattern(longdash_dot)), ///
legend(on order(1 "Output" 2 "Capital" 3 "Labor" 4 "Consumption") region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Index respect to actual scenario",margin(small) size(medlarge)) xline(0.1146, lwidth(vthin) lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lwidth(vthin) lpattern("---") lcolor(black)) yline(100, lpattern(dot) lwidth(thin) lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) ///
xlabel(0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50, format(%10.2f) labsize(medlarge))  ylabel(50(25)125,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black) lwidth(thin)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(50 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(50 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Macro_Aggregates.pdf", replace

********************************************************************************
// 3. Fiscal Policy Ratios
********************************************************************************
// Create Indexes
sort Agg_welf_change
qui gen Tr_Y_index = Tr_Y/Tr_Y[16]*100
qui gen G_Y_index = G_Y/G_Y[16]*100
qui gen T_Y_index = T_Y/T_Y[16]*100
qui gen T = T_Y*Y
qui gen Tr = Tr_Y*Y
qui gen G = G_Y*Y
qui gen Tr_index = Tr/Tr[16]*100
qui gen G_index = G/G[16]*100
qui gen T_index = T/T[16]*100
sort tau

grstyle set legend 2, inside nobox
twoway (line Tr_Y_index tau, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern("l")) ///
(line G_Y_index tau, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(shortdash_dot)) ///
(line T_Y_index tau, lcolor("255 190 93") msize(small) lwidth(medthick) lpattern(dash)), ///
legend(on order(1 "Transfers-to-Output" 2 "Gov. Expenditure-to-Output" 3 "Tax Revenue-to-Output") region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Index respect to actual scenario",margin(small)size(medlarge)) xline(0.1146, lwidth(vthin) lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lwidth(vthin) lpattern("---") lcolor(black)) yline(100, lpattern(dot) lwidth(thin) lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) ///
xlabel(0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50, format(%10.2f) labsize(medlarge)) ylabel(,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(50 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(50 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Fiscal_Ratios.pdf", replace

twoway (line Tr_index tau, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern("l")) ///
(line G_index tau, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(shortdash_dot)) ///
(line T_index tau, lcolor("255 190 93") msize(small)), ///
legend(on order(1 "Transfers" 2 "Gov. Expenditure" 3 "Tax Revenue") region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Index respect to actual scenario",margin(small)size(medlarge)) xline(0.1146, lwidth(vthin) lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lwidth(vthin) lpattern("---") lcolor(black)) yline(100, lpattern(dot) lwidth(thin) lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) ///
xlabel(0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50, format(%10.2f) labsize(medlarge)) ylabel(,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(70 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(70 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Fiscal_Aggregates.pdf", replace

********************************************************************************
// 4. Inequality on progressivity
********************************************************************************
grstyle set legend 2, inside nobox

// Create Indexes
sort Agg_welf_change
qui gen top1_wealth_index = p99p100_wealth/p99p100_wealth[16]*100
qui gen p95_99_wealth_index = p95p99_wealth/p95p99_wealth[16]*100
qui gen p90_95_wealth_index = p90p95_wealth/p90p95_wealth[16]*100
qui gen gini_wealth_index = gini_wealth/gini_wealth[16]*100
qui gen top1_income_index = p99p100_income/p99p100_income[16]*100
qui gen p95_99_income_index = p95p99_income/p95p99_income[16]*100
qui gen p90_95_income_index = p90p95_income/p90p95_income[16]*100
qui gen gini_income_index = gini_income/gini_income[16]*100
sort tau

twoway (line top1_wealth_index tau, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern(longdash_dot)) ///
(line p95_99_wealth_index tau, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(shortdash_dot)) ///
(line p90_95_wealth_index tau, lcolor("255 190 93") msize(small) lwidth(medthick) lpattern(dash)) ///
(line gini_wealth_index tau, lcolor("26 140 255") msize(small) lwidth(medthick) lpattern("l")), ///
legend(on order(1 "Top 1%" 2 "p95-p99" 3 "p90-p95" 4 "Gini Index") rows(2) region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Index respect to actual scenario",margin(small) size(medlarge)) xline(0.1146, lwidth(vthin) lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lwidth(vthin) lpattern("---") lcolor(black)) yline(100, lpattern(dot) lwidth(thin) lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) ///
xlabel(0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50, format(%10.2f) labsize(medlarge))  ylabel(25(25)125,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(30 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(30 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Ineq_wealth.pdf", replace

twoway (line top1_income_index tau, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern(longdash_dot)) ///
(line p95_99_income_index tau, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(shortdash_dot)) ///
(line p90_95_income_index tau, lcolor("255 190 93") msize(small) lwidth(medthick) lpattern(dash)) ///
(line gini_income_index tau, lcolor("26 140 255") msize(small) lwidth(medthick) lpattern("l")), ///
legend(on order(1 "Top 1%" 2 "p95-p99" 3 "p90-p95" 4 "Gini Index") rows(2) region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Index respect to actual scenario",margin(small) size(medlarge)) xline(0.1146, lwidth(vthin) lpattern("---") lcolor(black)) ///
xline(`=scalar(opt_tau)', lwidth(vthin) lpattern("---") lcolor(black)) yline(100, lpattern(dot) lwidth(thin) lcolor(black)) ///
xtitle("Progressivity level ({&tau})",margin(small) size(medlarge)) ///
xlabel(0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50, format(%10.2f) labsize(medlarge))  ylabel(25(25)125,labsize(medlarge) nogrid) ///
plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
text(30 `=scalar(0.1146+0.025)' "Actual", size(medlarge)) text(30 `=scalar(opt_tau+0.03)' "Optimal", size(medlarge))
graph export "$figurePath\Ineq_income.pdf", replace


********************************************************************************
// 5. Change in Average Tax Rates (model)
********************************************************************************
grstyle set legend 6, inside nobox

// Generate data, i.e. means of multiple income
set obs 200
gen x = .
local j = 0
forval i=.5(.1)16 {
	local j = `j'+1
	qui replace x = `i' in `j'
}
drop if x==.
// Baseline economy
gen be = 1 - .8924*(x)^(-.1146)
// Optimal economy
gen opt = 1 - 1.1832608*(x)^(-.23)

twoway (line be x, lcolor("218 108 122") msize(small) lwidth(medthick) lpattern("l")) ///
(line opt x, lcolor("0 64 129") msize(small) lwidth(medthick) lpattern(dash)), ///
legend(on order(1 "Actual ({&tau}=0.1146 & {&lambda}=0.8924)" 2 "Optimal ({&tau}=0.23 & {&lambda}=1.1833)") rows(2) region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
ytitle("Average tax rates",margin(small) size(medlarge)) xtitle("Multiples of mean household gross income",margin(small) size(medlarge)) ///
xlabel(1(1)16, nogrid format(%10.0f) labsize(medlarge)) ylabel(-0.50(.25).50, nogrid format(%10.2f)labsize(medlarge)) ///
yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white)
graph export "$figurePath\Avg_Tax_Rates_Optimal_Reform.pdf", replace

********************************************************************************
// 6. Average Tax Rates by Percentile of Gross Income
********************************************************************************
grstyle init
grstyle set graphsize 14cm 24cm
grstyle set legend 12, inside nobox

// Create intermediate dataset with change of avg. effective tax rates by percentiles
import excel "$tablePath\EffTaxRateChange.xlsx", sheet("Hoja1") clear
// Transpose the imported data to have it in Stata structure and rename variables
xpose, clear
rename (v1 v2 v5) (lb ub eff_avg_rate_change)
keep lb ub eff_avg_rate_change
save "$dataPath\Intermediate", replace
// Generate percentile group to later match microdata
gen match1_g = .
replace match1_g=1 if lb>=1 & ub<=10
replace match1_g=2 if lb>=12 & ub<=20
replace match1_g=3 if lb>=21 & ub<=30
replace match1_g=4 if lb>=31 & ub<=40
replace match1_g=5 if lb>=41 & ub<=50
replace match1_g=6 if lb>=51 & ub<=60
replace match1_g=7 if lb>=61 & ub<=70
replace match1_g=8 if lb>=71 & ub<=80
replace match1_g=9 if lb>=81 & ub<=90
replace match1_g=10 if lb>=90 & ub<=95
replace match1_g=11 if lb>=95 & ub<=99
replace match1_g=12 if lb>=99
drop if match1_g == .
drop lb ub
rename eff_avg_rate_change eff_avg_rate_change_1
save "$dataPath\Intermediate1_Change", replace
u "$dataPath\Intermediate", clear
erase "$dataPath\Intermediate.dta"
gen match2_g = .
replace match2_g=1 if lb==1 & ub==20
replace match2_g=2 if lb==21 & ub==40
replace match2_g=3 if lb==41 & ub==60
replace match2_g=4 if lb==61 & ub==80
replace match2_g=5 if lb==81 & ub>=100
drop if match2_g == .
drop lb ub
rename eff_avg_rate_change eff_avg_rate_change_2
save "$dataPath\Intermediate2_Change", replace

// Import microdata statistics
u "$dataPath\Income_Distribution_Statistics_2015_Households.dta", clear

// Keep Variables
keep Factor EDAD1 tipohogar laborg capitalg selfg grossinc avgtaxg 

// Define Restrictions
// restrictions on income
global restrict_g "grossinc>0 & laborg>=0 & capitalg>=0 & selfg>=0 & avgtaxg <= 0.56" 
// restrictions to create sub-samples
global nothing "_n>=1"
global restrict_sample $nothing

// Compute Mean Income
// Gross Income
sum grossinc [aw=Factor] if $restrict_g & $restrict_sample 
gen double meaninc = r(mean)
la var meaninc "Mean Income"
save "Intermediate.dta", replace
keep meaninc
rename meaninc meaninc_H2015
keep if _n == 1
save "$dataPath\meaninc_H2015", replace
use "Intermediate.dta", clear
erase "Intermediate.dta"

// Compute Normalized Income
gen double nincome = grossinc/meaninc if $restrict_g & $restrict_sample 
la var nincome "Normalized Income"

// Value Labels for Percentiles
// Quintiles
label define xaxis_lbl 20 "bottom 20%" 40 "20-40%" 60 "40-60%" 80 "60-80%" 100 "80-100%"
// Top and bottom 10%, 5% and 1%.
label define xaxis2_lbl 1 "bottom 1%" 5 "1-5%" 10 "5-10%" 90 "90-95%" 95 "95-99%" 99 "top 1%"

// Percentiles of Gross Income
sort grossinc
xtile pct_g = grossinc [pw=Factor] if $restrict_g & $restrict_sample , nq(100)
la var pct_g "Percentile of Gross Income"
tab pct_g [aw=Factor]

// Quintiles of Gross Incom
gen xaxis_g = .
la var xaxis_g "Quintiles of Gross Income"
replace xaxis_g=100 if pct_g<=100 
replace xaxis_g=80 if pct_g<=80 
replace xaxis_g=60 if pct_g<=60  
replace xaxis_g=40 if pct_g<=40 
replace xaxis_g=20 if pct_g<=20
label values xaxis_g xaxis_lbl
tab xaxis_g [aw=Factor]

// Top and Bottom 10%,5% and 1%
gen xaxis2_g=.
la var xaxis2_g "Gross Income Top and Bottom 10%,5% and 1%"
replace xaxis2_g=99 if pct_g<=100 & pct_g>99   
replace xaxis2_g=95 if pct_g<=99 & pct_g>95    
replace xaxis2_g=90 if pct_g<=95 & pct_g>90 
replace xaxis2_g=10 if pct_g<=10 & pct_g>5   
replace xaxis2_g=5 if pct_g<=5 & pct_g>1   
replace xaxis2_g=1 if pct_g<=1 
label values xaxis2_g xaxis2_lbl
tab xaxis2_g [aw=Factor]

// Generate percentile group to later match microdata
gen match1_g = .
replace match1_g=12 if pct_g<=100 & pct_g>99
replace match1_g=11 if pct_g<=99 & pct_g>95    
replace match1_g=10 if pct_g<=95 & pct_g>90
replace match1_g=9 if pct_g<=90
replace match1_g=8 if pct_g<=80
replace match1_g=7 if pct_g<=70
replace match1_g=6 if pct_g<=60
replace match1_g=5 if pct_g<=50
replace match1_g=4 if pct_g<=40
replace match1_g=3 if pct_g<=30
replace match1_g=2 if pct_g<=20
replace match1_g=1 if pct_g<=10
gen match2_g = .
replace match2_g=5 if pct_g<=100
replace match2_g=4 if pct_g<=80
replace match2_g=3 if pct_g<=60
replace match2_g=2 if pct_g<=40
replace match2_g=1 if pct_g<=20

// Macth with eff. avg. tac rate changes from actual to optimal scenario
merge m:1 match1_g using "$dataPath\Intermediate1_Change"
drop _merge
erase "$dataPath\Intermediate1_Change.dta"

merge m:1 match2_g using "$dataPath\Intermediate2_Change"
drop _merge
erase "$dataPath\Intermediate2_Change.dta"

// Generate Data for Figures and Figures
gen double avgtaxg_opt=.
replace avgtaxg_opt=avgtaxg*(1+eff_avg_rate_change_1)
replace avgtaxg_opt=. if xaxis_g == .
// replace avgtaxg_opt = 0 if avgtaxg_opt < 0
la var avgtaxg_opt "Average tax rate (wrt gross income) in the optimal scenario"

// Figure: Average Tax Rate by percentile
tabstat avgtaxg [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis_g) stat(mean) nototal format(%9.3f) save
matrix q = (r(Stat1)\r(Stat2)\r(Stat3)\r(Stat4)\r(Stat5))
tabstat avgtaxg [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis2_g) stat(mean) nototal format(%9.3f) save
matrix t = (r(Stat4)\r(Stat5)\r(Stat6))
matrix DATA_FIGURE = (q\t)
matrix rownames DATA_FIGURE = q1st q2nd q3rd q4th q5th t90-95 t95-99 t1

di "Average Tax Rates by Percentile"
matrix list DATA_FIGURE, format(%9.4f)

global relabel 1 "0-20" 2 "20-40" 3 "40-60" 4 "60-80" 5 "80-100" 6 "90-95" 7 "95-99" 8 "99-100"
svmat DATA_FIGURE, names(data_a)
// Figure: Average Tax Rate by percentile in the optimal scenario
tabstat avgtaxg_opt [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis_g) stat(mean) nototal format(%9.3f) save
matrix q = (r(Stat1)\r(Stat2)\r(Stat3)\r(Stat4)\r(Stat5))
tabstat avgtaxg_opt [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis2_g) stat(mean) nototal format(%9.3f) save
matrix t = (r(Stat4)\r(Stat5)\r(Stat6))
matrix DATA_FIGURE = (q\t)
matrix rownames DATA_FIGURE = q1st q2nd q3rd q4th q5th t90-95 t95-99 t1

di "Average Tax Rates by Percentile in the optimal scenario"
matrix list DATA_FIGURE, format(%9.4f)

global relabel 1 "0-20" 2 "20-40" 3 "40-60" 4 "60-80" 5 "80-100" 6 "90-95" 7 "95-99" 8 "99-100"
svmat DATA_FIGURE, names(data_b)
qui replace data_b1 = 0.00001 if data_b1 < 0

gen byte aux = _n if data_a1!=.
graph bar (asis) data_a1, over(aux, relabel($relabel )label(labsize(medlarge))) plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
bar(1,fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin)) blabel(total, format(%10.3f) size(medlarge)) ///
b1title("Percentiles of household gross income", margin(small)size(medlarge)) ///
yti("Average tax rates", margin(small) size(medlarge)) ylabel(0(.05).35,nogrid format(%10.2f) labsize(medlarge)) 
graph export "$figurePath\Avg_Tax_Rates_Income_Percentiles.pdf", as(pdf) replace

qui gen aux_a = aux - 0.2
qui gen aux_b = aux + 0.2
gen var_1= data_a1 + 0.03
gen var_2= data_b1 + 0.03
format data* %10.3f

twoway (bar data_a1 aux_a, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.37)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) ///
yti("Average tax rates", margin(small)size(medlarge)) ylabel(0(.1).40,nogrid format(%10.2f) labsize(medlarge)) ///
xti("Percentiles of household gross income", margin(small)size(medlarge)) xlabel($relabel , nogrid format(%10.1f) labsize(medlarge))) ///
(scatter var_1 aux_a, msym(none) mlab(data_a1) mlabpos(12) mlabcolor(black) mlabsize(medium) mlabangle(ver)) ///
(bar data_b1 aux_b, barwidth(.37)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(scatter var_2 aux_b, msym(none) mlab(data_b1) mlabpos(12) mlabcolor(black) mlabsize(medium) mlabangle(ver)), ///
legend(on order(1 "Actual" 3 "Optimal") rows(1) region(lstyle(none)) size(medlarge) region(lwidth(none))) 
graph export "$figurePath\Avg_Tax_Rates_Income_Percentiles_Optimal.pdf", as(pdf) replace
drop data* aux* var*
mat drop _all

*******************************************************************************
// 7. Who pays the income taxes?
********************************************************************************
replace avgtaxg=. if xaxis_g == .
gen double avg_tax = avgtaxg*grossinc
gen double avg_tax_opt = avgtaxg_opt*grossinc
// Figure: share of the total tax revenue
tabstat avg_tax [aw=Factor] if $restrict_g & $restrict_sample, stat(sum) format(%9.3f) save
mat tot = r(StatTotal)
tabstat avg_tax [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis_g) stat(sum) format(%9.3f) save
mat q = (r(Stat1)\r(Stat2)\r(Stat3)\r(Stat4)\r(Stat5))
tabstat avg_tax [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis2_g) stat(sum) format(%9.3f) save 
mat t = (r(Stat4)\r(Stat5)\r(Stat6))
matrix DATA_FIGURE = (q\t)
matrix rownames DATA_FIGURE = q1st q2nd q3rd q4th q5th t90-95 t95-99 t1
forval i=1/`=rowsof(DATA_FIGURE)' {
	mat DATA_FIGURE[`i',1] = DATA_FIGURE[`i',1] / tot[1,1]
}

di "Share of tax revenue"
matrix list DATA_FIGURE, format(%9.4f)
global relabel 1 "0-20" 2 "20-40" 3 "40-60" 4 "60-80" 5 "80-100" 6 "90-95" 7 "95-99" 8 "99-100"
svmat DATA_FIGURE, names(data_a)
// Figure: share of the total tax revenue in the optimal scenario
tabstat avg_tax_opt [aw=Factor] if $restrict_g & $restrict_sample, stat(sum) format(%9.3f) save
mat tot = r(StatTotal)
tabstat avg_tax_opt [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis_g) stat(sum) format(%9.3f) save
mat q = (r(Stat1)\r(Stat2)\r(Stat3)\r(Stat4)\r(Stat5))
tabstat avg_tax_opt [aw=Factor] if $restrict_g & $restrict_sample , by(xaxis2_g) stat(sum) format(%9.3f) save 
mat t = (r(Stat4)\r(Stat5)\r(Stat6))
matrix DATA_FIGURE = (q\t)
matrix rownames DATA_FIGURE = q1st q2nd q3rd q4th q5th t90-95 t95-99 t1
forval i=1/`=rowsof(DATA_FIGURE)' {
	mat DATA_FIGURE[`i',1] = DATA_FIGURE[`i',1] / tot[1,1]
}

di "Share of tax revenue in the optimal scenario"
matrix list DATA_FIGURE, format(%9.4f)
global relabel 1 "0-20" 2 "20-40" 3 "40-60" 4 "60-80" 5 "80-100" 6 "90-95" 7 "95-99" 8 "99-100"
svmat DATA_FIGURE, names(data_b)

gen byte aux = _n if data_a1!=. // For the axis
replace data_a1 = data_a1*100
replace data_b1 = data_b1*100
qui gen aux_a = aux - 0.2
qui gen aux_b = aux + 0.2
gen var_1= data_a1 + 7
gen var_2= data_b1 + 7
format data* %10.2f

grstyle set legend 10, inside nobox

twoway (bar data_a1 aux_a, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.37)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) ///
yti("Share of total income tax revenues (%)", margin(small)size(medlarge)) ylabel(0(20)100,nogrid format(%10.0f) labsize(medlarge)) ///
xti("Percentiles of household gross income", margin(small)size(medlarge)) xlabel($relabel , nogrid format(%10.1f) labsize(medlarge))) ///
(scatter var_1 aux_a, msym(none) mlab(data_a1) mlabpos(12) mlabcolor(black) mlabsize(medium) mlabangle(ver)) ///
(bar data_b1 aux_b, barwidth(.37)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(scatter var_2 aux_b, msym(none) mlab(data_b1) mlabpos(12) mlabcolor(black) mlabsize(medium) mlabangle(ver)), ///
legend(on order(1 "Actual" 3 "Optimal") rows(2) region(lstyle(none)) size(medlarge) region(lwidth(none))) 
graph export "$figurePath\Share_Tax_Revenue_Income_Percentiles_Optimal.pdf", as(pdf) replace
drop data* aux* var*
mat drop _all

*******************************************************************************
// 8. Fitness of the HSV function to the data (Average Tax Rate by Multiples of Income)
********************************************************************************
// Build matrix of income thresholds from which models were estimated.
save "Intermediate.dta", replace
use "$dataPath\Threshold_Estimates_H2015.dta", clear
mat THRESHOLDS = J(1,5,.)
local list1 log hsv power gs 
foreach i in `list1' {
	mkmat optthres_`i' in 1
	local col : list posof "`i'" in list1 	
	mat THRESHOLDS[1,`col'] = optthres_`i'[1,1]
	mat drop optthres_`i'
}
mat colnames THRESHOLDS = log hsv power gs 
mat rownames THRESHOLDS = thresholds
use "Intermediate.dta", clear
erase "Intermediate.dta"
// Multiples of income 
gen agin_scale = .
local min = 0.6
local max = 10
local step = 0.4
matrix AUX1 = `min'
forvalues i=`min'(`step')`max' {
	replace agin_scale = `i' if grossinc > `i'*meaninc
	if `i'>`min'{
		matrix AUX1 = (AUX1 \ `i' )
	}
}
replace agin_scale = . if grossinc > (`max'+`step')*meaninc
matrix colnames AUX1 = multiples
// Data of Average Tax Rates
tabstat avgtaxg [aw=Factor] if $restrict_g & $restrict_sample , by(agin_scale) stat(mean) nototal format(%25.4f) save
** Lines below avoid having to type r(Stat1) \ r(Stat2) \ ... \ r(Stat48) 
local n = int(((`max'-`min')/`step')+1)
matrix AUX2 = r(Stat1)
forvalues i = 2/`n'{
	matrix AUX2 = (AUX2 \ r(Stat`i')) 
}
matrix colnames AUX2 = data
// Log Model
estimates use "$dataPath\Parameters_Thresholds_H2015", number(1)
mat AUX3 = J(rowsof(AUX1),1,.)
local t = rowsof(AUX1)
forval i=1/`t'{
	if AUX1[`i',1] < THRESHOLDS[1,1] {
		mat AUX3[`i',1] = 0
	}
	else { 
		mat AUX3[`i',1] = _b[_cons] + _b[lnincome]*ln(AUX1[`i',1])
	}
}
mat colnames AUX3 = log
// Heathcote Model
estimates use "$dataPath\Parameters_Thresholds_H2015", number(2)
mat AUX4 = J(rowsof(AUX1),1,.)
forval i=1/`t'{
	if AUX1[`i',1] < THRESHOLDS[1,2] {
		mat AUX4[`i',1] = 0
	}
	else {
		mat AUX4[`i',1] = 1 - _b[lambda:_cons]*AUX1[`i',1]^(-_b[tau:_cons]) 
	}
}
mat colnames AUX4 = hsv
// Power Model
estimates use "$dataPath\Parameters_Thresholds_H2015", number(3)
mat AUX5 = J(rowsof(AUX1),1,.)
forval i=1/`t'{
	if AUX1[`i',1] < THRESHOLDS[1,3] {
		mat AUX5[`i',1] = 0
	}
	else {
		mat AUX5[`i',1] = _b[delta:_cons] + _b[gamma:_cons]*AUX1[`i',1]^_b[epsilon:_cons]
	}
}
mat colnames AUX5 = power
// Gouveia and Strauss
estimates use "$dataPath\Parameters_Thresholds_H2015", number(4)
sum meaninc [aw=Factor]
local m = r(mean)/1000
mat AUX6 = J(rowsof(AUX1),1,.)
forval i=1/`t'{
	if AUX1[`i',1] < THRESHOLDS[1,4] {
		mat AUX6[`i',1] = 0
	}
	else {
		mat AUX6[`i',1] = _b[b:_cons] * (1 - (_b[s:_cons]*(AUX1[`i',1]*`m')^_b[p:_cons]+1)^(-1/_b[p:_cons]))
	}
}
mat colnames AUX6 = gs
drop agin_scale
// Data for Figure 
matrix DATA_FIGURE = AUX1,AUX2,AUX3,AUX4,AUX5,AUX6
di "Average Tax Rate by Multiples of Income"
matrix list DATA_FIGURE, format(%9.4f)

// Figure 
grstyle set legend 12, inside nobox

svmat DATA_FIGURE, names(data)
twoway (line data2 data1, lwidth(medthick) plotregion(fcolor(white) lcolor(black) lwidth(medthin)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
lcolor(black) yti("Average tax rates",margin(small)size(medlarge)) ylabel(0(.1).40,nogrid format(%10.2f) labsize(medlarge)) ///
xti("Multiples of mean household gross income",margin(small) size(medlarge)) xlabel(1(1)10,nogrid format(%10.1f) labsize(medlarge)) ///
legend(lab(1 "Data")lab(2 "Estimated HSV tax function")order(1 2) region(lstyle(none)) size(medlarge) region(lwidth(none))) ///
legend(region(lcolor(white))rows(1)size(medlarge))) ///
(line data4 data1, lcolor("218 108 122") lwidth(medthick) lpattern(dash)) ///
// (connected data6 data1, lcolor("0 64 129") lwidth(thin) mcolor("0 64 129") msymbol(X)) 
// (connected data3 data1, lcolor(blue) lwidth(medthick) mcolor(blue) msymbol(triangle)) ///
// (connected data5 data1, lcolor(orange) lwidth(medthick) mcolor(orange) msymbol(plus)) 
graph export "$figurePath\Avg_Tax_Rates_Data_HSV.pdf", as(pdf) replace
drop data*
mat drop _all

*****************************************************
// 9. Average Tax Rates and Share of Tax Returns by Multiples of Income 
*****************************************************
// Build matrix of income thresholds from which models were estimated.
// Multiples of income 
gen agin_scale = .
local min = 0.5
local max = 10
local step = 0.5
matrix AUX1 = `min'
forvalues i=`min'(`step')`max' {
	replace agin_scale = `i' if grossinc > `i'*meaninc
	if `i'>`min'{
		matrix AUX1 = (AUX1 \ `i' )
	}
}
replace agin_scale = . if grossinc > (`max'+`step')*meaninc
matrix colnames AUX1 = multiples
// Data of Average Tax Rates
tabstat avgtaxg [aw=Factor] if $restrict_g & $restrict_sample , by(agin_scale) stat(mean) nototal format(%25.4f) save
// Lines below avoid having to type r(Stat1) \ r(Stat2) \ ... \ r(Stat48) 
local n = int(((`max'-`min')/`step')+1)
matrix AUX2 = r(Stat1)
forvalues i = 2/`n'{
	matrix AUX2 = (AUX2 \ r(Stat`i')) 
}
matrix colnames AUX2 = avgtax
// Share of tax returns
tabstat Factor if $restrict_g & $restrict_sample , by(agin_scale) stat(sum) nototal format(%25.4f) save
// Lines below avoid having to type r(Stat1) \ r(Stat2) \ ... \ r(Stat48) 
local n = int(((`max'-`min')/`step')+1)
matrix AUX3 = r(Stat1)
forvalues i = 2/`n'{
	matrix AUX3 = (AUX3 \ r(Stat`i')) 
}
tabstat Factor if $restrict_g & $restrict_sample , stat(sum) format(%25.4f) save
mat AUX = r(StatTotal)
forval i=1/`=rowsof(AUX3)' {
	mat AUX3[`i',1] = AUX3[`i',1] / AUX[1,1]
}
matrix colnames AUX3 = share
drop agin_scale
// Data for Figure 
matrix DATA_FIGURE = AUX1,AUX2,AUX3
matrix list DATA_FIGURE, format(%9.4f)

// Figure 
svmat DATA_FIGURE, names(data)
twoway bar data3 data1, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.37)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) ///
yti("Share of tax returns", margin(small)size(medlarge)) ylabel(0(.05).40,nogrid format(%10.2f) labsize(medlarge)) ///
xti("Multiples of mean household gross income", margin(small)size(medlarge)) xlabel(0(0.5)10,nogrid format(%10.1f) labsize(medlarge))
graph export "$figurePath\Share_Tax_Returns.pdf", as(pdf) replace
drop data*
mat drop _all

********************************************************************************
// 10. Welfare changes by household types
********************************************************************************
// Axis x and y just on wealth retiree plot
foreach j in "wealth" {
clear all
set more off
import excel "$tablePath\WelfareChange_HH_`j'.xlsx", sheet("Hoja1") clear
// Transpose the imported data to have it in Stata structure
xpose, clear
// Rename variables
foreach m of varlist v3-v11{
replace `m' = `m'*100
}
rename v1-v11 (pctile_lb pctile_ub all_s s_1 s_2 s_3 s_4 s_5 s_6 s_7 s_8)
gen byte aux = _n // For the axis
qui gen aux_1 = aux - 0.3
qui gen aux_2 = aux - 0.1
qui gen aux_3 = aux + 0.1
qui gen aux_4 = aux + 0.3
global relabel 1 "0-10" 2 "10-20" 3 "20-30" 4 "30-40" 5 "40-50" 6 "50-60" 7 "60-70" 8 "70-80" 9 "80-90" 10 "90-100"
// Retirees
twoway (bar s_5 aux_1, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.18)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti("% change in consumption", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti("Percentiles of household gross income", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large))) ///
(bar s_6 aux_2, barwidth(.18)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(bar s_7 aux_3, barwidth(.18)fcolor("255 190 93")lcolor("255 190 93")lwidth(medthin)) ///
(bar s_8 aux_4, barwidth(.18)fcolor("26 140 255")lcolor("26 140 255")lwidth(medthin)), ///
legend(on label(1 "s=5") label(2 "s=6") label(3 "s=7") label(4 "s=8") rows(1) region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_retirees.pdf", as(pdf) replace
}

// Axis x just on retiree plots (not retirees wealth)
foreach j in "income" "VF" {
clear all
set more off
import excel "$tablePath\WelfareChange_HH_`j'.xlsx", sheet("Hoja1") clear
// Transpose the imported data to have it in Stata structure
xpose, clear
// Rename variables
foreach m of varlist v3-v11{
replace `m' = `m'*100
}
rename v1-v11 (pctile_lb pctile_ub all_s s_1 s_2 s_3 s_4 s_5 s_6 s_7 s_8)
gen byte aux = _n // For the axis
qui gen aux_1 = aux - 0.3
qui gen aux_2 = aux - 0.1
qui gen aux_3 = aux + 0.1
qui gen aux_4 = aux + 0.3
global relabel 1 "0-10" 2 "10-20" 3 "20-30" 4 "30-40" 5 "40-50" 6 "50-60" 7 "60-70" 8 "70-80" 9 "80-90" 10 "90-100"
// Retirees
twoway (bar s_5 aux_1, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.18)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti(" ", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti("Percentiles of household gross income", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large))) ///
(bar s_6 aux_2, barwidth(.18)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(bar s_7 aux_3, barwidth(.18)fcolor("255 190 93")lcolor("255 190 93")lwidth(medthin)) ///
(bar s_8 aux_4, barwidth(.18)fcolor("26 140 255")lcolor("26 140 255")lwidth(medthin)), ///
legend(on label(1 "s=5") label(2 "s=6") label(3 "s=7") label(4 "s=8") rows(1) region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_retirees.pdf", as(pdf) replace
}

// Axis y just on wealth plots (not retirees wealth)
foreach j in "wealth" {
clear all
set more off
import excel "$tablePath\WelfareChange_HH_`j'.xlsx", sheet("Hoja1") clear
// Transpose the imported data to have it in Stata structure
xpose, clear
// Rename variables
foreach m of varlist v3-v11{
replace `m' = `m'*100
}
rename v1-v11 (pctile_lb pctile_ub all_s s_1 s_2 s_3 s_4 s_5 s_6 s_7 s_8)
gen byte aux = _n // For the axis
qui gen aux_1 = aux - 0.3
qui gen aux_2 = aux - 0.1
qui gen aux_3 = aux + 0.1
qui gen aux_4 = aux + 0.3
global relabel 1 "0-10" 2 "10-20" 3 "20-30" 4 "30-40" 5 "40-50" 6 "50-60" 7 "60-70" 8 "70-80" 9 "80-90" 10 "90-100"
// For any s
twoway bar all_s aux, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.67)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti("% change in consumption", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti(" ", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large)) ///
legend(on label(1 "For any state (s)") region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_all.pdf", as(pdf) replace
// Workers
twoway (bar s_1 aux_1, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.18)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti("% change in consumption", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti(" ", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large))) ///
(bar s_2 aux_2, barwidth(.18)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(bar s_3 aux_3, barwidth(.18)fcolor("255 190 93")lcolor("255 190 93")lwidth(medthin)) ///
(bar s_4 aux_4, barwidth(.18)fcolor("26 140 255")lcolor("26 140 255")lwidth(medthin)), ///
legend(on label(1 "s=1") label(2 "s=2") label(3 "s=3") label(4 "s=4") rows(1) region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_workers.pdf", as(pdf) replace
}

// Rest with no axis
foreach j in "income" "VF" {
clear all
set more off
import excel "$tablePath\WelfareChange_HH_`j'.xlsx", sheet("Hoja1") clear
// Transpose the imported data to have it in Stata structure
xpose, clear
// Rename variables
foreach m of varlist v3-v11{
replace `m' = `m'*100
}
rename v1-v11 (pctile_lb pctile_ub all_s s_1 s_2 s_3 s_4 s_5 s_6 s_7 s_8)
qui replace all_s = -15 if all_s<-15
qui replace s_4 = -15 if s_4<-15
gen byte aux = _n // For the axis
qui gen aux_1 = aux - 0.3
qui gen aux_2 = aux - 0.1
qui gen aux_3 = aux + 0.1
qui gen aux_4 = aux + 0.3
global relabel 1 "0-10" 2 "10-20" 3 "20-30" 4 "30-40" 5 "40-50" 6 "50-60" 7 "60-70" 8 "70-80" 9 "80-90" 10 "90-100"
// For any s
twoway bar all_s aux, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.67)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti(" ", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti(" ", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large)) ///
legend(on label(1 "For any state (s)") region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_all.pdf", as(pdf) replace
// Workers
twoway (bar s_1 aux_1, plotregion(fcolor(white) lcolor(black)) graphregion(fcolor(white) ifcolor(white) color(white) icolor(white)) bgcolor(white) ///
barwidth(.18)fcolor("218 108 122")lcolor("218 108 122")lwidth(medthin) yline(0, lpattern(dot) lwidth(thin) lcolor(black)) ///
yti(" ", margin(small)size(large)) ylabel(-40(10)30,nogrid format(%10.0f) labsize(large)) ///
xti(" ", margin(small)size(large)) xlabel($relabel , nogrid format(%10.1f) labsize(large))) ///
(bar s_2 aux_2, barwidth(.18)fcolor("0 64 129")lcolor("0 64 129")lwidth(medthin)) ///
(bar s_3 aux_3, barwidth(.18)fcolor("255 190 93")lcolor("255 190 93")lwidth(medthin)) ///
(bar s_4 aux_4, barwidth(.18)fcolor("26 140 255")lcolor("26 140 255")lwidth(medthin)), ///
legend(on label(1 "s=1") label(2 "s=2") label(3 "s=3") label(4 "s=4") rows(1) region(lstyle(none)) size(large) region(lwidth(none))) 
graph export "$figurePath\WelfareChange_`j'_workers.pdf", as(pdf) replace
}
