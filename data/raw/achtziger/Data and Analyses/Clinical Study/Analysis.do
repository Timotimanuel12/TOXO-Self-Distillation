*Change directory to where the RawData file is
cd "D:\Toxo Data to Upload\Clinical Study"

*import data
import delimited RawData.csv, delimiter(comma) bindquote(strict) varnames(1) case(preserve)  clear  

*Check the number of subj
tab subj
sum age female income
*Check the number of subj with negative RhD
tab subj if bloodpos==0
tab positive if bloodpos==0

* Lot transformed sum of the reaciton times from the reaction-time task of Havlivceketal01 and Novotnaetal08
sum logrt

*FMM classification depending on bloodtypes
fmm 2: regress logrt if bloodpos==0
estat lcmean 
estat lcprob  
predict clas* if bloodpos==0, classposteriorpr

*create variable of imputed positive based on the classification (take care of which group is slower).
*Be careful, this depends on the mean of RTs on the estimated classes it could be different for bloodpos==0 and ==1
gen slowerclass=1 if clas1<clas2   
replace slowerclass=0  if clas1>clas2 
 tab slowerclass if bloodpos==0, miss
*check if the 2 classifications are related
tab positive if bloodpos==0
tab slowerclass if bloodpos==0
tab slowerclass positive, chi2

 ttest logrt  if bloodpos==0, by(slowerclass)
 ranksum logrt  if bloodpos==0, by(slowerclass)

*Check if the 2 classificaitions are correlated
* Figure 2 (Right Panel)
preserve
 collapse  resultig clas2 bloodpos logrt, by(subj)
 sum resultig clas2 if bloodpos==0
pwcorr resultig clas2 if bloodpos==0, sig
gen logresultig = log(resultig)
sum logresultig
display log(3.0)
 #delimit ;
twoway (scatter logresultig clas2 if bloodpos==0, jitter(3)),
  title(" ") graphregion(color(white) fcolor(white) lcolor(white)) plotregion(fcolor(white) lcolor(white) ) 
 ytitle("Toxoplasma Gondii IgG Log(IU/ml)", size(large))
 xtitle("Posterior Probability of Being Infected", size(large))
 legend(off)   yscale(range(0 4))  xscale(range(0 1))
xline(0.5) yline(1.0986123) text(1.3 0.1  "Cutoff 3.0", color(red) size(large))
text(0.6 0.52 "Classified as Infected", color(red) orient(vertical)size(large))
;#delimit cr
restore

* Figure 2 (Left Panel)
preserve
gen Comp = 100
gen G =_n
 *Medical test sensitivity and specificity
gen sensitvmed = 98.25 if G==1
gen specificmed = 98.45 if G==2

gen sensitiv = 91.66 if G==1
* 11 out of 12 correctly
* 1 out of 12 false negative

gen specific = 97.01 if G==3
* false positive 2 our of 67 incorrectly

keep if G<4
*replace logrtneg = . if logrtneg<.707909 | logrtneg>5.244889
#delimit ;
twoway (bar  Comp G if G == 1 | G == 3, lcolor(black) fcolor(black))
	   (bar sensitiv G, lcolor(black) fcolor(green)) 
       (bar specific G, lcolor(black) fcolor(orange))
,
  graphregion(color(white) fcolor(white) lcolor(white)margin(0)) plotregion(fcolor(white) lcolor(white))
 ytitle("Correctly Classified (%)" " ", size(large))  yscale(range(0 120)) ylabel(0(20)100, grid)
legend(off) text(105 1 "11 of 12", size(4)) text(105 3 "65 of 67", size(4))  xlabel(1 "Sensitivity"  3 "Specificity", labsize(large) noticks) xtitle(" ")
;#delimit cr
restore

***sanity check the "positive" should be more likely to own cats
bysort slowerclass: sum cat if bloodpos==0 
ranksum cat if bloodpos==0, by(slowerclass)
bysort positive: sum cat   if bloodpos==0 
ranksum cat if bloodpos==0, by(positive)

*create variables for the forest plot: difference in means and CI (for variables which are in weird scales, normalize them)
* Owning or having owned a cat
preserve
keep if bloodpos==0 
ttest cat if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(cat) if positive==1
egen B = mean(cat) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Own or owned a cat"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save cat, replace
restore
*

*Employment status (Full+Part+Due Start Soon vs. Unemployed+Not in paid work)
tab withjob if bloodpos==0
bysort slowerclass: sum withjob if bloodpos==0
ranksum withjob if bloodpos==0, by(slowerclass)
bysort positive: sum withjob if bloodpos==0
ttest withjob if bloodpos==0, by(positive)
ranksum withjob if bloodpos==0, by(positive)
*

*Entrepreneurship/run your own business  (Yes/No)
bysort slowerclass: sum entrepreneur if bloodpos==0
ranksum entrepreneur if bloodpos==0, by(slowerclass)
bysort positive: sum entrepreneur if bloodpos==0
ranksum entrepreneur if bloodpos==0, by(positive)
ttest entrepreneur if bloodpos==0, by(positive)
*

*Have you ever considered a career move as a (paid) employee in a startup?  (Yes/No)
bysort slowerclass: sum Startup if bloodpos==0
ranksum Startup if bloodpos==0, by(slowerclass)
bysort positive: sum Startup if bloodpos==0
ranksum Startup if bloodpos==0, by(positive)
ttest Startup if bloodpos==0, by(positive)
*
 
*Neurodivergent (Yes/No)
bysort slowerclass: sum Neurodiv if bloodpos==0
ranksum Neurodiv if bloodpos==0, by(slowerclass)
bysort positive: sum Neurodiv if bloodpos==0
ranksum Neurodiv if bloodpos==0, by(positive)
ttest Neurodiv if bloodpos==0, by(positive)
*


*How many units of alcohol do you drink on average per week?
bysort slowerclass: sum alchool if bloodpos==0
ranksum alchool if bloodpos==0, by(slowerclass)
bysort positive: sum alchool if bloodpos==0
ranksum alchool if bloodpos==0, by(positive)
*
preserve
keep if bloodpos==0 
sum alchool, detail
replace alchool = alchool/r(max)
ttest alchool if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(alchool) if slowerclass==1
egen B = mean(alchool) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Nr. alcohol units /week"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff
save alchool, replace
restore
*


*Do you have – or have you had – a diagnosed, on-going mental health/illness/condition? (yes/no)
bysort slowerclass: sum mental if bloodpos==0
ranksum mental if bloodpos==0, by(slowerclass)
bysort positive: sum mental if bloodpos==0
ranksum mental if bloodpos==0, by(positive)
ttest mental if bloodpos==0, by(positive)
*


*Have you tried to access mental health support on the NHS in the last 12 months?[omitted “on the NHS” in the clinical study](yes/no)
bysort slowerclass: sum mental2 if bloodpos==0
ranksum mental2 if bloodpos==0, by(slowerclass)
bysort positive: sum mental2 if bloodpos==0
ranksum mental2 if bloodpos==0, by(positive)
ttest mental2 if bloodpos==0, by(positive)
*


*Do you experience depression? (yes/no)
bysort slowerclass: sum depression if bloodpos==0
ranksum depression if bloodpos==0, by(slowerclass)
bysort positive: sum depression if bloodpos==0
ranksum depression if bloodpos==0, by(positive)
*
preserve
keep if bloodpos==0 
ttest depression if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(depression) if slowerclass==1
egen B = mean(depression) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Depression"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save depression, replace
restore
*


*Do you experience anxiety?  (yes/no)
bysort slowerclass: sum anxiety if bloodpos==0
ranksum anxiety if bloodpos==0, by(slowerclass)
bysort positive: sum anxiety if bloodpos==0
ranksum anxiety if bloodpos==0, by(positive)
*
preserve
keep if bloodpos==0 
tab anxiety
ttest anxiety if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(anxiety) if positive==1
egen B = mean(anxiety) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Anxiety"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save anxiety, replace
restore
*

*Do you consider yourself to have attention deficit disorder (ADD)/attention deficit hyperactivity disorder (ADHD)?(yes/no)
bysort slowerclass: sum add if bloodpos==0 
ranksum add if bloodpos==0  , by(slowerclass)
bysort positive: sum add if bloodpos==0  
ranksum add if bloodpos==0   , by(positive)
ttest add if bloodpos==0, by(positive)
*


*What is your current smoking status? (current, recent, former vs. Never + <100 in my lifetime)
bysort slowerclass: sum smoke if bloodpos==0 
ranksum smoke  if bloodpos==0  , by(slowerclass)
bysort positive: sum smoke if bloodpos==0  
ranksum  smoke  if bloodpos==0   , by(positive)
ttest smoke if bloodpos==0, by(positive)


*What is your personal income per year (after tax) in GBP [or CHF]?
bysort slowerclass: sum income if bloodpos==0 
ranksum income if bloodpos==0  , by(slowerclass)
bysort positive: sum income if bloodpos==0  
ranksum income if bloodpos==0   , by(positive)
ttest income if bloodpos==0, by(positive)



*Economic Related Major (Management, Marketing, Commerce, Finance, Economics, Accounting vs. Everything else)
bysort slowerclass: sum econstudies if bloodpos==0 
ranksum econstudies if bloodpos==0  , by(slowerclass)
bysort positive: sum econstudies if bloodpos==0  
ranksum econstudies if bloodpos==0   , by(positive)
ttest econstudies if bloodpos==0, by(positive)


*Number of sexual partners
bysort slowerclass: sum number if bloodpos==0 
ranksum number if bloodpos==0  , by(slowerclass)
bysort positive: sum number if bloodpos==0  
ranksum number if bloodpos==0   , by(positive)
ttest number if bloodpos==0, by(positive)



*Lefthand (yes/no)
bysort slowerclass: sum lefthand if bloodpos==0 
ranksum lefthand if bloodpos==0  , by(slowerclass)
bysort positive: sum lefthand if bloodpos==0  
ranksum lefthand if bloodpos==0   , by(positive)
ttest lefthand if bloodpos==0, by(positive)


*Have investments in retirement (yes/no)
bysort slowerclass: sum retirement if bloodpos==0 
ranksum retirement if bloodpos==0  , by(slowerclass)
bysort positive: sum retirement if bloodpos==0  
ranksum retirement if bloodpos==0   , by(positive)
ttest retirement if bloodpos==0, by(positive)

*Have investments (yes/no)
bysort slowerclass: sum investment if bloodpos==0 
ranksum investment if bloodpos==0  , by(slowerclass)
bysort positive: sum investment if bloodpos==0  
ranksum investment if bloodpos==0   , by(positive)
ttest investment if bloodpos==0, by(positive)



* Risk attitude (Falketal18)
bysort slowerclass: sum risk if bloodpos==0 
ranksum risk if bloodpos==0  , by(slowerclass)
bysort positive: sum risk if bloodpos==0  
ranksum risk if bloodpos==0   , by(positive)
preserve
keep if bloodpos==0 
sum risk, detail
replace risk =risk/r(max)
ttest risk if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(risk) if slowerclass==1
egen B = mean(risk) if slowerclass==0
collapse A B se effectsize
gen diff = B-A
sum diff, detail
gen Paradigm = "Risk aversion"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save risk, replace
restore
*

*Number of car accidents
bysort slowerclass: sum accidents if bloodpos==0 
ranksum accidents if bloodpos==0  , by(slowerclass)
bysort positive: sum accidents if bloodpos==0  
ranksum accidents if bloodpos==0   , by(positive)
ttest accidents if bloodpos==0, by(positive)


*Number of accidents at work
bysort slowerclass: sum accidentswork if bloodpos==0 
ranksum accidentswork if bloodpos==0  , by(slowerclass)
bysort positive: sum accidentswork if bloodpos==0  
ranksum accidentswork if bloodpos==0   , by(positive)
ttest accidentswork if bloodpos==0, by(positive)


*CRT score
bysort slowerclass: sum crtscore if bloodpos==0 
ranksum crtscore if bloodpos==0  , by(slowerclass)
bysort positive: sum crtscore if bloodpos==0  
ranksum crtscore if bloodpos==0   , by(positive)
preserve
keep if bloodpos==0 
sum crtscore, detail
replace crtscore = crtscore/r(max)
ttest crtscore if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(crtscore) if slowerclass==1
egen B = mean(crtscore) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "CRT score"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff diff05 diff95
save crtscore, replace
restore
*


*Stress-scale score Cohen et al. (1994):
bysort slowerclass: sum stress if bloodpos==0 
ranksum stress if bloodpos==0  , by(slowerclass)
bysort positive: sum stress if bloodpos==0  
ranksum stress if bloodpos==0   , by(positive)
preserve
keep if bloodpos==0 
sum stress, detai
replace stress = stress/r(max)
ttest stress if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(stress) if positive==1
egen B = mean(stress) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Stress"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save stress, replace
restore
*

*Self-control score (TangneyBaumeisterBoone04)
bysort slowerclass: sum selfcontrol  if bloodpos==0 
ranksum selfcontrol  if bloodpos==0  , by(slowerclass)
bysort positive: sum selfcontrol  if bloodpos==0  
ranksum selfcontrol  if bloodpos==0   , by(positive)
preserve
keep if bloodpos==0 
sum selfcontrol
replace selfcontrol = selfcontrol/r(max)
ttest selfcontrol if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(selfcontrol) if positive==1
egen B = mean(selfcontrol) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Self control"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save selfcontrol, replace
restore
*

* Big 5 (Gosling, S. D., P. J. Rentfrow, and W. B. Swann Jr., 2003)
*Extraversion
bysort slowerclass: sum extraversion  if bloodpos==0 
ranksum extraversion  if bloodpos==0  , by(slowerclass)
bysort positive: sum extraversion  if bloodpos==0  
ranksum extraversion  if bloodpos==0   , by(positive)
ttest extraversion  if bloodpos==0, by(positive)


*Agreeable
bysort slowerclass: sum  agreeableness  if bloodpos==0 
ranksum  agreeableness  if bloodpos==0  , by(slowerclass)
bysort positive: sum  agreeableness  if bloodpos==0  
ranksum  agreeableness  if bloodpos==0   , by(positive)
ttest agreeableness  if bloodpos==0, by(positive)


*Conscient
bysort slowerclass: sum conscientiousness  if bloodpos==0 
ranksum  conscientiousness if bloodpos==0  , by(slowerclass)
bysort positive: sum  conscientiousness if bloodpos==0  
ranksum conscientiousness  if bloodpos==0   , by(positive)
ttest conscientiousness  if bloodpos==0, by(positive)


*Stability
bysort slowerclass: sum stability  if bloodpos==0 
ranksum  stability if bloodpos==0  , by(slowerclass)
bysort positive: sum  stability if bloodpos==0  
ranksum stability  if bloodpos==0   , by(positive)
preserve
keep if bloodpos==0 
sum stability, detail
replace stability = stability/r(max)
ttest stability if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(stability) if positive==1
egen B = mean(stability) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Neuroticism"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save stability, replace
restore
*

*Openess
bysort slowerclass: sum openness if bloodpos==0 
ranksum openness if bloodpos==0  , by(slowerclass)
bysort positive: sum openness if bloodpos==0  
ranksum openness  if bloodpos==0   , by(positive)
ttest openness  if bloodpos==0, by(positive)


*Sensationseeking (Cooketal15)
bysort slowerclass: sum sensationseeking if bloodpos==0 
ranksum sensationseeking if bloodpos==0  , by(slowerclass)
bysort positive: sum sensationseeking if bloodpos==0  
ranksum sensationseeking if bloodpos==0   , by(positive)
*
preserve
keep if bloodpos==0 
sum sensationseeking, detail
replace sensationseeking= sensationseeking/r(max)
ttest sensationseeking if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(sensationseeking) if slowerclass==1
egen B = mean(sensationseeking) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Sensation seeking"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save sensationseeking, replace
restore
*

*Patience (Falketal18)
bysort slowerclass: sum Patience if bloodpos==0 
ranksum Patience if bloodpos==0  , by(slowerclass)
bysort positive: sum Patience if bloodpos==0  
ranksum Patience if bloodpos==0   , by(positive)
*
preserve
keep if bloodpos==0 
sum Patience, detail
replace Patience= Patience/r(max)
ttest Patience if bloodpos==0, by(positive)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(Patience) if positive==1
egen B = mean(Patience) if positive==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Patience"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save Patience, replace
restore
*

* Merge toghether signle tests to make a "forest plot" type graph
use depression, clear
append using anxiety
append using stress
append using sensationseeking 
append using cat

append using alchool
append using selfcontrol
append using Patience
append using risk
append using stability
append using crtscore

meta set diff diff05 diff95, civarlevel(95)  studylabel(Paradigm) 

* Figure 5
*manually reorder variables based on exposition in the paper
sum diff diff05 diff95
  metan diff diff05 diff95 ,dp(3) nostats  nointeger nooverall scheme(s1color) nowt null(0) lcols(Paradigm)  graphregion(color(white) fcolor(white)  lcolor(white)) plotregion(fcolor(white) lcolor(white) )   nobox  xscale(range(-40 60)) xlabels(-40(40)80, labsize(small)) force title("Toxoplasmosis-infected - Toxoplasmosis-free") xtitle("Difference between Means (in %)", size(large))
	*trick the command in order to get results without having CIs
replace effectsize = abs(effectsize)
gen effectsize05 = effectsize-0.0001 
gen effectsize95 = effectsize+0.0001

sum effectsize
  metan effectsize effectsize05  effectsize95 ,dp(3) nostats  nointeger nooverall scheme(s1color) nowt null(0) lcols(Paradigm)  graphregion(color(white) fcolor(white)  lcolor(white)) plotregion(fcolor(white) lcolor(white) )   nobox xscale(range(0 1.5)) xlabels(0(0.5)1.5, labsize(small)) force title("Toxoplasmosis-infected - Toxoplasmosis-free") xtitle("Effect Size", size(large))
  
	