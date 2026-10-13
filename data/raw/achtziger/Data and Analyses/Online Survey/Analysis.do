*Change directory to where the RawData file is
cd "D:\Toxo Data to Upload\Online Survey"

*import data
import delimited RawData.csv, delimiter(comma) bindquote(strict) varnames(1) case(preserve)  clear  

*Check the number of subj
tab subj
sum age female income
sum age ,detail

* Sum of reaciton time from the reaction-time task of Havlivceketal01 and Novotnaetal08
*summarize response times (meanrt2 is without outliers)
sum meanrt meanrt2

*use log rt
gen logrt2= log(meanrt2)

*FMM classification depending on bloodtypes
fmm 2: regress logrt2 if bloodpos==1
*get weights
estat lcmean  
estat lcprob  
*get probability each participant is belonging to each class
predict pr* if bloodpos==1, classposteriorpr

fmm 2: regress logrt2 if bloodpos==0
estat lcmean 
estat lcprob  
predict clas* if bloodpos==0, classposteriorpr

*create variable of imputed afflicted based on the classification (take care of which group is slower).
*Be careful, this depends on the mean of RTs on the estimated classes it could be different for bloodpos==1 and ==0!!
gen slowerclass=1 if pr1<pr2 & bloodpos==1
replace slowerclass=0  if pr1>pr2 & bloodpos==1

replace slowerclass=1 if clas1<clas2 & bloodpos==0
replace slowerclass=0  if clas1>clas2 & bloodpos==0 
 tab slowerclass if bloodpos==0 

*check size of the groups
tab slowerclass if bloodpos==1
tab slowerclass if  bloodpos==0 
tab bloodpos

 ttest logrt  if bloodpos==0, by(slowerclass)
 ranksum logrt  if bloodpos==0, by(slowerclass)

***sanity check the "afflicted" should be more likely to own cats
bysort slowerclass: sum cat if bloodpos==0 
ranksum cat if bloodpos==0, by(slowerclass)
bysort slowerclass: sum cat if bloodpos==1 
ranksum cat if bloodpos==1, by(slowerclass)

*create variables for the forest plot: difference in means and CI (for variables which are in weird scales, normalize them)

* Owning or having owned a cat
preserve
keep if bloodpos==0 
ttest cat if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(cat) if slowerclass==1
egen B = mean(cat) if slowerclass==0
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
bysort slowerclass: sum withjob if bloodpos==0
ranksum withjob if bloodpos==0, by(slowerclass)
bysort slowerclass: sum withjob if bloodpos==1
ranksum withjob if bloodpos==1, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest withjob if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(withjob) if slowerclass==1
egen B = mean(withjob) if slowerclass==0
collapse A B se effectsize
gen diff = B-A
sum diff, detail
gen Paradigm = "Unemployed"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save withjob, replace
restore
*


*Entrepreneurship/run your own business  (Yes/No)
bysort slowerclass: sum entrepreneur if bloodpos==1
ranksum entrepreneur if bloodpos==1, by(slowerclass)
bysort slowerclass: sum entrepreneur if bloodpos==0
ranksum entrepreneur if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest entrepreneur if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(entrepreneur) if slowerclass==1
egen B = mean(entrepreneur) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Entrepreneur initiative"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save entrepreneur, replace
restore
*

*Have you ever considered a career move as a (paid) employee in a startup?  (Yes/No)
bysort slowerclass: sum Startup if bloodpos==1
ranksum Startup if bloodpos==1, by(slowerclass)
bysort slowerclass: sum Startup if bloodpos==0
ranksum Startup if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest Startup if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(Startup) if slowerclass==1
egen B = mean(Startup) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Startup exp."
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save Startup, replace
restore
*


*Went to University (Yes/No)
bysort slowerclass: sum university if bloodpos==1
ranksum university if bloodpos==1, by(slowerclass)
bysort slowerclass: sum university if bloodpos==0
ranksum university if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest university if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(university) if slowerclass==1
egen B = mean(university) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "University exp."
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save university, replace
restore
*

 
*Neurodivergent (Yes/No)
bysort slowerclass: sum Neurodiv if bloodpos==1
ranksum Neurodiv if bloodpos==1, by(slowerclass)
bysort slowerclass: sum Neurodiv if bloodpos==0
ranksum Neurodiv if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest Neurodiv if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(Neurodiv) if slowerclass==1
egen B = mean(Neurodiv) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Neurodivergent"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save Neurodiv, replace
restore
*


*How many units of alcohol do you drink on average per week?
bysort slowerclass: sum alchool if bloodpos==1
ranksum alchool if bloodpos==1, by(slowerclass)
bysort slowerclass: sum alchool if bloodpos==0
ranksum alchool if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
sum alchool, detail
replace alchool = alchool/r(max)
ttest alchool if bloodpos==0, by(slowerclass)
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
bysort slowerclass: sum mental if bloodpos==1
ranksum mental if bloodpos==1, by(slowerclass)
bysort slowerclass: sum mental if bloodpos==0
ranksum mental if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest mental if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(mental) if slowerclass==1
egen B = mean(mental) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Mental disorders (12 m.)"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save mental, replace
restore
*

*Have you tried to access mental health support on the NHS in the last 12 months?[omitted “on the NHS” in the clinical study](yes/no)
bysort slowerclass: sum mental2 if bloodpos==1
ranksum mental2 if bloodpos==1, by(slowerclass)
bysort slowerclass: sum mental2 if bloodpos==0
ranksum mental2 if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest mental2 if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(mental2) if slowerclass==1
egen B = mean(mental2) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Mental disorders ongoing"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save mental2, replace
restore
*

*Do you experience depression? (yes/no)
bysort slowerclass: sum depression if bloodpos==1
ranksum depression if bloodpos==1, by(slowerclass)
bysort slowerclass: sum depression if bloodpos==0
ranksum depression if bloodpos==0, by(slowerclass)
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
bysort slowerclass: sum anxiety if bloodpos==1
ranksum anxiety if bloodpos==1, by(slowerclass)
bysort slowerclass: sum anxiety if bloodpos==0
ranksum anxiety if bloodpos==0, by(slowerclass)
*
preserve
keep if bloodpos==0 
tab anxiety
ttest anxiety if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(anxiety) if slowerclass==1
egen B = mean(anxiety) if slowerclass==0
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
bysort slowerclass: sum add if bloodpos==1 
ranksum add if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum add if bloodpos==0  
ranksum add if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest add if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(add) if slowerclass==1
egen B = mean(add) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "ADHD"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save add, replace
restore
*

*What is your current smoking status? (current, recent, former vs. Never + <100 in my lifetime)
bysort slowerclass: sum smoke if bloodpos==1 
ranksum smoke  if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum smoke if bloodpos==0  
ranksum  smoke  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest smoke if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(smoke) if slowerclass==1
egen B = mean(smoke) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Smoker"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save smoke, replace
restore
*

*What is your personal income per year (after tax) in GBP [or CHF]?
bysort slowerclass: sum income if bloodpos==1 
ranksum income if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum income if bloodpos==0  
ranksum income if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum income, detail
replace income = (income/r(max))*100
ttest income if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(income) if slowerclass==1
egen B = mean(income) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Income in % (avg.diff. £2.5k)"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
save income, replace
restore
*



*Economic Related Major (Management, Marketing, Commerce, Finance, Economics, Accounting vs. Everything else)
tab slowerclass  if bloodpos==0
tab slowerclass  if bloodpos==0  & econstudies ==0 
tab econstudies  if bloodpos==0
tab econstudies  if bloodpos==0 & university==0
tab slowerclass  if bloodpos==0 & econstudies==1
tab slowerclass  if bloodpos==0 & econstudies==0

prtest slowerclass  if bloodpos==0, by(econstudies)

bysort slowerclass: sum econstudies if bloodpos==1 
ranksum econstudies if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum econstudies if bloodpos==0  
ranksum econstudies if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest econstudies if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(econstudies) if slowerclass==1
egen B = mean(econstudies) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Econ. major"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save econstudies, replace
restore
*

*Where would you put yourself on the socioeconomic ladder?
*Think of a ladder as representing where people stand in society. At the top of the ladder are the people who are best off—those who have the most money, most education and the best jobs. At the bottom are the people who are worst off—who have the least money, least education and the worst jobs or no job. The higher up you are on this ladder, the closer you are to people at the very top and the lower you are, the closer you are to the bottom. Choose the number whose position best represents where you would be on this ladder. 
bysort slowerclass: sum socioecon if bloodpos==1 
ranksum socioecon if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum socioecon if bloodpos==0  
ranksum socioecon if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum socioecon, detail
replace socioecon = socioecon/r(max)
ttest socioecon if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(socioecon) if slowerclass==1
egen B = mean(socioecon) if slowerclass==0
collapse A B se effectsize
gen diff = B-A
sum diff, detail
gen Paradigm = "SocioEcon status"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save socioecon, replace
restore
*

*Number of sexual partners
bysort slowerclass: sum number if bloodpos==1 
ranksum number if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum number if bloodpos==0  
ranksum number if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum number, detail
replace number = number/r(max)
ttest number if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(number) if slowerclass==1
egen B = mean(number) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Nr. sex partners"
gen diff05 = diff - (1.26 * se)
gen diff95 = diff + (1.26 * se)
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff diff05 diff95
save number, replace
restore
*



*Lefthand (yes/no)
bysort slowerclass: sum lefthand if bloodpos==1 
ranksum lefthand if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum lefthand if bloodpos==0  
ranksum lefthand if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest lefthand if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(lefthand) if slowerclass==1
egen B = mean(lefthand) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Lefthanded"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save lefthand, replace
restore
*

*Have investments in retirement (yes/no)
bysort slowerclass: sum retirement if bloodpos==1 
ranksum retirement if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum retirement if bloodpos==0  
ranksum retirement if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest retirement if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(retirement) if slowerclass==1
egen B = mean(retirement) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Invested for retirement"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save retirement, replace
restore
*

*Have investments (yes/no)
bysort slowerclass: sum investment if bloodpos==1 
ranksum investment if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum investment if bloodpos==0  
ranksum investment if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest investment if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(investment) if slowerclass==1
egen B = mean(investment) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Has investments"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save investment, replace
restore
*


*Access to universal credit
bysort slowerclass: sum univcredit if bloodpos==1 
ranksum univcredit if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum univcredit if bloodpos==0  
ranksum univcredit if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
ttest univcredit if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(univcredit) if slowerclass==1
egen B = mean(univcredit) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Access univ. credit"
gen diff05 = diff - (1.26 * se)
gen diff95 = diff + (1.26 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save univcredit, replace
restore
*

* Risk attitude (Falketal18)
bysort slowerclass: sum risk if bloodpos==1 
ranksum risk if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum risk if bloodpos==0  
ranksum risk if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum risk, detail
replace risk =risk/r(max)
ttest risk if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(risk) if slowerclass==1
egen B = mean(risk) if slowerclass==0
collapse A B se effectsize
gen diff = B-A
sum diff, detail
gen Paradigm = "Risk aversion"
gen diff05 = diff - (1.26 * se)
gen diff95 = diff + (1.26 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save risk, replace
restore
*

*Number of car accidents
bysort slowerclass: sum accidents if bloodpos==1 
ranksum accidents if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum accidents if bloodpos==0  
ranksum accidents if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum accidents, detail
replace accidents = accidents/r(max)
ttest accidents if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(accidents) if slowerclass==1
egen B = mean(accidents) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Nr. accidents"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)

replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff diff05 diff95
save accidents, replace
restore
*

*Number of accidents at work
bysort slowerclass: sum accidentswork if bloodpos==1 
ranksum accidentswork if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum accidentswork if bloodpos==0  
ranksum accidentswork if bloodpos==0   , by(slowerclass)

preserve
keep if bloodpos==0 
sum accidentswork , detail
*replace accidentswork= accidentswork  / r(max)
ttest accidentswork if bloodpos==0, by(slowerclass)
sum 
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(accidentswork) if slowerclass==1
egen B = mean(accidentswork) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Nr. accidents at work"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff diff05 diff95
save accidentswork, replace
restore
*


*CRT score
bysort slowerclass: sum crtscore if bloodpos==1 
ranksum crtscore if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum crtscore if bloodpos==0  
ranksum crtscore if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum crtscore, detail
replace crtscore = crtscore/r(max)
ttest crtscore if bloodpos==0, by(slowerclass)
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
bysort slowerclass: sum stress if bloodpos==1 
ranksum stress if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum stress if bloodpos==0  
ranksum stress if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum stress, detai
replace stress = stress/r(max)
ttest stress if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(stress) if slowerclass==1
egen B = mean(stress) if slowerclass==0
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
bysort slowerclass: sum selfcontrol  if bloodpos==1 
ranksum selfcontrol  if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum selfcontrol  if bloodpos==0  
ranksum selfcontrol  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum selfcontrol
replace selfcontrol = selfcontrol/r(max)
ttest selfcontrol if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(selfcontrol) if slowerclass==1
egen B = mean(selfcontrol) if slowerclass==0
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
bysort slowerclass: sum extraversion  if bloodpos==1 
ranksum extraversion  if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum extraversion  if bloodpos==0  
ranksum extraversion  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum extraversion, detail
*replace extraversion = extraversion/r(max)
ttest extraversion if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(extraversion) if slowerclass==1
egen B = mean(extraversion) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Extraversion"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save extraversion, replace
restore
*

*Agreeable
bysort slowerclass: sum  agreeableness  if bloodpos==1 
ranksum  agreeableness  if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum  agreeableness  if bloodpos==0  
ranksum  agreeableness  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum agreeableness, detail
*replace agreeableness = agreeableness/r(max)
ttest agreeableness if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(agreeableness) if slowerclass==1
egen B = mean(agreeableness) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Agreeableness"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save agreeableness, replace
restore
*


*Conscient
bysort slowerclass: sum conscientiousness  if bloodpos==1 
ranksum  conscientiousness if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum  conscientiousness if bloodpos==0  
ranksum conscientiousness  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum conscientiousness, detail
*replace conscientiousness = conscientiousness/r(max)
ttest conscientiousness if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(conscientiousness) if slowerclass==1
egen B = mean(conscientiousness) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Conscientiousness"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save conscientiousness, replace
restore
*


*Stability
bysort slowerclass: sum stability  if bloodpos==1 
ranksum  stability if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum  stability if bloodpos==0  
ranksum stability  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum stability, detail
replace stability = stability/r(max)
ttest stability if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(stability) if slowerclass==1
egen B = mean(stability) if slowerclass==0
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
bysort slowerclass: sum openness if bloodpos==1 
ranksum openness if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum openness if bloodpos==0  
ranksum openness  if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum openness, detail
*replace openness = openness/r(max)
ttest openness if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(openness) if slowerclass==1
egen B = mean(openness) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Openness"
gen diff05 = diff - (1.96 * se)
gen diff95 = diff + (1.96 * se)
sum diff diff05 diff95
replace diff = diff*100
replace diff05 = diff05*100
replace diff95 = diff95*100
sum diff, detail
save openness, replace
restore
*


*Sensationseeking (Cooketal15)
bysort slowerclass: sum sensationseeking if bloodpos==1 
ranksum sensationseeking if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum sensationseeking if bloodpos==0  
ranksum sensationseeking if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum sensationseeking, detail
replace sensationseeking= sensationseeking/r(max)
ttest sensationseeking if bloodpos==0, by(slowerclass)
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
bysort slowerclass: sum Patience if bloodpos==1 
ranksum Patience if bloodpos==1  , by(slowerclass)
bysort slowerclass: sum Patience if bloodpos==0  
ranksum Patience if bloodpos==0   , by(slowerclass)
*
preserve
keep if bloodpos==0 
sum Patience, detail
replace Patience= Patience/r(max)
ttest Patience if bloodpos==0, by(slowerclass)
gen effectsize =  (r(mu_1)-r(mu_2))/r(sd)
gen se = r(se)
egen A = mean(Patience) if slowerclass==1
egen B = mean(Patience) if slowerclass==0
collapse A B se effectsize
gen diff = A-B
sum diff, detail
gen Paradigm = "Patience"
gen diff05 = diff - (1.26 * se)
gen diff95 = diff + (1.26 * se)
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

append using accidents
append using number


append using alchool
append using selfcontrol
append using Patience
append using risk
append using stability
append using crtscore

append using smoke
append using mental
append using mental2

append using Startup
append using entrepreneur
append using withjob

append using income
append using socioecon
append using univcredit
append using econstudies

meta set diff diff05 diff95, civarlevel(95)  studylabel(Paradigm) 

* Figure 4
*manually reorder variables based on exposition in the paper
sum diff diff05 diff95
  metan diff diff05 diff95 , dp(3) nostats  nointeger nooverall scheme(s1color) nowt null(0) lcols(Paradigm)  graphregion(color(white) fcolor(white)  lcolor(white)) plotregion(fcolor(white) lcolor(white) )   nobox  xscale(range(-20 20)) xlabels(-20(10)20, labsize(small)) force title("Toxoplasmosis-infected - Toxoplasmosis-free") xtitle("Difference between Means (in %)", size(large))
 
		*trick the command in order to get results without having CIs	
replace effectsize = abs(effectsize)
gen effectsize05 = effectsize-0.0001 
gen effectsize95 = effectsize+0.0001

sum effectsize
  metan effectsize effectsize05  effectsize95 ,dp(3) nostats  nointeger nooverall scheme(s1color) nowt null(0) lcols(Paradigm)  graphregion(color(white) fcolor(white)  lcolor(white)) plotregion(fcolor(white) lcolor(white) )   nobox   xscale(range(0 0.5)) xlabels(0(0.1)0.5, labsize(small)) force xtitle("Effect Size", size(large)) title("Toxoplasmosis-infected - Toxoplasmosis-free")
		
	