# Overview of a case-control study
This retrospective study starts with the final health outcome and then traces back to identify potential causes that could have led to this outcome

# Important historical case-control studies
- 1950's: [cigarette smoking and lung cancer](https://www.thelancet.com/pdfs/journals/lanonc/PIIS1470204509704012.pdf)
- 1970's: post-memopausal estrogens and endometrial cancer
- 1980's: Aspirin and Reyes syndrome; Tampon use and toxic shock syndrome; AIDS and sexual practices
- 1990's: vaccine effectiveness; diet and cancer

# Features of a case-control study
- An observational study
- A retrospective study that depends on recalling and recording events that happened before the selection of cases and controls
- Use odds ratio to measure the association between the exposure and disease

# Basic idea
- We hypothesize that the exposure to a factor increases an individual's risk of developing the disease
- We select the individuals with the disease as cases, and the individuals without the disease who are comparable to the cases as controls
- We compare the odds of exposure in cases and controls to derive the odds ratio

# Selection of cases
- Define cases: 1) homogeneity of disease entity; 2) strict diagnostic criteria; 3) prevalent or incident cases
- Source of cases: hospitals; disease registries; medical practices; general population
- We typically take all the cases from the given source at the specified time

# Selection of controls - selection bias
- The contols are from the **same source population** as the cases, such that they would have been identified and included as cases had they developed the disease
- Define controls: 1) free of the disease in question; 2) at risk of the disease; 3) other eligibility requirements as needed
- Source of controls: hospitals; medical practices; population rosters; workforce; coomunity organizations; **relatives of cases**; friends or neighbors of cases
- We typicall sample among potential controls
- Who is the best control?
  - If cases are a random sample of all cases in the population, controls should be a random sample of all non-cases in the population sampled at the same time
  - If study cases are not a random sample of all cases, it is unlikely that a random sample of the population of non-cases will make up a good control population; in other words, controls must be selected so as to **mirror the same biases** that entered into the selection of cases
  - **Comparability** is more important than representativeness in the selection of controls; the control must be at risk of developing the disease; the control should **resemble the case in all respects except for the presence of disease**
  - ! Do not use convenient controls; use controls who are well defined

# Source of data - exposure - recall bias
- Interviews with cases and controls
- Interviews with surrogates (relatives)
- Medical records
- Employment records
- Biologic specimens: **genetic characteristics**; cholesterol; infection; markers of environmental exposures
- Determine the accuracy of the exposure information since it is retrospective

# Study base
- It is composed of a population at risk of exposure over a period of risk of exposure
- Both cases and controls should emerge from the same study base
- If cases are selected exclusively from hospitalized patients, controls must also be selected from hospitalized patients
- If cases must have gone through a certain ascertainment process (e.g., mammogram-detected screening), controls must have also
- If cases must have reached a certain age before they can become cases, so must controls
- If the exposure of interest is cumulative over time, the controls and cases must each have the same opportunity to be exposed to that exposure

# Issues in matching controls
- Identify the pool from which controls may come. This pool is likely to reflect the way cases and controls were ascertained (e.g., hospital, screening test, telephone survey)
- Control selection is typically through matching
  - Define matching variables (e.g., age) and criteria (e.g., control must be within the same 5 year age group) in advance
- Controls can be individually matched or frequency matched
  - Individual matching: select one or more controls who meet the matching criteria for each case
  - Frequency matching: select a population of controls such that the overall characteristics of the group match the overall characteristics of the cases. e.g., if 20% of cases are female, 20% of controls are female, as well.
- Avoid over-matching: match only on factors known to be causes of the disease
- Obtain power: match more than one control per case. Generally, the number of controls should be <= 5, since there is no further gain of power above 4 controls per case
- Obtain generalizability: match more than one type of control

# Issues with matching
- It is difficult to find a good control if there are matching criteria on many variables
- Once the variables are matched, we cannot explore the association of the disease with these matched variables

# Measure the strength of association: odds ratio

| Exposure      | Have disease   | No disease     |
|---------------|----------------|----------------|
| Yes           | A              | B              |
| No            | C              | D              |

- Odds of exposure in the cases = A/C
- Odds of exposure in the controls = B/D
- Odds ratio (OR) = (A/C)/(B/D) = AD/BC
- OR > 1: positive association (increased odds); 0 < OR < 1: negative association (decreased odds); OR = 1.0: no association

# Interpreting the OR
- What does the OR of age (1.02) mean in a multivariable logistic regression model?
  ```
  Patients with one-year older are 2% more likely to have the disease, after controlling for other variables. 
  ```
- What does the OR of age (0.94) mean in a multivariable logistic regression model?
  ```
  Patients with one-year older are 6% less likely to have the disease, after controlling for other variables.
  ```
- We cannot calculate the incidence rate directly in a case-control study, because we do not have the population at risk. Instead, we calculate the odds of having the exposure in the cases and controls
- The OR approximates relative risk (RR) when 1) controls are representative of the general population; 2) cases are representative of all diseased; 3) the **disease is rare**

# Advantages of case-control studies
- Well suited to the study of **rare diseases** or those with long latency
- Allows study of **multiple potential causes** of a disease
- Exisitng records can occasionally be used
- Relatively efficient to organize and conduct
- Relatively inexpensive compared to cohort and clinical trial
- Relatively short time commitment from respondents compared to cohort and clinical trial
- Investigators can obtain answers quickly

# Disadvantages of case-control studies
- Limitations in recall: sometimes, validation of exposure information is difficult or impossible
- Control of confounding variables may not be complete
- Restricted to a single outcome
- Difficult to study **rare exposures**

# Compare case-control and cohort studies
- Case-control study
  - Only estimate the odds ratio
  - Potentially weaker causal investigation
  - Inexpensive
  - Short-term study
  - Can be powerful with small sample of cases
  - Efficient design for rare disease
  - Good for multiple exposures
  - Suffer recall bias
  - Less likely to suffer loss of follow-up
  - The results may not be generalizable
  - Cannot examine the natural course of disease, survival
- Cohort study
  - Can calculate the incidence rate, risk, and relative risk
  - Potentially greater causal investigation
  - Expensive
  - Long-term study
  - Large sample size required
  - Efficient design for rare exposure
  - Good for multiple outcomes
  - Less likely to suffer recall bias
  - More likely to suffer loss of follow-up
  - The results may be generalizable
  - Can examine the natural course of disease, survival
