clear
use Brain_analysis.dta


 * cox:inactive - ref
local covariates1 age_accel sex  //model 1
local covariates2 age_accel sex TDI education new_ethnicity //model 2
local covariates3 age_accel sex i.smoke i.alcohol_frequency diet_score sleep_score  //model 3
local covariates4 age_accel sex BMI hypertension diabetes other_cancer  //model 4
local covariates5 age_accel sex TDI education new_ethnicity i.smoke i.alcohol_frequency diet_score sleep_score BMI hypertension diabetes other_cancer  //model 5

local nn "HES_dementia HES_stroke HES_PD HES_depressive HES_anxiety HES_bd"
local mm "wk_warrior_mvpa"

foreach v of varlist `nn'{
    foreach w of varlist `mm'{
        stset follow_brain_health, failure(`v'==1)
		stcox i.`w' `covariates1'
        est store `w'_1  //model 1
		
        stcox i.`w' `covariates2'
        est store `w'_2  //model 2
		
		stcox i.`w' `covariates3'
        est store `w'_3  //model 3
		
        stcox i.`w' `covariates4'
        est store `w'_4  //model 4
		
		stcox i.`w' `covariates5'
        est store `w'_5  //model 5
		asdoc estat phtest, detail, save(PH_model5) title("`v'_`w'") append 
        
        esttab `w'_1 `w'_2 `w'_3 `w'_4 `w'_5 using Brain_related_inactive_ref.rtf, title("`v'_`w'") eform ci(2) wide nostar b(%9.2f) keep(*.`w') append  //HR with CI
        
        esttab `w'_1 `w'_2 `w'_3 `w'_4 `w'_5 using Brain_related_inactive_ref.rtf, p wide keep(*.`w') append  //P value
    }
}

 * Subgroup analysis: by age and sex 
 ** age in all adjusted model
recode age_accel (65/max = 1) (min/65 = 0), prefix(sub_)  //generate age group by 65 
 
local covariates5 sex TDI education new_ethnicity i.smoke i.alcohol_frequency diet_score sleep_score BMI hypertension diabetes other_cancer  //all adjusted except for age

local nn "HES_dementia HES_stroke HES_PD HES_depressive HES_anxiety HES_bd"
local mm "wk_warrior_mvpa"

foreach v of varlist `nn'{
    foreach w of varlist `mm'{
		stset follow_brain_health, failure(`v'==1)  //data set
		
		stcox i.`w' `covariates5' if sub_age_accel == 0
        est store `w'_65low  //  < 65 years
		
		stcox i.`w' `covariates5' if sub_age_accel == 1
        est store `w'_65high  //  >= 65 years
        
        esttab `w'_65low `w'_65high using Subgroup_age.rtf, title("`v'_`w'") eform ci(2) wide nostar b(%9.2f) keep(*.`w') append  //HR and CI
        
        esttab `w'_65low `w'_65high using Subgroup_age.rtf, p wide keep(*.`w') append  //P value
		
		stcox i.`w'##sub_age_accel `covariates5'
		asdoc testparm `w'#sub_age_accel, save(p_sub_age) title("`v'_`w'") append  //interaction test
    }
}

 ** sex in all adjusted model
local covariates5 age_accel TDI education new_ethnicity i.smoke i.alcohol_frequency diet_score sleep_score BMI hypertension diabetes other_cancer  //all adjusted except for sex

local nn "HES_dementia HES_stroke HES_PD HES_depressive HES_anxiety HES_bd"
local mm "wk_warrior_mvpa"

foreach v of varlist `nn'{
    foreach w of varlist `mm'{
		stset follow_brain_health, failure(`v'==1)  //data set
		
		stcox i.`w' `covariates5' if sex == 0
        est store `w'_female  //  female
		
		stcox i.`w' `covariates5' if sex == 1
        est store `w'_male  //  male
        
        esttab `w'_female `w'_male using Subgroup_sex.rtf, title("`v'_`w'") eform ci(2) wide nostar b(%9.2f) keep(*.`w') append  //HR and CI
        
        esttab `w'_female `w'_male using Subgroup_sex.rtf, p wide keep(*.`w') append  //P value
		
		stcox i.`w'##sex `covariates5'
		asdoc testparm `w'#sex, save(p_sub_sex) title("`v'_`w'") append  //interaction test
    }
} 


 * cox:inactive - 25% 50% 75% threshold
local covariates5 age_accel sex TDI education new_ethnicity i.smoke i.alcohol_frequency diet_score sleep_score BMI hypertension diabetes other_cancer  //model 5

local nn "HES_dementia HES_stroke HES_PD HES_depressive HES_anxiety HES_bd"
local mm "warrior_mvpa_25 warrior_mvpa_50 warrior_mvpa_75"

foreach v of varlist `nn'{
    foreach w of varlist `mm'{
        stset follow_brain_health, failure(`v'==1)
		
		stcox i.`w' `covariates5'
        est store `w'_5  //model 5
        
        esttab `w'_5 using Brain_related_different_threshold.rtf, title("`v'_`w'") eform ci(2) wide nostar b(%9.2f) keep(*.`w') append  //HR and CI
        
        esttab `w'_5 using Brain_related_different_threshold.rtf, p wide keep(*.`w') append  //P value
    }
}
