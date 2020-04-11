*-----------------------------------------------------------------------------
*! v 2.4   	10apr2020				  			by JPA     		
*	lnsd: fixed
*	mz 	: multiple poverty lines
*   mmu	: multiple mean values
* v 2.3.1   08apr2020				  			by JPA     		
*   add SD was an option when estimating groupped data
*	Remove PW since it is not supported by SUMARIZE
*   Type 1 grouped data: P=Cumulative proportion of population, L=Cumulative 
*		proportion of income held by that proportion of the population
*   Type 2 grouped data: Q=Proportion of population, R=Proportion of incometype 
*   Type 5 grouped data: W=Percentage of the population in a given interval of 
*		incomes, X=The mean income of that interval.
*   Type 6 grouped data: W=Percentage of the population in a given interval of 
*		incomes, X=The max income of that interval.
*   Unit record data: Percentage of the population with same income level, 
*		The income level.
*		improve the layout
* v 2.2   06apr2020				  				by JPA     		
*   dependencies checks run quietly
*   apoverty and ainequal added to the dependencies check
* v 2.1   05apr2020				  				by JPA     		
*   changed ado name from grouppov to groupdata
* v 2.0   02apr2020				  				by JPA     		
*   changes made to use this method to estimate learning poverty 
* 	add support to aweight
*   replace wtile2 by alorenz
*   add microdata value as benchmark
* v 1.1   14jan2014				  				by SM and JPA
*   change ado name from povcal to grouppov
*   technical note on Global Poverty Estimation: Theoratical and Empirical 
*   Validity of Parametric Lorenz Curve Estiamtes and Revisitng Non-parametric 
*   techniques. (January, 2014), for discussions on the World Bank Global 
*   Poverty Monitoring Working Group.
* v 1.0   02fev2012				  				by SM and JPA 			
*   povcal.ado created by Joao Pedro Azevedo (JPA) and Shabana Mitra (SM)
*-----------------------------------------------------------------------------

program define groupdata, rclass

    version 8.0

    syntax varlist(numeric min=1 max=1)         ///
                 [in] [if]                      ///
                 [fweight aweight pweight]      ///
                 ,                              ///
                         z(real)                ///
                    [							///
						 Mu(real -99)          	///
                         GROUPed                ///
                         BINs(real -99)         ///
                         REGress				///
						 BENCHmark				///
						 NOFIGures				///
						 UNITRECord				///
						 type(string) 			///
						 NOElasticities			///
						 NOLorenz				///
						 NOChecks				///
						 min(string)			///
						 max(string)			///
						 sd(real -99)			///
						 debug					///
					]

quietly {

preserve 
		
*-----------------------------------------------------------------------------
* 	Temp names 
*-----------------------------------------------------------------------------
		
		tempname A  gq cofb cof  gqg cofbg tmp
	
		tempvar  temp touse rnd lninc lnmpce mpce pp pw L p y1 y2 a b c  x1 x2  Lg pg yg ag bg cg yg2 x1g x2g  type2 model var value

*-----------------------------------------------------------------------------
* 	Weights
*-----------------------------------------------------------------------------
		
		* keep original weights
		local wtg2 = "`weight'"
		local weight2 = "`weight'"
		local exp2 = subinstr("`exp'","=","",.)

		* set-up weights when it is not available
		if ("`weight'" == "") {
			tempvar wtg
			gen `wtg' = 1
			loc weight "fw"
			loc exp    "=`wtg'"
			local pop "`wtg'"
			local `mpce' = `pop'
		}
		else {
			local pop =subinstr("`exp'","=","",.)
			local `mpce' = `pop'
}
		
*-----------------------------------------------------------------------------
* Check options
*-----------------------------------------------------------------------------

	if (`mu' != -99) & ("`benchmark'" != "") {
		noi di ""
        di as err "Option benchmark only works with original mean."
        exit 198
		noi di ""
	}
	  
	if (`mu' == -99) & ("`type'" != "") {
		noi di ""
        di as err "Estimates based on group data require the user to provide the mean value of the distribution"
        exit 198
		noi di ""
	}
	if (strmatch(" 1 2 5 6","*`type'*") == 0) {
		noi di ""
        di as err "Please select a valid data type. see help."
		exit 198
		noi di ""
		noi di as text "Type 1 grouped data: " as res "P=Cumulative proportion of population, L=Cumulative proportion of income held by that proportion of the population" 
		noi di as text "Type 2 grouped data: " as res "Q=Proportion of population, R=Proportion of incometype" 			
		noi di as text "Type 5 grouped data: " as res "W=Percentage of the population in a given interval of incomes, X=The mean income of that interval"
		noi di as text "Type 6 grouped data: " as res "W=Percentage of the population in a given interval of incomes, X=The max income of that interval"
		noi di ""
		
	}

	if ("`type'" == "1") {
					
		if ("`wtg2'" == "") {
			noi di as err "Type 1 only accepts accept AW."
			exit 198
		}

		if strmatch("pweight","*`wtg2'*") == 1 { 
			noi di as err "Type 1 only accepts accept AW."
			exit 198
		}
			
		if strmatch("fweight","*`wtg2'*") == 1 { 
			noi di as err "Type 1 only accepts accept AW."
			exit 198
		}
	}
		
	if ("`type'" == "2") {

		if ("`wtg2'" == "") {
			noi di as err "Type 2 only accepts accept AW."
			exit 198
		}

		if strmatch("pweight","*`wtg2'*") == 1 { 
			noi di as err "Type 2 only accepts accept AW."
			exit 198
		}
				
		if strmatch("fweight","*`wtg2'*") == 1 { 
			noi di as err "Type 2 only accepts accept AW."
			exit 198
		}
	
	}

	if ("`type'" == "5") {

		if (substr(trim("`wtg2'"),1,2) == "aw") { 
			noi di as err "Type 5 does not accept AW weights. Please use either PW, FW or no weights."
			exit 198
		}
		if strmatch("pweight","*`wtg2'*") == 1 { 
			local weight2 = "aweight"
		}

	}

	if ("`type'" == "6") {

		if (substr(trim("`wtg2'"),1,2) == "aw") { 
			noi di as err "Type 6 does not accept AW weights. Please use either PW, FW or no weights."
			exit 198
		}
		if strmatch("pweight","*`wtg2'*") == 1 { 
			local weight2 = "aweight"
		}

	}	
	
*-----------------------------------------------------------------------------
* Download and install required user written ado's
*-----------------------------------------------------------------------------
* Fill this list will all user-written commands this project requires

		  local user_commands groupfunction alorenz which_version apoverty ainequal estout

* Loop over all the commands to test if they are already installed, if not, then install
		  qui foreach command of local user_commands {
			cap which `command'
			if _rc == 111 { 
				ssc install `command'
			}
			else {
				which_version groupfunction 
				if  (`s(version)' < 2.0) {
					ado update groupfunction , update
				}
				which_version alorenz
				if  (`s(version)' < 3.1) {
					ado update alorenz , update
				}
			}
		  }

*-----------------------------------------------------------------------------
* 	Display Options 
*-----------------------------------------------------------------------------

	* show regression outputs
		if ("`regress'" != "") {
			loc noireg "noi "
		}
		else {
			loc noireg ""
		}
	
	* does not show lorenz
		if ("`nolorenz'" != "") {
			loc noilor ""
		}
		else {
			loc noilor "noi"
		}
	
	* does not show elasticities
		if ("`noelasticities'" != "") {
			loc noelast ""
		}
		else {
			loc noelast "noi"
		}
		
	* does not show checks
		if ("`nochecks'" != "") {
			loc nocheck ""
		}
		else {
			loc nocheck "noi"
		}

	* debug
		if ("`debug'" == "") {
			loc noidebug ""
		}
		else {
			loc noidebug "noi"
		}

*-----------------------------------------------------------------------------
* 	Filters
*-----------------------------------------------------------------------------

		tokenize `varlist'
		local inc `1'

		mark `touse' `if' `in' [`weight'`exp']
		
		markout `touse' `varlist'

*-----------------------------------------------------------------------------
* 	Data sort
*-----------------------------------------------------------------------------

		set seed 1234568
		
        gen double `rnd' = uniform()		if `touse'
        gen `lninc' 	= ln(`inc') 			if `touse'
		gen `lnmpce' 	= ln(`inc') 			if `touse'
        sort `inc' `rnd'

*-----------------------------------------------------------------------------
* 	Mean values
*-----------------------------------------------------------------------------

	* generate mean values if unit record data is provided
	if ("`type'" == "") {
		
	    if (`mu' == -99) {
            * generate mean and stadard deviation for unit record data
			sum `inc' [`weight2'`exp']			if `touse'
            local mu = `r(mean)'

            sum `lnmpce' [`weight2'`exp']		if `touse'
            local lnmu = r(mean)
            local lnsd = r(sd)
        }
	
        if (`mu') != -99 {
			* use the mean provided as an option 
			sum `lnmpce' [`weight2'`exp']		if `touse'
            local lnmu = ln(`mu')
            local lnsd = r(sd)
        }
	}
	
	if ("`type'" != "") {
            
			local lnmu = ln(`mu')
            
			if (`sd' == -99) {
				sum `lnmpce' [`weight2'`exp']	if `touse'
				local lnmu = ln(`mu')
				local lnsd = r(sd)
*				local lnsd = ln(.5)
			}
			
			else {
				local lnsd = ln(`sd')
			}
			
			keep if `touse'

			local N = _N
			
			if (`N' < 24) {
				set obs 24
			}
	}
	

*-----------------------------------------------------------------------------
* 	Unit Record 
*-----------------------------------------------------------------------------
		
		if ("`benchmark'" == "benchmark") {
			
			* unit record poverty estimates
			apoverty `inc' [`weight'`exp'] 	if `touse', line(`z')  fgt3  pgr
			local afgt0 = r(head_1) 
			local afgt1 = r(pogapr_1) 
			local afgt2 = r(fogto3_1) 
			
			* unit record inequality estimates
			ainequal `inc' [`weight'`exp']  if `touse'
			local agini = r(gini_1) 

			* create row labels for output matrix
			local rownames_unitrecord " fgt0 fgt1 fgt2 gini "
		
		}

*-----------------------------------------------------------------------------
* 	Unit Record
*-----------------------------------------------------------------------------

		
        qui if ("`unitrecord'" == "unitrecord") {
			
			noi di ""
			noi di ""
			noi di "Estimation using unit record information."
		
			
            ************************************
            ** cumulative distribution
            ************************************

            egen double `pw' = pc(`inc')			if `touse', prop
            egen double `pp' = pc(`pop')			if `touse', prop

            gen double 	`L' = `pw'					if `touse'
            replace 	`L' = `pw'+`L'[_n-1] in 2/l	if `touse'

            gen double 	`p' = `pp'					if `touse'
            replace 	`p' = `pp'+`p'[_n-1] in 2/l	if `touse'

            ************************************
            ** generate variables (GQ Lorenz Curve)
            ************************************

            gen double `y1' = `L'*(1-`L')			if `touse'
            gen double `a' 	= ((`p'^2)-`L')			if `touse'
            gen double `b' 	= `L'*(`p'-1)			if `touse'
            gen double `c' 	= (`p'-`L')				if `touse'

            ************************************
            ** generate variables Beta Lorenz Curve
            ************************************

            gen double `y2'	=	ln(`p'-`L')			if `touse'
            gen double `x1'	=	ln(`p')				if `touse'
            gen double `x2'	=	ln(1-`p')			if `touse'

            local last = _N-1

            ************************************
			** Plot Figure 
            ************************************

			if ("`nofigures'" == "") {
			    
				local mustr = strofreal(`mu',"%9.2f")
				local intercept00 = _N + 1
				replace `L' = 0 in `intercept00'
				replace `p' = 0 in `intercept00'
				
				graph twoway lowess `L' `p'		if `touse', 						///
					ytitle("Lorenz") xtitle("Population (Cumulative)") 				///
					note("mean: `mustr' [`bins' bins]") name(lorenz, replace)
								
				kdensity `inc' 					if `touse', 						///
					xline(`z') xtitle("`inc'") name(pdf, replace)

				graph twoway lowess `inc' `p'	if `touse', 						///
					yline(`z') ytitle("`inc'") xtitle("Population (Cumulative)") 	///
					note("mean: `mustr' [`bins' bins]") name("pen", replace)
	
			}

            ************************************
            ** Estimation: GQ Lorenz Curve
            ************************************

			label var `y1' 	"`inc'" 
			label var `a'  	"A"
			label var `b' 	"B"
			label var `c'	"C"

            qui reg `y1' `a' `b' `c' in 1/`last' if `touse', noconstant
            est store coefgq
            mat `gq' = e(b)
            mat `cof' = e(b)

            ************************************
            ** Estimation: Beta Lorenz Curve
            ************************************

			label var `y2' 	"`inc'"
			label var `x1'	"B"
			label var `x2'	"C"
			
            qui reg `y2' `x1' `x2' in 1/`last' if `touse'
            est store coefbeta
            mat `cofb' = e(b)

        }
		
*-----------------------------------------------------------------------------
* 	Group data (provided)
*-----------------------------------------------------------------------------

		
        qui if ("`type'" != "") {
			
			noi di ""
			noi di ""
			noi di "Estimation using provided distribution (Type `type')"

			if ("`type'" == "1") {
				noi di ""
				noi di "Type 1 grouped data: P=Cumulative proportion of population, L=Cumulative proportion of income held by that proportion of the population" 
			}
			if ("`type'" == "2") {
				noi di ""
				noi di "Type 2 grouped data: Q=Proportion of population, R=Proportion of incometype" 			
			}
			if ("`type'" == "5") {
				noi di ""
				noi di "Type 5 grouped data: W=Percentage of the population in a given interval of incomes, X=The mean income of that interval"
			}
			if ("`type'" == "6") {
				noi di "Type 6 grouped data: W=Percentage of the population in a given interval of incomes, X=The max income of that interval"
			}

			************************************
            ** Type 1 grouped data: P=Cumulative proportion of population, L=Cumulative proportion of income held by that proportion of the population
            ************************************
			
			if ("`type'" == "1") {
	
				sum `inc'									if `touse'
				local bins = r(N)
				local last = `bins'-1
				
				
				if (substr(trim("`wtg2'"),1,2) == "aw") { 
					
					gen `Lg' = `inc'/100					if `touse'
					gen `pg' = `exp2'/100					if `touse'
				
				}
			}
						
			* noi list `pg' mean_score_lp `inc' `delta' `inc2' `Lg' in 1/`bins'
			
			 
			************************************
            ** Type 2 grouped data: Q=Proportion of population, R=Proportion of incometype
            ************************************
			
			if ("`type'" == "2") {
					
				sum `inc'										if `touse'
				local bins = r(N)
				local last = `bins'-1
				
				if (substr(trim("`wtg2'"),1,2) == "aw") { 

					gen `Lg' = `inc'/100						if `touse'
					replace `Lg' = `Lg'[_n]+`Lg'[_n-1] in 2/l	if `touse'	
					
					gen `pg' = `exp2'/100						if `touse'
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l	if `touse'	
					
				}
			}
						
			* noi list `pg' mean_score_lp `inc' `delta' `inc2' `Lg' in 1/`bins'
			 
            ************************************
            ** Type 5 grouped data: W=Percentage of the population in a given interval of incomes, X=The mean income of that interval
            ************************************
			
			if ("`type'" == "5") {

			*set trace on

				sum `inc'											if `touse'
				local bins = r(N)
				local last = `bins'-1
				
				if ("`wtg2'" == "") {

					gen double 	`pg' 	= 	1/`bins'				if `touse'
					replace 	`pg' = `pg'[_n]+`pg'[_n-1] in 2/l	if `touse'	

					sum `inc'										if `touse'
					local sumL = r(sum)							
					gen double 	`Lg' = `inc'/`sumL'					if `touse'
					replace `Lg' = `Lg'[_n]+`Lg'[_n-1] in 2/l		if `touse'
				
					`noidebug' list `pg' `inc' `Lg' `LLLLL' `PPPPP'

				}
			
			*set trace off

				
				
				if (substr(trim("`wtg2'"),1,2) == "pw") { 

					tempvar LLLLL PPPPP
				
					gen 	`pg' = `exp2'						
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l	
					
					gen double `PPPPP' = `exp2'*`inc'*100000
					sum `PPPPP' 						
					local sumL = r(sum)
					gen double 	`LLLLL' = `PPPPP'/`sumL'
					replace 	`LLLLL' = `LLLLL'[_n]+`LLLLL'[_n-1] in 2/l	
					
					gen double `Lg' = `LLLLL'
					
					`noidebug' list `pg' `inc' `Lg' `LLLLL' `PPPPP'

				}
				
				
				if (substr(trim("`wtg2'"),1,2) == "fw") { 
					
					sum `exp2'									if `touse'
					local sumP = r(sum)
					gen double 	`pg' = `exp2'/`sumP'			if `touse'
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l	if `touse'

					sum `inc' 	[`weight'`exp']					if `touse'
					local sumL = r(sum)
					gen doulbe 	`Lg' = `inc'/`sumL'				if `touse'
					replace `Lg' = `Lg'[_n]+`Lg'[_n-1] in 2/l	if `touse'
				
				}
				
			}
						
            ************************************
            ** Type 6 grouped data: W=Percentage of the population in a given interval of incomes, X=The max income of that interval
            ************************************
			
			tempvar inc2 delta
			
			if ("`type'" == "6") {

				sum `inc'															if `touse'
				local bins = r(N)
				local bins = `bins'+1
				local last = `bins'-1

				noi di "min: " `min'
				noi di "max: " `max'
				noi di "bins: " `bins'
					
				if ("`wtg2'" == "") {

					gen double 	`pg' 	= 	1/`bins'								if `touse'

					gen 	double `delta' = .
					replace `delta' = (`inc'[_n]-`min')			/2 	in 1			if `touse'
					replace `delta' = (`inc'[_n]-`inc'[_n-1])	/2 	in 2/`last'		if `touse'
					replace `delta' = (`max'	-`inc'[_n-1])	/2 	in `bins'		if `touse'
					
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l						if `touse'
					
					gen double `inc2' = .											if `touse'
					replace `inc2' = `inc' - `delta'				in 1/`last'		if `touse'
					replace `inc2' = `max' - `delta'				in `bins'		if `touse'
					
					sum `inc2'														if `touse'
					local sumL = r(sum)
					gen double 	`Lg' = `inc2'/`sumL'								if `touse'
					replace `Lg' = `Lg'[_n]+`Lg'[_n-1] in 2/l						if `touse'
					
					`noidebug' list `pg' `inc' `delta' `inc2' `Lg' in 1/`bins'
				
				}

				if (substr(trim("`wtg2'"),1,2) == "pw") { 
	
					tempvar LLLLL PPPPP
	
					gen double 	`pg' 	= 	`exp2'									if `touse'

					gen 	double `delta' = .										if `touse'
					replace `delta' = (`inc'[_n]-`min')			/2 	in 1			if `touse'
					replace `delta' = (`inc'[_n]-`inc'[_n-1])	/2 	in 2/`last'		if `touse'
					replace `delta' = (`max'	-`inc'[_n-1])	/2 	in `bins'		if `touse'
					
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l						if `touse'
					
					gen double `inc2' = .											if `touse'
					replace `inc2' = `inc' - `delta'				in 1/`last'		if `touse'
					replace `inc2' = `max' - `delta'				in `bins'		if `touse'

					
					gen double `PPPPP' = `exp2'*`inc2'*100000
					sum `PPPPP' 						
					local sumL = r(sum)
					gen double 	`LLLLL' = `PPPPP'/`sumL'
					replace 	`LLLLL' = `LLLLL'[_n]+`LLLLL'[_n-1] in 2/`bins'	
					
					gen double `Lg' = `LLLLL'
	
					`noidebug' list `pg' `inc' `delta' `inc2' `Lg' in 1/`bins'
					
				}
				
				if (substr(trim("`wtg2'"),1,2) == "fw") { 
	
					sum `exp2'														if `touse'
					local sumP = r(sum)
					gen double 	`pg' = `exp2'/`sumP'								if `touse'

					gen 	double `delta' = .										if `touse'
					replace `delta' = (`inc'[_n]-`min')			/2 	in 1			if `touse'
					replace `delta' = (`inc'[_n]-`inc'[_n-1])	/2 	in 2/`last'		if `touse'
					replace `delta' = (`max'	-`inc'[_n-1])	/2 	in `bins'		if `touse'
					
					replace `pg' = `pg'[_n]+`pg'[_n-1] in 2/l						if `touse'
					
					gen double `inc2' = .											if `touse'
					replace `inc2' = `inc' - `delta'				in 1/`last'		if `touse'
					replace `inc2' = `max' - `delta'				in `bins'		if `touse'
					
					sum `inc2'		[`weight'`exp']									if `touse'
					local sumL = r(sum)
					gen double 	`Lg' = `inc2'/`sumL'								if `touse'
					replace `Lg' = `Lg'[_n]+`Lg'[_n-1] in 2/l						if `touse'
					
					`noidebug' list `pg' `inc' `delta' `inc2' `Lg' in 1/`bins'
					
				}

			}
			
            ************************************
            ** cumulative distribution
            ************************************
/*
            egen double `pw' = pc(`inc')			if `touse', prop
            egen double `pp' = pc(`pop')			if `touse', prop

            gen double `L' = `pw'					if `touse'
            replace `L' = `pw'+`L'[_n-1] in 2/l		if `touse'

            gen double `p' = `pp'					if `touse'
            replace `p' = `pp'+`p'[_n-1] in 2/l		if `touse'
*/

			sum `pg' 
			local s = r(sum)
			if (`s'>99) {
				gen double `p' = `pg'/100					
			} 
			else {
				gen double `p' = `pg'					
			}

			sum `Lg' 
			local s = r(sum)
			if (`s'>99) {
				gen double `L' = `Lg'/100					
			} 
			else {
				gen double `L' = `Lg'					
			}

			
            ************************************
            ** generate variables (GQ Lorenz Curve)
            ************************************

            gen double `y1' = `L'*(1-`L')			
            gen double `a' = ((`p'^2)-`L')			
            gen double `b' = `L'*(`p'-1)			
            gen double `c' = (`p'-`L')				

            ************************************
            ** generate variables Beta Lorenz Curve
            ************************************

            gen double `y2'=ln(`p'-`L')				
            gen double `x1'=ln(`p')					
            gen double `x2'=ln(1-`p')				

            ************************************
			** Plot Figure 
            ************************************

			if ("`nofigures'" == "") {
			    
				local mustr = strofreal(`mu',"%9.2f")
				local intercept00 = `bins' + 1
				replace `L' = 0 in `intercept00'
				replace `p' = 0 in `intercept00'
				
				graph twoway lowess `L' `p'		, 						///
					ytitle("Lorenz") xtitle("Population (Cumulative)") 				///
					note("mean: `mustr' [`bins' bins]") name(lorenz, replace)
								
				kdensity `inc' 					, 						///
					xline(`z') xtitle("`inc'") name(pdf, replace)

				graph twoway lowess `inc' `p'	, 						///
					yline(`z') ytitle("`inc'") xtitle("Population (Cumulative)") 	///
					note("mean: `mustr' [`bins' bins]") name("pen", replace)
	
			}

			
			`noidebug' list `Lg' `L' `pg' `p'  `y1' `a' `b' `c'  `y2' `x1' `x2' if `y1' !=.
			
            ************************************
            ** Estimation: GQ Lorenz Curve
            ************************************

			label var `y1' 	"`inc'" 
			label var `a'  	"A"
			label var `b' 	"B"
			label var `c'	"C"

            qui reg `y1' `a' `b' `c' in 1/`last' , noconstant
            est store coefgq
            mat `gq' = e(b)
            mat `cof' = e(b)

            ************************************
            ** Estimation: Beta Lorenz Curve
            ************************************

			label var `y2' 	"`inc'"
			label var `x1'	"B"
			label var `x2'	"C"
			
            qui reg `y2' `x1' `x2' in 1/`last' 
            est store coefbeta
            mat `cofb' = e(b)

        }

*-----------------------------------------------------------------------------
* 	Group data (estimated)
*-----------------------------------------------------------------------------
		
        qui if ("`grouped'" == "grouped") {
		
			noi di ""
			noi di "Estimation using grouped data..."

            ** cumulative distribution (grouped data)
			
			if (`bins'!= 0) {
				
				if (`bins' == -99) {
					local bins = 20
					noi di ""
					noi di "... creating groupped data with `bins' bins."
					noi di ""
					alorenz `inc' [`weight'`exp']	if `touse', points(`bins') 
				}
				else {
					noi di ""
					noi di "... creating groupped data with `bins' bins."
					noi di ""
					alorenz `inc' [`weight'`exp']	if `touse', points(`bins') 
				}
			}
			
			mat `A' = r(lorenz1)
			
			svmat double `A'
			
			return matrix data = `A'

            ** generate variables (GQ Lorenz Curve) (grouped data)

            gen double `Lg' = `A'2/100
            gen double `pg' = `A'4/100

            gen double `yg' = `Lg'*(1-`Lg')
            gen double `ag' = ((`pg'^2)-`Lg')
            gen double `bg' = `Lg'*(`pg'-1)
            gen double `cg' = (`pg'-`Lg')

            ** generate variables Beta Lorenz Curve (Grouped data)

            gen double `yg2' = ln(`pg'-`Lg')
            gen double `x1g' = ln(`pg')
            gen double `x2g' = ln(1-`pg')

            local lastg = `bins'-1

            ************************************
			** Plot Figure 
            ************************************

			if ("`nofigures'" == "") {
				
				local mustr = strofreal(`mu',"%9.2f")
				local intercept00 = `bins'+1
				replace `Lg' = 0 in `intercept00'
				replace `pg' = 0 in `intercept00'
				
				graph twoway lowess `Lg' `pg', 										///
					ytitle("Lorenz") xtitle("Population (Cumulative)") 				///
					note("mean: `mustr' [`bins' bins]") name(lorenz, replace)
								
				kdensity `A'6, 														///
					xline(`z') xtitle("`inc'") name(pdf, replace)

				graph twoway lowess `A'3 `pg', 										///
					yline(`z') ytitle("`inc'") xtitle("Population (Cumulative)") 	///
					note("mean: `mustr' [`bins' bins]") name("pen", replace)
	
			}

            ************************************
            ** Estimation: GQ Lorenz Curve (grouped data)
            ************************************
			
			label variable `yg'     "`inc'"
			label variable `ag'		"A" 
			label variable `bg' 	"B"
			label variable `cg' 	"C"
						
            qui reg `yg' `ag' `bg' `cg' in 1/`lastg', noconstant
            est store coefgqg
            mat `gqg' = e(b)

            ************************************
            ** Estimation: Beta Lorenz Curve (Grouped data)
            ************************************

			label variable `yg2' 	"`inc'"
			label variable `x1g'  	"B"
			label variable `x2g'	"C"
			
            qui reg `yg2' `x1g' `x2g'  in 1/`lastg'
            est store coefbetag
            mat `cofbg' = e(b)

        }
		
        /**************************************
        ** Test
        **************************************

        save tmp0 , replace

        preserve
        use y a b c using tmp0 , clear
        save tmp1, replace
        use yg ag bg cg using tmp0 , clear
        rename yg y
        rename ag a
        rename bg b
        rename cg c
        gen group = 2
        drop if y == .
        save tmp2, replace
        use tmp1, clear
        gen group = 1
        append using tmp2
        tab group

        ** GQ Lorenz Curve

        regress y a b c  if group == 1, noconstant
        est store gq

        ** GQ Lorenz Curve (grouped data)

        regress y a b c   if group == 2, noconstant
        est store gqg

        suest gq gqg
        test [gq_mean = gqg_mean]

        restore
*/

*-----------------------------------------------------------------------------
* Choice of the Lorenz curve                   
*-----------------------------------------------------------------------------
/*
        ******test stat for GQ Lorenz******

        estimates restore gqg
        predict lhatg if _est_gqg==1
        gen Lhatg=lhatg
        drop lhatg
        if _est_gqg==1 {
            local k_t=(pg<`H')
        }
        else {
        local k_t=.
        }
        if `k_t'==1{
        local testGQ=(`Lhatg'-yg)^2
        }
        else {
        local testGQ=.
        }
        sum `testGQ'
        local GQ_TT=r(sum)

        *******test stat for Beta Lorenz**************
        estimates restore gqg
        predict lhatgBeta if  _est_blcg==1
        local LhatgBeta=lhatg
        drop lhatgBeta
         if  _est_blcg==1 {
        local k_tBeta=(pg<`hcrb')
        }
        else {
        local k_tBeta=.
        }
        if `k_tBeta'==1 {
        local testBeta=(`LhatgBeta'-yg)^2
        }
        else {
        local testBeta=.
        }

        sum `testBeta'
        local Beta_TT=r(sum)

        local test=(`GQ_TT'<`Beta_TT')

        if `test==1{
            noi di as res "GQ has a lower statistic"
        }
        else {
            noi di as res "Beta has a lower statistic"
        }
*/

        **************************************
        ** call sub routine
        **************************************

	_subroutine_groupdata.ado
	
		
restore
	
end
