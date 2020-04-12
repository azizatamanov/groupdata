*-----------------------------------------------------------------------------
*! v 2.4   	10apr2020				  			by JPA     		
*	lnsd: fixed
*-----------------------------------------------------------------------------

program define _subroutine_groupdata, rclass

    version 8.0

    syntax          							///
                 ,                              ///
                         z(real)                ///
						 Mu(real -99)          	///
						 lnsd(real)				///
						 lnmu(real)				///
						 a(real)				///
						 b(real)				///
						 c(real)				///
 						 aatheta(real) 			///
						 aagama(real)			///
						 aagama2(real)			///
						 aadelta(real)			///
						 aadelta2(real)			///
                    [							///
                         BINs(real -99)         ///
						 BENCHmark				///
						 type(string) 			///
						 NOElasticities			///
						 NOLorenz				///
						 NOChecks				///
						 min(string)			///
						 max(string)			///
						 sd(real -99)			///
					]


*-----------------------------------------------------------------------------
* 	Temp names 
*-----------------------------------------------------------------------------
		
		tempname A  gq cofb cof  gqg cofbg tmp
	
		tempvar  temp touse rnd lninc lnmpce mpce pp pw L p y1 y2 a b c  x1 x2  Lg pg yg ag bg cg yg2 x1g x2g  type2 model var value

*-----------------------------------------------------------------------------
* 	Table 2 (Datt, 1998)             
*-----------------------------------------------------------------------------
* 		GQ Lorenz Curve 
*-----------------------------------------------------------------------------

        if ("`grouped'" != "") {
            local a = `gqg'[1,1]
            local b = `gqg'[1,2]
            local c = `gqg'[1,3]
        }
        else {
            local a = `gq'[1,1]
            local b = `gq'[1,2]
            local c = `gq'[1,3]
        }

        local e     = -(`a'+`b'+`c'+1)
        local m     = `b'*`b' - (4*`a')
        local n     = (2*`b'*`e') - (4*`c')
        local r     = sqrt((`n'*`n') - (4*`m'*(`e'*`e')))
        local s1    = (`r'-`n')/(2*`m')
        local s2    = -(`r'+`n')/(2*`m')

        local H = -(1/(2*`m'))*(`n'+`r'*(`b'+2*`z'/`mu')*(1/sqrt((`b'+2*`z'/`mu')*(`b'+2*`z'/`mu')-`m')))
        local lH = -(1/2)*(`b'*`H' + `e' + sqrt(`m'*`H'*`H' + `n'*`H' + `e'*`e'))
        local PG = `H'-(`mu'/`z')*`lH'
        local SPG = 2*`PG' - `H' - ((`mu'/`z')*(`mu'/`z')) * (`a'*`H' + `b'*`lH' - (`r'/16) * ln((1-`H'/`s1')/(1-`H'/`s2')))

        /*** Second derivative of Lorenz curve*/
        local ldph = ((`r'*`r')/8)*((`m'*(`H'*`H')+ `n'*`H' + `e'*`e')^(-3/2))

        /*** Gini */

        #delim ;
        if `m'<0 { ;	
        	local gini_tt = (`e'/2)- `n'*(`b'+2)/(4*`m') +
        	((`r'^2)/(8*`m'*sqrt(-`m')))* 	(asin((2*`m'+`n')/`r') - asin(`n'/`r')) ;	
        };
        else {;
        	local gini_tt = (`e'/2)- `n'*(`b'+2)/(4*`m') - ((`r'^2)/(8*`m'*sqrt(`m')))
        	*ln(abs(((2*`m')+`n'+(2*(sqrt(`m'))*(`a'+`c'-1)))/(`n'-(2*`e'*sqrt(`m')))));
        };
        #delim cr

        /*** Gini */

        local gini_ln = (2*normal(`lnsd'/sqrt(2))) - 1
        local lnsd_tt = sqrt(2)*invnormal((`gini_tt'+1)/2)
        local dirsigma = normal((1/`lnsd')*ln(`z'/`mu')+(`lnsd'/2))
        local ginisigma = normal((1/`lnsd_tt')*ln(`z'/`mu')+(`lnsd_tt'/2))
        local pglnd = `dirsigma' *((`z' - (normal((1/`lnsd')*ln(`z'/`mu')-(`lnsd'/2))* `mu' )/`dirsigma')/`z')
        local pglng = `ginisigma' * ((`z' - (normal((1/`lnsd_tt')*ln(`z'/`mu')-(`lnsd_tt'/2))* `mu' )/`ginisigma')/`z')
        local spglnd = `pglnd'*`pglnd'/`dirsigma'
        local spglng = `pglng'*`pglng'/`ginisigma'

        if "`gini_G'" !=""{
        	local gini_gg=`gini_G'/100
        	local lnsd_gg = sqrt(2)*invnormal((`gini_gg'+1)/2)
        	local HcGg = normal((1/`lnsd_gg')*ln(`z'/`mu')+(`lnsd_gg'/2))
        	local PgGg = (`z' - (normal((1/`lnsd_gg')*ln(`z'/`mu')-(`lnsd_gg'/2))* `mu' )/`HcGg')/`z'
        	local SPgGg = `PgGg'*`PgGg'
        }

        if "`nsmean'" != "" {   		
        	local tem3 = (1/`lnsd')*(ln(`z'/`nmu'))+(`lnsd'/2)  		
        	local tem4 = (1/`lnsd1')*(ln(`z'/`nmu'))+(`lnsd1'/2)  		
        	local dis12 = normal(`tem3')  		
        	local dis13 = normal(`tem4')  		
        }  	

        if "`npovline'" != "" & "`nsmean'" != "" {  		
        	local tem5 = (1/`lnsd')*(ln(`n'*`z'/`nmu'))+(lnsd/2)  		
        	local tem6 = (1/`lnsd1')*(ln(`n'*`z'/`nmu'))+(lnsd1/2)  		
        	local dis14 = normal(`tem5')  		
        	local dis15 = normal(`tem6')  		
        }

        /*** Elasticities QG Lorenz */

        local elhmu         = -(`z'/(`mu'*`H'*`ldph'))
        local elhgini       = (1-(`z'/`mu'))/ (`H'*`ldph')
        local elpgmu        = 1-(`H'/`PG')
        local elpggini      = 1+(((`mu'/`z')-1)* (`H'/`PG'))
        local elspgmu       = 2*(1-`PG'/`SPG')
        local elspggini     = 2*(1+((`mu'/`z')-1)*(`PG'/`SPG'))

*-----------------------------------------------------------------------------
* 		Beta Lorenz Curve               
*-----------------------------------------------------------------------------

        if ("`grouped'" != "") {
            local aatheta   =   exp(`cofbg'[1,3])
            local aagama    =   `cofbg'[1,1]
            local aagama2   =   2*`aagama'
            local aadelta   =   `cofbg'[1,2]
            local aadelta2  =   2*`aadelta'
        }
        else {
            local aatheta   =   exp(`cofb'[1,3])
            local aagama    =   `cofb'[1,1]
            local aagama2   =   2*`aagama'
            local aadelta   =   `cofb'[1,2]
            local aadelta2  =   2*`aadelta'
        }


        /*** Poverty */

        		local Xx=.00001*(1-`z'/`mu')
        		local hcrb0 = .01
        		local j = 1
        		local ff=`aatheta'*`hcrb0'^`aagama'*(1-`hcrb0')^`aadelta'*((`aagama'/`hcrb0')-(`aadelta'/(1-`hcrb0')))-1+`z'/`mu'
        		while ( `Xx' <`ff' | `ff' < -`Xx') & `j'<51 {
        			local i=1
        			local hcrb1 = 0
        			local hcrb2 = `hcrb0'
        			local ff1=`aatheta'*`hcrb2'^`aagama'*(1-`hcrb2')^`aadelta'*((`aagama'/`hcrb2')-(`aadelta'/(1-`hcrb2')))-1+`z'/`mu'

        			while `hcrb2'-`hcrb1'>.0001 & `i'<500{
        				local ff2=`aatheta'*`hcrb2'^`aagama'*(1-`hcrb2')^`aadelta'*(`aagama'*(`aagama'-1)/*
        				*//`hcrb2'^2-2*`aagama'*`aadelta'/(`hcrb2'*(1-`hcrb2'))+`aadelta'*(`aadelta'-1)/(1-`hcrb2')^2)
        				local hcrb1=`hcrb2'
        				local hcrb2 = `hcrb2' - (`ff1'/`ff2')
        				local ff1=`aatheta'*`hcrb2'^`aagama'*(1-`hcrb2')^`aadelta'*((`aagama'/`hcrb2')-(`aadelta'/(1-`hcrb2')))-1+`z'/`mu'
        				local i =`i'+1
        			}
        			local ff = `ff1'
        			local hcrb0 = `hcrb0' + .01
        			local j = `j'+1
        		}
        		local ff =`ff'
        		local hcrb = `hcrb2'
        		if abs(`ff')> `Xx'{
        			local hcrb4=.0001
        			local hcrb3 = .0001
        			local j=1
        			local fff=`aatheta'*`hcrb3'^`aagama'*(1-`hcrb3')^`aadelta'*((`aagama'/`hcrb3')-(`aadelta'/(1-`hcrb3')))-1+`z'/`mu'
        			while abs(`fff')>`Xx' & `j'<10000{
        				local hcrb3 = `hcrb3'+.0001
        				local j = `j'+ 1
        				local fff1=`aatheta'*`hcrb3'^`aagama'*(1-`hcrb3')^`aadelta'*((`aagama'/`hcrb3')-(`aadelta'/(1-`hcrb3')))-1+`z'/`mu'
        				if abs(`fff1')>abs(`fff'){
        					local fff = `fff'
        				}
        			else{
        				local fff =`fff1'
        				local hcrb4=`hcrb3'
        			}
        		}
        		local hcrb=`hcrb4'
        		local ff =`ff1'
        	}

        	local LhBeta= `hcrb' - `aatheta'* `hcrb'^`aagama'*(1- `hcrb')^`aadelta'
        	local muDiz = `mu'/`z'

            /*** Poverty Gap (Beta Lorenz) */
        	local PgBeta = `hcrb' - `muDiz'*`LhBeta'

            /*** Poverty Gap Saqured (Beta Lorenz) */
        	local ibita1 = (ibeta(`aagama2'-1,`aadelta2'+1,`hcrb'))*exp(lngamma(`aagama2'-1))*exp(lngamma(`aadelta2'/*
        	*/+1))/exp(lngamma(`aagama2'+`aadelta2'))
        	local ibita2 = (ibeta(`aagama2',`aadelta2',`hcrb'))*exp(lngamma(`aagama2'))*exp(lngamma(`aadelta2'/*
        	*/))/exp(lngamma(`aagama2'+`aadelta2'))
        	local ibita3 = (ibeta(`aagama2'+1,`aadelta2'-1,`hcrb'))*exp(lngamma(`aagama2'+1))*exp(lngamma(`aadelta2'/*
        	*/-1))/exp(lngamma(`aagama2'+`aadelta2'))
        	local FgtBeta = (1-`muDiz')*(2*`PgBeta'-(1-`muDiz')*`hcrb')+ `aatheta'*`aatheta'*`muDiz'*`muDiz'*((`aagama'*`aagama'*`ibita1')-/*
        	*/2*`aagama'*`aadelta'*`ibita2' + `aadelta'*`aadelta'*`ibita3')

            /*** Gini (Beta Lorenz) */
        	local GiniBeta = 2*`aatheta'*exp(lngamma(1+`aagama'))*exp(lngamma(1+`aadelta'))/exp(lngamma/*
        	*/(2+`aagama'+`aadelta'))

        	local ldpBeta = `aatheta'*`hcrb'^`aagama'*(1-`hcrb')^`aadelta'*((`aagama'*(1-`aagama')/`hcrb'*`hcrb')+(2*`aagama'*`aadelta'/*
        	*//(`hcrb'*(1-`hcrb')))+(`aadelta'*(1-`aadelta')/((1-`hcrb')*(1-`hcrb'))))

            /*** Elasticities (Beta Lorenz) */
        	local elhmub       = -(`z'/(`mu'*`hcrb'*`ldpBeta'))
        	local elhginib     = (1-(`z'/`mu'))/ (`hcrb'*`ldpBeta')
        	local elpgmub      = 1-(`hcrb'/`PgBeta')
        	local elpgginib    = 1+(((`mu'/`z')-1)* (`hcrb'/`PgBeta'))
        	local elspgmub     = 2*(1-`PgBeta'/`FgtBeta')
        	local elspgginib   = 2*(1+((`mu'/`z')-1)*(`PgBeta'/`FgtBeta'))

			
			
*-----------------------------------------------------------------------------
 * 	Output
*-----------------------------------------------------------------------------

        /*** Display results */

        gen `type2'  = .
        gen `model' = .
        gen `var'   = .
        gen `value' = .

        replace `var'   = _n in 1/24

        * main results type = 1
		replace `type2'  = 1     in 1/8
        
		*elasticities lines (type=2 and type=3)
		replace `type2'  = 2     	in  9
        replace `type2'  = 3     	in  10
        replace `type2'  = 2     	in  11
        replace `type2'  = 3     	in  12
        replace `type2'  = 2     	in  13
        replace `type2'  = 3     	in  14
        replace `type2'  = 2     	in  15
        replace `type2'  = 3     	in  16
        replace `type2'  = 2     	in  17
        replace `type2'  = 3     	in  18
        replace `type2'  = 2     	in  19
        replace `type2'  = 3     	in  20

		* model types 	
        replace `model' = 	1   	in 1/4
        replace `model' = 	1     	in 9/14
        replace `model' = 	2     	in 5/8
        replace `model' = 	2     	in 15/20
		
		
		* main results values
        replace `value' = `H'*100               in  1
        replace `value' = `PG'*100              in  2
        replace `value' = `SPG'*100             in  3
        replace `value' = `gini_ln'             in  4
        replace `value' = `hcrb'*100            in  5
        replace `value' = `PgBeta'*100          in  6
        replace `value' = `FgtBeta'*100         in  7
        replace `value' = `GiniBeta'            in  8
        
		* elasticities
		replace `value' = `elhmu'               in  9
        replace `value' = `elhgini'             in  10
        replace `value' = `elpgmu'              in  11
        replace `value' = `elpggini'            in  12
        replace `value' = `elspgmu'             in  13
        replace `value' = `elspggini'           in  14
        replace `value' = `elhmub'              in  15
        replace `value' = `elhginib'            in  16
        replace `value' = `elpgmub'             in  17
        replace `value' = `elpgginib'           in  18
        replace `value' = `elspgmub'            in  19
        replace `value' = `elspgginib'          in  20

		* output label
		replace `var' = 1       in  1
        replace `var' = 2       in  2
        replace `var' = 3       in  3
        replace `var' = 4       in  4
        replace `var' = 1       in  5
        replace `var' = 2       in  6
        replace `var' = 3       in  7
        replace `var' = 4       in  8
        replace `var' = 1       in  9
        replace `var' = 1       in  10
        replace `var' = 2       in  11
        replace `var' = 2       in  12
        replace `var' = 3       in  13
        replace `var' = 3       in  14
        replace `var' = 1       in  15
        replace `var' = 1       in  16
        replace `var' = 2       in  17
        replace `var' = 2       in  18
        replace `var' = 3       in  19
        replace `var' = 3       in  20

		if ("`benchmark'" == "benchmark") {
			* add apoverty and ainequal to main results (type=1) 
			replace `type2'	=	1 		in 21
			replace `type2'	=	1 		in 22
			replace `type2'	=	1 		in 23
			replace `type2'	=	1 		in 24
			* add apoverty and ainequal model (model=0)
			replace `model'	=	0		in 21/24
			* add apoverty and ainequal values
			replace `value' = `afgt0'            in  21
			replace `value' = `afgt1'            in  22
			replace `value' = `afgt2'            in  23
			replace `value' = `agini'            in  24
			* add output label
			replace `var' = 1       in  21
			replace `var' = 2       in  22
			replace `var' = 3       in  23
			replace `var' = 4       in  24
		}
		
        label define var 1 "FGT(0)", add modify
        label define var 2 "FGT(1)", add modify
        label define var 3 "FGT(2)", add modify
        label define var 4 "Gini", add modify

        label define model 0 "Unit Record", add modify
        label define model 1 "QG Lorenz Curve", add modify
        label define model 2 "Beta Lorenz Curve", add modify

        label define type 1 "Estimated Value", add modify
        label define type 2 "with respect to the Mean", add modify
        label define type 3 "with respect to the Gini", add modify

*         gen type = `type'
*         gen model = `model'
*         gen var     = `var'
*         gen value   = `value'

        label values `model' model
        label values `type2' type
        label values `var' var
		
		label var `model' Model
		label var `type2'  Type
		label var `var'   Indicator

		
*-----------------------------------------------------------------------------
* Display Lorenz
*-----------------------------------------------------------------------------
		
		label var `pg' p
		label var `Lg' Lorenz

		format `pg' %16.2f
		format `Lg' %16.3f  
	
		`noilor' di 			""
		`noilor' di as text 	"{hline 15}    Distribution    {hline 15}"
		`noilor' di as text 	_col(5) "i "    _col(15) "P"   _col(40) "L" 
		`noilor' di as text 	"{hline 50}"
		
		forvalues l = 1(1)`bins' {
			local P = `pg' in `l'
			local L = `Lg' in `l'
			`noilor' di as text _col(5) "`l'"  as res  _col(15) %5.4f `P'   _col(40) %5.4f `L'
		}
		
		`noilor' di as text 	"{hline 50}"
	
	
*-----------------------------------------------------------------------------
* Display Regression results 
*-----------------------------------------------------------------------------
		 
	if ("`grouped'" == "grouped") {
		 
        `noireg' di ""
		`noireg' di ""
        `noireg' di as text "Estimation: " as res "GQ Lorenz Curve (grouped data)"
		`noireg' estout  coefgqg, cells("b(star fmt(%9.3f)) se t p")                ///
			  stats(r2_a F rmse mss rss N, fmt(%9.3f %9.0g) labels("Adj. R-squared"))      ///
			  legend label  

		`noireg' di ""
        `noireg' di ""
        `noireg' di as text "Estimation: " as res "Beta Lorenz Curve (Grouped data)"
		`noireg' estout coefbetag, cells("b(star fmt(%9.3f)) se t p")                ///
			  stats(r2_a F rmse mss rss N, fmt(%9.3f %9.0g) labels("Adj. R-squared" F-sta RMSE MSS RSS Obs))      ///
			  legend label  varlabels(_cons A)  

	}

	if ("`grouped'" == "") {

		`noireg' di ""
        `noireg' di ""
        `noireg' di as text "Estimation: " as res "GQ Lorenz Curve"
		`noireg' estout  coefgq, cells("b(star fmt(%9.3f)) se t p")                ///
			  stats(r2_a F rmse mss rss N, fmt(%9.3f %9.0g) labels("Adj. R-squared" F-sta RMSE MSS RSS Obs))      ///
			  legend label    

		`noireg' di ""
        `noireg' di ""
        `noireg' di as text "Estimation: " as res "Beta Lorenz Curve"
		`noi2reg' estout coefbeta, cells("b(star fmt(%9.3f)) se t p")                ///
			  stats(r2_a F rmse mss rss N, fmt(%9.3f %9.0g) labels("Adj. R-squared" F-sta RMSE MSS RSS Obs))      ///
			  legend label  varlabels(_cons A)  

	}
	
*-----------------------------------------------------------------------------
* Display POverty and Inequality Results 
*-----------------------------------------------------------------------------
		
        noi di ""
        noi di ""
        noi di "Estimated Poverty and Inequality Measures:"
        noi tabdisp `var' `model' if `var' != . & `type2' == 1, cell(`value')
        noi di "Mean `inc': " as res %16.2f `mu'

*-----------------------------------------------------------------------------
* Display Elasticities 
*-----------------------------------------------------------------------------
	

		`noelast' di ""
        `noelast' di ""
        `noelast' di "Estimated Elasticities:"
        `noelast' tabdisp `var' `model' `type2' if `var' != . & `type2' != 1 & `value' != . , cell(`value')

	
*-----------------------------------------------------------------------------
*  Checking for consistency of lorenz curve estimation (section 4)
*-----------------------------------------------------------------------------

		noi di as text "Estimation Validity"

        ***********************
        /* GQ Lorenz Curve */
        ***********************
        quietly {

            `nocheck' di ""
            `nocheck' di ""
            `nocheck' di as text "Checking for consistency of lorenz curve estimation: " as res "GQ Lorenz Curve"

            /** Condition 1 */
            if (`e' < 0) {
                `nocheck' di as text "L(0;pi)=0: " as res  "OK"
                local ccheck1 = 1
            }
            else {
                `nocheck' di as text "L(0;pi)=0: " as err "FAIL"
                local ccheck1 = 0
            }

            /** Condition 2 */
            local t = (`a'+`c')
            if (`t' >= 1) {
                `nocheck' di as text "L(1;pi)=1: " as res "OK (value=" %9.4f `t' ")"
                local ccheck2 = 1
            }
            else {
                noi`nocheck'di as text "L(1;pi)=1: " as err "FAIL (value=" %9.4f `t' ")"
                local ccheck2 = 0
            }

            /** Condition 3 */
            if (`c' >= 0) {
                `nocheck' di as text "L'(0+;pi)>=0: " as res  "OK"
                local ccheck3 = 1
            }
            else {
                `nocheck' di as text "L'(0+;pi)>=0: " as err "FAIL"
                local ccheck3 = 0
            }


            /** Condition 4 */

            if ( `m' < 0 | (( 0 < `m' <(`n'^2/(4*`e'^2)))	& `n' >= 0) | ((0 < `m' < (-`n'/2)) & (`m' < (`n'^2 /(4*`e'^2))))) {
                `nocheck' di as text "L''(p;pi)>=0 for p within (0,1): " as res  "OK"
                local ccheck4 = 1
            }
            else {
                `nocheck' di as text "L''(p;pi)>=0 for p within (0,1): " as err "FAIL"
                local ccheck4 = 0
            }

        }

        ***********************
        /* Beta Lorenz curve */
        ***********************

        `nocheck' di ""
        `nocheck' di as text "Checking for consistency of lorenz curve estimation: " as res "Beta Lorenz curve"

        /** Condition 1 */
        * automatically satisfied by the functional form

        /** Condition 2 */
        * automatically satisfied by the functional form

        /** Condition 3 */
        	* We check the validity of the Beta Lorenz curve
        	local check1 = 1- `aatheta'*.001^`aagama'*.999^`aadelta'*(`aagama'/.001-`aadelta'/.999)

        /** Condition 4 */

        	local check2 = 0
        	local i=.01
        	while `i'<1{
        		local chk = `aatheta'*`i'^`aagama'*(1-`i')^`aadelta'*((`aagama'*(1-`aagama')/`i'^2)+(2*`aagama'*`aadelta'/*
        		*//(`i'*(1-`i')))+(`aadelta'*(1-`aadelta')/(1-`i')^2))
        		if `chk'<0{
        			local check2=1
        		}
        		else{
        		}
        		local i=`i'+.01
        	}

        `nocheck' di as text "L(0;pi)=0: " as res "OK (automatically satisfied by the functional form)"

        `nocheck' di as text "L(1;pi)=1: " as res "OK (automatically satisfied by the functional form)"

        if `check1'>=0  {
            `nocheck' di as text "L'(0+;pi)>=0: " as res  "OK"
			local bcheck3 = 1
        }
        else {
            `nocheck' di as text "L'(0+;pi)>=0: " as err "FAIL "
			local bcheck3 = 0
        }

        if `check2'==0 {
            `nocheck' di as text "L''(p;pi)>=0 for p within (0,1): " as res  "OK"
			local bcheck4 = 1
        }
        else {
            `nocheck' di as text "L''(p;pi)>=0 for p within (0,1): " as err "FAIL"
			local bcheck4 = 0
        }

		`nocheck' di ""
		`nocheck' di ""
		
				
*-----------------------------------------------------------------------------
* Store results 
*-----------------------------------------------------------------------------

		mkmat  `var' `model' `type2' `value' if `value' != ., matrix(`tmp') 

		matrix colnames `tmp' = indicator model type value
		
		mat check = `tmp'
		
		matrix rownames `tmp' = H PG SPG gini_ln hcrb PgBeta	FgtBeta	GiniBeta  elhmu	 elhgini	elpgmu	elpggini	elspgmu	elspggini	elhmub	elhginib	elpgmub	elpgginib	elspgmub	elspgginib `rownames_unitrecord'  
		
		return matrix results = `tmp'
		
        return scalar Hgq   		= `H'*100
        return scalar PGgq  		= `PG'*100
        return scalar SPGgq 		= `SPG'*100
        return scalar GINIgq  		= `gini_ln'
        return scalar Hb    		= `hcrb'*100
        return scalar PGb   		= `PgBeta'*100
        return scalar SPGb  		= `FgtBeta'*100
        return scalar GINIb 		= `GiniBeta'
        return scalar elhmu       	= `elhmu'
        return scalar elhgini     	= `elhgini'
        return scalar elpgmu      	= `elpgmu'
        return scalar elpggini    	= `elpggini'
        return scalar elspgmu     	= `elspgmu'
        return scalar elspggini   	= `elspggini'
        return scalar elhmub      	= `elhmub'
        return scalar elhginib    	= `elhginib'
        return scalar elpgmub     	= `elpgmub'
        return scalar elpgginib   	= `elpgginib'
        return scalar elspgmub    	= `elspgmub'
        return scalar elspgginib  	= `elspgginib'
	if ("`nochecks'" != "") {
        return scalar check1b   	= 1
        return scalar check2b   	= 1
        return scalar check3b   	= `bcheck3'
        return scalar check4b   	= `bcheck4'
        return scalar check1gq  	= `ccheck1'
        return scalar check2gq  	= `ccheck2'
        return scalar check3gq  	= `ccheck3'
        return scalar check4gq  	= `ccheck4'
        return scalar t         	= `t'
	}
        return scalar mu        	= `mu'
		
}

	return add

	cap: drop yg ag bg cg yg2 x1g x2g

	cap: drop y1 a b c y2 x1 x2
	
end