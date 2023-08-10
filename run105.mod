;;----------------------------------------------------------------------------;;
;; MORU NONMEM Pyronaridine -- PD linear model
;; 
;;----------------------------------------------------------------------------;;
;; Kuangyi Zhang
;; Date: 2023-07-01
;;----------------------------------------------------------------------------;;

$PROBLEM Pyronaridine -- PD linear model
$INPUT 		ID 
			DATE = DROP 						
			TIME
			IPK1 = DROP
			IPK
			HR
			UQT
			RR
			QTCF
			DV
			MDV
			EVID
			AMT
			OCC  				;=ADMIN						
			CMT
			IF1
			IMTT
			ICL
			IV2
			IV3
			IV4
			IQ1
			IQ2
			PD
			WT

$DATA 		ECG.csv
			IGNORE = @ 					; Ignore non-numerical values

$SUBROUTINES ADVAN5 TRANS1

$MODEL 		COMP = (absorp)					; Absorption compartment
			COMP = (central)			 		; Central compartment
			COMP = (peri1) 					; Peripheral compartment 1
			COMP = (peri2) 					; Peripheral compartment 2
			COMP = (trans1) 					; Transit compartment 1
			COMP = (trans2) 					; Transit compartment 2
			COMP = (trans3) 					; Transit compartment 3
			COMP = (trans4) 					; Transit compartment 4

$PK			

			F1 = IF1
			MTT = IMTT			; Mean transit time
			CL = ICL 			; Individual clearance
			V2 = IV2			; Individual volumn for central compartment
			V3 = IV3 			; Individual volumn for peripheral compartment 1
			Q1 = IQ1
			V4 = IV4 			; Individual volumn for peripheral compartment 2
			Q2 = IQ2

			TVSLOPE = THETA(1)
			SLOPE = TVSLOPE + ETA(1)

			TVBASE = THETA(2)
			BASE = TVBASE + ETA(2)


; INTERCEPT

			KTR = 5 / MTT
			K20 = CL / V2
			K23 = Q1 / V2
			K32 = Q1 / V3
			K24 = Q2 / V2
			K42 = Q2 / V4
			S2 = V2 * 1000
			K15 = KTR
			K56 = KTR
			K67 = KTR
			K78 = KTR
			K82 = KTR

			IF(NEWIND.LE.1) TDOS = 0
 
  			IF(EVID.EQ.1) THEN
    			TDOS = TIME
  			ENDIF
 
  			TAD = TIME-TDOS


$ERROR 		
			CONC = A(2) / S2
			IPRED = SLOPE * CONC + BASE			; add intercept
			Y = IPRED + EPS(1)


$THETA 		(-10,0.001) 								; SLOPE
$THETA 		(-5,-2) 								; BASE

$OMEGA 		0 FIX 								; SLOPE
$OMEGA 		0.3 								; BASE

$SIGMA 		0.05

$ESTIMATION
		     MAXEVAL=9999                                   	; Number of evaluations
		     PRINT=5                                        	; Print every 5th iteration
		     METHOD=1                                       	; Estimation method: 1=First Order Conditional (FOCE)
		     INTER                                          	; Estimation with interaction between ETA's and EPS's

$COV 		PRINT = E
			UNCONDITIONAL

$TABLE ID TIME IPK QTCF TAD EVID DV MDV CMT OCC WT CONC IPRED SLOPE CWRES NOPRINT ONEHEADER FILE=mytab105