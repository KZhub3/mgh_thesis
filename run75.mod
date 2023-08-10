;;----------------------------------------------------------------------------;;
;; MORU NONMEM Pyronaridine -- drug-drug interaction V3
;; 
;;----------------------------------------------------------------------------;;
;; Kuangyi Zhang
;; Date: 2023-05-16
;;----------------------------------------------------------------------------;;

$PROBLEM Pyronaridine -- drug-drug interaction V3
$INPUT 		ID 
			DATE = DROP 						
			TIME 						
			EVID 						; Dose = 1, observation = 0
			MDV 						; Missing data value
			NDDV 						; Non-transformed DV
			DV 							; Logarithm DV
			AMT 						
			CMT 						; Dose = cmt(1), observation = cmt(2)
			ADMIN 						; visit 3 = 3, visit 4 = 4
			PQFLAG 						; PA + PQ = 1, PA = 0
			SEX 						; Male = 0, Female = 1
			HT
			WT

$DATA 		PYN.csv
			IGNORE = @ 					; Ignore non-numerical values

$SUBROUTINES ADVAN5 TRANS1

$MODEL 		COMP = (1)					; Absorption compartment
			COMP = (2)			 		; Central compartment
			COMP = (3) 					; Peripheral compartment 1
			COMP = (4) 					; Peripheral compartment 2
			COMP = (5) 					; Transit compartment 1
			COMP = (6) 					; Transit compartment 2
			COMP = (7) 					; Transit compartment 3
			COMP = (8) 					; Transit compartment 4

$PK			
			IF(PQFLAG.EQ.1) COV1 = (1 + THETA(9))  				; PA + PQ
			IF(PQFLAG.EQ.0) COV1 = 1	

			F1 = THETA(1) * EXP(ETA(1)) 
			MTT = THETA(2) * EXP(ETA(2)) 												; Mean transit time
			CL = THETA(3) * (WT/58.5)**0.75 * EXP(ETA(3))   						; Individual clearance
			V2 = THETA(4) * (WT/58.5)**1 * EXP(ETA(4)) 									; Individual volumn for central compartment
			V3 = THETA(5) * (WT/58.5)**1 * EXP(ETA(5)) * COV1									; Individual volumn for peripheral compartment 1
			Q1 = THETA(6) * (WT/58.5)**0.75 * EXP(ETA(6)) 	
			V4 = THETA(7) * (WT/58.5)**1 * EXP(ETA(7)) 									; Individual volumn for peripheral compartment 2
			Q2 = THETA(8) * (WT/58.5)**0.75 * EXP(ETA(8)) 

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

$ERROR 		IPRED = A(2) / S2
			IF(IPRED.GT.0) 	IPRED = LOG(IPRED)
			Y = IPRED + EPS(1)

			W = SQRT(SIGMA(1,1))
			IRES = DV - IPRED
			IWRES = IRES / W

$THETA 		1 FIX 								; F1 
$THETA 		(0,0.5) 							; MTT
$THETA 		(0,50)						 		; CL
$THETA 		(0,5000) 					 		; V2
$THETA 		(0,5000) 							; V3
$THETA 		(0,100) 							; Q1
$THETA 		(0,5000) 							; V4
$THETA 		(0,100) 							; Q2
$THETA 		(0,0.1) 								; COV1


$OMEGA 		0.01 								; F1
$OMEGA 		0.01 								; MTT
$OMEGA 		0.01								; CL
$OMEGA 		0.01 								; V2
$OMEGA 		0.01 								; V3
$OMEGA 		0.01 								; Q1
$OMEGA 		0.01 								; V4
$OMEGA 		0.01 								; Q2


$SIGMA 		0.05

$ESTIMATION
		     MAXEVAL=9999                                   	; Number of evaluations
		     PRINT=5                                        	; Print every 5th iteration
		     METHOD=1                                       	; Estimation method: 1=First Order Conditional (FOCE)
		     INTER                                          	; Estimation with interaction between ETA's and EPS's

;$COV 		PRINT = E
;			UNCONDITIONAL

$TABLE ID TIME TAD EVID NDDV DV MDV CMT ADMIN PQFLAG SEX HT WT IPRED CWRES IWRES NOPRINT ONEHEADER FILE=mytab75