$PROBLEM Pyronaridine -- 2 compartment model + 4 transit (F1)
$INPUT 		ID 
			DATE = DROP 						
			TIME 						
			EVID 						; Dose = 1, observation = 0
			MDV 							; Missing data value
			NDDV 						; Non-transformed DV
			DV 							; Logarithm DV
			AMT 						
			CMT 							; Dose = cmt(1), observation = cmt(2)
			ADMIN 						; visit 3 = 3, visit 4 = 4
			PQFLAG 						; PA + PQ = 1, PA = 0
			SEX 							; Male = 0, Female = 1
			HT
			WT

$DATA 		PYN.csv
			IGNORE = @ 					; Ignore non-numerical values

$SUBROUTINES ADVAN5 TRANS1

$MODEL 		COMP = (1)					; Absorption compartment
			COMP = (2)			 		; Central compartment
			COMP = (3) 					; Peripheral compartment
			COMP = (4) 					; Transit compartment 1
			COMP = (5) 					; Transit compartment 2
			COMP = (6) 					; Transit compartment 3
			COMP = (7) 					; Transit compartment 4

$PK			
			F1 = THETA(1) * EXP(ETA(1))
			MTT = THETA(2) * EXP(ETA(2)) 		; Mean transit time
			CL = THETA(3) * EXP(ETA(3)) 		; Individual clearance
			V2 = THETA(4) * EXP(ETA(4)) 		; Individual volume for central compartment
			V3 = THETA(5) * EXP(ETA(5)) 		; Individual volume for peripheral compartment
			Q  = THETA(6) * EXP(ETA(6))	

			KTR = 5 / MTT

			K20 = CL / V2
			K23 = Q / V2
			K32 = Q / V3
			S2 = V2 * 1000
			K14 = KTR
			K45 = KTR
			K56 = KTR
			K67 = KTR
			K72 = KTR

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

$THETA 		1 FIX 						; F1
$THETA 		(0,5) 						; MTT
$THETA 		(0,50)						; CL
$THETA 		(0,1000) 					 	; V2
$THETA 		(0,10000) 				 	; V3
$THETA 		(0,100) 						; Q

$OMEGA 		0.01 							; F1
$OMEGA 		0.01 							; MTT
$OMEGA 		0.01							; CL
$OMEGA 		0.01 							; V2
$OMEGA 		0.01 							; V3
$OMEGA 		0.01							; Q

$SIGMA 		0.02

$ESTIMATION
     MAXEVAL=9999                                   	; Number of evaluations
     PRINT=5                                        	; Print every 5th iteration
     METHOD=1                                       	; Estimation method: 1=First Order Conditional (FOCE)
     INTER                                          	; Estimation with interaction between ETA's and EPS's

;$COV 		PRINT = E
;			UNCONDITIONAL

$TABLE ID TIME TAD EVID NDDV DV MDV CMT ADMIN PQFLAG SEX HT WT IPRED IWRES NOPRINT ONEHEADER FILE=mytab44
