$PROBLEM Pyronaridine -- 3 cmp + 4 transit (IOV, F1+MTT, fixed iiv Q1)
$INPUT 		ID 
			DATE = DROP 						
			TIME 						
			EVID 						; Dose = 1, observation = 0
			MDV 						; Missing data value
			NDDV 						; Non-transformed DV
			DV 							; Logarithm DV
			AMT 						
			CMT 						; Dose = cmt(1), observation = cmt(2)
			OCC 						; visit 3 = 3, visit 4 = 4
			PQFLAG 						; PA + PQ = 1, PA = 0
			SEX 						; Male = 0, Female = 1
			HT = DROP
			WT
			AGE
			WBC
			HB
			HCT
			PLT
			NEUT
			LYMPH
			MONO
			EO
			BASO
			METHM
			PROTN
			GLUC
			BUN
			CREATE
			DIRB
			TOTB
			ALBM
			ALP
			AST
			ALT
			NA
			K
			CK

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
  			; Occasions
  			OC1 = 0
  			OC2 = 0
  			IF(OCC.EQ.3) OC1 = 1
  			IF(OCC.EQ.4) OC2 = 1

  			IOV = ETA(9) * OC1 + ETA(10) * OC2
  			IOV2 = ETA(11) * OC1 + ETA(12) * OC2

			TVF1 = THETA(1)
			F1 = TVF1 * EXP(ETA(1) + IOV2)

			TVMTT = THETA(2)
			MTT = TVMTT * EXP(ETA(2) + IOV) 				; Mean transit time

			TVCL = THETA(3)
			CL = TVCL * (WT/58.5)**0.75 * EXP(ETA(3))  		; Individual clearance

			TVV2 = THETA(4)
			V2 = TVV2 * (WT/58.5)**1 * EXP(ETA(4)) 			; Individual volume for central compartment

			TVV3 = THETA(5)
			V3 = TVV3 * (WT/58.5)**1 * EXP(ETA(5)) 			; Individual volume for peripheral compartment 1

			TVQ1 = THETA(6)
			Q1 = TVQ1 * (WT/58.5)**0.75 * EXP(ETA(6))

			TVV4 = THETA(7) 	
			V4 = TVV4 * (WT/58.5)**1 * EXP(ETA(7)) 			; Individual volume for peripheral compartment 2

			TVQ2 = THETA(8)
			Q2 = TVQ2 * (WT/58.5)**0.75 * EXP(ETA(8)) 

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
$THETA 		(0,0.9) 							; MTT
$THETA 		(0,40)						 		; CL
$THETA 		(0,1500) 					 		; V2
$THETA 		(0,5000) 							; V3
$THETA 		(0,60) 							; Q1
$THETA 		(0,10000) 							; V4
$THETA 		(0,30) 							; Q2


$OMEGA 		0.04 								; F1
$OMEGA 		0 FIX 								; MTT
$OMEGA 		0 FIX								; CL
$OMEGA 		0.02 								; V2
$OMEGA 		0 FIX 								; V3
$OMEGA 		0 FIX 								; Q1
$OMEGA 		0 FIX 								; V4
$OMEGA 		0.07 								; Q2
$OMEGA 		BLOCK(1) 
0.1  ; IOV_MTT
$OMEGA		BLOCK(1) SAME
$OMEGA		BLOCK(1) 
0.1  ; IOV_F
$OMEGA		BLOCK(1) SAME


$SIGMA 		0.05

$ESTIMATION
		     MAXEVAL=9999                                   	; Number of evaluations
		     PRINT=5                                        	; Print every 5th iteration
		     METHOD=1                                       	; Estimation method: 1=First Order Conditional (FOCE)
		     INTER                                          	; Estimation with interaction between ETA's and EPS's

;$COV 		PRINT = E
;			UNCONDITIONAL

$TABLE ID TIME TAD EVID NDDV DV MDV CMT OCC PQFLAG WT IPRED CWRES IWRES F1 MTT IOV NOPRINT ONEHEADER FILE=mytab87
$TABLE ID TIME F1 MTT CL V2 V3 V4 Q1 Q2 NOPRINT ONEHEADER FILE=patab87
