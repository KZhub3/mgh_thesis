$PROBLEM Pyronaridine -- final model (calculate secondary parameters)
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

$ABBREVIATED COMRES=2 					; For Cmax and Tmax

$SUBROUTINES ADVAN6 TOL=6; TRANS=1

$MODEL 		COMP = (absorp)					; Absorption compartment
			COMP = (central)			 	; Central compartment
			COMP = (peri1) 					; Peripheral compartment 1
			COMP = (peri2) 					; Peripheral compartment 2
			COMP = (trans1) 				; Transit compartment 1
			COMP = (trans2) 				; Transit compartment 2
			COMP = (trans3) 				; Transit compartment 3
			COMP = (trans4) 				; Transit compartment 4

			COMP = (auc) 					; AUC

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

			IF(NEWIND.LE.1) THEN
				TDOS = 0
				COM(1) = -1 			;Cmax
				COM(2) = -1 			;Tmax
			ENDIF

  			IF(EVID.EQ.1) THEN
    			TDOS = TIME
  			ENDIF
 
  			TAD = TIME-TDOS

  			A0 = K20*K32*K42 									; For half-lives
			A1 = K20*K42+K32*K42+K32*K24+K20*K32+K42*K23
			A2 = K23+K32+K24+K42+K20
			 
			P1  = A1-(A2*A2/3)
			QP  = (2*A2*A2*A2/27)-(A1*A2/3)+A0
			RT1 = SQRT(-1*P1*P1*P1/27)
			PHI = ACOS((-1*QP/2)/RT1)/3
			RT2 = 2*EXP(LOG(RT1)/3)
			ROOT1 = -1*(COS(PHI)*RT2-A2/3)
			ROOT2 = -1*(COS(PHI+2*3.141593/3)*RT2-A2/3)
			ROOT3 = -1*(COS(PHI+4*3.141593/3)*RT2-A2/3)
			 
			HLA = LOG(2)/ROOT2
			HLB = LOG(2)/ROOT3
			HLG = LOG(2)/ROOT1
			 
			HL = 0 												; HL would be the longest half-life among HLA, HLB and HLC
			IF (ROOT1.LT.ROOT2.AND.ROOT1.LT.ROOT3) HL=HLG
			IF (ROOT2.LT.ROOT1.AND.ROOT2.LT.ROOT3) HL=HLA
			IF (ROOT3.LT.ROOT1.AND.ROOT3.LT.ROOT2) HL=HLB

$DES
			DADT(1) = -KTR*A(1)
			DADT(2) = KTR*A(8) + K32*A(3) + K42*A(4) - K23*A(2) - K24*A(2) - K20*A(2)
			DADT(3) = K23*A(2) - K32*A(3)
			DADT(4) = K24*A(2) - K42*A(4)
			DADT(5) = KTR*A(1) - KTR*A(5)
			DADT(6) = KTR*A(5) - KTR*A(6)
			DADT(7) = KTR*A(6) - KTR*A(7)
			DADT(8) = KTR*A(7) - KTR*A(8)

			DADT(9) = A(2)

			AUC = A(9) / S2 				; Calculate AUC

			CT = A(2) / S2
			IF(CT.GT.COM(1)) THEN 			; For Cmax and Tmax
				COM(1) = CT
				COM(2) = TAD
			ENDIF


$ERROR 		IPRED = A(2) / S2
			IF(IPRED.GT.0) 	IPRED = LOG(IPRED)
			Y = IPRED + EPS(1)

			W = SQRT(SIGMA(1,1))
			IRES = DV - IPRED
			IWRES = IRES / W

			CMAX = COM(1)
			TMAX = COM(2)

$THETA 		1 FIX 								; F1 
$THETA 		(0,0.9) 							; MTT
$THETA 		(0,40)						 		; CL
$THETA 		(0,1500) 					 		; V2
$THETA 		(0,5000) 							; V3
$THETA 		(0,60) 								; Q1
$THETA 		(0,10000) 							; V4
$THETA 		(0,30) 								; Q2


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

$TABLE ID TIME TAD EVID NDDV DV MDV CMT OCC PQFLAG WT IPRED CWRES IWRES F1 MTT IOV AUC CMAX TMAX HLA HLB HLG HL NOPRINT ONEHEADER FILE=mytab100
$TABLE ID TIME F1 MTT CL V2 V3 V4 Q1 Q2 NOPRINT ONEHEADER FILE=patab100
