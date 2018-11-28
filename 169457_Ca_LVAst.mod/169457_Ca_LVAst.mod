:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Ca_LVAst
	USEION ca READ eca WRITE ica
	RANGE gCa_LVAstbar, gCa_LVAst, ica, offma, offmt, offha, offht, sloma, slomt, sloha, sloht, taummin, taumdiff, tauhmin, tauhdiff
	GLOBAL eca
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_LVAstbar = 0.00001 (S/cm2)
	offma = -40.0 (mV)
	offmt = -35.0 (mV)
	offha = -90.0 (mV)
	offht = -50.0 (mV)
	sloma = 6.0 (mV)
	slomt = 5.0 (mV)
	sloha = 6.4 (mV)
	sloht = 7.0 (mV)
	taummin = 5.0 (ms)
	taumdiff = 20.0 (ms)
	tauhmin = 20.0 (ms)
	tauhdiff = 50.0 (ms)
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_LVAst = gCa_LVAstbar*m*m*h
	ica = gCa_LVAst*(v-eca)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		mInf = 1.0000/(1+ exp((offma-v)/sloma))
		mTau = (taummin + taumdiff/(1+exp(-(offmt-v)/slomt)))/qt
		hInf = 1.0000/(1+ exp(-(offha-v)/sloha))
		hTau = (tauhmin + tauhdiff/(1+exp(-(offht-v)/sloht)))/qt
	UNITSON
}
