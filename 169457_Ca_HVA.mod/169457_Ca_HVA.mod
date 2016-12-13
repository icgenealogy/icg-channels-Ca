:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica, offma, offmb, offha, offhb, sloma, slomb, sloha, slohb, tauma, taumb, tauha, tauhb
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_HVAbar = 0.00001 (S/cm2) 
        offma = -27 (mV)
        offmb = -75 (mV)
        offha = -13 (mV)
        offhb = -15 (mV)
        sloma = 3.8 (mV)
        slomb = 17 (mV)
        sloha = 50 (mV)
        slohb = 28 (mV)
	tauma = 18.1818 (ms)
	taumb = 1.06383 (ms)
	tauha = 2188.18 (ms)
	tauhb = 153.846 (ms)
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_HVA	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_HVA = gCa_HVAbar*m*m*h
	ica = gCa_HVA*(v-eca)
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
	UNITSOFF
        if((v == offma) ){        
            v = v+0.0001
        }
		mAlpha =  (offma-v)/tauma/(exp((offma-v)/sloma) - 1)        
		mBeta  =  exp((offmb-v)/slomb)/taumb
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta)
		hAlpha =  exp((offha-v)/sloha)/tauha
		hBeta  =  1.0/tauhb/(exp((offhb-v)/slohb)+1)
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}
