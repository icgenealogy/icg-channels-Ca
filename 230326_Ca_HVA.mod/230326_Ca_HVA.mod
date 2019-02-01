:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca, cai, cao WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica
	GLOBAL use_ghk
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_HVAbar = 0.00001 (S/cm2) 
	use_ghk = 0

}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	cai (mM)
	cao (mM)
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

UNITSOFF
FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}
FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
UNITSON

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_HVA = gCa_HVAbar*m*m*h
	if (use_ghk == 0) {
		ica = gCa_HVA*(v-eca)
	}
	if (use_ghk == 1) {
		ica = gCa_HVA*(ghk(v,cai,cao)-106)
	}
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
        if((v == -27) ){        
            v = v+0.0001
        }
		mAlpha =  (0.055*(-27-v))/(exp((-27-v)/3.8) - 1)        
		mBeta  =  (0.94*exp((-75-v)/17))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta)
		hAlpha =  (0.000457*exp((-13-v)/50))
		hBeta  =  (0.0065/(exp((-v-15)/28)+1))
		hInf = hAlpha/(hAlpha + hBeta)
		hTau = 1/(hAlpha + hBeta)
	UNITSON
}




