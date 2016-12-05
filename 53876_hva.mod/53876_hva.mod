COMMENT
This file, hva.mod, implements the high voltage activated calcium current gCa(HVA)
Quadroni and Knopfel 1994 table 1
ENDCOMMENT

NEURON {
	SUFFIX hva
	:	NONSPECIFIC_CURRENT i
	USEION ca READ eca WRITE ica
	RANGE gbar
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 1110e-6	(S/cm2) < 0, 1e9 >
	:Erev = 80 (mV)
}

ASSIGNED {
        eca (mV)
	ica (mA/cm2)
	:i (mA/cm2)
	v (mV)
	g (S/cm2)
	sinf
	rinf
	tau_s (ms)
	tau_r (ms)
}

STATE {	s r }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * s*s * r
	ica = g * (v - eca)
	:i = ica	: supplied for diagnostic graphing
}

INITIAL {
	: assume that v has been constant for a long time
	s = alphas(v)/(alphas(v) + betas(v))
	r = alphar(v)/(alphar(v) + betar(v))
}

DERIVATIVE states {
	rates(v)
	s' = (sinf - s)/tau_s
	r' = (rinf - r)/tau_r
}

LOCAL alpha_r, alpha_s	: stores value to save a couple of function calls

FUNCTION alphas(Vm (mV)) (/ms) {
	UNITSOFF
	alphas = 2.0 /(1 + exp(-(Vm + 2.0)* 0.054))
	UNITSON
}

FUNCTION betas(Vm (mV)) (/ms) {
	UNITSOFF
	betas =  -0.08 * (Vm + 15.9) / (1 - exp( (Vm + 15.9)*0.2 ))
	UNITSON
}

FUNCTION taus(Vm (mV)) (/ms) {
	UNITSOFF
	taus = 1.0 / (alpha_s + betas(Vm))	: taus only called from rates
	UNITSON
}

FUNCTION alphar(Vm (mV)) (/ms) {
	UNITSOFF
	alphar = 0.01 * exp( -(Vm + 60)/20 )
	if (alphar > 0.01) {
		alphar = 0.01
	}
	UNITSON
}

FUNCTION betar(Vm (mV)) (/ms) {
	UNITSOFF
	betar =  0.01 - alphar(Vm)
	UNITSON
}

FUNCTION taur(Vm (mV)) (/ms) {
	UNITSOFF
	taur = 1.0 / (alpha_r + betar(Vm))	: taur only called from rates
	UNITSON
}

::::: special warning - if any of the above rate functions are desired to be
::::: called from hoc, the rates function below needs to be called first each
::::: time to set the alpha_r, and alpha_s variables.

PROCEDURE rates(Vm (mV)) {
	alpha_s = alphas(Vm)
	tau_s = taus(Vm)
	sinf = alphas(Vm) * tau_s

	alpha_r = alphar(Vm)
	tau_r = taur(Vm)
	rinf = alpha_r * tau_r
}
