COMMENT
Described by Gruber et al. 2003, which they based on Bargas et al. 1994

Gruber, A.J., Solla, S.A., Surmeier, D.J., and Houk, J.C.
Modulation of striatal single units by expected reward: 
a spiny neuron model displaying dopamine-induced bistability.
J. Neurophysiol. 90:1095-1114, 2003.

Bargas, J., Howe, A., Eberwine, J., and Surmeier, D.J.
Cellular and molecular characterization of Ca2+ currents in acutely isolated, 
adult rat neostriatal neurons.
J. Neurosci. 14:6667-6686, 1994.

Unlike the formulation used by Gruber et al., 
which assumed instantaneous activation, 
this implementation assumes a constant activation time constant 
that is relatively fast compared to the time scale of the model 
(100-1000 ms).
ENDCOMMENT

NEURON {
	SUFFIX caL
	USEION ca READ cai, cao WRITE ica
	RANGE Pbar, P, i
	GLOBAL minf, mtau
	POINTER mu : hoc level DAsyn[i].msg--see dasyn.mod
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	vh = -35	(mV)	: half activation
	ve = 6.1	(mV)	: slope
	mtauconst = 0.1	(ms)	: m activates much faster than 100-1000 ms
: Pbar value in Table 1 of Gruber et al.
: produces a current that is 10 times too small
:	Pbar = 4.2	(nanometer/s)	<0,1e9>
	Pbar = 42	(nanometer/s)	<0,1e9>
	cao = 2		(mM)
	cai = 10e-6	(mM)
}

ASSIGNED {
	celsius	(degC)
	v	(mV)
	i	(uA/cm2)	: for consistency with their usage of uA/cm2
	ica	(mA/cm2)
	minf	(1)
	mtau	(ms)
	P	(nanometer/s)
	zFRT	(1/volt)
	zVFRT	(1)
	ghk	(coulomb/liter)
	mu	(1)
}


STATE {
	m
}


INITIAL {
	zFRT = (1000)*2*FARADAY/(R*(273+celsius))
	rates(v)
	m = minf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	P = mu*Pbar*m
	zVFRT = (0.001)*zFRT*v
	ghk = 2*FARADAY*(cai - cao*exp(-zVFRT))*gtrap(zVFRT)
	i = (1e-4)*P*ghk
	ica = (1e-3)*i
}


DERIVATIVE states { 
	rates(v)
	m' = (minf-m)/mtau
}


: rates() computes rate and other constants at present v
: call once from hoc to initialize inf at resting v
PROCEDURE rates(v(mV)) {
UNITSOFF
	: "m" gca activation
	mtau = mtauconst
	minf = 1/(1 + exp(-(v - vh)/ve))
}
UNITSON

: traps for 0 in denominator of ghk
FUNCTION gtrap(x) {
	if (fabs(x) < 1e-6) {
		gtrap = 1 + x/2
	} else {
		gtrap = x/(1 - exp(-x))
	}
}
