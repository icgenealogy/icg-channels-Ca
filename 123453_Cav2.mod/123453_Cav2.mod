TITLE Cav2.1 (P-type) calcium channel

COMMENT

NEURON implementation of a Cav2.1 calcium channel
Kinetical scheme: Hodgkin-Huxley (m), no inactivation

Model includes a calculation of the gating current

Modified from Khaliq et al., J. Neurosci. 23(2003)4899

Reference: Akemann et al. (2009) 96:3959-3976

Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Date of Implementation: April 2007
Contact: akemann@brain.riken.jp

ENDCOMMENT

NEURON {
	SUFFIX Cav2
	USEION ca READ cai, cao WRITE ica
	NONSPECIFIC_CURRENT i
	RANGE pbar, ica, i, igate, nc
	GLOBAL minf, taum
	GLOBAL gateCurrent, punit 
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {
	e0 = 1.60217646e-19 (coulombs)
	q10 = 2.7
	F = 9.6485e4 (coulombs)
	R = 8.3145 (joule/kelvin)

	cv = 19 (mV)
	ck = 5.5 (mV)

	zm = 4.6244 (1)				: gating charge
}

PARAMETER {
	gateCurrent = 0 (1)			: gating currents ON = 1 OFF = 0
	
	pbar = 3e-5 (cm/s)
	punit = 3.290e-13 (cm3/s)		: unitary calcium permeability
	
	monovalConc = 140 (mM)
	monovalPerm = 0 (1)
}

ASSIGNED {
	celsius (degC)
	v (mV)
	cai (mM)
	cao (mM)

	ica (mA/cm2)
	i (mA/cm2)
	igate (mA/cm2)
	
	nc (1/cm2)				: membrane density of channel

      minf (1) 
	taum (ms)
	T (kelvin)
	E (volt)
	zeta (1)
	qt (1)
}

STATE { m }

INITIAL {
	nc = pbar / punit
	qt = q10^((celsius-22 (degC))/10 (degC))
	T = kelvinfkt( celsius )
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = (1e3) * pbar * m * ghk(v, cai, cao, 2)
	igate = nc * (1e6) * e0 * zm * mgateFlip()

	if (gateCurrent != 0) { 
		i = igate
	}
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
}

FUNCTION ghk( v (mV), ci (mM), co (mM), z )  (coulombs/cm3) { 
	E = (1e-3) * v
      zeta = (z*F*E)/(R*T)	
	
	: ci = ci + (monovalPerm) * (monovalConc) :Monovalent permeability

	if ( fabs(1-exp(-zeta)) < 1e-6 ) {
	ghk = (1e-6) * (z*F) * (ci - co*exp(-zeta)) * (1 + zeta/2)
	} else {
	ghk = (1e-6) * (z*zeta*F) * (ci - co*exp(-zeta)) / (1-exp(-zeta))
	}
}

PROCEDURE rates( v (mV) ) {
	minf = 1 / ( 1 + exp(-(v+cv)/ck) )
	taum = (1e3) * taumfkt(v)/qt
}

FUNCTION taumfkt( v (mV) ) (s) {
	UNITSOFF
	if ( v > -50 ) {
	taumfkt = 0.000191 + 0.00376 * exp(-((v+41.9)/27.8)^2)
	} else {
	taumfkt = 0.00026367 + 0.1278 * exp(0.10327*v)
	}
	UNITSON
}

FUNCTION kelvinfkt( t (degC) )  (kelvin) {
	UNITSOFF
	kelvinfkt = 273.19 + t
	UNITSON
}

FUNCTION mgateFlip() (1/ms) {
	mgateFlip = (minf-m)/taum
}
