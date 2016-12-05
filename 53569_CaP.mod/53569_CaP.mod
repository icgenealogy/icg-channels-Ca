TITLE persistent calcium current

COMMENT
12/1/2005 NTC Made compatible with adaptive integration
Unused stuff removed
ENDCOMMENT

:   Modified by Steven Prescott based on current described below
:   Prescott and De Koninck. 2005. J Neurosci 25: 4743-4754
:   low threshold but persistent calcium current, 
:   acts synergistically with persistent sodium current to prolong subthreshold depolarization
:   in tonic spinal lamina I neurons
:   The suffix CaP stands for persistent calcium.  It does not necessarily represent a P-type current
:
:   original current described below...
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by standard equations (NOT GHK)
:   using a m2h format, according to the voltage-clamp data
:   (whole cell patch clamp) of Huguenard & Prince, J Neurosci.
:   12: 3804-3817, 1992.  The model was introduced in Destexhe et al.
:   J. Neurophysiology 72: 803-818, 1994.
:
:    - Kinetics adapted to fit the T-channel of reticular neuron
:    - Q10 changed to 5 and 3
:    - Time constant tau_h fitted from experimental data
:    - shift parameter for screening charge
:
:   ACTIVATION FUNCTIONS FROM EXPERIMENTS (NO CORRECTION)
:
:   Reversal potential taken from Nernst Equation
:
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:
:   Modifications by Arthur Houweling for use in MyFirstNEURON

NEURON {
	SUFFIX CaP
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, shift
	RANGE delta
	RANGE ica
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v			(mV)
	celsius		(degC)
	gcabar = .00002 	(mho/cm2)
	shift	= 2 		(mV)		: screening charge for Ca_o = 2 mM
	cai			(mM)		
	cao			(mM)
	delta = 60		(mV)		: increases tau_m
}

STATE {
	m 
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	phi_m
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

 	m'= (m_inf-m) / tau_m
}

UNITSOFF
INITIAL {
:
:   Activation functions and kinetics were obtained from
:   Huguenard & Prince, and were at 23-25 deg.
:   Transformation to 36 deg assuming Q10 of 5 and 3 for m and h
:   (as in Coulter et al., J Physiol 414: 587, 1989)
:
	phi_m = 5.0 ^ ((celsius-24)/10)
	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) { 
:
:   Time constants were obtained from J. Huguenard
:
	m_inf = 1.0/(1+exp(-(v+shift+50)/(7.4)))
	tau_m = ( 3+delta+1.0 /(exp((v+shift+25)/10)+exp(-(v+shift+100)/15)))/phi_m
}
UNITSON
