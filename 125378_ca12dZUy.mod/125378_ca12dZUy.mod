COMMENT

T-type Ca2+ channel
	large scaled Morkov model for T-type ca2+ current
		C0<->C1<->C2<->C3<->C4<->O
		I0<->I1<->I2<->I3<->I4<->Io
		C0<->I0,C1<->I1,C2<->I2,C3<->I4,O<->Io

Reference:
	Serrano JR, Perez-Reyes E,& Jones SW, 1999, J Gen Physiol
Written by:
	Yuki HAYASHIDA, 2001 Sep-Oct
   	

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ca12dZUy
	USEION ca READ cai, cao WRITE ica
	RANGE C0, C1, C2, C3, C4, O, I0, I1, I2, I3, I4, Io
	RANGE kv, kv0, k_v, k_v0, ko, k_o, ki, k_i
	RANGE f, h, p
	RANGE vkv, vk_v
	GLOBAL vmin, vmax, vshift
}

UNITS {
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
} 

PARAMETER {
	p    = 7e-6  	(cm/s)		: max permeability
	v 		(mV)

	: max rates
	kv0 = 3.7	(/ms)
	k_v0 = 18e-3	(/ms)
	ko = 20		(/ms)
	k_o = 1.65	(/ms)
	
	ki = 62e-3	(/ms)
	k_i = 0.13e-3	(/ms)

	f = 0.224
	h = 0.125

	: voltage dependence
	vkv = 25.5	(mV)
	vk_v = -15.3	(mV)

	vshift = 0	(mV)

	celsius		(degC)

	vmin = -200	(mV)
	vmax = 200	(mV)
} 


ASSIGNED {
	ica 		(mA/cm2)
	cao		(mM)
	cai		(mM)
	kv		(1/ms)
	k_v		(1/ms)
}


STATE { C0 C1 C2 C3 C4 O I0 I1 I2 I3 I4 Io }

INITIAL { 
	:at Vh=-62mV
	C0 = 0.0076057097
	C1 = 0.0095551048
	C2 = 0.0045015537
	C3 = 0.00094255543
	C4 = 7.4008676e-05
	O = 0.00089707485
	I0 = 7.9627334e-05
	I1 = 0.0035727272
	I2 = 0.060113056
	I3 = 0.44952644
	I4 = 0.035296445
	Io = (1-C0-C1-C2-C3-C4-O-I0-I1-I2-I3-I4)
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse
	ica = O * p * ghk(v,cai,cao)
} 


KINETIC kstates {
	~ C0 <-> C1	(4*kv,k_v)
	~ C1 <-> C2	(3*kv,2*k_v)	     
	~ C2 <-> C3	(2*kv,3*k_v)
	~ C3 <-> C4	(kv,4*k_v)
	~ C4 <-> O	(ko,k_o)
	~ I0 <-> I1	(4*kv/f,h*k_v)
	~ I1 <-> I2	(3*kv/f,2*h*k_v)
	~ I2 <-> I3	(2*kv/f,3*h*k_v)
	~ I3 <-> I4	(kv,4*k_v)
	~ I4 <-> Io	(ko,k_o)
	~ C0 <-> I0	(f*f*f*ki,k_i/h/h/h)
	~ C1 <-> I1	(f*f*ki,k_i/h/h)
	~ C2 <-> I2	(f*ki,k_i/h)
	~ C3 <-> I3	(ki,k_i)
	~ C4 <-> I4	(ki,k_i)
	~ O <-> Io	(ki,k_i)

	CONSERVE C0+C1+C2+C3+C4+O+I0+I1+I2+I3+I4+Io = 1
}	
	
PROCEDURE rates(v(mV)) {
	TABLE kv, k_v
	DEPEND kv0, k_v0, vkv, vk_v
	FROM vmin TO vmax WITH 401

	kv = kv0 * exp((v-vshift)/vkv)
	k_v = k_v0 * exp((v-vshift)/vk_v)
}

: Special gear for calculating the Ca2+ reversal potential
: via Goldman-Hodgkin-Katz eqn.
: [Ca2+]o "cao" and [Ca2+]i "cai" are assumed to be set elsewhere
FUNCTION ghk(v(mV), ci(mM), co(mM)) (0.001 coul/cm3) {
	LOCAL z

	z = (0.001)*2*F*v/(R*(celsius+273.15))
	ghk = (0.001)*2*F*(ci*efun(-z) - co*efun(z))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
