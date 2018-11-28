TITLE Ca R-type channel with high threshold for activation

: HVA calcium channels are inserted in the spine head
: Activation and inactivation parameters taken from
: Foehring RC, Mermelstein PG, Song W, Ulrich S and Surmeier DJ
: Unique properities of R-type calcium currents in neucortical and neostriatal neurons
: J Neurophysiol (2000) 84: 2225 - 2236
:
: written by Lei Tian on 04/11/06 

NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica
    RANGE gcabar, m, h, g, p, eca
	RANGE inf, fac, tau, k
	GLOBAL irtype
	:EXTERNAL Area_canmda
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {	: parameters that can be entered when function is called in cell-setup
    v               (mV)
    celsius = 30	(degC)
	dt              (ms)
    gcabar = 0.351  (mho/cm2) : initialized conductance 
	:eca = 10		(mV)      : Ca++ reversal potential was choosen to best fit the GHK between -40 and -10 mV	

	Area = 1.11e-8 (cm2)           
	k = 1e-06		(mA/nA)

        }  

STATE {	m h }               

ASSIGNED {
        eca (mV)                  
	ica             (mA/cm2)
    inf[2]
	fac[2]
	tau[2]
	irtype
	g                       :R_type channel total conductance
	p
	
}

BREAKPOINT {
	SOLVE states
	ica = gcabar*m*m*m*h*(v - eca)
	irtype= -gcabar*m*m*m*h*(v - eca)
	g = gcabar*m*m*m*h*Area*1e6	:[uS]
	p = m*m*m*h
	}

INITIAL {
	:Area = Area_canmda
    m = 0                               : initial activation parameter value
	h = 0.5                             : initial inactivation parameter value
	states()
	ica = gcabar*m*m*m*h*(v - eca)      : initial Ca++ current value
    irtype=-gcabar*m*m*m*h*(v - eca) 	: the ca current through R_type channel
	g = gcabar*m*m*m*h*Area*1e6 		:[uS]
	p = m*m*m*h
	}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	}	

PROCEDURE states() {                    : exact when v held constant
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION varss(v, i) {
	if (i==0) {
           varss = 1 / (1 + exp((v+14)/(-6.7)))	: Ca activation
	}
	else if (i==1) {    
        varss = 1/ (1 + exp((v+65)/(11.8)))     : Ca inactivation
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = 3.6		: activation variable time constant 
        }
	else if (i==1) {
           vartau = 200		: inactivation variable time constant 
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b :rest = -70
:	TABLE inf, fac DEPEND dt, celsius FROM -100 TO 100 WITH 200
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}


