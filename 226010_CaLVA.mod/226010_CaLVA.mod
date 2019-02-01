TITLE LVA calcium current (CaLVA-current) 
 
COMMENT
written for NEURON by Antonios Dougalis, 23 Feb 2015, London, UK
based on voltage clamp data from Dougalis et al., 2017 J Comput Neurosci 
ENDCOMMENT
 
UNITS {
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
 }
 
NEURON {
        SUFFIX caLVA
        USEION ca READ eca WRITE ica
        RANGE gcaLVAbar, icaLVA, ica
        RANGE dLVAinf, fLVAinf
		RANGE dLVAtau, fLVAtau 
		RANGE vhalfAct,slopeAct,vhalfInact,slopeInact
		RANGE vhalfTact,slopeTact,vhalfTInact,slopeTInact
		
}
 

 
PARAMETER {
        v (mV)
        dt (ms)
        gcaLVAbar =  0.00004 (S/cm2)
        eca = 120 (mV)
		vhalfAct = -57.5 (mV)
		vhalfInact = -83.0 (mV)
		slopeAct = 6.5
		slopeInact = -6.1
		vhalfTact = -68.8 (mV)
		slopeTact = -5.0
		vhalfTInact = -24.6 (mV)
		slopeTInact = -8.6
        
}
 
STATE {
        dLVA 
		fLVA
}
 
ASSIGNED {
        ica (mA/cm2)
        icaLVA (mA/cm2)
        dLVAinf (1) 
		fLVAinf (1)
		dLVAtau (ms)
		fLVAtau (ms)
		
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        icaLVA = gcaLVAbar*dLVA*fLVA*(v - eca)
        ica = icaLVA
}
 
UNITSOFF
 
INITIAL {
        dLVA = dLVAinf
        fLVA = fLVAinf
		
}

DERIVATIVE states { 
LOCAL dLVAinf,dLVAtau,fLVAinf,fLVAtau
        dLVAinf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
		fLVAinf = 1/(1 + exp(-(v - vhalfInact)/slopeInact))
        dLVAtau = 4.45/(1 + exp(-(v - vhalfTact)/slopeTact)) + 0.95
		fLVAtau = 5/(1 + exp(-(v - vhalfTInact)/slopeTInact)) + 15
        dLVA' = (dLVAinf-dLVA)/dLVAtau
        fLVA' = (fLVAinf-fLVA)/fLVAtau
		
}
 
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}

 
UNITSON

