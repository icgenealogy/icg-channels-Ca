TITLE HVA calcium currents (CaHVA-current) 
 
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
        SUFFIX caHVA
        USEION ca READ eca WRITE ica
        RANGE gcaHVAbar, icaHVA, ica 
		RANGE dHVAinf,fHVAinf
		RANGE dHVAtau, fHVAtau
		RANGE vhalfAct,slopeAct,vhalfInact,slopeInact
		RANGE vhalfTact,slopeTact,vhalfTinact,slopeTinact
}
 

PARAMETER {
        v (mV)
        dt (ms)
        gcaHVAbar =  0.00004 (S/cm2)
        eca = 120 (mV)
		vhalfAct = -22(mV)
		vhalfInact = -40.0 (mV)
		slopeAct = 5.0
		slopeInact = -7.0
		vhalfTact = -40.0(mV)
		vhalfTinact = -39 (mV)
		slopeTact = -3
		slopeTinact = -2.6	
}
 
STATE {
        dHVA fHVA  
}
 
ASSIGNED {
        ica (mA/cm2)
		icaHVA (mA/cm2)
        dHVAinf 
		fHVAinf 
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        icaHVA = gcaHVAbar*dHVA*fHVA*(v - eca)
		ica = icaHVA
}
 
UNITSOFF
 
INITIAL {
        dHVA = dHVAinf
        fHVA = fHVAinf
}

DERIVATIVE states {  
LOCAL dHVAinf,dHVAtau,fHVAinf,fHVAtau
        dHVAinf = 1/(1 + exp(-(v - vhalfAct)/slopeAct))
		fHVAinf = 1/(1 + exp(-(v - vhalfInact)/slopeInact))
        dHVAtau = 1/(1 + exp(-(v - vhalfTact)/slopeTact)) + 0.78
		fHVAtau = 77/(1 + exp(-(v - vhalfTinact)/slopeTinact)) + 17
        dHVA' = (dHVAinf-dHVA)/dHVAtau
        fHVA' = (fHVAinf-fHVA)/fHVAtau

}
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}

UNITSON

