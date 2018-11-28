: Rod  Photoreceptor Ca and Calcium  channel
: Ref. Kourenny and  Liu 2002   ABME 30 : 1196-1203
: Modification 2004-02-07
NEURON 
{
	SUFFIX Ca
	
	:USEION Ca WRITE iCa VALENCE 2
        USEION ca READ eca WRITE ica
        RANGE gCabar,VhalfCam,SCam
        RANGE VhalfCah,SCah
        RANGE aomCa,bomCa
        RANGE gammaohCa,deltaohCa


}

UNITS
{
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millimho)
	(mol)= (1)
	(M)  = (mol/liter)
	(uM) = (micro M)
}

PARAMETER
{
       
       : Calcium channel 
       gCabar = 2 (mS/cm2) <0,1e9> :different from ABME paper
       :eCa =  40 (mV)
       aomCa = 50  (/s)  : changed from 3.10/s, 20/s
       bomCa = 50  (/s)
       gammaohCa = 1 (/s)
       deltaohCa =1 (/s)  
 
       VhalfCam=-20.0 (mV)
       VhalfCah=10 (mV)
       SCam =6.0      (mV) 
       
       SCah =9        (mV)   
     
}


STATE
{

	mCa
	hCa
	
}

ASSIGNED
{
        ica (mA/cm2)
        eca (mV)
	gCa (mho/cm2)
    
	v (mV)
	
	:iCa (mA/cm2)

	infmCa
	taumCa  (ms) 
	


	infhCa
	tauhCa (ms)



}

INITIAL
{
	rate(v)
	mCa = infmCa
	hCa = infhCa

}




BREAKPOINT
{
	SOLVE states METHOD cnexp
	gCa = (0.001)*gCabar*mCa*hCa
	: g is in unit of S/cm2 ,i is in unit of mA/cm2 and v is in mV
	
	ica = gCa*(v - eca)
	: the current is in the unit of mA/cm2
	
	
}

DERIVATIVE states
{
	rate(v)
	mCa' = (infmCa - mCa)/taumCa
	hCa'= (infhCa-hCa)/tauhCa


}




FUNCTION alphamCa(v(mV))(/ms)
{ 
	alphamCa = 0.001*aomCa*exp( (v - VhalfCam)/(2*SCam)   )
}

FUNCTION betamCa(v(mV))(/ms)
{ 
	betamCa = 0.001*bomCa*exp( - ( v-VhalfCam)/(2*SCam) )
}
FUNCTION gammahCa(v(mV))(/ms)
{ 
	gammahCa = 0.001*gammaohCa*exp( (v - VhalfCah)/(2*SCah))
}

FUNCTION deltahCa(v(mV))(/ms)
{ 
	deltahCa = 0.001*deltaohCa*exp( - ( v-VhalfCah)/(2*SCah) )
}


PROCEDURE rate(v (mV))
{
        LOCAL a, b,c, d


	a = alphamCa(v)
	b = betamCa(v)
	taumCa = 1/(a + b)
	infmCa = a/(a + b)
	
	c = gammahCa(v)
	d = deltahCa(v)
	tauhCa = 1/(c + d)
	infhCa = d/(c + d)

}

