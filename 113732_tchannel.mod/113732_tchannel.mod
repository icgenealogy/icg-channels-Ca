INDEPENDENT {t FROM 0 TO 1 WITH 1 (MS)}

NEURON {
     SUFFIX tchannel
     USEION ca READ cai, cao WRITE ica
     RANGE pbar
     RANGE tt_inf, u_inf
     RANGE tau_tt, tau_u
     RANGE tt_exp, u_exp
     RANGE itt, itmax

     RANGE pcabar,qt, vshift, GHKa, GHKb,GHKc, GHKo, it

}

UNITS {
     (molar)=(1/liter)
     (mM) = (millimolar)
     (mA) = (milliamp)
     (mV) = (millivolt)
     (nA) = (nanoamp)

     FARADAY = (faraday)(coulombs)
     R = (k-mole)(joule/degC)
}

PARAMETER {
     celsius (degC)
     pbar = .000017 (cm/s) <0,1e9>
     pcabar = .0002 (cm/s) <0,1e9>
     vshift (mV)
     v (mV)
     dt (ms)
     qt=1
     cai (mM)
     cao (mM)
     GHKo
     GHKa
     GHKb
     GHKc

}

STATE {
     tt u
  

}

INITIAL {
     tt = .0003
     u = .9

}

ASSIGNED {
     ica (nA/cm2)
     itt (nA/cm2)
     itmax (nA/cm2)
     it (nA/cm2)

     tt_inf u_inf
     tau_tt tau_u
     tt_exp u_exp

}

BREAKPOINT {
     SOLVE states
    



     
     itt = (tt*tt*tt*u*pbar*ghk(v,cai, cao)) 
     itmax = (pbar*ghk(v,cai, cao))/10^6 
     ica = itt



:    GHKa = (4*(FARADAY^2)*(v*.001))*(1/(R*(celsius+273.18)))

:    GHKb = ((cai*(.001))-(cao*(.001)))*exp((-2*FARADAY*(v*.001))/(R*(celsius+273.18)))
     
:    GHKc = 1/(1-exp((-2*FARADAY*(v*.001))/(R*(celsius+273.18))))
     
     
:    GHKo = GHKa*GHKb*GHKc

:    it=tt*tt*tt*u*pbar*GHKo

:    ica = it
}

FUNCTION ghk(v(mV), ci(mM), co(mM))(.001 coul/cm3){
     LOCAL z, eci, eco
     z = 1e-3*2*FARADAY*v/(R*(celsius+273.18))
      
     eco = co*efun(z)
     eci = ci*efun(-z)
     ghk = (.001)*2*FARADAY*(eci-eco)
     
}
FUNCTION efun(z) {
     if(fabs(z)<1e-4) {
     efun = 1-z/2
}else{
     efun = z/(exp(z)-1)
     }
}

PROCEDURE states() { :exact when v held constant
     evaluate_fct(v)
     tt = tt + tt_exp*(tt_inf-tt)
     u = u + u_exp*(u_inf-u)

VERBATIM
return 0;
ENDVERBATIM

}

UNITSOFF

PROCEDURE evaluate_fct(v(mV)) {LOCAL a,b,y,d


:LVA calcium channel
:        qt=1
qt=2.9^((celsius-18)/10)
vshift = (15-cao)*(-.757575)


: activation
tt_inf =( 1.0/(1+exp(((-56.4+vshift)-(v))/13.2)))^3
	if(v<-60+vshift)
	{
        tau_tt = (44.2 + .8014(v+vshift) +.0049(v+vshift)^2+((9.7*(10^-6))*( v+vshift)^3))*(1/qt)
	}
	else{
        tau_tt = (.5187 + .01550(v+vshift) +.0013(v+vshift)^2+((-.4064*10^-6)*(v+vshift)^3))*(1/qt)
        }
:inactivation - u
 u_inf = 1.0/(1+exp(((-68.2+vshift)-(v))/-7.63))
	if(v<=-55+vshift)
	{
	tau_u =((33.08 +(346.6-33.08)/(1+exp((-82.2-(v+vshift))/-3.362))))*(1/qt)       
	}
	else{
	tau_u = (((.01149*exp(-.1348*(v-vshift)))+(31.42*exp(-.0007572*(v-vshift)))+(-13.38)))*(1/qt)
        }

:states vars to infinity

	tt_exp = 1-exp(-dt/tau_tt)
	u_exp = 1-exp(-dt/tau_u)
}
UNITSON
