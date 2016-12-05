NEURON {
    SUFFIX calRT03 }
NEURON {
    USEION ca READ eca WRITE ica }
ASSIGNED {
    ica
    eca (mV)
}

PARAMETER {
	:erev 		= 125    (mV)
	gmax 		= 1.0    (mho/cm^2)
        vrest           = 0

	maflag 		= 2
	malphaA 	= 1.6
	malphaB		= -13.889
	malphaV0	= 5.
	mbflag 		= 3
	mbetaA 		= 0.02
	mbetaB		= 5.
	mbetaV0		= -8.9
	exptemp		= 37
	mq10		= 1
	mexp 		= 2

	haflag 		= 0
	halphaA 	= 0
	halphaB		= 0
	halphaV0	= 0
	hbflag 		= 0
	hbetaA 		= 0
	hbetaB		= 0
	hbetaV0		= 0
	hq10		= 3
	hexp 		= 0
}

INCLUDE "custom_code/inc_files/64229_geneval_cvode.inc"

PROCEDURE iassign () { i=g*(v-eca) ica=i }

