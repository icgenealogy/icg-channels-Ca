TITLE CaL channel

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX CaL
    USEION ca READ eca WRITE ica
    RANGE gbar
}

PARAMETER {
    gbar = 0.0025 (mho/cm2)
} 

ASSIGNED {
    v (mV)
    eca (mV)
    ica (mA/cm2)
}

BREAKPOINT {
    LOCAL den, m
    den = 1 + exp((v+20)/-9)
    m = 1/den
    ica  = gbar*m*m*(v-eca)
}



