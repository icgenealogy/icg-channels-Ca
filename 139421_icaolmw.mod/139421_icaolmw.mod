: $Id: icaolmw.mod,v 1.2 2010/11/30 16:44:13 samn Exp $ 
COMMENT

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
//
// Copyright 2007, The University Of Pennsylvania
// 	School of Engineering & Applied Science.
//   All rights reserved.
//   For research use only; commercial use prohibited.
//   Distribution without permission of Maciej T. Lazarewicz not permitted.
//   mlazarew@seas.upenn.edu
//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENDCOMMENT

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
}

NEURON {
  SUFFIX ICaolmw
  USEION ca READ eca WRITE ica
  RANGE gca
}

PARAMETER {
  gca = 1    (mS/cm2)
  :eca = 120  (mV)
}
    
ASSIGNED { 
  eca (mV)
  ica (mA/cm2)    
  v   (mV)
}

PROCEDURE iassign () { ica = (1e-3) * gca * mcainf(v)^2 * (v-eca) }

INITIAL {
  iassign()
}

BREAKPOINT { iassign() }

FUNCTION mcainf(v(mV)) { mcainf = fun2(v, -20, 1,  -9)*1(ms) }

INCLUDE "custom_code/inc_files/139421_aux_fun.inc"
