!***********************************************************************
!
!    MODULE:       EL_EOS.INC
!    TYPE:         INCLUDE FILE
!    AUTHOR:       F. DOUGLAS SWESTY
!    DATE:         2/12/91
!
!
!    CALL LINE:    INCLUDE 'EL_EOS.INC'
!
!***********************************************************************
!
!
!                              Electron density and concentration
      DOUBLE PRECISION NSUBE, QSUBE
!
!                              Temporary coefficients
      DOUBLE PRECISION ACOEF, BCOEF, CCOEF, DBDT
!
!                              Electron chemical potential
      DOUBLE PRECISION MUSUBE
!
!                              Electron pressure, internal energy,
!                              Helmholtz free energy, and entropy
      DOUBLE PRECISION EPRESS, EU, FSUBE, ES
!
!                              Photon pressure, internal energy,
!                              Helmholtz free energy, and entropy
      DOUBLE PRECISION PPRESS, PS, PU, PF
!
!
      DOUBLE PRECISION DEMUDT, DEMUDN, DEMUDY
      DOUBLE PRECISION DEPDT, DEPDN, DEPDY
      DOUBLE PRECISION DESDT, DESDN, DESDY
      DOUBLE PRECISION DEUDT, DEUDN, DEUDY
!
!
      DOUBLE PRECISION DPPDN, DPPDT, DPPDY
      DOUBLE PRECISION DPSDN, DPSDT, DPSDY
      DOUBLE PRECISION DPUDN, DPUDT, DPUDY
!
!
!                              This common block contains all the
!                              electron and phton variables needed
!                              by the outside world
      COMMON /ELECTR/ NSUBE, MUSUBE, 
     1  EPRESS, EU, FSUBE, ES, PPRESS, PS, PU, PF,
     2  DEMUDT, DEMUDN, DEMUDY, DEPDT, DEPDN, DEPDY,
     3  DESDT, DESDN, DESDY, DEUDT, DEUDN, DEUDY,
     4  DPPDN, DPPDT, DPPDY, DPSDN, DPSDT, DPSDY,
     5  DPUDN, DPUDT, DPUDY
!

!
!