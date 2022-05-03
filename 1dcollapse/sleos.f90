!***********************************************************************
!                                                                       
!    FILE:         INVEOS                                               
!    MODULE:       INVEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         5/23/90                                              
!                  Bug fixed on (5/24/90)                               
!                                                                       
!                                                                       
!    CALL LINE:    CALL INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,XPREV
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  T_OLD = INITIAL GUESS AT THE TEMPERATURE             
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                  IFLAG = 1 --> INPVAR IS TEMPERATURE                  
!                          2 --> INPVAR IS INTERNAL ENERGY              
!                          3 --> INPVAR IS ENTROPY (NOT IMPLEM)         
!                                                                       
!    OUTPUTS       EOSFLG = 1 --> "NO NUCLEI" EOS                       
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS FOR DENSITIES ABOVE NUCLEAR  
!                  XPREV = PREVIOUS VALUE OF X                          
!                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY            
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M1D.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!
      SUBROUTINE INVEOS(INPVAR,T_OLD,YE,BRYDNS,IFLAG,EOSFLG,            &
     &                  FORFLG,SF,XPREV,P_PREV)                         
!                                                                       
!                                                                       
!                                                                       
                                                                        
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
!                         Local variables                               
      INCLUDE 'eos_m4a.inc' 
!                                                                       
      DOUBLE PRECISION INP_V, INP_VO, INP_VN, UFTN, DUFTN, DT 
      DOUBLE PRECISION T_OLD, T_NEW, T_TEMP, T_LB, T_UB, PERDIF 
      DOUBLE PRECISION eps, prec, var, t_old0, al, fa, bl, fb, cl,      &
     &                 fc, d, e, tm, tol1, s, p, ql, r                  
      INTEGER LOOP, SF, NEW_F 
      INTEGER itmax, IFLAG0 
!                                                                       
      data itmax/300/ , eps/1.d-20/ 
      data prec/1.d-14/ 
!                                                                       
      RSFLAG = 1 
!                         Input is the temperature; call the EOS        
!                         normally and then return                      
      IF(IFLAG.EQ.1) THEN 
        CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,               &
     & XPREV,P_PREV)                                                    
        T_OLD = INPVAR(1) 
        RETURN 
      ENDIF 
!                                                                       
!--the calling variable is either internal energy or entropy            
!                                                                       
      IFLAG0=IFLAG 
      var=INPVAR(1) 
      t_old0=T_OLD 
!                                                                       
!--find lower bound in T                                                
!                                                                       
      INPVAR(1)=0.5*t_old0 
      al=INPVAR(1) 
      IFLAG=1 
      CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,                 &
     &             XPREV,P_PREV)                                        
      IF(IFLAG0.eq.2)THEN 
         fa=UTOT - var 
      ELSE 
         fa=STOT - var 
      END IF 
!                                                                       
!--find upper bound in T                                                
!                                                                       
      INPVAR(1)=2.*t_old0 
      bl=INPVAR(1) 
      IFLAG=1 
      CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,                 &
     &             XPREV,P_PREV)                                        
      IF(IFLAG0.eq.2)THEN 
         fb=UTOT - var 
      ELSE 
         fb=STOT - var 
      END IF 
!                                                                       
!--root has been bracketed, iterate to find solution                    
!                                                                       
   10 continue 
      fc=fb 
      do i=1,itmax 
!                                                                       
!--rename a,b,c and adjust bounding interval                            
!                                                                       
         if(fb*fc.gt.0.) then 
            cl=al 
            fc=fa 
            d=bl-al 
            e=d 
         end if 
         if(dabs(fc).lt.dabs(fb)) then 
            al=bl 
            bl=cl 
            cl=al 
            fa=fb 
            fb=fc 
            fc=fa 
         end if 
!                                                                       
!--check for convergence                                                
!                                                                       
         tm=0.5*(cl-bl) 
         tol1=2.*eps*dabs(bl) 
         if(dabs(tm).lt.tol1.or.dabs(fb/var).lt.prec) then 
            T_OLD=bl 
            return 
         end if 
!                                                                       
!--attempt inverse quadratic interpolation                              
!                                                                       
         if(dabs(e).ge.tol1.and.dabs(fa).gt.dabs(fb)) then 
            s=fb/fa 
            if(al.eq.cl) then 
               p=2.*tm*s 
               ql=1.-s 
            else 
               ql=fa/fc 
               r=fb/fc 
               p=s*(2.*tm*ql*(ql-r)-(bl-al)*(r-1.)) 
               ql=(ql-1.)*(r-1.)*(s-1.) 
            end if 
!                                                                       
!--check wether in bound                                                
!                                                                       
            if(p.gt.0.) ql=-ql 
            p=dabs(p) 
!                                                                       
!--accept or refuse interpolation                                       
!                                                                       
            if(2.*p.lt.dmin1(3.*tm*ql-dabs(tol1*ql),dabs(e*ql))) then 
!                                                                       
!--accept interpolation                                                 
!                                                                       
               e=d 
               d=p/ql 
            else 
!                                                                       
!--interpolation failed use bisection                                   
!                                                                       
               d=tm 
               e=d 
            end if 
!                                                                       
!--bound decreasing to slowly use bisection                             
!                                                                       
         else 
            d=tm 
            e=d 
         end if 
!                                                                       
!--move last guess to a                                                 
!                                                                       
         al=bl 
         fa=fb 
!                                                                       
!--evalue new trial point                                               
!                                                                       
         if(dabs(d).gt.tol1) then 
            bl=bl + d 
         else 
            bl=bl + dsign(tol1,tm) 
         end if 
         INPVAR(1)=bl 
         IFLAG=1 
         CALL EOS_M4A(INPVAR,YE,BRYDNS,1,EOSFLG,FORFLG,SF,              &
     &                XPREV,P_PREV)                                     
         IF(IFLAG0.eq.2)THEN 
            fb=UTOT - var 
         ELSE 
            fb=STOT - var 
         END IF 
!                                                                       
      enddo 
!                                                                       
!--did not converge write(iprint,error message                          
!                                                                       
      write(*,*)'iterations did not converge!' 
!                                                                       
      return 
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EOS4B.FOR                                            
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EOS_M4B                                              
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/15/90 Modified from model 4A to include the       
!                  phase boundary cutoffs and Maxwell construction      
!                  boundaries.                                          
!                  7/13/90 Modified from model 1-d to include Maxwell   
!                  construction                                         
!                  5/25/90  MODEL 1D   (RELEASE # 1.2)                  
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL EOS_M4A(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,    
!                  XPREV,P_PREV)                                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                  IFLAG = 1 --> INPVAR IS TEMPERATURE                  
!                          2 --> INPVAR IS INTERNAL ENERGY (NOT IMPLEM) 
!                          3 --> INPVAR IS ENTROPY (NOT IMPLEMENTED)    
!                  FFLAG = "FORCING FLAG"  0 --> NO FORCING             
!                                          1 --> FORCE A PARTICULAR     
!                                                SCHEME TO BE USED      
!                                                                       
!                                                                       
!    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B           
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS (includes alpha's)           
!                  XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON     
!                          FIRST CALL)                                  
!                  P_PREV = PREVIOUS VALUE OF PROTON DENSITY (MUST BE   
!                          SUPPLIED ON FIRST CALL)                      
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EOS_M4A(INPVAR,YE,BRYDNS,IFLAG,EOSFLG,FFLAG,SSFLAG,    &
     &                   XPREV,P_PREV)                                  
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION OUTVAR(4) 
!                                                                       
!                                                                       
!                       This include file contains all variable         
!                       declarations.  NOTE:: no implicit typing        
!                       scheme is followed in this code; if you         
!                       have any doubt as to a variables type CHECK     
!                       IT!!!!.  Also note that ALL variables are       
!                       declared explicitly.                            
!                                                                       
      INCLUDE 'eos_m4a.inc' 
!                                                                       
!                                                                       
!                         Set the "switch" flag to zero                 
      SWTFLG = 0 
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options should go through 
!                         INVEOS untill further notice)                 
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                         If the "forcing" flag is set then skip        
!                         the EOS determination logic and go straight   
!                         to the EOS determined by EOSFLG               
      IF(FFLAG.EQ.1) THEN 
        GOTO 10 
      ELSE 
!                         Otherwise let the EOS logic module determine  
!                         the correct EOS to use                        
        CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG) 
      ENDIF 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                        Try NUCEOS first and if not successfull        
!                        then try bulk EOS                              
   10 CONTINUE 
      IF(EOSFLG.EQ.1) THEN 
!                                                                       
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR) 
          OUTVAR(1) = INPVAR(1) 
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
          IF(SSFLAG.NE.1) THEN 
            OUTVAR(2) = 0.155 
            OUTVAR(3) = -15.0 
            OUTVAR(4) = -20.0 
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
          ENDIF 
!                                                                       
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN 
!                    Set EOS flag to full scheme                        
          EOSFLG = 2 
!                                                                       
!                    Else if fraction of nuclei is less than the minimum
!                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE 
          IF(FFLAG.NE.1) THEN 
!                                                                       
            CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
            IF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN 
              EOSFLG = 1 
              WRITE(*,*) 'A2 failed at try = ',T,BRYDNS,YE 
              GOTO 999 
            ENDIF 
!                                                                       
!                    Set nuclei to bulk EOS                             
            EOSFLG = 3 
!                    Save value of proton fraction                      
            P_PREV = YE*BRYDNS 
!                                                                       
            GOTO 999 
!                                                                       
          ELSE 
            IF(NF_FLG.EQ.1)                                             &
     &          WRITE(*,*) 'NUC failed at t,rho = ',t,brydns            
            GOTO 999 
          ENDIF 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!                                                                       
!          End of NUCEOS--BULK EOS calculations                         
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                            CALCULATE FULL EOS (INCLUDING NUCLEI)      
      IF(EOSFLG.EQ.2) THEN 
!                                                                       
!                    Call the nuclear EOS                               
        CALL NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
        IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
!ccc          WRITE(*,*) ' EOS_M4A:: r.i.gs.'                           
          CALL RESET(INPVAR,YE,BRYDNS,OUTVAR) 
          OUTVAR(1) = INPVAR(1) 
          CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
          IF(SSFLAG.NE.1) THEN 
            OUTVAR(2) = 0.155 
            OUTVAR(3) = -15.0 
            OUTVAR(4) = -20.0 
            CALL NUCEOS(OUTVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
          ENDIF 
!                                                                       
!                                                                       
!                                                                       
          IF(SSFLAG.NE.1) THEN 
!ccc            WRITE(*,*) '     r.i.gs. failure @ try: ',inpvar        
            GOTO 999 
          ELSE 
            INPVAR(2) = OUTVAR(2) 
            INPVAR(3) = OUTVAR(3) 
            INPVAR(4) = OUTVAR(4) 
          ENDIF 
!                    Otherwise quit and return                          
        ELSEIF((SSFLAG.NE.1).AND.(FFLAG.EQ.1)) THEN 
          GOTO 999 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    If fraction of heavies is greater than the minimum 
!                    parameter, then this EOS is OK                     
        IF((XH.GT.HEAVCT).AND.(SSFLAG.EQ.1)) THEN 
!                    Set EOS flag to full scheme                        
          EOSFLG = 2 
!                                                                       
!                    Else if fraction of nuclei is less than the minimum
!                    or if NUCEOS was unsuccessful use the no nuclei EOS
        ELSE 
!                    If the forcing flag is not set                     
          IF(FFLAG.NE.1) THEN 
!                    Set nuclei to no nuclei EOS                        
            EOSFLG = 3 
!                    Set flag to indicate switch is being made          
            SWTFLG = 1 
!                                                                       
            WRITE(*,*) ' NUCEOS failed at try =',t,brydns,ye 
            WRITE(*,*) ' where it shouldnt have; Bulk EOS was used' 
            WRITE(*,*) ' IV = ',INPVAR 
            WRITE(*,*) ' ' 
!                                                                       
!                    Branch to bulk EOS                                 
            GOTO 50 
!                                                                       
!                    Otherwise since forcing flag is set then declare   
!                    a failure and return                               
          ELSE 
!                      If the failure message flag is set then announce 
!                      the failure                                      
            IF(NF_FLG.EQ.1)                                             &
     &          WRITE(*,*) 'NUC failed at t,r = ',t,brydns              
            GOTO 999 
          ENDIF 
        ENDIF 
!                                                                       
      ENDIF 
!                              END OF FULL EOS CALULATIONS              
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                              CALCULATE BULK EOS                       
   50 CONTINUE 
      IF(EOSFLG.EQ.3) THEN 
!                                                                       
        CALL ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
        IF((SSFLAG.EQ.0).AND.(FFLAG.EQ.1).AND.(NF_FLG.EQ.1)) THEN 
          WRITE(*,*) 'A1 failed at t,rho = ',t,brydns 
          GOTO 999 
        ENDIF 
!                           If this EOS was used as a result of the     
!                           nuclear EOS failing then set the            
!                           success flag to indicate a warning          
        IF(SWTFLG.EQ.1) THEN 
          SSFLAG = 2 
        ENDIF 
!                                                                       
!                           Save the value of the proton fraction       
        P_PREV = YE*BRYDNS 
!                                                                       
        GOTO 999 
!                                                                       
      ENDIF 
!                END OF BULK EOS CALCULATIONS                           
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                              CALCULATE VIA MAXWELL CONSTRUCTION       
      IF(EOSFLG.EQ.4) THEN 
!                                                                       
        CALL MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                 Save the value of the proton fraction                 
        P_PREV = YE*BRYDNS 
!                                                                       
!                 If Maxwell EOS failed then announce the failure       
        IF(SSFLAG.NE.1) THEN 
          WRITE(*,*) ' MAXWEL failed at try = ' 
          WRITE(*,*) T,BRYDNS,YE 
        ENDIF 
!                                                                       
          GOTO 999 
!                                                                       
      ENDIF 
!                END OF MAXWELL CONSTRUCTION CALCULATIONS               
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         NUCEOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       NUCEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/13/90 Modified from model 1-d                      
!                                                                       
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!    CALL LINE:    CALL NUCEOS(INPVAR,YE,BRYDNS,X_PREV,SSFLAG)          
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:      XPREV = PREVIOUS VALUE OF X (MUST BE SUPPLIED ON     
!                  FIRST CALL)                                          
!                  SSFLAG = SUCCESS FLAG 0 --> FAILURE                  
!                                        1 --> SUCCESS                  
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE NUCEOS(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                       Function type declarations                      
!                                                                       
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFI, FHALFO 
      double precision fhalf 
!                                                                       
      DOUBLE PRECISION ZNG, ZPG 
      INTEGER TCFLAG, ftflag 
!                                                                       
      INTEGER KKI,LLI 
      DOUBLE PRECISION RESULT(5), R_CHECK(5) 
      double precision a_tmp(5,5) 
      DOUBLE PRECISION NI_MIN 
!                                                                       
      integer cflag, schflg 
      double precision dtst1, dtst2 
      double precision break, dnsi, dtmp8 
      double precision dtmp1,dtmp2,dtmp3,dtmp4,dtmp5,dtmp6,dtmp7 
!c      double precision tbsph, tbph, tbnh, tbspl, tbpl, tbnl           
!c      double precision dbspdx, dbpdx, dbndx, dbspdu, dbpdu, dbndu     
!c      double precision tsgl, tsgh, thl, thh, dsgdx, dhfdx, ds2dx,dzdx 
!c      double precision dpt1dx, dpt2dx                                 
!                                                                       
!                                                                       
!                         Set the scheme flag to zero                   
      SCHFLG = 0 
!                                                                       
!                                                                       
    5 CONTINUE 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options are not implemente
!                         in this version)                              
      T = INPVAR(1) 
      NSUBI = INPVAR(2) 
      ETA_PO = INPVAR(3) 
      ETA_NO = INPVAR(4) 
!                                                                       
!                                                                       
!                         Calc the quantum concentration of nucleons    
      NQ = 2.36D-4*T**1.5 
!                                                                       
!                         Calc the Fermi integral coefficent            
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
      KQ = ((T/UQ)**2.5)/(2.0*PI**2) 
!                                                                       
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2)) 
!                                                                       
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ) 
!                                                                       
      IF(ETA_PO.GE.ETAMAX) ETA_PO = ETAMAX-0.1 
      IF(ETA_NO.GE.ETAMAX) ETA_NO = ETAMAX-0.1 
      NI_MIN = DMAX1(4.5D-2,BRYDNS) 
      IF(NSUBI.LT.NI_MIN) NSUBI = NI_MIN+1.0D-3 
!                                                                       
      TCFLAG = 0 
!                                                                       
      cflag = 0 
!                                                                       
      NEWFLG = 1 
!                                                                       
!                    Start Newton-Raphson iteration here                
!                                                                       
!                                                                       
      DO 30 I=1,MAXIT,1 
!                                                                       
        IT_NUM = I 
!                       Set the "Negative" flag                         
        NGFLAG = 0 
!                                                                       
!                                                                       
!                                                                       
        NNOUT = MQ*F_1_2(ETA_NO)/(2.0*PI**2) 
        NPOUT = MQ*F_1_2(ETA_PO)/(2.0*PI**2) 
!                                                                       
        NOUT = NPOUT+NNOUT 
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
!                                                                       
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
!                                                                       
        MUALFA = 2.0*MUN_O+2.0*MUP_O+BALPHA-BPROUT*V_ALFA 
!                                                                       
        IF(ABS(MUALFA/T).LT.30.0) THEN 
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
        ELSEIF((MUALFA/T).LT.-30.0) THEN 
          ALFDNS = 0.0 
        ELSE 
          ALFDNS = 8.0*NQ*DEXP(3.0D1) 
        ENDIF 
!                                                                       
!                                                                       
!                   These statements take out the alfas if the          
!                   alpha particle enable flag is not set               
        IF(ALFLAG.NE.1) THEN 
          ALFDNS = 0.0 
          MUALFA = -300.0 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                                                                       
        BPRALF = ALFDNS*T 
!                                                                       
!---------------------------------------------------                    
!                                                                       
!                                                                       
!             Calculate fraction of space occupied by nuclei            
        U_NUC = (BRYDNS-EXALFA*NOUT-4.0*ALFDNS)/                        &
     &        (NSUBI-EXALFA*NOUT-4.0*ALFDNS)                            
!                                                                       
!                                                                       
!            Is volume occupied by nuclei within acceptable limits?     
!c        IF((U_NUC.LT.0.0).OR.((U_NUC-1.0).GT.-1.0E-20)) THEN          
        IF((U_NUC.LT.0.0).OR.(U_NUC.GT.0.996)) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!            Volume exclusion factor due to nuclei                      
        EXCLU = 1.0-U_NUC 
!                                                                       
!                                                                       
!            If calculated nucleon and alfa densities are larger        
!            than the baryon density then reduce the eta's              
        IF((EXCLU*EXALFA*NOUT+EXCLU*4.0*ALFDNS).GT.BRYDNS) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!            Calculate the internal (inside nuclei) proton fraction     
!                                                                       
        X = (BRYDNS*YE-(1.0-U_NUC)*(EXALFA*NPOUT+2.0*ALFDNS))/          &
     &    (U_NUC*NSUBI)                                                 
        COMPX = 1.0-X 
!                                                                       
!                                                                       
!            Is X within reasonable (but not necessarily correct)       
!            limits? (YE may not be the lower bound on X !!!)           
!ccc        X_MIN = DMAX1(1.0D-2,(YE-0.05))                             
        X_MIN = DMAX1(1.0D-2,(0.8*YE)) 
!c        x_min = 0.95*ye                                               
        IF((X.LT.X_MIN).OR.(X.GT.0.6)) THEN 
          NGFLAG = 1 
          GOTO 29 
        ENDIF 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                     Calculate critical temperature & its X derivative 
        TSC_12 = 87.76*((COMP/375.0)**0.5)*((0.155/NSUBS)**OVR3) 
!                                                                       
!ccdebug      tsc_12 = 1.0d8                                            
!                                                                       
        TSUBC = TSC_12*X*COMPX 
        DTCDX = TSC_12*(1.0-2.0*X) 
        DTCDXX = -2.0*TSC_12 
!                                                                       
!c        tsubc = tsc_12*0.25                                           
!c        dtcdx = 0.0                                                   
!c        dtcdxx = 0.0                                                  
!                                                                       
!                                                                       
!C        TSUBC = 80.0*X*COMPX                                          
!C        DTCDX = 80.0*(1.0-2.0*X)                                      
!C        DTCDXX = -160.0                                               
!                                                                       
!                     If the X is such that T is greater than the       
!                     critical temperature then fix NSUBI so that       
!                     it lies in the bounds of acceptable parameter     
!                     space                                             
        ftflag = 0 
        IF((T.GT.TSUBC).AND.(SCHFLG.EQ.0)) THEN 
!                       If this is an initial guess, then lower         
!                       NSUBI untill we get a good X                    
          IF(NEWFLG.EQ.1) THEN 
!c        write(*,*) ' nuc exceeded Tc'                                 
!c        write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc                   
            ZNG = 2.0*(PI**2)*BRYDNS*(1.0-0.1*YE)/(1.5*MQ) 
            ZPG = 2.0*(PI**2)*BRYDNS*0.1*YE/(1.5*MQ) 
            IF(TCFLAG.NE.1) THEN 
              TCFLAG = 1 
              ETA_PO = FINV12(ZPG)-0.0 
              ETA_NO = FINV12(ZNG)-0.0 
            ELSE 
              ETA_PO = ETA_PO-2.0/T 
              ETA_NO = ETA_NO-2.0/T 
              NSUBI = DMAX1(0.9*NSUBI,5.1D-2) 
            ENDIF 
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '1',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
            GOTO 30 
          ELSE 
!                       Otherwise go back and cut the stepsize in       
!                       half since it was obviously too big             
            NGFLAG = 1 
            GOTO 29 
          ENDIF 
        ELSEIF((T.GT.TSUBC).AND.(SCHFLG.EQ.1)) THEN 
          ftflag = 1 
          tsubc = 80.0*(0.25+0.5*ye)*(0.75-0.25*ye) 
!                                                                       
        ENDIF 
!                                                                       
!                                                                       
        R_0 = (0.75/(PI*NSUBS))**OVR3 
        Q = (384.0*PI*(R_0**2)*SIG_S/SYM_S)-16.0 
!                                                                       
!                        Calculate surface functions of the internal    
!                        (nuclear) proton fraction, X                   
        SIGMA = 1.0/(Q+1.0/(X**3)+1.0/(COMPX**3)) 
        OVRX4 = (1.0/X**4)-(1.0/COMPX**4) 
        DSIGDX = 3.0*(SIGMA**2)*OVRX4 
        SIGSGP = DSIGDX/SIGMA 
        SIGSG2 = 18.0*(SIGMA**2)*OVRX4**2-12.0*SIGMA*((1.0/X**5)+       &
     &  (1.0/COMPX**5))                                                 
!                                                                       
!                        If T is less than critical temp then           
        IF(T.LT.TSUBC) THEN 
!                        Calculate the surface energy temperature factor
!                        and its X and T derivatives                    
          H = 1.0-2.0*(T/TSUBC)**2+(T/TSUBC)**4 
          HPRIM = -4.0*T/(TSUBC**2)+4.0*((T/TSUBC)**3)/TSUBC 
          HPPRIM = -4.0/(TSUBC**2)+12.0*(T**2)/(TSUBC**4) 
          DHDX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDX 
          DHDXX = 4.0*(T**2/TSUBC**3-T**4/TSUBC**5)*DTCDXX+             &
     &    4.0*(-3.0*T**2/TSUBC**4+5.0*T**4/TSUBC**6)*(DTCDX**2)         
          HX = DHDX/H 
          DHDTDX = 8.0*(T/TSUBC**3-2.0*(T**3)/TSUBC**5)*DTCDX 
!                                                                       
!                                                                       
!                        X independent version of TZERO                 
!          TZERO = 0.25*TSC_12                                          
!          DTZDX = 0.0                                                  
!          DTZDXX = 0.0                                                 
!                        X dependent version of TZERO                   
!          TZERO = TSUBC                                                
!          DTZDX = DTCDX                                                
!          DTZDXX = DTCDXX                                              
!                                                                       
!                                                                       
!                                                                       
!                        Coulomb liquid correction factors and their    
!                        derivatives                                    
!          W = 1-(T/TZERO)**2                                           
!          DWDX = 2.0*(T**2)*DTZDX/(TZERO**3)                           
!          DWDT = -2.0*T/(TZERO**2)                                     
!          DWDTDX = 4.0*T*DTZDX/(TZERO**3)                              
!          DWDXDX = 2.0*(T**2)*                                         
!     1    (DTZDXX/(TZERO**3)-3.0*(DTZDX**2)/(TZERO**4))                
!          DWDTDT = -2.0/(TZERO**2)                                     
!                                                                       
          w = 1.0 
          dwdt = 0.0 
          dwdx = 0.0 
          dwdtdx = 0.0 
          dwdxdx = 0.0 
          dwdtdt = 0.0 
!                                                                       
!                                                                       
!                                                                       
!                        Calc lattice factor & derivatives & products   
!                                                                       
          EXCLU = 1.0-U_NUC 
          COMPU = 1.0-U_NUC 
!                                                                       
          DU = DMAX1(1.0D-15, (1.0-1.5*W*U_NUC**OVR3+0.5*U_NUC)) 
          DMU = DMAX1(1.0D-15,(1.0-1.5*W*(1.0-U_NUC+1.0E-20)**OVR3+     &
     &    0.5*(1.0-U_NUC)))                                             
!                                                                       
          DUP = -0.5*W*U_NUC**M2OVR3+0.5 
          DMUP =-0.5*W*(1.0-U_NUC+1.0E-20)**M2OVR3+0.5 
          DUPP = OVR3*W*((U_NUC+1.0D-20)**M5OVR3) 
          DMUPP = OVR3*W*((1.0-U_NUC)+1.0E-20)**M5OVR3 
!                                                                       
!                Derivatives w.r.t. T                                   
!                                                                       
          DUT = -1.5*DWDT*U_NUC**OVR3 
          DMUT = -1.5*DWDT*(1.0-U_NUC+1.0E-20)**OVR3 
          DUPT = -0.5*DWDT*U_NUC**M2OVR3 
          DMUPT = -0.5*DWDT*(1.0-U_NUC+1.0E-20)**M2OVR3 
!                                                                       
!                Derivatives w.r.t. X                                   
!                                                                       
          DUX = -1.5*DWDX*U_NUC**OVR3 
          DMUX = -1.5*DWDX*(1.0-U_NUC+1.0E-20)**OVR3 
          DUPX = -0.5*DWDX*U_NUC**M2OVR3 
          DMUPX = -0.5*DWDX*(1.0-U_NUC+1.0E-20)**M2OVR3 
!                                                                       
!                Second derivatives w.r.t. X                            
!                                                                       
          DUXX = -1.5*DWDXDX*U_NUC**OVR3 
          DMUXX = -1.5*DWDXDX*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                Second derivatives w.r.t. T                            
!                                                                       
          DUTT = -1.5*DWDTDT*U_NUC**OVR3 
          DMUTT = -1.5*DWDTDT*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                Second derivatives w.r.t. X & T                        
!                                                                       
          DUXT = -1.5*DWDTDX*U_NUC**OVR3 
          DMUXT = -1.5*DWDTDX*(1.0-U_NUC+1.0E-20)**OVR3 
!                                                                       
!                                                                       
          TMP1 = (U_NUC**2)+(COMPU**2)+0.6*(U_NUC*COMPU)**2 
          TMP1P = 4.0*U_NUC-2.0+                                        &
     &    2.0*0.6*(U_NUC*COMPU**2-COMPU*U_NUC**2)                       
          TMP1PP = 4.0+2.0*0.6*(COMPU**2-4.0*U_NUC*COMPU+U_NUC**2) 
!                                                                       
          TMP2 = COMPU*(DU**OVR3) 
          TMP2P = -1.0*DU**OVR3+OVR3*COMPU*(DU**M2OVR3)*DUP 
          TMP2PP = -OVR23*(DU**M2OVR3)*DUP-OVR29*COMPU*                 &
     &    (DU**M5OVR3)*DUP**2+OVR3*COMPU*(DU**M2OVR3)*DUPP              
!                                                                       
          TMP2T = OVR3*COMPU*(DU**M2OVR3)*DUT 
          TMP2X = OVR3*COMPU*(DU**M2OVR3)*DUX 
          TMP2XX = OVR3*COMPU*(DU**M2OVR3)*DUXX+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUX**2)                   
          TMP2TT = OVR3*COMPU*(DU**M2OVR3)*DUTT+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*(DUT**2)                   
          TMP2XT = OVR3*COMPU*(DU**M2OVR3)*DUXT+                        &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUX*DUT                    
          TMP2PT = -OVR3*(DU**M2OVR3)*DUT+                              &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUT+                   &
     &        OVR3*COMPU*(DU**M2OVR3)*DUPT                              
          TMP2PX = -OVR3*(DU**M2OVR3)*DUX+                              &
     &        M2OVR3*OVR3*COMPU*(DU**M5OVR3)*DUP*DUX+                   &
     &        OVR3*COMPU*(DU**M2OVR3)*DUPX                              
!                                                                       
!                                                                       
!                                                                       
          TMP3 = U_NUC*(DMU**OVR3) 
          TMP3P = (DMU**OVR3)-OVR3*U_NUC*(DMU**M2OVR3)*DMUP 
          TMP3PP = -OVR23*(DMU**M2OVR3)*DMUP-OVR29*U_NUC*               &
     &    (DMU**M5OVR3)*(DMUP**2)+OVR3*U_NUC*(DMU**M2OVR3)*DMUPP        
!                                                                       
          TMP3T = OVR3*U_NUC*(DMU**M2OVR3)*DMUT 
          TMP3X = OVR3*U_NUC*(DMU**M2OVR3)*DMUX 
          TMP3XX = OVR3*U_NUC*(DMU**M2OVR3)*DMUXX+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUX**2)                 
          TMP3TT = OVR3*U_NUC*(DMU**M2OVR3)*DMUTT+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*(DMUT**2)                 
          TMP3XT = OVR3*U_NUC*(DMU**M2OVR3)*DMUXT+                      &
     &        M2OVR3*OVR3*U_NUC*(DMU**M5OVR3)*DMUX*DMUT                 
          TMP3PT = OVR3*(DMU**M2OVR3)*DMUT-OVR3*M2OVR3*U_NUC*           &
     &      (DMU**M5OVR3)*DMUP*DMUT-OVR3*U_NUC*(DMU**M2OVR3)*DMUPT      
                                                                        
          TMP3PX = OVR3*(DMU**M2OVR3)*DMUX-OVR3*M2OVR3*U_NUC*           &
     &      (DMU**M5OVR3)*DMUP*DMUX-OVR3*U_NUC*(DMU**M2OVR3)*DMUPX      
!                                                                       
!                                                                       
!                 Combination D function                                
!                                                                       
          SCRDU = U_NUC*COMPU*(TMP2+TMP3)/TMP1 
          SCRDUT = U_NUC*COMPU*(TMP2T+TMP3T)/TMP1 
          SCRDUX = U_NUC*COMPU*(TMP2X+TMP3X)/TMP1 
          SCRDXX = U_NUC*COMPU*(TMP2XX+TMP3XX)/TMP1 
          SCRDTT = U_NUC*COMPU*(TMP2TT+TMP3TT)/TMP1 
          SCRDXT = U_NUC*COMPU*(TMP2XT+TMP3XT)/TMP1 
!                                                                       
          SCRD = SCRDU/U_NUC 
          SCRDT = SCRDUT/U_NUC 
          SCRDX = SCRDUX/U_NUC 
!                                                                       
          SCRD2 = SCRDU/COMPU 
          SCRD2T = SCRDUT/COMPU 
          SCRD2X = SCRDUX/COMPU 
!                                                                       
          SCRDUP = SCRD-SCRD2+U_NUC*COMPU*                              &
     &    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)                
!                                                                       
          SCRDPT = SCRDT-SCRD2T+U_NUC*COMPU*                            &
     &    ((TMP2PT+TMP3PT)/TMP1-(TMP2T+TMP3T)*TMP1P/TMP1**2)            
!                                                                       
          SCRDPX = SCRDX-SCRD2X+U_NUC*COMPU*                            &
     &    ((TMP2PX+TMP3PX)/TMP1-(TMP2X+TMP3X)*TMP1P/TMP1**2)            
!                                                                       
          SCRDPP = (SCRDUP-SCRD)/U_NUC-(SCRD2+SCRDUP)/COMPU+            &
     &    (1.0-2.0*U_NUC)*                                              &
     &    ((TMP2P+TMP3P)/TMP1-(TMP2+TMP3)*TMP1P/TMP1**2)+U_NUC*COMPU*   &
     &    ((TMP2PP+TMP3PP)/TMP1-2.0*(TMP2P+TMP3P)*TMP1P/TMP1**2-        &
     &    (TMP2+TMP3)*TMP1PP/TMP1**2+                                   &
     &    2.0*(TMP2+TMP3)*(TMP1P**2)/TMP1**3)                           
!                                                                       
!                                                                       
!                                                                       
!           bubble D function                                           
!bub          scrdu = (1.0-u_nuc)*dmu**ovr3                             
!bub          scrd = scrdu/u_nuc                                        
!bub          scrd2 = dmu**ovr3                                         
!bub          scrdup = -1.0*dmu**ovr3-                                  
!bub     1    ovr3*(1.0-u_nuc)*dmup*dmu**m2ovr3                         
!bub          scrdpp = ovr23*dmup*dmu**m2ovr3-ovr29*(1.0-u_nuc)*        
!bub     1    dmu**m5ovr3*dmup**2+ovr3*(1.0-u_nuc)*dmu**m2ovr3*dmupp    
!                                                                       
!                                                                       
!           nuclei D function                                           
!nuc          scrdu = u_nuc*du**ovr3                                    
!nuc          scrd = du**ovr3                                           
!nuc          scrd2 = scrdu/(1.0-u_nuc)                                 
!nuc          scrdup = du**ovr3+ovr3*u_nuc*dup*du**m2ovr3               
!nuc          scrdpp = ovr23*dup*du**m2ovr3-ovr29*u_nuc*                
!nuc     1    (du**m5ovr3)*(dup**2)+ovr3*u_nuc*(du**m2ovr3)*dupp        
!                                                                       
!                                                                       
!                                                                       
          ZETA_0 = CSSCAL*6.035204*(SIG_S*(16.0+Q))**OVR23 
!                                                                       
!                        Surface energy coefficent                      
          ZETA = ZETA_0*(H*SIGMA*X*NSUBI)**OVR23 
!                                                                       
!                        Derivative of Zeta w.r.t. X                    
          DZDT = OVR23*ZETA*HPRIM/H 
!                                                                       
!                        Derivative of Zeta w.r.t. X                    
          DZDX = OVR23*ZETA*(DHDX/H+SIGSGP+1.0/X) 
!                                                                       
!                        Derivative of Zeta w.r.t. NSUBI                
          DZDNI = OVR23*ZETA/NSUBI 
!                                                                       
!                                                                       
!                                                                       
!                        Nuclear radius                                 
          RSUBN = 9.0*H*SIGMA*SIG_0*U_NUC*(1.0-U_NUC)/                  &
     &    (2.0*ZETA*SCRDU)                                              
!                                                                       
!                        Nuclear volume                                 
          VSUBN = 4.0*PI*(RSUBN**3)/3.0 
!                                                                       
!                        Atomic number                                  
          A = NSUBI*VSUBN 
!                                                                       
!                        Now calc surface, Coulomb free energies        
!                                                                       
          FSUBSC = ZETA*SCRDU/BRYDNS 
          FSUBS = OVR23*ZETA*SCRDU/BRYDNS 
          FSUBC = OVR3*ZETA*SCRDU/BRYDNS 
!                                                                       
!                                                                       
!                                                                       
!                   Translational chemical potential                    
          MUSUBT = TRSCAL*                                              &
     &        T*DLOG((1.0-U_NUC)*(U_NUC*NSUBI)/(NQ*AZERO**2.5))         
!                                                                       
!                   Derivative of trans. chem. potential w.r.t. T       
          DMUTDT = TRSCAL*(MUSUBT/T-1.5) 
!                                                                       
!                   Translational free energy per baryon                
          FTRANS = TRSCAL*H*(MUSUBT-T)/AZERO 
!                                                                       
!                            if T is above the critical temperature     
        ELSE 
          A = 0.0 
          RSUBN = 0.0 
          VSUBN = 0.0 
          FSUBS = 0.0 
          FSUBC = 0.0 
          FTRANS = 0.0 
        ENDIF 
!                            Calc ratio of NSUBI to NSUBS               
        NRATIO = NSUBI/NSUBS 
!                                                                       
!                                                                       
        VNI = 2.0*AA*NSUBI+4.0*BB*X*NSUBI+CC*(1.0+DD)*NSUBI**DD 
!                                                                       
        VPI = 2.0*AA*NSUBI+4.0*BB*(1.0-X)*NSUBI+                        &
     &    CC*(1.0+DD)*NSUBI**DD+DELTAM                                  
!                                                                       
!---------------------------------------------------                    
!                                                                       
        ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ 
!                                                                       
        ZPI = 2.0*(PI**2)*NSUBI*X/MQ 
!                                                                       
        ETA_NI = FINV12(ZNI) 
!                                                                       
        ETA_PI = FINV12(ZPI) 
!                                                                       
        MUN_I = T*ETA_NI+VNI 
!                                                                       
        MUP_I = T*ETA_PI+VPI 
!                                                                       
        F32_NI = F_3_2(ETA_NI) 
!                                                                       
        F32_PI = F_3_2(ETA_PI) 
!                                                                       
        PSUBI = LQ*(F32_NI+F32_PI)+                                     &
     &    (NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+DD*CC*NSUBI**(1.0+DD)        
!                                                                       
!                                                                       
        BN = OVR23*ZETA*SCRD*(SIGSGP+HX+1.5*SCRDUX/SCRDU)*X/NSUBI-      &
     &  TRSCAL*(1.0-U_NUC)*(MUSUBT*(H-X*DHDX)/AZERO+X*DHDX*T/AZERO)     
!                                                                       
        BP = -OVR23*ZETA*SCRD*                                          &
     & ((SIGSGP+HX+1.5*SCRDUX/SCRDU)*COMPX+1.0/X)/NSUBI-                &
     & TRSCAL*(1.0-U_NUC)*                                              &
     & (MUSUBT*(H+DHDX*COMPX)/AZERO-DHDX*T*COMPX/AZERO)                 
!                                                                       
        BSUBP = ZETA*SCRDUP-OVR23*ZETA*SCRD-                            &
     &        TRSCAL*U_NUC*NSUBI*H*MUSUBT/AZERO                         
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!c        GPI = 2.0*FHALFI(ETA_PI)                                      
!c        GPI = 2.0*FHALFI(2.0*(pi**2)*x*nsubi/mq)                      
!c        GNI = 2.0*FHALFI(ETA_NI)                                      
!c        GNI = 2.0*FHALFI(2.0*(pi**2)*(1.0-x)*nsubi/mq)                
!                                                                       
!c        GPO = 2.0*FHALFO(ETA_PO)                                      
!c        GNO = 2.0*FHALFO(ETA_NO)                                      
!                                                                       
!                                                                       
        GPO = 2.0*FHALF(ETA_PO) 
        GNO = 2.0*FHALF(ETA_NO) 
        GPI = 2.0*FHALF(ETA_PI) 
        GNI = 2.0*FHALF(ETA_NI) 
!                                                                       
!                  Derivatives of inside potentials                     
!                                                                       
        DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
        DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
        DVNIDP = DVPIDN 
        DVNIDN = DVPIDP 
!                                                                       
!                  Derivatives of outside potentials                    
!                                                                       
        DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
        DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0))) 
        DVNODP = DVPODN 
        DVNODN = DVPODP 
!                                                                       
!                  Derivatives of inside K.E. densities                 
!                                                                       
        MSSCON = 3.0*MASSN/((HBAR*C)**2) 
        DTPIDP = MSSCON*T*GPI 
        DTPIDN = 0.0 
        DTNIDP = 0.0 
        DTNIDN = MSSCON*T*GNI 
!                                                                       
!                  Derivatives of outside K.E. densities                
!                                                                       
        DTPODP = MSSCON*T*GPO 
        DTPODN = 0.0 
        DTNODP = 0.0 
        DTNODN = MSSCON*T*GNO 
!                                                                       
!                                                                       
!                  Derivatives of inside chem. potentials               
!                                                                       
        DMPIDP = T*GPI/(X*NSUBI)+DVPIDP 
        DMPIDN = DVPIDN 
        DMNIDP = DVNIDP 
        DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN 
!                                                                       
!                  Derivatives of outside chem. potentials              
!                                                                       
        DMPODP = T+DVPODP*NPOUT/GPO 
        DMPODN = DVPODN*NNOUT/GNO 
        DMNODP = DVNODP*NPOUT/GPO 
        DMNODN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of inside pressure                       
!                                                                       
        DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP 
        DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN 
!                                                                       
!                  Derivatives of outside pressure                      
!                                                                       
        DPODP = NPOUT*DMPODP+NNOUT*DMNODP 
        DPODN = NPOUT*DMPODN+NNOUT*DMNODN 
!                                                                       
!                  Derivatives of alpha pressure                        
!                                                                       
        DPADP = ALFDNS*                                                 &
     &  ( (2.0-NPOUT*V_ALFA)*DMPODP+(2.0-NNOUT*V_ALFA)*DMNODP )         
        DPADN = ALFDNS*                                                 &
     &  ( (2.0-NPOUT*V_ALFA)*DMPODN+(2.0-NNOUT*V_ALFA)*DMNODN )         
!                                                                       
!                                                                       
        N1 = NSUBI-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
        N2 = NSUBI*X-EXALFA*NPOUT-2.0*ALFDNS 
!                                                                       
!                  Derivatives of U                                     
!                                                                       
        DUDPO = -EXCLU*(EXALFA*NPOUT/GPO+                               &
     &           (4.0-NOUT*V_ALFA)*DPADP/T)/N1                          
        DUDNO = -EXCLU*(EXALFA*NNOUT/GNO+                               &
     &           (4.0-NOUT*V_ALFA)*DPADN/T)/N1                          
        DUDNI = -U_NUC/N1 
!                                                                       
!                  Derivatives of X                                     
!                                                                       
        DXDPO = -(N2*DUDPO+EXCLU*(EXALFA*NPOUT/GPO+                     &
     &           (2.0-NPOUT*V_ALFA)*DPADP/T))/(U_NUC*NSUBI)             
        DXDNO = -(N2*DUDNO+EXCLU*(2.0-NPOUT*V_ALFA)*DPADN/T)/           &
     &           (U_NUC*NSUBI)                                          
        DXDNI = (N2-X*N1)/(NSUBI*N1) 
!                                                                       
!                  Derivatives of B's w.r.t. NSUBI                      
!                                                                       
        DB1DNI = TRSCAL*( -U_NUC*H*(MUSUBT+T)/AZERO )+                  &
     &      OVR23*ZETA*(SCRDUP-OVR23*SCRD)/NSUBI                        
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
!                                                                       
        DB2DNI = -2.0*ZETA*SCRD*TMP4/(9.0*NSUBI**2)-                    &
     &  TRSCAL*( (COMPU*T/(AZERO*NSUBI))*(H+COMPX*DHDX) )               
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU) 
        DB3DNI = -2.0*ZETA*SCRD*X*TMP4/(9.0*NSUBI**2)-                  &
     &          TRSCAL*( ((COMPU*T)/(AZERO*NSUBI))*(H-X*DHDX) )         
!                                                                       
!                                                                       
!                                                                       
!                  Derivatives of B's w.r.t. X                          
!                                                                       
        DB1DX = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*(SIGSGP+DHDX/H+1.0/X)+   &
     &  OVR23*ZETA*(SCRDPX-OVR23*SCRDX)-                                &
     &  TRSCAL*( U_NUC*NSUBI*DHDX*MUSUBT/AZERO )                        
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
!                                                                       
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+(X**(-2))+(X-1.0)*        &
     &  (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+1.5*SCRDXX/SCRDU-       &
     &  1.5*(SCRDUX/SCRDU)**2)                                          
!                                                                       
!                                                                       
        DB2DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+      &
     &      OVR23*ZETA*SCRD*TMP5/NSUBI-TRSCAL*(                         &
     &      COMPU*(DHDX*MUSUBT+(DHDXX*(1.0-X)-DHDX)*(MUSUBT-T))/AZERO)  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X 
!                                                                       
        TMP5 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU+X*                        &
     &         (SIGSG2-(SIGSGP**2)-(DHDX/H)**2+DHDXX/H+                 &
     &       1.5*SCRDXX/SCRDU-1.5*(SCRDUX/SCRDU)**2)                    
!                                                                       
        DB3DX = OVR23*(ZETA*SCRDUX+SCRDU*DZDX)*TMP4/(U_NUC*NSUBI)+      &
     &      OVR23*ZETA*SCRD*TMP5/NSUBI-                                 &
     &      TRSCAL*( COMPU*(DHDX*T-X*DHDXX*(MUSUBT-T))/AZERO )          
!                                                                       
!                                                                       
!                                                                       
!                  Derivatives of B's w.r.t. U_NUC                      
!                                                                       
        DB1DU = ZETA*(SCRDPP-OVR23*SCRDUP/U_NUC+OVR23*SCRD/U_NUC)-      &
     &  TRSCAL*( NSUBI*H*(MUSUBT+T*(1.0-2.0*U_NUC)/(1.0-U_NUC))/AZERO ) 
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
        TMP5 = (X-1.0)*1.5*(SCRDPX/SCRDU-SCRDUX*SCRDUP/SCRDU**2) 
        DB2DU = (OVR23*ZETA*SCRD/NSUBI)*TMP4*(SCRDUP/SCRDU-1.0/U_NUC)+  &
     &    OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+                          &
     &    TRSCAL*( (H*MUSUBT+DHDX*COMPX*(MUSUBT-T))/AZERO-              &
     &    (T*(1.0-2.0*U_NUC)/U_NUC)*(H+DHDX*COMPX)/AZERO )              
!                                                                       
        TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*X 
        TMP5 = X*1.5*(SCRDPX/SCRDU-SCRDUP*SCRDUX/SCRDU**2) 
        DB3DU = OVR23*ZETA*SCRD*TMP4*(U_NUC*SCRDUP/SCRDU-1.0)/          &
     & (U_NUC*NSUBI)+OVR23*ZETA*SCRDU*TMP5/(U_NUC*NSUBI)+               &
     &  TRSCAL*( (H*MUSUBT-X*DHDX*(MUSUBT-T))/AZERO-                    &
     & T*(1.0-2.0*U_NUC)*(H-X*DHDX)/(AZERO*U_NUC) )                     
!                                                                       
!                                                                       
!                      A1 derivatives                                   
!                                                                       
        DA1ID1 = X*DPIDP+(1.0-X)*DPIDN+NSUBI*(DPIDP-DPIDN)*DXDNI 
        DA1ID2 = NSUBI*(DPIDP-DPIDN)*DXDPO 
        DA1ID3 = NSUBI*(DPIDP-DPIDN)*DXDNO 
!                                                                       
        DA1OD1 = 0.0 
        DA1OD2 = DPODP+DPADP 
        DA1OD3 = DPODN+DPADN 
!                                                                       
        DB1D1 = DB1DNI+DB1DX*DXDNI+DB1DU*DUDNI 
        DB1D2 = DB1DX*DXDPO+DB1DU*DUDPO 
        DB1D3 = DB1DX*DXDNO+DB1DU*DUDNO 
!                                                                       
        DA1D1 = DA1ID1-DB1D1-DA1OD1 
        DA1D2 = DA1ID2-DB1D2-DA1OD2 
        DA1D3 = DA1ID3-DB1D3-DA1OD3 
!                                                                       
!                      A3 derivatives                                   
!                                                                       
        DA3ID1 = X*DMNIDP+(1.0-X)*DMNIDN+NSUBI*(DMNIDP-DMNIDN)*DXDNI 
        DA3ID2 = NSUBI*(DMNIDP-DMNIDN)*DXDPO 
        DA3ID3 = NSUBI*(DMNIDP-DMNIDN)*DXDNO 
!                                                                       
        DA3OD1 = 0.0 
        DA3OD2 = DMNODP 
        DA3OD3 = DMNODN 
!                                                                       
        DB3D1 = DB3DNI+DB3DX*DXDNI+DB3DU*DUDNI 
        DB3D2 = DB3DX*DXDPO+DB3DU*DUDPO 
        DB3D3 = DB3DX*DXDNO+DB3DU*DUDNO 
!                                                                       
        DA3D1 = DA3ID1-DB3D1-DA3OD1 
        DA3D2 = DA3ID2-DB3D2-DA3OD2 
        DA3D3 = DA3ID3-DB3D3-DA3OD3 
!                                                                       
!                      A2 derivatives                                   
!                                                                       
        DA2ID1 = X*DMPIDP+(1.0-X)*DMPIDN+NSUBI*(DMPIDP-DMPIDN)*DXDNI 
        DA2ID2 = NSUBI*(DMPIDP-DMPIDN)*DXDPO 
        DA2ID3 = NSUBI*(DMPIDP-DMPIDN)*DXDNO 
!                                                                       
        DA2OD1 = 0.0 
        DA2OD2 = DMPODP 
        DA2OD3 = DMPODN 
!                                                                       
        DB2D1 = DB2DNI+DB2DX*DXDNI+DB2DU*DUDNI 
        DB2D2 = DB2DX*DXDPO+DB2DU*DUDPO 
        DB2D3 = DB2DX*DXDNO+DB2DU*DUDNO 
!                                                                       
        DA2D1 = DA2ID1-DB2D1-DA2OD1 
        DA2D2 = DA2ID2-DB2D2-DA2OD2 
        DA2D3 = DA2ID3-DB2D3-DA2OD3 
!                                                                       
!                                                                       
!                      Eta derivatives                                  
!                                                                       
        DNDETN = NNOUT/GNO 
        DPDETP = NPOUT/GPO 
!                                                                       
        DA1DN = DA1D1 
        DA1ETP = DA1D2 
        DA1ETN = DA1D3 
!                                                                       
        DA2DN = DA2D1 
        DA2ETP = DA2D2 
        DA2ETN = DA2D3 
!                                                                       
        DA3DN = DA3D1 
        DA3ETP = DA3D2 
        DA3ETN = DA3D3 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        A1 = PSUBI-BSUBP-BPROUT-BPRALF 
        A2 = MUP_I-BP-MUP_O 
        A3 = MUN_I-BN-MUN_O 
!                                                                       
!                                                                       
!                          Unset the "new" flag                         
        NEWFLG = 0 
!                                                                       
        DETERM = DA1DN*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)-                   &
     &           DA1ETP*(DA2DN*DA3ETN-DA2ETN*DA3DN)+                    &
     &           DA1ETN*(DA2DN*DA3ETP-DA2ETP*DA3DN)                     
!                                                                       
        DNSUBI = -1.0*(A1*(DA2ETP*DA3ETN-DA2ETN*DA3ETP)+                &
     &           A2*(DA3ETP*DA1ETN-DA1ETP*DA3ETN)+                      &
     &           A3*(DA1ETP*DA2ETN-DA1ETN*DA2ETP))/DETERM               
!                                                                       
!                                                                       
        DETAP = -1.0*(A1*(DA2ETN*DA3DN-DA2DN*DA3ETN)+                   &
     &          A2*(DA1DN*DA3ETN-DA1ETN*DA3DN)+                         &
     &          A3*(DA1ETN*DA2DN-DA1DN*DA2ETN))/DETERM                  
!                                                                       
!                                                                       
        DETAN = -1.0*(A1*(DA2DN*DA3ETP-DA2ETP*DA3DN)+                   &
     &          A2*(DA1ETP*DA3DN-DA1DN*DA3ETP)+                         &
     &          A3*(DA1DN*DA2ETP-DA1ETP*DA2DN))/DETERM                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                        Check the step size in NSUBI                   
        IF(ABS(DNSUBI/NSUBI).GT.0.04) THEN 
          DNSUBI = 0.04*DNSUBI*NSUBI/ABS(DNSUBI) 
        ENDIF 
   26   CONTINUE 
        NSUBIN = NSUBI+DNSUBI 
        IF((NSUBIN.LT.DMAX1(4.5D-2,BRYDNS)).OR.(NSUBIN.GT.0.25)) THEN 
          DNSUBI = 0.5*DNSUBI 
          GOTO 26 
        ENDIF 
!                                                                       
!                        Check the step size in ETA_PO                  
        IF(ABS(DETAP).GT.4.0) THEN 
          DETAP = 4.0*DETAP/ABS(DETAP) 
        ENDIF 
   27   CONTINUE 
        NETAP = ETA_PO+DETAP 
        IF((NETAP.LT.-5000.0).OR.(NETAP.GT.ETAMAX)) THEN 
          DETAP = 0.5*DETAP 
          GOTO 27 
        ENDIF 
!                                                                       
!                        Check the step size in ETA_NO                  
        IF(ABS(DETAN).GT.4.0) THEN 
          DETAN = 4.0*DETAN/ABS(DETAN) 
        ENDIF 
   28   CONTINUE 
        NETAN = ETA_NO+DETAN 
        IF((NETAN.LT.-5000.0).OR.(NETAN.GT.ETAMAX)) THEN 
          DETAN = 0.5*DETAN 
          GOTO 28 
        ENDIF 
!                                                                       
!                                                                       
!                        Update the variables                           
!cc        if(i.lt.30) write(*,1205) i,nsubi,eta_no,eta_po,x,u_nuc      
 1205   format(i3,1p9e21.14) 
!                                                                       
        NSUBI = NSUBI+DNSUBI 
        ETA_PO = ETA_PO+DETAP 
        ETA_NO = ETA_NO+DETAN 
!                                                                       
!                                                                       
!                                                                       
!                        If the required tolarences have been met       
!                        break out of the loop                          
        IF((ABS(DNSUBI).LT.NSIACC).AND.(ABS(DETAP).LT.PRTACC)           &
     &    .AND.(ABS(DETAN).LT.NUTACC) ) THEN                            
          GOTO 40 
        ELSE 
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '2',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
          GOTO 30 
        ENDIF 
!                                                                       
!                                                                       
   29   CONTINUE 
        IF(NEWFLG.NE.1) THEN 
          cflag = cflag+1 
          DNSUBI = 0.5*DNSUBI 
          NSUBI = NSUBI-DNSUBI 
          DETAP = 0.5*DETAP 
          ETA_PO = ETA_PO-DETAP 
          DETAN = 0.5*DETAN 
          ETA_NO = ETA_NO-DETAN 
          IF(DBFLAG.EQ.1) THEN 
            WRITE(*,2000) '3',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
          ENDIF 
          GOTO 30 
        ELSE 
          NSUBI = NSUBS 
!c          ETA_PO = ETA_PO-0.5/T                                       
!c          ETA_NO = ETA_NO-0.5/T                                       
          ETA_PO = ETA_PO-2.0/T 
          ETA_NO = ETA_NO-2.0/T 
        ENDIF 
!                                                                       
!                                                                       
      IF(DBFLAG.EQ.1) THEN 
        WRITE(*,2000) '4',i,NSUBI,ETA_PO,ETA_NO,DNSUBI 
      ENDIF 
 2000   FORMAT(t2,a,1x,i3,1x,f8.5,3(1X,G13.5)) 
!                                                                       
!                                                                       
   30 END DO 
!                                                                       
!            If scheme 1 has failed try scheme 2                        
      if(schflg.eq.0) then 
        schflg = 1 
        goto 5 
      endif 
!                                                                       
!                                                                       
      SSFLAG = 0 
      GOTO 999 
!                                                                       
!                    Branch label to break out of DO 30 iteration       
   40 CONTINUE 
!                                                                       
!                                                                       
!                    The following logic determines whether this was    
!                    the correct scheme to use, and if not then which   
!                    one should be used                                 
!                                                                       
      if(ftflag.ne.0) then 
        ssflag = 4 
        goto 999 
      endif 
!                                                                       
!                    If calculated critical temperature is less than T, 
!                    then switch to the scheme with no nuclei           
      IF(T.GE.TSUBC) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                    If fraction of nuclei present is zero and no switch
!                    has been made then switch to the no nuclei scheme  
      IF(U_NUC.LE.0.0) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ELSEIF(U_NUC.GT.1.0) THEN 
!                    Set flag to indicate FAILURE                       
        SSFLAG = 0 
        GOTO 999 
      ELSE 
!                    Set flag to indicate success                       
        SSFLAG = 1 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                    If eqns aren't really zeroed then fail             
!                                                                       
!                                                                       
      IF( (ABS(A1).GT.1.0D-5).OR.(ABS(A2).GT.1.0D-5).OR.                &
     &    (ABS(A3).GT.1.0D-5) ) THEN                                    
        SSFLAG = 0 
!c        WRITE(*,*) ' NUCEOS: False convg; A = ',A1,A2,A3              
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF(NSUBI.LT.0.05) THEN 
        WRITE(*,*) 'NUCEOS:: <<WARNING>> NSUBI GETTING CLOSE TO LB' 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
      ZNI = 2.0*(PI**2)*NSUBI*(1.0-X)/MQ 
!                                                                       
      ZPI = 2.0*(PI**2)*NSUBI*X/MQ 
!                                                                       
      ETA_NI = FINV12(ZNI) 
!                                                                       
      ETA_PI = FINV12(ZPI) 
!                                                                       
      MUN_I = T*ETA_NI+VNI 
!                                                                       
      MUP_I = T*ETA_PI+VPI 
!                                                                       
      F32_NI = F_3_2(ETA_NI) 
!                                                                       
      F32_PI = F_3_2(ETA_PI) 
!                                                                       
      EXCLU = 1.0-U_NUC 
      EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                                                                       
!                                                                       
!                    Calculate particle fractions                       
!                                                                       
      XALFA = 4.0*EXCLU*ALFDNS/BRYDNS 
      XNUT = NNOUT*EXCLU*EXALFA/BRYDNS 
      XPROT = NPOUT*EXCLU*EXALFA/BRYDNS 
      XH = 1.0-XPROT-XNUT-XALFA 
      XHCHK = U_NUC*NSUBI/BRYDNS 
!                                                                       
      IF((XH.LT.HEAVCT).OR.(XHCHK.LT.HEAVCT)) THEN 
!                    Set flag to indicate switch is being made          
        SSFLAG = 0 
!c        write(*,*) ' xh,xhchk = ',xh,xhchk                            
        GOTO 999 
      ENDIF 
!                                                                       
      IF((XALFA.LT.0.0).OR.(XH.LT.0.0).OR.                              &
     &   (XNUT.LT.0.0).OR.(XPROT.LT.0.0)) THEN                          
        SSFLAG = 0 
        write(*,*) ' Xs hnpa = ',xh,xnut,xprot,xalfa 
        GOTO 999 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Baryons                                            
!                                                                       
!                                                                       
      MUPROT = MUP_O 
      MUN = MUN_O 
      MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
      IF(ABS((XH-XHCHK)/XHCHK).GT.1.0D-4) THEN 
        SSFLAG = 0 
        GOTO 999 
!CC        WRITE(*,*) ' INCONSISTENCEY IN XH AT',T,BRYDNS,YE,XH,XHCHK   
      ENDIF 
!                                                                       
      NUCDNS = BRYDNS*XH 
!                                                                       
      TAU_PO = KQ*F32_PO 
      TAU_PI = KQ*F32_PI 
!                                                                       
      TAU_NO = KQ*F32_NO 
      TAU_NI = KQ*F32_NI 
!                                                                       
      IF(NOUT.GT.0.0) XOUT = NPOUT/NOUT 
!                                                                       
!                                                                       
!                    Calculate internal energy of outside nucleons,     
!                    alpha particles, and nuclei (per baryon)           
!                                                                       
      BUOUT = (EXCLU*EXALFA/BRYDNS)*( UQ*(TAU_PO+TAU_NO)+EIFLAG*        &
     &    ( (NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+                            &
     &    CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )                             
!                                                                       
      BUNUC = XH*( ( UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*                     &
     & (AA+4.0*BB*X*(1.0-X))+CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM )/       &
     & NSUBI)+FSUBSC*(1.0-T*(SCRDUT/SCRDU+OVR23*HPRIM/H))+              &
     & TRSCAL*                                                          &
     & (1.0-U_NUC)*XH*(FTRANS*(1.0-T*HPRIM/H)-H*(MUSUBT-2.5*T)/AZERO)   
!                                                                       
!                                                                       
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA) 
!                                                                       
      BU = BUOUT+BUALFA+BUNUC 
!                                                                       
!                                                                       
      BSOUT = (EXCLU*EXALFA/BRYDNS)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)- &
     & NNOUT*ETA_NO-NPOUT*ETA_PO )                                      
!                                                                       
!                                                                       
!                    Calculate entropy of alpha particles (per baryon)  
      BSALFA = -0.25*XALFA*(MUALFA/T-2.5) 
!                                                                       
!                                                                       
      BSNUC = XH*( (5.0*UQ/(3.0*T))*(TAU_NI+TAU_PI)-                    &
     & NSUBI*(1.0-X)*ETA_NI-NSUBI*X*ETA_PI )/NSUBI-                     &
     & FSUBSC*(SCRDUT/SCRDU+OVR23*HPRIM/H)-                             &
     & XH*TRSCAL*(1.0-U_NUC)*                                           &
     & ((FTRANS*HPRIM/H)+H*(MUSUBT/T-2.5)/AZERO)                        
!                                                                       
!                    Calculate total baryon entropy (per baryon)        
      BS = BSOUT+BSNUC+BSALFA 
!                                                                       
!                    Calculate free energy of outside nucleons (per bary
      BFOUT = BUOUT-T*BSOUT 
!                                                                       
!                    Calculate free energy of alpha particles (per baryo
      BFALFA = BUALFA-T*BSALFA 
!                                                                       
!                    Calculate free energy of nuclei (per baryon)       
      BFNUC = BUNUC-T*BSNUC 
!                                                                       
!                    Calculate total baryon free energy (per baryon)    
      BFTOT = BFOUT+BFNUC+BFALFA 
!                                                                       
!                    Calculate pressure due to nuclei                   
      BPRNUC = -ZETA*(SCRDU-U_NUC*SCRDUP)+                              &
     & TRSCAL*U_NUC*NSUBI*H*((1.0-U_NUC)*T-U_NUC*MUSUBT)/AZERO          
!                                                                       
!                                                                       
!                    Calculate total baryon pressure                    
      BPRESS = BPROUT+BPRALF+BPRNUC 
!                                                                       
!                                                                       
!                    Leptons & Photons                                  
!                                                                       
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
!                                                                       
!                    Total pressure and eng/ent per baryon              
!                                                                       
      FBARY = BFTOT+FSUBE 
      PBARY = BPRESS+EPRESS 
      MUBARY = YE*MUPROT+(1.0-YE)*MUN 
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN 
!                                                                       
      FTOT = BFTOT+FSUBE+PF 
      UTOT = BU+EU+PU 
      STOT = BS+ES+PS 
      PTOT = BPRESS+EPRESS+PPRESS 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
!                Derivatives of thermodynamic variables                 
!-----------------------------------------------------------------------
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of exterior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and ETA's)    !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                  Derivatives of exterior potentials                   
!                  w.r.t. particle densities                            
      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVNODP = DVPODN 
      DVNODN = DVPODP 
!                                                                       
!                                                                       
!                  Derviatives of exterior chem. pot. w.r.t. ETA's      
!                  (at fixed T)                                         
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of pressure potential w.r.t.             
!                  particle densities                                   
      DV_DPO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
      DV_DNO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
!                                                                       
!                  Derivatives of pressure potential w.r.t. ETA's       
!                  (at fixed T)                                         
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
!                  Derivatives of outside pressure w.r.t. ETA's         
!                  (at fixed T)                                         
      DPODEP = NPOUT*T+DV_DEP 
      DPODEN = NNOUT*T+DV_DEN 
!                                                                       
!                  Derivatives of alpha density w.r.t. ETA's            
!                  (at fixed T)                                         
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T 
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T 
!                                                                       
!                  Derivatives of alpha pressure w.r.t. ETA's           
!                  (at fixed T)                                         
      DPADEP = T*DNADEP 
      DPADEN = T*DNADEN 
!                                                                       
!                  Derivatives of particle densities w.r.t. T           
!                  (at fixed ETA's)                                     
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
!                  Derivatives of exterior chem. pot. w.r.t. T          
!                  (at fixed ETA's)                                     
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
!                  (at fixed ETA's)                                     
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
!                  (at fixed ETA's)                                     
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                  Derivative of alpha chem. pot. w.r.t. T              
!                  (at fixed ETA's)                                     
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
!                  (at fixed ETA's)                                     
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                  Derivative of alpha particle pressure w.r.t. T       
!                  (at fixed ETA's)                                     
      DPADT = ALFDNS+T*DNADT 
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of interior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and density)  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                   Derivatives of kinetic energy densities w.r.t. T    
!                   (holding the number densities (X & NSUBI) fixed)    
      DTPIDT =2.5*TAU_PI/T-2.25*X*NSUBI*GPI/UQ 
      DTNIDT =2.5*TAU_NI/T-2.25*(1.0-X)*NSUBI*GNI/UQ 
!                                                                       
!                   Derivatives of pressures w.r.t. T                   
!                   (holding the number densities (X & NSUBI) fixed)    
      DPIDT = OVR23*UQ*(DTPIDT+DTNIDT) 
!                                                                       
!                   Derivatives of interior chem. pot. w.r.t. T         
!                   (holding the number densities (X & NSUBI) fixed)    
      DMPIDT = ETA_PI-1.5*GPI 
      DMNIDT = ETA_NI-1.5*GNI 
!                                                                       
!                                                                       
!                  Derivatives of inside potentials w.r.t.              
!                  interior proton and neutron densities                
!                  (at fixed T)                                         
      DVPIDP = 2.0*AA+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
      DVPIDN = 2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NSUBI**(DD-1.0)) 
      DVNIDP = DVPIDN 
      DVNIDN = DVPIDP 
!                                                                       
!                   Derivatives of interior chemical potentials         
!                   w.r.t. interior neutron and proton densities        
!                  (at fixed T)                                         
      DMPIDP = T*GPI/(X*NSUBI)+DVPIDP 
      DMPIDN = DVPIDN 
      DMNIDP = DVNIDP 
      DMNIDN = T*GNI/((1.0-X)*NSUBI)+DVNIDN 
!                                                                       
!                   Derivatives of interior pressure                    
!                   w.r.t. interior neutron and proton densities        
!                  (at fixed T)                                         
      DPIDP = X*NSUBI*DMPIDP+(1.0-X)*NSUBI*DMNIDP 
      DPIDN = X*NSUBI*DMPIDN+(1.0-X)*NSUBI*DMNIDN 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of "B" terms    !                  
!                 !      from the chemical and       !                  
!                 !      pressure equilibrium        !                  
!                 !      equations                   !                  
!                 !                                  !                  
!                 !      (w.r.t. Temperature )       !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!             Derivative of term from pressure equilibrium eqn.         
!                                                                       
      DB1DT = OVR23*ZETA*(SCRDUP-OVR23*SCRD)*HPRIM/H+                   &
     &    ZETA*(SCRDPT-OVR23*SCRDT)-                                    &
     &    TRSCAL*U_NUC*NSUBI*(HPRIM*MUSUBT+H*DMUTDT)/AZERO              
!                                                                       
!                                                                       
!             Derivative of term from proton equilibrium eqn.           
!                                                                       
      TMP4 = (SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU)*(X-1.0)-1.0/X 
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+                                  &
     & 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2                      
!                                                                       
      DB2DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*TMP4+                   &
     &    OVR23*ZETA*SCRDT*TMP4/NSUBI+                                  &
     &    OVR23*(ZETA*SCRD/NSUBI)*(X-1.0)*TMP5-                         &
     &    TRSCAL*EXCLU*(DMUTDT*(H+DHDX*(1.0-X))+MUSUBT*                 &
     &    (HPRIM+DHDTDX*(1.0-X))-DHDX*(1.0-X)-T*DHDX*(1.0-X))/AZERO     
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!             Derivative of term from neutron equilibrium eqn.          
!                                                                       
      TMP4 = SIGSGP+DHDX/H+1.5*SCRDUX/SCRDU 
      TMP5 = DHDTDX/H-DHDX*HPRIM/H**2+                                  &
     & 1.5*SCRDXT/SCRDU-1.5*SCRDUX*SCRDUT/SCRDU**2                      
      DB3DT = OVR49*(ZETA*SCRD*HPRIM/(H*NSUBI))*X*TMP4+                 &
     &        OVR23*(ZETA*SCRDT/NSUBI)*X*TMP4+                          &
     &        OVR23*(ZETA*SCRD/NSUBI)*X*TMP5-                           &
     &        TRSCAL*EXCLU*(HPRIM*MUSUBT+H*DMUTDT-X*DHDTDX*(MUSUBT-T)-  &
     &        X*DHDX*(DMUTDT-1.0))/AZERO                                
!                                                                       
!                                                                       
!                                                                       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of constraint   !                  
!                 !      and equilibrium equations   !                  
!                 !      with respect to the five    !                  
!                 !      compositional variables     !                  
!                 !      (U,x,n_i,eta_po,eta_no)     !                  
!                 !      and the three independent   !                  
!                 !      variables                   !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 1 (Baryon conservation)                       
!                                                                       
      DFDOM(1,1) = NOUT*EXALFA+4.0*ALFDNS-NSUBI 
!                                                                       
      DFDOM(1,2) = 0.0 
!                                                                       
      DFDOM(1,3) = -U_NUC 
!                                                                       
      DFDOM(1,4) = -EXCLU*EXALFA*NPOUT/GPO+                             &
     &             V_ALFA*DNADEP*EXCLU*NOUT-4.0*EXCLU*DNADEP            
!                                                                       
      DFDOM(1,5) = -EXCLU*EXALFA*NNOUT/GNO+                             &
     &             V_ALFA*DNADEN*EXCLU*NOUT-4.0*EXCLU*DNADEN            
!                                                                       
!                                                                       
!                                                                       
      DFDL_1(1) = -1.0 
!                                                                       
      DFDL_2(1) = EXCLU*EXALFA*(DNPODT+DNNODT)-EXCLU*V_ALFA*NOUT*DNADT+ &
     &     4.0*EXCLU*DNADT                                              
!                                                                       
      DFDL_3(1) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 2 (Charge conservation)                       
!                                                                       
      DFDOM(2,1) = EXALFA*NPOUT+2.0*ALFDNS-X*NSUBI 
!                                                                       
      DFDOM(2,2) = -U_NUC*NSUBI 
!                                                                       
      DFDOM(2,3) = -X*U_NUC 
!                                                                       
      DFDOM(2,4) = -EXCLU*EXALFA*NPOUT/GPO+                             &
     &     V_ALFA*EXCLU*NPOUT*DNADEP-2.0*EXCLU*DNADEP                   
!                                                                       
      DFDOM(2,5) = V_ALFA*EXCLU*NPOUT*DNADEN-2.0*EXCLU*DNADEN 
!                                                                       
!                                                                       
!                                                                       
      DFDL_1(2) = -1.0*YE 
!                                                                       
      DFDL_2(2) = EXCLU*EXALFA*DNPODT-V_ALFA*EXCLU*NPOUT*DNADT+         &
     &     2.0*EXCLU*DNADT                                              
!                                                                       
      DFDL_3(2) = -1.0*BRYDNS 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 3 (Proton chemical equilibrium)               
!                                                                       
      DFDOM(3,1) = -DB2DU 
!                                                                       
      DFDOM(3,2) = NSUBI*(DMPIDP-DMPIDN)-DB2DX 
!                                                                       
      DFDOM(3,3) = (1.0-X)*DMPIDN+X*DMPIDP-DB2DNI 
!                                                                       
      DFDOM(3,4) = -DMPDEP 
!                                                                       
      DFDOM(3,5) = -DMPDEN 
!                                                                       
      DFDL_1(3) = 0.0 
      DFDL_2(3) = -1.0*(DMPIDT-DMPODT-DB2DT) 
      DFDL_3(3) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 4 (Neutron chemical equilibrium)              
!                                                                       
      DFDOM(4,1) = -DB3DU 
!                                                                       
      DFDOM(4,2) = NSUBI*(DMNIDP-DMNIDN)-DB3DX 
!                                                                       
      DFDOM(4,3) = (1.0-X)*DMNIDN+X*DMNIDP-DB3DNI 
!                                                                       
      DFDOM(4,4) = -DMNDEP 
!                                                                       
      DFDOM(4,5) = -DMNDEN 
!                                                                       
      DFDL_1(4) = 0.0 
      DFDL_2(4) = -1.0*(DMNIDT-DMNODT-DB3DT) 
      DFDL_3(4) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 5 (Pressure equilibrium)                      
!                                                                       
      DFDOM(5,1) = -DB1DU 
!                                                                       
      DFDOM(5,2) = NSUBI*(DPIDP-DPIDN)-DB1DX 
!                                                                       
      DFDOM(5,3) = (1.0-X)*DPIDN+X*DPIDP-DB1DNI 
      ncomp = dfdom(5,3) 
!                                                                       
      DFDOM(5,4) = -DPODEP-DPADEP 
!                                                                       
      DFDOM(5,5) = -DPODEN-DPADEN 
!                                                                       
      DFDL_1(5) = 0.0 
      DFDL_2(5) = -1.0*(DPIDT-DPODT-DPADT-DB1DT) 
      DFDL_3(5) = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!      write(*,*) ' '                                                   
!c      write(*,*) ' '                                                  
!c      write(*,7010) db1dx,db2dx,db3dx                                 
!c      write(*,7010) db1du,db2du,db3du                                 
!c      write(*,7010) db1dni,db2dni,db3dni                              
!c      write(*,7010) db1dt,db2dt,db3dt                                 
 7010 format(3(1x,g13.6)) 
!      write(*,7000) x,u_nuc,nsubi,eta_po,eta_no                        
!      write(*,*) 'eta_i ',eta_pi,eta_ni                                
!      write(*,*) ' as ',a1,a2,a3                                       
!      write(*,*) ' '                                                   
!      write(*,7000) (dfdom(1,i),i=1,5,1)                               
!      write(*,7000) (dfdom(2,i),i=1,5,1)                               
!      write(*,7000) (dfdom(3,i),i=1,5,1)                               
!      write(*,7000) (dfdom(4,i),i=1,5,1)                               
!      write(*,7000) (dfdom(5,i),i=1,5,1)                               
!      write(*,*) ' '                                                   
!c      write(*,*) ' dna: ',dnadpo,dnadno                               
!c      write(*,*) ' dt: ',dmpidt,dmpodt                                
!      write(*,7000) (dfdl_1(i),i=1,5,1)                                
 7000 format(5(1x,g13.6)) 
!                                                                       
!      pause                                                            
!                    Invert the DFDOM matrix                            
!                                                                       
      CALL MATINV(DFDOM,DFDOMI,5) 
!  IMSL subroutine call to invert the matrix                            
!CC      CALL DLINRG(5,DFDOM,5,DFDOMI,5)                                
!                                                                       
!c      call matmul(dfdom,dfdomi,a_tmp,5,5,5)                           
!                                                                       
!c      write(*,*) ' '                                                  
!c      write(*,7000) (a_tmp(1,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(2,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(3,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(4,i),i=1,5,1)                              
!c      write(*,7000) (a_tmp(5,i),i=1,5,1)                              
!                                                                       
!c      write(*,7000) (dfdomi(1,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(2,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(3,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(4,i),i=1,5,1)                             
!c      write(*,7000) (dfdomi(5,i),i=1,5,1)                             
!c      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '      
!c      write(*,7000) (dfdl_2(i),i=1,5,1)                               
!c      write(*,*) ' >>>>>>>>>>>>>>>    dfdl_2 <<<<<<<<<<<<<<<<< '      
!                                                                       
!      DO 800 LLI=1,5,1                                                 
!        R_CHECK(LLI) = 0.0                                             
!        DO 801 KKI=1,5,1                                               
!          R_CHECK(LLI) = R_CHECK(LLI)+DFDOM(LLI,KKI)*RESULT(KKI)       
! 801    CONTINUE                                                       
!        r_check(lli) = r_check(lli)-dfdl_1(lli)                        
! 800  CONTINUE                                                         
!      write(*,*) ' >>>>>>>>>>>>>>>    R check <<<<<<<<<<<<<<<<< '      
!      write(*,7000) (r_check(i),i=1,5,1)                               
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                    Multiply the DFDL_1 vector by the DFDOMI matrix    
!                    to get the density derivatives                     
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_1,RESULT,5) 
!                                                                       
      DU_DN = RESULT(1) 
      DX_DN = RESULT(2) 
      DNI_DN = RESULT(3) 
      DEP_DN = RESULT(4) 
      DEN_DN = RESULT(5) 
!                                                                       
!                                                                       
!                    Multiply the DFDL_2 vector by the DFDOMI matrix    
!                    to get the Temperature derivatives                 
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_2,RESULT,5) 
!                                                                       
      DU_DT = RESULT(1) 
      DX_DT = RESULT(2) 
      DNI_DT = RESULT(3) 
      DEP_DT = RESULT(4) 
      DEN_DT = RESULT(5) 
!                                                                       
!                    Multiply the DFDL_3 vector by the DFDOMI matrix    
!                    to get the Ye derivatives                          
!                                                                       
      CALL MV_MUL(DFDOMI,DFDL_3,RESULT,5) 
!                                                                       
      DU_DY = RESULT(1) 
      DX_DY = RESULT(2) 
      DNI_DY = RESULT(3) 
      DEP_DY = RESULT(4) 
      DEN_DY = RESULT(5) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of finite size  !                  
!                 !      terms in the internal       !                  
!                 !      energy and entropy          !                  
!                 !      densities w.r.t. to U,X,n_i !                  
!                 !      and T.  These are used in   !                  
!                 !      calculating the derivatives !                  
!                 !      w.r.t. the independant vars !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                        Free energy Surface & Coulomb terms            
!                                  (Densities)                          
!                                                                       
      F_SC = ZETA*SCRDU 
!                                                                       
      DFSCDU = ZETA*SCRDUP 
!                                                                       
      DFSCDX = ZETA*SCRDUX+SCRDU*DZDX 
!                                                                       
      DFSCDN = SCRDU*DZDNI 
!                                                                       
      DFSCDT = ZETA*SCRDUT+SCRDU*DZDT 
!                                                                       
!                                                                       
!                        Free energy translational terms                
!                                  (Densities)                          
      FTR = U_NUC*EXCLU*NSUBI*FTRANS 
!                                                                       
      DFTRDT = FTR*(HPRIM/H+1.0/T)-                                     &
     &    1.5*TRSCAL*U_NUC*EXCLU*NSUBI*H/AZERO                          
!                                                                       
      DFTRDX = FTR*DHDX/H 
!                                                                       
      DFTRDU = FTR/U_NUC-FTR/EXCLU+                                     &
     &    TRSCAL*NSUBI*H*(1.0-2.0*U_NUC)/AZERO                          
!                                                                       
      DFTRDN = FTR/NSUBI+TRSCAL*U_NUC*EXCLU*H*T/AZERO 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!                                                                       
!                        Internal energy Surface & Coulomb terms        
!                                  (Densities)                          
!                                                                       
      TMP4 = 1.0-T*SCRDUT/SCRDU-OVR23*T*HPRIM/H 
!                                                                       
      E_SC = F_SC*TMP4 
!                                                                       
      DESCDU = DFSCDU*TMP4+                                             &
     &    F_SC*(T*SCRDUT*SCRDUP/SCRDU**2-T*SCRDPT/SCRDU)                
!                                                                       
      DESCDX = DFSCDX*TMP4+                                             &
     &    F_SC*(T*SCRDUT*SCRDUX/SCRDU**2-T*SCRDXT/SCRDU+                &
     &    OVR23*T*HPRIM*DHDX/H**2-OVR23*T*DHDTDX/H)                     
!                                                                       
      DESCDN = DFSCDN*TMP4 
!                                                                       
      DESCDT = DFSCDT*TMP4+F_SC*                                        &
     &   (T*(SCRDUT**2)/SCRDU**2-SCRDUT/SCRDU-T*SCRDTT/SCRDU+           &
     &    OVR23*T*(HPRIM**2)/H**2-OVR23*HPRIM/H-OVR23*T*HPPRIM/H)       
!                                                                       
!                        Internal energy translational terms            
!                                  (Densities)                          
!                                                                       
      TMP4 = 1.5*H*T/AZERO-T*HPRIM*(MUSUBT-T)/AZERO 
!                                                                       
      E_TR = TRSCAL*EXCLU*BRYDNS*XH*TMP4 
!                                                                       
      DETRDU = TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4-                      &
     &    NSUBI*(T**2)*HPRIM*(1.0-2.0*U_NUC)/AZERO)                     
!                                                                       
      DETRDX = TRSCAL*BRYDNS*XH*EXCLU*                                  &
     &    (1.5*T*DHDX/AZERO-T*(MUSUBT-T)*DHDTDX/AZERO)                  
!                                                                       
      DETRDN = TRSCAL*(U_NUC*EXCLU*TMP4-                                &
     &    BRYDNS*XH*EXCLU*(T**2)*HPRIM/(NSUBI*AZERO))                   
!                                                                       
      DETRDT = TRSCAL*BRYDNS*XH*EXCLU*                                  &
     &    (1.5*(H+T*HPRIM)/AZERO-(HPRIM+T*HPPRIM)*(MUSUBT-T)/AZERO-     &
     &    T*HPRIM*(MUSUBT/T-2.5)/AZERO )                                
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!                                                                       
!                        Entropy Surface & Coulomb terms                
!                                  (Densities)                          
!                                                                       
      S_SC = (E_SC-F_SC)/T 
!                                                                       
      DSSCDU = (DESCDU-DFSCDU)/T 
!                                                                       
      DSSCDX = (DESCDX-DFSCDX)/T 
!                                                                       
      DSSCDN = (DESCDN-DFSCDN)/T 
!                                                                       
      DSSCDT = (DESCDT-DFSCDT)/T-(E_SC-F_SC)/T**2 
!                                                                       
!                        Entropy translational terms                    
!                                  (Densities)                          
!                                                                       
      TMP4 = MUSUBT*(HPRIM+H/T)/AZERO-(T*HPRIM+2.5*H)/AZERO 
!                                                                       
      S_TR = -TRSCAL*BRYDNS*XH*EXCLU*TMP4 
!                                                                       
      DSTRDU = -TRSCAL*(NSUBI*(1.0-2.0*U_NUC)*TMP4+                     &
     &    NSUBI*T*(1.0-2.0*U_NUC)*(HPRIM+H/T)/AZERO)                    
!                                                                       
      DSTRDX = -TRSCAL*BRYDNS*XH*EXCLU*                                 &
     &    (MUSUBT*(DHDTDX+DHDX/T)/AZERO-                                &
     &    (T*DHDTDX+2.5*DHDX)/AZERO)                                    
!                                                                       
      DSTRDN = -TRSCAL*                                                 &
     &    (U_NUC*EXCLU*TMP4+U_NUC*EXCLU*T*(HPRIM+H/T)/AZERO)            
!                                                                       
      DSTRDT = -(BRYDNS*XH*EXCLU*((MUSUBT/T-1.5)*(HPRIM+H/T)/AZERO+     &
     &    MUSUBT*(HPPRIM+HPRIM/T-H/T**2)/AZERO-                         &
     &    (3.5*HPRIM+T*HPPRIM)/AZERO ))*TRSCAL                          
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of interior bulk !                 
!                 !      terms in the internal        !                 
!                 !      energy and entropy           !                 
!                 !      densities w.r.t. to U,X,n_i  !                 
!                 !      and T.  These are used in    !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
      S_NUC =(OVR53*UQ/T)*(TAU_NI+TAU_PI)-                              &
     &    NSUBI*((1.0-X)*ETA_NI+X*ETA_PI)                               
!                                                                       
      E_NUC = UQ*(TAU_PI+TAU_NI)+(NSUBI**2)*(AA+4.0*BB*X*(1.0-X))+      &
     &    CC*NSUBI**(1.0+DD)+X*NSUBI*DELTAM                             
!                                                                       
!                                                                       
!                    Interior particle densties                         
      NPI = X*NSUBI 
      NNI = (1.0-X)*NSUBI 
!                                                                       
      DTPIDT = 2.5*TAU_PI/T-2.25*NPI*GPI/UQ 
      DTNIDT = 2.5*TAU_NI/T-2.25*NNI*GNI/UQ 
!                                                                       
!               Derivative of interior entropy density w.r.t. T         
      DSIDT = UQ*(DTPIDT+DTNIDT)/T 
!                                                                       
!               Derivative of interior internal energy density w.r.t. T 
      DEIDT = T*DSIDT 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Derivatives of eta's w.r.t. X and NSUBI            
      DETPDX = GPI/X 
      DETNDX = -GNI/(1.0-X) 
      DETPDN = GPI/NSUBI 
      DETNDN = GNI/NSUBI 
!                                                                       
!                    Derivatives of Tau's w.r.t. X and NSUBI            
      DTPIDX = 1.5*T*NPI*DETPDX/UQ 
      DTNIDX = 1.5*T*NNI*DETNDX/UQ 
      DTPDNI = 1.5*T*NPI*DETPDN/UQ 
      DTNDNI = 1.5*T*NNI*DETNDN/UQ 
!                                                                       
!                                                                       
!                                                                       
!           Derivative of interior entropy density w.r.t. X             
      DSIDX = OVR53*UQ*(DTPIDX+DTNIDX)/T-NSUBI*(ETA_PI-ETA_NI)-         &
     &    NSUBI*((1.0-X)*DETNDX+X*DETPDX)                               
!                                                                       
!           Derivative of interior internal energy density w.r.t. X     
      DEIDX = UQ*(DTPIDX+DTNIDX)+                                       &
     &    (NSUBI**2)*4.0*BB*(1.0-2.0*X)+NSUBI*DELTAM                    
!                                                                       
!                                                                       
!           Derivative of interior entropy density w.r.t. NSUBI         
      DSIDN = OVR53*UQ*(DTPDNI+DTNDNI)/T-((1.0-X)*ETA_NI+X*ETA_PI)-     &
     &    NSUBI*((1.0-X)*DETNDN+X*DETPDN)                               
!                                                                       
!                                                                       
!           Derivative of interior internal energy density w.r.t. NSUBI 
      DEIDN = UQ*(DTPDNI+DTNDNI)+2.0*NSUBI*(AA+4.0*BB*X*(1.0-X))+       &
     &    CC*(1.0+DD)*(NSUBI**DD)+X*DELTAM                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of exterior bulk !                 
!                 !      nucleon internal energy &    !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO 
!                                                                       
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*                                &
     &((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)  
!                                                                       
!                   Derivative of exterior entropy density w.r.t. T     
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-                   &
     &     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T                            
!                                                                       
      DEODT = T*DSODT 
!                                                                       
!                    Derivatives of exterior particle densities w.r.t.  
!                    Temperature (ETA's fixed)                          
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                                                                       
      DNPDEP = NPOUT/GPO 
      DNNDEN = NNOUT/GNO 
!                                                                       
      DTPDEP = 1.5*T*NPOUT/UQ 
      DTNDEN = 1.5*T*NNOUT/UQ 
!                                                                       
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP 
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN 
!                                                                       
!                                                                       
!                    Exterior particle potentials                       
      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD ) 
      VPOUT = EIFLAG*                                                   &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)        
!                                                                       
!                                                                       
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP 
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN 
!                                                                       
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of alpha         !                 
!                 !      particle internal energy &   !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      S_ALFA = ALFDNS*(2.5-MUALFA/T) 
!                                                                       
!                                                                       
      E_ALFA = ALFDNS*(1.5*T-BALPHA) 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                                                                       
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                                                                       
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2 
!                                                                       
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS 
!                                                                       
!                                                                       
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP 
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN 
!                                                                       
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP 
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN 
!                                                                       
      DNADEP = ALFDNS*DMADEP/T 
      DNADEN = ALFDNS*DMADEN/T 
!                                                                       
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T 
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T 
!                                                                       
      DEADEP = DNADEP*(1.5*T-BALPHA) 
      DEADEN = DNADEN*(1.5*T-BALPHA) 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      S_DENS = U_NUC*S_NUC+EXCLU*EXALFA*S_OUT+EXCLU*S_ALFA+S_SC+S_TR 
!                                                                       
      E_DENS = U_NUC*E_NUC+EXCLU*EXALFA*E_OUT+EXCLU*E_ALFA+E_SC+E_TR 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !      Temperature Derivatives     !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT 
!                                                                       
!                                                                       
      DBSDT = (DU_DT*S_NUC-                                             &
     &    DU_DT*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DT*S_OUT                  &
     &    -DU_DT*S_ALFA+                                                &
     &    U_NUC*(DSIDT+DSIDX*DX_DT+DSIDN*DNI_DT)+                       &
     &    EXCLU*EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+             &
     &    EXCLU*(DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT)+                    &
     &    DSSCDT+DSSCDU*DU_DT+DSSCDX*DX_DT+DSSCDN*DNI_DT+               &
     &    DSTRDT+DSTRDU*DU_DT+DSTRDX*DX_DT+DSTRDN*DNI_DT)/BRYDNS        
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDT = T*DBSDT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+                 &
     &    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)                 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDT = BRYDNS*(DBMUDT-DBFDT) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !       Density Derivatives        !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN 
!                                                                       
!                                                                       
      DBSDN = (DU_DN*S_NUC-                                             &
     &    DU_DN*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DN*S_OUT-DU_DN*S_ALFA+    &
     &    U_NUC*(DSIDX*DX_DN+DSIDN*DNI_DN)+                             &
     &    EXCLU*EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+                   &
     &    EXCLU*(DSADEP*DEP_DN+DSADEN*DEN_DN)+                          &
     &    DSSCDU*DU_DN+DSSCDX*DX_DN+DSSCDN*DNI_DN+                      &
     &    DSTRDU*DU_DN+DSTRDX*DX_DN+DSTRDN*DNI_DN)/BRYDNS-              &
     &    S_DENS/BRYDNS**2                                              
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDN = (DU_DN*E_NUC-                                             &
     &    DU_DN*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DN*E_OUT-DU_DN*E_ALFA+    &
     &    U_NUC*(DEIDX*DX_DN+DEIDN*DNI_DN)+                             &
     &    EXCLU*EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+                   &
     &    EXCLU*(DEADEP*DEP_DN+DEADEN*DEN_DN)+                          &
     &    DESCDU*DU_DN+DESCDX*DX_DN+DESCDN*DNI_DN+                      &
     &    DETRDU*DU_DN+DETRDX*DX_DN+DETRDN*DNI_DN)/BRYDNS-              &
     &    E_DENS/BRYDNS**2                                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDN = DBUDN-T*DBSDN 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+                        &
     &    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)                        
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !         Ye Derivatives           !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY 
!                                                                       
!                                                                       
      DBSDY = (DU_DY*S_NUC-                                             &
     &    DU_DY*EXALFA*S_OUT-EXCLU*V_ALFA*DNA_DY*S_OUT-DU_DY*S_ALFA+    &
     &    U_NUC*(DSIDX*DX_DY+DSIDN*DNI_DY)+                             &
     &    EXCLU*EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+                   &
     &    EXCLU*(DSADEP*DEP_DY+DSADEN*DEN_DY)+                          &
     &    DSSCDU*DU_DY+DSSCDX*DX_DY+DSSCDN*DNI_DY+                      &
     &    DSTRDU*DU_DY+DSTRDX*DX_DY+DSTRDN*DNI_DY)/BRYDNS               
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDY = (DU_DY*E_NUC-                                             &
     &    DU_DY*EXALFA*E_OUT-EXCLU*V_ALFA*DNA_DY*E_OUT-DU_DY*E_ALFA+    &
     &    U_NUC*(DEIDX*DX_DY+DEIDN*DNI_DY)+                             &
     &    EXCLU*EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+                   &
     &    EXCLU*(DEADEP*DEP_DY+DEADEN*DEN_DY)+                          &
     &    DESCDU*DU_DY+DESCDX*DX_DY+DESCDN*DNI_DY+                      &
     &    DETRDU*DU_DY+DETRDX*DX_DY+DETRDN*DNI_DY)/BRYDNS               
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDY = DBUDY-T*DBSDY 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+                 &
     &    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN                    
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDY = BRYDNS*(DBMUDY-DBFDY) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                End of derivatives of thermodynamic variables          
!-----------------------------------------------------------------------
!                                                                       
!                  Total derivatives                                    
!                  (Baryons+Electrons+Photons)                          
!                                                                       
      DUDT = DBUDT+DEUDT+DPUDT 
      DUDN = DBUDN+DEUDN+DPUDN 
      DUDY = DBUDY+DEUDY+DPUDY 
!                                                                       
!                                                                       
      DSDT = DBSDT+DESDT+DPSDT 
      DSDN = DBSDN+DESDN+DPSDN 
      DSDY = DBSDY+DESDY+DPSDY 
!                                                                       
!                                                                       
      DPDT = DBPDT+DEPDT+DPPDT 
      DPDN = DBPDN+DEPDN+DPPDN 
      DPDY = DBPDY+DEPDY+DPPDY 
!                                                                       
!                                                                       
      DMUDT = DBMUDT+YE*DEMUDT 
      DMUDN = DBMUDN+YE*DEMUDN 
      DMUDY = DBMUDY+YE*DEMUDY 
!                                                                       
!                Calculate the adiabatic index                          
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
!                Set the value of XPREV to X for use the next           
!                time through                                           
!                                                                       
      XPREV = X 
!                                                                       
!                Save the value of the proton density to be used        
!                by the "no nuclei" scheme on the next call             
      P_PREV = NPOUT 
!                                                                       
!                                                                       
!                Return the three internal compositional variables      
      INPVAR(2) = NSUBI 
      INPVAR(3) = ETA_PO 
      INPVAR(4) = ETA_NO 
!                                                                       
!                                                                       
!                                                                       
!                Rejoice for this routine is finished!!!!!!!            
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         ALFEOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       ALFEOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         8/30/90 Modified from model 4-A                      
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU                   
!                            FSWESTY@SBAST3.SUNYSB.EDU                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL ALFEOS(INPVAR,YE,BRYDNS)                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:                                                           
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE ALFEOS(INPVAR,YE,BRYDNS,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                       Function type declarations                      
!                                                                       
      DOUBLE PRECISION F_1_2, F_3_2, FINV12, FHALFO, fhalf 
!                                                                       
!                                                                       
!                                                                       
!                         Ratio of baryon density to saturation density 
      Y = BRYDNS/NSUBS 
!                                                                       
!                                                                       
!                         Set T equal to the input variable (the entropy
!                         and internal energy options are not implemente
!                         in this version)                              
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                         Calc the quantum concentration of nucleons    
      NQ = 2.36D-4*T**1.5 
!                                                                       
!                         Calc the Fermi integral coefficent            
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
      KQ = ((T/UQ)**2.5)/(2.0*PI**2) 
!                                                                       
      LQ = UQ*(MQ**OVR53)/(3.0*(PI**2)) 
!                                                                       
      ETAMAX = 0.95*FINV12(2.0*(PI**2)*BRYDNS/MQ) 
!                                                                       
!                                                                       
!                                                                       
!                              Set the proton density to its old value  
      NPOUT = P_PREV 
!                                                                       
      IF(BRYDNS.GT.(0.98*2.0/(YE*V_ALFA))) THEN 
        NPOUT = YE*BRYDNS 
        NNOUT = (1.0-YE)*BRYDNS 
        NOUT = BRYDNS 
!                                                                       
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
!                                                                       
        ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
        ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
        ETA_NO = FINV12(ZNO) 
!                                                                       
        ETA_PO = FINV12(ZPO) 
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        TAU_NO = KQ*F32_NO 
        TAU_PO = KQ*F32_PO 
!                                                                       
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
        MUN = MUN_O 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
        MUPROT = MUP_O 
!                                                                       
!                              Calculate diff. of chem. potentials      
        MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
!                              Calculate the alpha particle             
!                              chemical potential                       
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
        ALFDNS = 0.0 
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
      ELSE 
!                                                                       
!                              Calculate the neutron density            
        NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+         &
     &            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/                   &
     &            (2.0-BRYDNS*YE*V_ALFA)                                
!                                                                       
!                              Calculate density of outside nucleons    
        NOUT = NPOUT+NNOUT 
!                                                                       
        VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD) 
!                                                                       
        VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                       &
     &    CC*(1.0+DD)*NOUT**DD+DELTAM)                                  
!                                                                       
!                                                                       
        ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
        ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
        ETA_NO = FINV12(ZNO) 
!                                                                       
        ETA_PO = FINV12(ZPO) 
!                                                                       
        F32_NO = F_3_2(ETA_NO) 
        F32_PO = F_3_2(ETA_PO) 
!                                                                       
        TAU_NO = KQ*F32_NO 
        TAU_PO = KQ*F32_PO 
!                                                                       
!                                                                       
        BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                            &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )      
!                                                                       
!                                                                       
        MUN_O = T*ETA_NO+VNOUT 
        MUN = MUN_O 
!                                                                       
        MUP_O = T*ETA_PO+VPOUT 
        MUPROT = MUP_O 
!                                                                       
!                              Calculate diff. of chem. potentials      
        MUHAT = MUN-MUPROT 
!                                                                       
!                                                                       
!                              Calculate the alpha particle             
!                              chemical potential                       
        MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
!                              Calculate density of alpha particles     
!                                                                       
        IF(ABS(MUALFA/T).LT.30.0) THEN 
          ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
        ELSEIF((MUALFA/T).LT.-30.0) THEN 
          ALFDNS = 0.0 
        ELSE 
          ALFDNS = 8.0*NQ*DEXP(3.0D1) 
        ENDIF 
!                                                                       
!                                                                       
        EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                              Calculate "non-zeroness" of baryon       
!                              conservation equation and save the       
!                              value to be used in the finite           
!                              difference approximation of DGDPRT       
        GOLD = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
        PRTOLD = NPOUT 
!                                                                       
!                              Take a small step to get derivative      
        NPOUT = NPOUT+0.001*BRYDNS 
!                                                                       
        DO 11 I=1,30,1 
!                                                                       
!                              Calculate the neutron density            
          NNOUT = 2.0*BRYDNS*(1.0-2.0*YE)/(2.0-BRYDNS*YE*V_ALFA)+       &
     &            NPOUT*(2.0-(1.0-YE)*BRYDNS*V_ALFA)/                   &
     &            (2.0-BRYDNS*YE*V_ALFA)                                
!                                                                       
!                              Calculate density of outside nucleons    
          NOUT = NPOUT+NNOUT 
!                                                                       
          VNOUT = EIFLAG*                                               &
     &      (2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD)             
!                                                                       
          VPOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NNOUT+                     &
     &      CC*(1.0+DD)*NOUT**DD+DELTAM)                                
!                                                                       
!                                                                       
          ZNO = 2.0*(PI**2)*NNOUT/MQ 
!                                                                       
          ZPO = 2.0*(PI**2)*NPOUT/MQ 
!                                                                       
          ETA_NO = FINV12(ZNO) 
!                                                                       
          ETA_PO = FINV12(ZPO) 
!                                                                       
          F32_NO = F_3_2(ETA_NO) 
!                                                                       
          F32_PO = F_3_2(ETA_PO) 
!                                                                       
          TAU_NO = KQ*F32_NO 
          TAU_PO = KQ*F32_PO 
!                                                                       
          BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*                           &
     &      (AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)) )   
!                                                                       
          MUN_O = T*ETA_NO+VNOUT 
          MUN = MUN_O 
!                                                                       
          MUP_O = T*ETA_PO+VPOUT 
          MUPROT = MUP_O 
!                                                                       
!                              Calc difference of potentials            
          MUHAT = MUN-MUPROT 
!                                                                       
!                              Calc alpha particle chemical potentials  
          MUALFA = 2.0*MUN+2.0*MUPROT+BALPHA-BPROUT*V_ALFA 
!                                                                       
!                              Calc alpha particle density              
!                                                                       
          IF(ABS(MUALFA/T).LT.30.0) THEN 
            ALFDNS = 8.0*NQ*DEXP(MUALFA/T) 
          ELSEIF((MUALFA/T).LT.-30.0) THEN 
            ALFDNS = 0.0 
          ELSE 
            ALFDNS = 8.0*NQ*DEXP(3.0D1) 
          ENDIF 
!                                                                       
!                                                                       
          EXALFA = 1.0-ALFDNS*V_ALFA 
!                                                                       
!                              Calc "non-zeroness" of baryon cons. eqn. 
          G = BRYDNS-EXALFA*(NNOUT+NPOUT)-4.0*ALFDNS 
!                                                                       
!                              Calculate derivative of baryon conservati
!                              equation w.r.t. proton density by finite 
!                              diference approximation                  
          DGDPRT = (G-GOLD)/(NPOUT-PRTOLD) 
!                                                                       
!                              Calculate new Newton-Raphson step        
          DPRT = G/DGDPRT 
!                                                                       
!                              Save old value of proton density & G     
          PRTOLD = NPOUT 
          GOLD = G 
!                                                                       
!                                                                       
   13     CONTINUE 
!                                                                       
!                              Potential "new" value of proton density  
          PRTNEW = NPOUT-DPRT 
!                                                                       
!                              If new proton density is less than the   
!                              baryon density and greater than zero     
!                              then update the proton density           
          IF(PRTNEW*(BRYDNS-PRTNEW).GT.0.0) THEN 
            NPOUT = NPOUT-DPRT 
!                              Else cut the step size in half and try ag
          ELSE 
            DPRT = DPRT*0.5 
            GOTO 13 
          ENDIF 
!                                                                       
!                              If step size is small enough break out of
!                              the DO 11 loop, otherwise continue       
          IF(ABS(DPRT/NPOUT).LT.10E-11) GOTO 12 
   11   CONTINUE 
!                                                                       
!      write(*,*) 'A failed to converge; switching to F' ! take out late
        SSFLAG = 0 
        GOTO 999 
!                                                                       
!                                                                       
   12   CONTINUE 
!                                                                       
      ENDIF 
!                              Set the success flag                     
      SSFLAG = 1 
!                                                                       
!                                                                       
!                              Calc outside nucleon density             
      NOUT = NNOUT+NPOUT 
!                                                                       
!                              Calc outside nucleon fraction            
      XOUT = NPOUT/NOUT 
!                                                                       
!                              Calculate particle fractions             
      XALFA = 4.0*ALFDNS/BRYDNS 
      XPROT = EXALFA*NPOUT/BRYDNS 
      XNUT = EXALFA*NNOUT/BRYDNS 
      XH = 0.0 
!                                                                       
!                              Baryons                                  
!                                                                       
      F32_NO = F_3_2(ETA_NO) 
!                                                                       
      F32_PO = F_3_2(ETA_PO) 
!                                                                       
      TAU_PO = KQ*F32_PO 
!                                                                       
      TAU_NO = KQ*F32_NO 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                    Calculate internal energy of outside nucleons      
      BUOUT = (XNUT+XPROT)*( UQ*(TAU_PO+TAU_NO)+                        &
     &    EIFLAG*((NOUT**2)*AA+                                         &
     &   4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM) )/NOUT      
!                                                                       
!                                                                       
!                                Calc alfa particle internal energy     
      BUALFA = 0.25*XALFA*(1.5*T-BALPHA) 
!                                                                       
!                                Set nuclei internal energy to zero     
      BUNUC = 0.0 
!                                Calculate total baryon internal energy 
!                                (per baryon)                           
      BU = BUOUT+BUALFA+BUNUC 
!                                                                       
!                                                                       
!                                Calc entropy of outside nucleons       
      BSOUT = (XNUT+XPROT)*( (5.0*UQ/(3.0*T))*(TAU_NO+TAU_PO)-          &
     &   NNOUT*ETA_NO-NPOUT*ETA_PO )/NOUT                               
!                                                                       
!                                Calc alpha particle entropy            
      BSALFA = 0.25*XALFA*(2.5-MUALFA/T) 
!                                                                       
!                                Set nuclei entropy to zero             
      BSNUC = 0.0 
!                                                                       
!                                Calc total baryon entropy (per baryon) 
      BS = BSOUT+BSALFA+BSNUC 
!                                                                       
!                                                                       
!                                                                       
!                                Calc outside free energy               
      BFOUT = BUOUT-T*BSOUT 
!                                Calc alpha particle free energy        
      BFALFA = BUALFA-T*BSALFA 
!                                Set nuclei free energy to zero         
      BFNUC = BUNUC-T*BSNUC 
!                                Calc total baryon free energy (per nucl
      BFTOT = BFOUT+BFALFA+BFNUC 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                Calc outside pressure                  
      BPROUT = LQ*(F32_PO+F32_NO)+EIFLAG*(                              &
     &    AA*(NOUT**2)+4.0*BB*NPOUT*NNOUT+DD*CC*(NOUT**(1.0+DD)))       
!                                                                       
!                                Calc alfa particle pressure            
      BPRALF = ALFDNS*T 
!                                                                       
!                                Set nuclei pressure to zero            
      BPRNUC = 0.0 
!                                                                       
!                                Calc total baryon pressure             
      BPRESS = BPROUT+BPRALF+BPRNUC 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                           Leptons & Photons                           
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
!                                                                       
!                           Total pressure and eng/ent per baryon       
!                                                                       
      FBARY = BFTOT+FSUBE 
      PBARY = BPRESS+EPRESS 
      MUBARY = YE*MUPROT+(1.0-YE)*MUN 
      MU_MAT = YE*(MUPROT+MUSUBE)+(1.0-YE)*MUN 
!                                                                       
      FTOT = BFTOT+FSUBE+PF 
      UTOT = BU+EU+PU 
      STOT = BS+ES+PS 
      PTOT = BPRESS+EPRESS+PPRESS 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
!                Derivatives of thermodynamic variables                 
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!c      GPO = 2.0*FHALFO(ETA_PO)                                        
!c      GNO = 2.0*FHALFO(ETA_NO)                                        
!                                                                       
!                                                                       
      GPO = 2.0*FHALF(ETA_PO) 
      GNO = 2.0*FHALF(ETA_NO) 
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of exterior     !                  
!                 !      quantities                  !                  
!                 !      (w.r.t. Temp. and ETA's)    !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                  Derivatives of exterior potentials                   
!                  w.r.t. particle densities                            
      DVPODP = EIFLAG*(2.0*AA+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVPODN = EIFLAG*(2.0*AA+4.0*BB+DD*(1.0+DD)*CC*(NOUT**(DD-1.0)) ) 
      DVNODP = DVPODN 
      DVNODN = DVPODP 
!                                                                       
!                                                                       
!                  Derviatives of exterior chem. pot. w.r.t. ETA's      
!                  (at fixed T)                                         
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                  Derivatives of pressure potential w.r.t.             
!                  particle densities                                   
      DV_DPO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
      DV_DNO = EIFLAG*                                                  &
     &    (2.0*AA*NOUT+4.0*BB*NPOUT+CC*DD*(1.0+DD)*(NOUT**DD) )         
!                                                                       
!                  Derivatives of pressure potential w.r.t. ETA's       
!                  (at fixed T)                                         
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
!                  Derivatives of outside pressure w.r.t. ETA's         
!                  (at fixed T)                                         
      DPODEP = NPOUT*T+DV_DEP 
      DPODEN = NNOUT*T+DV_DEN 
!                                                                       
!                  Derivatives of alpha density w.r.t. ETA's            
!                  (at fixed T)                                         
      DNADEP = ALFDNS*(2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP)/T 
      DNADEN = ALFDNS*(2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN)/T 
!                                                                       
!                                                                       
!                  Derivatives of particle densities w.r.t. T           
!                  (at fixed ETA's)                                     
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
!                  Derivatives of exterior chem. pot. w.r.t. T          
!                  (at fixed ETA's)                                     
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
!                  (at fixed ETA's)                                     
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
!                  (at fixed ETA's)                                     
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                  Derivative of alpha chem. pot. w.r.t. T              
!                  (at fixed ETA's)                                     
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
!                  (at fixed ETA's)                                     
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                  Derivative of alpha particle pressure w.r.t. T       
!                  (at fixed ETA's)                                     
      DPADT = ALFDNS+T*DNADT 
!                                                                       
!                                                                       
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !      Derivatives of constraint   !                  
!                 !      and equilibrium equations   !                  
!                 !      with respect to the five    !                  
!                 !      compositional variables     !                  
!                 !      (U,x,n_i,eta_po,eta_no)     !                  
!                 !      and the three independent   !                  
!                 !      variables                   !                  
!                 !      (Baryon density, T, and Ye) !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 1 (Baryon conservation)                       
!                                                                       
!                                                                       
!                                                                       
      DG1DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NOUT-4.0)*DNADEP 
!                                                                       
      DG1DO2 = -EXALFA*NNOUT/GNO+(V_ALFA*NOUT-4.0)*DNADEN 
!                                                                       
!                                                                       
      DG1DL1 = 1.0 
!                                                                       
      DG1DL2 = -EXALFA*(DNNODT+DNPODT)+(V_ALFA*NOUT-4.0)*DNADT 
!                                                                       
      DG1DL3 = 0.0 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                Equation 2 (Charge conservation)                       
!                                                                       
!                                                                       
      DG2DO1 = -EXALFA*NPOUT/GPO+(V_ALFA*NPOUT-2.0)*DNADEP 
!                                                                       
      DG2DO2 = (V_ALFA*NPOUT-2.0)*DNADEN 
!                                                                       
!                                                                       
      DG2DL1 = YE 
!                                                                       
      DG2DL2 = -EXALFA*DNPODT+(V_ALFA*NPOUT-2.0)*DNADT 
!                                                                       
      DG2DL3 = BRYDNS 
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      DET_GT = DG1DO1*DG2DO2-DG1DO2*DG2DO1 
!                                                                       
!                                                                       
      DEP_DN = (DG1DO2*DG2DL1-DG2DO2*DG1DL1)/DET_GT 
      DEN_DN = (DG2DO1*DG1DL1-DG1DO1*DG2DL1)/DET_GT 
!                                                                       
!                                                                       
      DEP_DT = (DG1DO2*DG2DL2-DG2DO2*DG1DL2)/DET_GT 
      DEN_DT = (DG2DO1*DG1DL2-DG1DO1*DG2DL2)/DET_GT 
!                                                                       
!                                                                       
      DEP_DY = (DG1DO2*DG2DL3-DG2DO2*DG1DL3)/DET_GT 
      DEN_DY = (DG2DO1*DG1DL3-DG1DO1*DG2DL3)/DET_GT 
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of exterior bulk !                 
!                 !      nucleon internal energy &    !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
      S_OUT =(OVR53*UQ/T)*(TAU_NO+TAU_PO)-NNOUT*ETA_NO-NPOUT*ETA_PO 
!                                                                       
      E_OUT = UQ*(TAU_PO+TAU_NO)+EIFLAG*                                &
     &((NOUT**2)*AA+4.0*BB*NPOUT*NNOUT+CC*NOUT**(1.0+DD)+NPOUT*DELTAM)  
!                                                                       
!                   Derivative of exterior entropy density w.r.t. T     
      DSODT =  OVR53*UQ*(1.5*(TAU_PO+TAU_NO)/(T**2))-                   &
     &     1.5*(NPOUT*ETA_PO+NNOUT*ETA_NO)/T                            
!                                                                       
      DEODT = T*DSODT 
!                                                                       
!                    Derivatives of exterior particle densities w.r.t.  
!                    Temperature (ETA's fixed)                          
      DNPODT = 1.5*NPOUT/T 
      DNNODT = 1.5*NNOUT/T 
!                                                                       
      DMPODT = ETA_PO+DVPODP*DNPODT+DVPODN*DNNODT 
      DMNODT = ETA_NO+DVNODP*DNPODT+DVNODN*DNNODT 
!                                                                       
!                                                                       
      DNPDEP = NPOUT/GPO 
      DNNDEN = NNOUT/GNO 
!                                                                       
      DTPDEP = 1.5*T*NPOUT/UQ 
      DTNDEN = 1.5*T*NNOUT/UQ 
!                                                                       
      DSODEP = (OVR53*UQ/T)*DTPDEP-NPOUT-ETA_PO*DNPDEP 
      DSODEN = (OVR53*UQ/T)*DTNDEN-NNOUT-ETA_NO*DNNDEN 
!                                                                       
!                                                                       
!                    Exterior particle potentials                       
      VNOUT = EIFLAG*(2.0*AA*NOUT+4.0*BB*NPOUT+CC*(1.0+DD)*NOUT**DD ) 
      VPOUT = EIFLAG*                                                   &
     &    (2.0*AA*NOUT+4.0*BB*NNOUT+CC*(1.0+DD)*NOUT**DD+DELTAM)        
!                                                                       
!                                                                       
      DEODEP = UQ*DTPDEP+VPOUT*DNPDEP 
      DEODEN = UQ*DTNDEN+VNOUT*DNNDEN 
!                                                                       
      DMPDEP = T+DVPODP*NPOUT/GPO 
      DMPDEN = DVPODN*NNOUT/GNO 
      DMNDEP = DVNODP*NPOUT/GPO 
      DMNDEN = T+DVNODN*NNOUT/GNO 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                 -------------------------------------                 
!                 !      Derivatives of alpha         !                 
!                 !      particle internal energy &   !                 
!                 !      entropy densities and the    !                 
!                 !      chem. pot.  w.r.t. to eta_p, !                 
!                 !      ate_n & T. These are used in !                 
!                 !      calculating the derivatives  !                 
!                 !      w.r.t. the independant vars  !                 
!                 !      (Baryon density, T, and Ye)  !                 
!                 !                                   !                 
!                 -------------------------------------                 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      S_ALFA = ALFDNS*(2.5-MUALFA/T) 
!                                                                       
!                                                                       
      E_ALFA = ALFDNS*(1.5*T-BALPHA) 
!                                                                       
!                  Derivative of pressure potential w.r.t. T            
      DV_DT = DV_DPO*DNPODT+DV_DNO*DNNODT 
!                                                                       
!                  Derivative of outside pressure w.r.t. T              
      DPODT = OVR23*UQ*2.5*(TAU_PO+TAU_NO)/T+DV_DT 
!                                                                       
!                                                                       
      DMUADT = 2.0*DMPODT+2.0*DMNODT-V_ALFA*DPODT 
!                                                                       
!                  Derivative of alpha particle density w.r.t. T        
      DNADT = 1.5*ALFDNS/T-ALFDNS*MUALFA/(T**2)+ALFDNS*DMUADT/T 
!                                                                       
!                                                                       
      DSADT = DNADT*(2.5-MUALFA/T)-ALFDNS*DMUADT/T+ALFDNS*MUALFA/T**2 
!                                                                       
      DEADT = DNADT*(1.5*T-BALPHA)+1.5*ALFDNS 
!                                                                       
!                                                                       
      DV_DEP = DV_DPO*NPOUT/GPO 
      DV_DEN = DV_DNO*NNOUT/GNO 
!                                                                       
      DPODEP = OVR23*UQ*DTPDEP+DV_DEP 
      DPODEN = OVR23*UQ*DTNDEN+DV_DEN 
!                                                                       
      DMADEP = 2.0*DMPDEP+2.0*DMNDEP-V_ALFA*DPODEP 
      DMADEN = 2.0*DMPDEN+2.0*DMNDEN-V_ALFA*DPODEN 
!                                                                       
      DNADEP = ALFDNS*DMADEP/T 
      DNADEN = ALFDNS*DMADEN/T 
!                                                                       
      DSADEP = DNADEP*(2.5-MUALFA/T)-ALFDNS*DMADEP/T 
      DSADEN = DNADEN*(2.5-MUALFA/T)-ALFDNS*DMADEN/T 
!                                                                       
      DEADEP = DNADEP*(1.5*T-BALPHA) 
      DEADEN = DNADEN*(1.5*T-BALPHA) 
!                                                                       
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\                                         
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
      S_DENS = EXALFA*S_OUT+S_ALFA 
!                                                                       
      E_DENS = EXALFA*E_OUT+E_ALFA 
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !      Temperature Derivatives     !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
      DNA_DT = DNADT+DNADEP*DEP_DT+DNADEN*DEN_DT 
!                                                                       
!                                                                       
      DBSDT = (-V_ALFA*DNA_DT*S_OUT+                                    &
     &    EXALFA*(DSODT+DSODEP*DEP_DT+DSODEN*DEN_DT)+                   &
     &    (DSADT+DSADEP*DEP_DT+DSADEN*DEN_DT) )/BRYDNS                  
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBUDT = T*DBSDT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDT = DBUDT-S_DENS/BRYDNS-T*DBSDT 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDT = YE*(DMPODT+DMPDEP*DEP_DT+DMPDEN*DEN_DT)+                 &
     &    (1.0-YE)*(DMNODT+DMNDEP*DEP_DT+DMNDEN*DEN_DT)                 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDT = BRYDNS*(DBMUDT-DBFDT) 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !       Density Derivatives        !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
      DNA_DN = DNADEP*DEP_DN+DNADEN*DEN_DN 
!                                                                       
!                                                                       
      DBSDN = (-V_ALFA*DNA_DN*S_OUT+                                    &
     &    EXALFA*(DSODEP*DEP_DN+DSODEN*DEN_DN)+                         &
     &   (DSADEP*DEP_DN+DSADEN*DEN_DN) )/BRYDNS-S_DENS/BRYDNS**2        
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBUDN = (-V_ALFA*DNA_DN*E_OUT+                                    &
     &    EXALFA*(DEODEP*DEP_DN+DEODEN*DEN_DN)+                         &
     &   (DEADEP*DEP_DN+DEADEN*DEN_DN) )/BRYDNS-E_DENS/BRYDNS**2        
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDN = DBUDN-T*DBSDN 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDN = YE*(DMPDEP*DEP_DN+DMPDEN*DEN_DN)+                        &
     &    (1.0-YE)*(DMNDEP*DEP_DN+DMNDEN*DEN_DN)                        
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDN = BRYDNS*(DBMUDN-DBFDN)+MUBARY-BFTOT 
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!                 ------------------------------------                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !         Ye Derivatives           !                  
!                 !                                  !                  
!                 !                                  !                  
!                 !                                  !                  
!                 ------------------------------------                  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DNA_DY = DNADEP*DEP_DY+DNADEN*DEN_DY 
!                                                                       
!                                                                       
      DBSDY = (-V_ALFA*DNA_DY*S_OUT+                                    &
     &    EXALFA*(DSODEP*DEP_DY+DSODEN*DEN_DY)+                         &
     &   (DSADEP*DEP_DY+DSADEN*DEN_DY) )/BRYDNS                         
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBUDY = (-V_ALFA*DNA_DY*E_OUT+                                    &
     &    EXALFA*(DEODEP*DEP_DY+DEODEN*DEN_DY)+                         &
     &   (DEADEP*DEP_DY+DEADEN*DEN_DY) )/BRYDNS                         
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
!                                                                       
      DBFDY = DBUDY-T*DBSDY 
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBMUDY = YE*(DMPDEP*DEP_DY+DMPDEN*DEN_DY)+MUPROT+                 &
     &    (1.0-YE)*(DMNDEP*DEP_DY+DMNDEN*DEN_DY)-MUN                    
!                                                                       
!~~~~~~~~~~~~~~~~~~                                                     
!                                                                       
      DBPDY = BRYDNS*(DBMUDY-DBFDY) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                End of derivatives of thermodynamic variables          
!-----------------------------------------------------------------------
!                                                                       
!                  Total derivatives                                    
!                  (Baryons+Electrons+Photons)                          
!                                                                       
      DUDT = DBUDT+DEUDT+DPUDT 
      DUDN = DBUDN+DEUDN+DPUDN 
      DUDY = DBUDY+DEUDY+DPUDY 
!                                                                       
!                                                                       
      DSDT = DBSDT+DESDT+DPSDT 
      DSDN = DBSDN+DESDN+DPSDN 
      DSDY = DBSDY+DESDY+DPSDY 
!                                                                       
!                                                                       
      DPDT = DBPDT+DEPDT+DPPDT 
      DPDN = DBPDN+DEPDN+DPPDN 
      DPDY = DBPDY+DEPDY+DPPDY 
!                                                                       
!                                                                       
      DMUDT = DBMUDT+YE*DEMUDT 
      DMUDN = DBMUDN+YE*DEMUDN 
      DMUDY = DBMUDY+YE*DEMUDY 
!                                                                       
!                                                                       
!                  Adiabatic index                                      
      GAM_S = BRYDNS*DPDN/PTOT+T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
      INPVAR(2) = NSUBS 
      INPVAR(3) = ETA_PO 
      INPVAR(4) = ETA_NO 
!                                                                       
!                                                                       
!                           Approximate the nuclear density             
      NSUBI = NSUBS 
!                                                                       
!                           Use 0.45 as the nuclear proton fraction     
      X = 0.45 
!                                                                       
!                           Save the proton number density for use      
!                           as the initial guess on next call           
      P_PREV = NPOUT 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         MAXWEL.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MAXWEL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/13/90                                              
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU                   
!                                                                       
!                                                                       
!    CALL LINE:    CALL MAXWEL(INPVAR,YE,BRYDNS)                        
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:                                                           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC                                        
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE MAXWEL(INPVAR,YE,BRYDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION OUTVAR(4) 
      DOUBLE PRECISION DPN_DT, DPN_DN, DPN_DY 
      DOUBLE PRECISION DSN_DT, DSN_DN, DSN_DY 
      DOUBLE PRECISION DSB_DT, DSB_DN, DSB_DY 
      DOUBLE PRECISION DMU_DT, DMU_DN, DMU_DY 
      DOUBLE PRECISION DPHADT, DPHADY, DELDNS 
      DOUBLE PRECISION N_XH, N_XA, N_XN, N_XP, B_XA, B_XN, B_XP 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'el_eos.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
!                                                                       
!                   Set the temperature                                 
      T = INPVAR(1) 
!                                                                       
!                                                                       
!                   Calculate and save chem. pot. and thermodynamic     
!                   quantaties from low end of two phase region         
      CALL NUCEOS(INPVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                                                                       
!                    If the nuclear EOS failed and the reset flag is set
!                    then reset the initial guesses and try again       
      IF((SSFLAG.NE.1).AND.(RSFLAG.EQ.1)) THEN 
        CALL RESET(INPVAR,YE,LOWDNS,OUTVAR) 
        OUTVAR(1) = INPVAR(1) 
        CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
!                                                                       
!                                                                       
!                    Make a last ditch effort at convergence            
        IF(SSFLAG.NE.1) THEN 
          OUTVAR(2) = 0.155 
          OUTVAR(3) = -15.0 
          OUTVAR(4) = -20.0 
          CALL NUCEOS(OUTVAR,YE,LOWDNS,XPREV,P_PREV,SSFLAG) 
        ELSE 
          INPVAR(2) = OUTVAR(2) 
          INPVAR(3) = OUTVAR(3) 
          INPVAR(4) = OUTVAR(4) 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      PRLOW = PTOT-PPRESS 
      S_LOW = STOT-PS 
      F_LOW = FTOT-PF 
      MUTLOW = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE) 
      MUELOW = MUSUBE 
      MUHLOW = MUHAT 
!                                                                       
      DPN_DT = DPDT 
      DPN_DN = DPDN 
      DPN_DY = DPDY 
!                                                                       
      DMU_DT = DMUDT 
      DMU_DN = DMUDN 
      DMU_DY = DMUDY 
!                                                                       
      DSN_DT = DSDT-DPSDT 
      DSN_DN = DSDN 
      DSN_DY = DSDY 
!                                                                       
      N_XH = XH 
      N_XA = XALFA 
      N_XP = XPROT 
      N_XN = XNUT 
!                                                                       
!                                                                       
      IF(SSFLAG.NE.1) THEN 
        WRITE(*,*) 'MAXWEL:  Nuclear EOS failed at try:' 
        WRITE(*,*) T,LOWDNS,YE 
        WRITE(*,*) INPVAR 
        GOTO 999 
      ENDIF 
!                   Calculate and save chem. pot. and thermodynamic     
!                   quantaties from high end of two phase region        
      CALL ALFEOS(INPVAR,YE,HIDNS,P_PREV,SSFLAG) 
!                                                                       
      PRHI = PTOT-PPRESS 
      S_HI = STOT-PS 
      F_HI = FTOT-PF 
      MUTHI = (1.0-YE)*MUN+YE*(MUPROT+MUSUBE) 
      MUEHI = MUSUBE 
      MUHHI = MUHAT 
!                                                                       
!                                                                       
      DSB_DT = DSDT-DPSDT 
      DSB_DN = DSDN 
      DSB_DY = DSDY 
!                                                                       
!                                                                       
      B_XA = XALFA 
      B_XP = XPROT 
      B_XN = XNUT 
!                                                                       
!                                                                       
      IF(SSFLAG.NE.1) THEN 
        WRITE(*,*) 'MAXWEL:  Alfa EOS failed at try:' 
        WRITE(*,*) T,HIDNS,YE 
        WRITE(*,*) INPVAR 
        GOTO 999 
      ENDIF 
!                                                                       
!                   Calculate "average" chem. pot. and pressure         
!                   in order to avoid numerical problems                
      MUTILD = (MUTLOW+MUTHI)/2.0 
      PRTILD = (PRLOW+PRHI)/2.0 
!                                                                       
!                   Calculate phase fraction                            
      PHASEF = (BRYDNS-LOWDNS)/(HIDNS-LOWDNS) 
!                                                                       
!                                                                       
!                   Electron number density                             
      NSUBE = BRYDNS*YE 
!                                                                       
!                   Call electron EOS to determine the                  
!                   electron chemical potential                         
      CALL EL_EOS(T,YE,BRYDNS) 
!                                                                       
!                                                                       
      MUHAT = MUSUBE+(1.0-PHASEF)*(MUHLOW-MUELOW)+PHASEF*(MUHHI-MUEHI) 
!                                                                       
      MUN = MUTILD+YE*(MUHAT-MUSUBE) 
!                                                                       
      MUPROT = MUN-MUHAT 
!                                                                       
!                   Calculate thermodynamic quantities                  
!                                                                       
      STOT = ((1.0-PHASEF)*S_LOW*LOWDNS+PHASEF*S_HI*HIDNS)/BRYDNS+PS 
!                                                                       
      FTOT = (LOWDNS*F_LOW+MUTILD*(BRYDNS-LOWDNS))/BRYDNS+PF 
!                                                                       
      UTOT = FTOT+T*STOT+PU 
!                                                                       
      PTOT = PRTILD+PPRESS 
!                                                                       
!                                                                       
      XH = (1.0-PHASEF)*N_XH 
      XALFA = (1.0-PHASEF)*N_XA 
      XNUT = (1.0-PHASEF)*N_XN 
      XPROT = (1.0-PHASEF)*N_XP 
      XALFA2 = PHASEF*B_XA 
      XNUT2 = PHASEF*B_XN 
      XPROT2 = PHASEF*B_XP 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DELDNS = HIDNS-LOWDNS 
!                                                                       
!                                                                       
      DPHADT = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DT-           &
     &    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DT                            
!                                                                       
      DPDT = DPN_DT+DPN_DN*DNL_DT 
      DMUDT = DMU_DT+DMU_DN*DNL_DT 
      DSDT = (1.0-PHASEF)*LOWDNS*(DSN_DT+DSN_DN*DNL_DT)/BRYDNS+         &
     & (1.0-PHASEF)*S_LOW*DNL_DT/BRYDNS-LOWDNS*S_LOW*DPHADT/BRYDNS+     &
     &    (DPHADT*S_HI*HIDNS+PHASEF*DNH_DT*S_HI+                        &
     &    PHASEF*HIDNS*(DSB_DT+DSB_DN*DNH_DT))/BRYDNS+DPSDT             
      DUDT = DMUDT-DPDT/BRYDNS+STOT+T*DSDT 
!                                                                       
!                                                                       
      DPDN = 0.0 
      DMUDN = 0.0 
      DSDN = -DPDT/BRYDNS**2 
      DUDN = (LOWDNS*(MUTILD-FTOT)/BRYDNS**2)+T*DSDN 
!                                                                       
!                                                                       
      DPHADY = ((BRYDNS-LOWDNS)/DELDNS**2-1.0/DELDNS)*DNL_DY-           &
     &    ((BRYDNS-LOWDNS)/DELDNS**2)*DNH_DY                            
!                                                                       
      DPDY = DPN_DY+DPN_DN*DNL_DY 
      DMUDY = DMU_DY+DMU_DN*DNL_DY 
      DSDY = (1.0-PHASEF)*LOWDNS*(DSN_DY+DSN_DN*DNL_DY)/BRYDNS+         &
     & (1.0-PHASEF)*S_LOW*DNL_DY/BRYDNS-LOWDNS*S_LOW*DPHADY/BRYDNS+     &
     &    (DPHADY*S_HI*HIDNS+PHASEF*DNH_DY*S_HI+                        &
     &    PHASEF*HIDNS*(DSB_DY+DSB_DN*DNH_DY))/BRYDNS                   
      DUDY = DMUDY-DPDY/BRYDNS+T*DSDY 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!             Adiabatic index                                           
!             (Note that the first term vanishes in this expression)    
      GAM_S = T*(DPDT**2)/(BRYDNS*PTOT*DUDT) 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EOSLOG.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EOSLOG                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/15/90                                             
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL EOSLOG(INPVAR,YE,BRYDNS,EOSFLG)                 
!                                                                       
!                                                                       
!    INPUTS:       INPVAR = TEMP, INTERNAL ENG, OR ENTROPY              
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!                                                                       
!                                                                       
!    OUTPUTS:      EOSFLG = 1 --> Not implemented in model 4B           
!                           2 --> GENERAL EOS                           
!                           3 --> BULK EOS (includes alpha's)           
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC, MAXWEL.INC                            
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EOSLOG(INPVAR,YE,BRYDNS,EOSFLG) 
!                                                                       
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                       This include file contains all variable         
!                       declarations.  NOTE:: no implicit typing        
!                       scheme is followed in this code; if you         
!                       have any doubt as to a variables type CHECK     
!                       IT!!!!.  Also note that ALL variables are       
!                       declared explicitly.                            
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
      DOUBLE PRECISION NLOW, NHI, N_CUT, TEMP_1, TEMP_2, T_BNDY 
!                                                                       
      DOUBLE PRECISION LMM, LMP, LPM, LPP 
      DOUBLE PRECISION DNDY1, DNDY2 
!                                                                       
!                                                                       
!                                                                       
   10 CONTINUE 
!                                                                       
!                                                                       
!                         Set T equal to the input variable (any calls  
!                         with entropy or internal energy should go     
!                         the the EOS_M4B subroutine)                   
!                                                                       
      T = INPVAR(1) 
!                                                                       
!-----------------------------------------------------------------------
!         code to figure out the boundaries from the tables             
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      IF(YE.GT.Y_HI) THEN 
!                         Ye is too large for EOS                       
!                                                                       
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time' 
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_HI,' instead' 
        YE = Y_HI 
        GOTO 10 
!                                                                       
      ELSEIF(YE.GT.Y_LOW) THEN 
!                         Calculate high and low boundary densities     
!                         for the Maxwell construction                  
!                                                                       
!----------------------------------------------------------             
!           Calc Ye index                                               
!----------------------------------------------------------             
!                                                                       
        YFRAC = (YE-Y_LOW)/(Y_HI-Y_LOW) 
        J_MXWL = INT(YFRAC*(NUMYE-1))+1 
        DELT_Y = (Y_HI-Y_LOW)/DBLE(NUMYE-1) 
!                                                                       
        YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
        YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
!                                                                       
!                                                                       
        IF((YE.GE.YMINUS).AND.(YE.LE.YPLUS)) THEN 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
        ELSEIF(YE.GT.YPLUS) THEN 
          J_MXWL = J_MXWL+1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ELSE 
          J_MXWL = J_MXWL-1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ENDIF 
!                                                                       
!                                                                       
        IF(J_MXWL.GT.(NUMYE-1)) THEN 
          J_MXWL = NUMYE-1 
          J_BD = J_MXWL 
          J_BNDY = J_MXWL 
          YMINUS = Y_LOW+DBLE(J_MXWL-1)*DELT_Y 
          YPLUS = Y_LOW+DBLE(J_MXWL)*DELT_Y 
        ENDIF 
!                                                                       
!                                                                       
        YINTRP = (YE-YMINUS)/(YPLUS-YMINUS) 
!                                                                       
!----------------------------------------------------------             
!           Calc T index                                                
!----------------------------------------------------------             
!                                                                       
!                                                                       
        TFRAC = (T-T_LOW)/(T_HI-T_LOW) 
        I_MXWL = INT(TFRAC*(NUMTMP-1))+1 
        DELT_T = (T_HI-T_LOW)/DBLE(NUMTMP-1) 
!                                                                       
        TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
        TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
!                                                                       
!                                                                       
        IF((T.GT.TMINUS).AND.(T.LE.TPLUS)) THEN 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ELSEIF(T.GT.TPLUS) THEN 
          I_MXWL = I_MXWL+1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ELSE 
          I_MXWL = I_MXWL-1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ENDIF 
!                                                                       
!                                                                       
        IF(I_MXWL.GT.(NUMTMP-1)) THEN 
          I_MXWL = NUMTMP-1 
          TMINUS = T_LOW+DBLE(I_MXWL-1)*DELT_T 
          TPLUS = T_LOW+DBLE(I_MXWL)*DELT_T 
        ENDIF 
!                                                                       
!                                                                       
        TINTRP = (T-TMINUS)/(TPLUS-TMINUS) 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                Find the temperature and density at the top of the     
!                Maxwel construction                                    
!                                                                       
!C      T_MXWL = YINTRP*(T_H(J_MXWL+1)-T_H(J_MXWL))+T_H(J_MXWL)         
!C      D_MXWL = YINTRP*(D_H(J_MXWL+1)-D_H(J_MXWL))+D_H(J_MXWL)         
        T_MXWL = DMIN1(T_H(J_MXWL+1),T_H(J_MXWL)) 
        IF(T_H(J_MXWL+1).GT.T_H(J_MXWL)) THEN 
          D_MXWL = D_H(J_MXWL) 
        ELSE 
          D_MXWL = D_H(J_MXWL+1) 
        ENDIF 
!                                                                       
!                                                                       
!                                                                       
!--------------------------------------------------------------------   
!            Interpolate to get Maxwell construction densities          
!--------------------------------------------------------------------   
!                                                                       
!                                                                       
!                                                                       
        DNS_1 = YINTRP*(BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))+ &
     &               BRYLOW(I_MXWL,J_MXWL)                              
        DNS_2 = YINTRP*                                                 &
     &        (BRYLOW(I_MXWL+1,J_MXWL+1)-BRYLOW(I_MXWL+1,J_MXWL))+      &
     &               BRYLOW(I_MXWL+1,J_MXWL)                            
!                                                                       
        LOWDNS = TINTRP*(DNS_2-DNS_1)+DNS_1 
!                                                                       
!                Derivative of lower density w.r.t. T                   
        DNL_DT = (DNS_2-DNS_1)/DELT_T 
!                                                                       
        DNDY1 = (BRYLOW(I_MXWL,J_MXWL+1)-BRYLOW(I_MXWL,J_MXWL))/DELT_Y 
        DNDY2 = (BRYLOW(I_MXWL+1,J_MXWL+1)-                             &
     &      BRYLOW(I_MXWL+1,J_MXWL))/DELT_Y                             
        DNL_DY = TINTRP*(DNDY2-DNDY1)+DNDY1 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
        IF(YE.GT.Y_CUT) THEN 
!                                                                       
          DNS_1 = YINTRP*                                               &
     &        (BRYHI(I_MXWL,J_MXWL+1)-BRYHI(I_MXWL,J_MXWL))+            &
     &        BRYHI(I_MXWL,J_MXWL)                                      
          DNS_2 = YINTRP*                                               &
     &        (BRYHI(I_MXWL+1,J_MXWL+1)-BRYHI(I_MXWL+1,J_MXWL))+        &
     &               BRYHI(I_MXWL+1,J_MXWL)                             
!                                                                       
          HIDNS = TINTRP*(DNS_2-DNS_1)+DNS_1 
!                                                                       
!                Derivative of higher density w.r.t. T                  
          DNH_DT = (DNS_2-DNS_1)/DELT_T 
!                                                                       
!                                                                       
        DNDY1 = (BRYHI(I_MXWL,J_MXWL+1)-                                &
     &      BRYHI(I_MXWL,J_MXWL))/DELT_Y                                
        DNDY2 = (BRYHI(I_MXWL+1,J_MXWL+1)-                              &
     &      BRYHI(I_MXWL+1,J_MXWL))/DELT_Y                              
        DNH_DY = TINTRP*(DNDY2-DNDY1)+DNDY1 
!                                                                       
!                                                                       
        ELSE 
          HIDNS = LOWDNS 
        ENDIF 
!                                                                       
!                                                                       
!--------------------------------------------------------------------   
!--------------------------------------------------------------------   
!                                                                       
!                       Ye is too low                                   
      ELSE 
        WRITE(*,*) ' EOSLOG:: Cant do Ye = ',YE, 'at this time' 
        stop 
        WRITE(*,*) ' EOSLOG:: assuming YE =',Y_LOW,' instead' 
        YE = Y_LOW 
        GOTO 10 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      DLTLN1 = (LNCUT-LNLOW)/DBLE(NUMLOW-1) 
      DLTLN2 = (LNHI-LNCUT)/DBLE(NUMHI-1) 
!                                                                       
!                                                                       
      NLOW = 10.0**LNLOW 
      NHI = 10.0**LNHI 
      N_CUT = 10.0**LNCUT 
      LOGBRY = DLOG10(BRYDNS) 
      LOGBCH = LOGBRY 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!----------------------------------------------------------             
!           Calc T index                                                
!----------------------------------------------------------             
!                                                                       
!                                                                       
      IF(LOGBRY.GE.LNHI) THEN 
        I_BD = NBPNTS 
        I_BNDY = NBPNTS 
        T_BNDY = YINTRP*                                                &
     &           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+       &
     &            LBOUND(I_BNDY,J_BNDY)                                 
        GOTO 70 
      ELSEIF((LOGBRY.LT.LNHI).AND.(LOGBRY.GT.LNCUT)) THEN 
!                                                                       
        I_BD = INT((LOGBRY-LNCUT)/DLTLN2)+NUMLOW 
        LNMINS = LNCUT+DBLE(I_BD-NUMLOW)*DLTLN2 
        LNPLUS = LNCUT+DBLE(I_BD-NUMLOW+1)*DLTLN2 
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN 
          I_BNDY = I_BD 
        ELSEIF(LOGBCH.GT.LNPLUS) THEN 
          I_BD = I_BD+1 
          I_BNDY = I_BD 
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
        ELSE 
          I_BD = I_BD-1 
          I_BNDY = I_BD 
          LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
          LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
        ENDIF 
!                                                                       
      ELSEIF((LOGBRY.LE.LNCUT).AND.(LOGBRY.GT.LNLOW)) THEN 
!                                                                       
        I_BD = INT((LOGBRY-LNLOW)/DLTLN1)+1 
        LNMINS = LNLOW+DBLE(I_BD-1)*DLTLN1 
        LNPLUS = LNLOW+DBLE(I_BD)*DLTLN1 
        IF((LOGBCH.LE.LNPLUS).AND.(LOGBCH.GE.LNMINS)) THEN 
          I_BNDY = I_BD 
        ELSEIF(LOGBCH.GT.LNPLUS) THEN 
          I_BD = I_BD+1 
          I_BNDY = I_BD 
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1 
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1 
        ELSE 
          I_BD = I_BD-1 
          I_BNDY = I_BD 
          LNMINS = LNLOW+DBLE(I_BNDY-1)*DLTLN1 
          LNPLUS = LNLOW+DBLE(I_BNDY)*DLTLN1 
        ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      IF(I_BNDY.GT.(NBPNTS-1)) THEN 
        I_BD = NBPNTS-1 
        I_BNDY = I_BD 
        LNMINS = LNCUT+DBLE(I_BNDY-NUMLOW)*DLTLN2 
        LNPLUS = LNCUT+DBLE(I_BNDY-NUMLOW+1)*DLTLN2 
      ENDIF 
!                                                                       
!                                                                       
!                                                                       
      LMM = LBOUND(I_BNDY,J_BNDY) 
      LPM = LBOUND(I_BNDY+1,J_BNDY) 
      LMP = LBOUND(I_BNDY,J_BNDY+1) 
      LPP = LBOUND(I_BNDY+1,J_BNDY+1) 
!                                                                       
      LNFRAC = (LOGBCH-LNMINS)/(LNPLUS-LNMINS) 
!                                                                       
!                Interpolate in Ye first                                
!                                                                       
      TEMP_1 = YINTRP*                                                  &
     &           (LBOUND(I_BNDY,J_BNDY+1)-LBOUND(I_BNDY,J_BNDY))+       &
     &            LBOUND(I_BNDY,J_BNDY)                                 
      TEMP_2 = YINTRP*                                                  &
     &        (LBOUND(I_BNDY+1,J_BNDY+1)-LBOUND(I_BNDY+1,J_BNDY))+      &
     &               LBOUND(I_BNDY+1,J_BNDY)                            
!                                                                       
!                Interpolate in density between the two Ye              
!                interpolated values                                    
!                                                                       
      T_BNDY = LNFRAC*(TEMP_2-TEMP_1)+TEMP_1 
!                                                                       
!                                                                       
!----------------------------------------------------------             
!----------------------------------------------------------             
!                                                                       
   70 CONTINUE 
!                                                                       
      TCHK_B = 1.01*T_BNDY 
      TCHK_N = 0.95*T_BNDY 
!                                                                       
      IF((LMM.GE.LPM).OR.(LMP.GT.LPP)) THEN 
        TCHK_N = DMAX1(0.0D0,DMIN1(0.95*TCHK_N,T_BNDY-3.0)) 
      ENDIF 
!                                                                       
!-----------------------------------------------------------------------
!               EOS Logic                                               
!-----------------------------------------------------------------------
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                     If T is below the maximum for maxwel construction 
      IF(T.LT.T_MXWL) THEN 
!                       If rho is greater than the upper max. con.      
!                       density the use the bulk EOS                    
        IF(BRYDNS.GT.HIDNS) THEN 
          EOSFLG = 3 
!                       Else if rho is greater than the lower max. con. 
!                       density then                                    
        ELSEIF(BRYDNS.GT.LOWDNS) THEN 
!                         If Ye is large enough to have a signifigant   
!                         max con then use the maxwell con. EOS         
          IF(YE.GT.Y_CUT) THEN 
            EOSFLG = 4 
!                         Otherwise use the bulk EOS                    
          ELSE 
            EOSFLG = 3 
          ENDIF 
!                                                                       
!                       If density is greater than the minimum          
!                       Maxwell con. density, then we know that we are  
!                       in the Nuclear EOS density                      
        ELSEIF(BRYDNS.GT.D_MXWL) THEN 
          EOSFLG = 2 
!                                                                       
!                                                                       
!                       Otherwise check the Boundary table              
        ELSE 
!                                                                       
!                         If T is well below the phase boundary curve   
!                         then use the nuclear EOS                      
          IF(T.LT.TCHK_N) THEN 
            EOSFLG = 2 
!                         Otherwise if T is near the boundary, first    
!                         try the nuclear EOS and if not successfull    
!                         then use the bulk EOS                         
          ELSEIF(T.LT.TCHK_B) THEN 
            EOSFLG = 1 
          ELSE 
!                         Otherwise T is well above the boundary so     
!                         use the bulk EOS                              
            EOSFLG = 3 
          ENDIF 
        ENDIF 
!                                                                       
!                     Otherwise T is above the maximum for a maxwell    
!                     construction                                      
      ELSE 
!                       If density is greater than that at the top of   
!                       the maxwell construction then use the bulk EOS  
        IF(BRYDNS.GT.D_MXWL) THEN 
          EOSFLG = 3 
!                                                                       
!                       Otherwise density is below the maxwell con.     
        ELSE 
!                                                                       
!                         If T is well below the phase boundary curve   
!                         then use the nuclear EOS                      
          IF(T.LT.TCHK_N) THEN 
            EOSFLG = 2 
!                                                                       
!                         Otherwise if T is near the phase boundary     
!                         curve then try the nuclear EOS and if not     
!                         successfull then use the bulk EOS             
          ELSEIF(T.LT.TCHK_B) THEN 
            EOSFLG = 1 
!                                                                       
!                         Otherwise T is well above the phase boundary  
!                         curve so use the bulk EOS                     
          ELSE 
            EOSFLG = 3 
          ENDIF 
        ENDIF 
      ENDIF 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                         Done with EOS logic so return EOSFLG          
!-----------------------------------------------------------------------
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         RESET.FOR                                            
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       RESET                                                
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         12/21/90                                             
!                                                                       
!                  Please report any problems to me at:                 
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@ASTRO.SUNYSB.EDU or                
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!                                                                       
!    CALL LINE:    CALL RESET(INPVAR,YE,BRYDNS,OUTVAR)                  
!                                                                       
!                                                                       
!    INPUTS:       INPVAR = TEMP, NSUBI, ETA_PO, ETA_NO                 
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!                                                                       
!                                                                       
!    OUTPUTS:      OUTVAR = ARRAY OF LENGTH 4 CONTAINING RESET VALUES   
!                  FOR THE INITIAL GUESSES                              
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES: NONE                                                
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE RESET(INPVAR,YE,BRYDNS,OUTVAR) 
!                                                                       
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
!                      Subroutine parameters                            
!                                                                       
      DOUBLE PRECISION INPVAR(4), OUTVAR(4), YE, BRYDNS 
!                                                                       
!                                                                       
!                      Local variables                                  
!                                                                       
      DOUBLE PRECISION ZPG, ZNG, ETA_PG, ETA_NG, PI, UQ, MQ, T, EFRAC 
!                                                                       
!                      Functions                                        
!                                                                       
      DOUBLE PRECISION FINV12 
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
      T = INPVAR(1) 
!                                                                       
!                                                                       
      PI = 3.1415927 
      UQ = 20.721 
!                                                                       
      MQ = (T/UQ)**1.5 
!                                                                       
!                                                                       
      EFRAC = 0.5*YE 
!                                                                       
      ZNG = 2.0*(PI**2)*BRYDNS*(1.0-EFRAC)/MQ 
!                                                                       
      ZPG = 2.0*(PI**2)*BRYDNS*EFRAC/MQ 
!                                                                       
      ETA_NG = FINV12(ZNG) 
!                                                                       
      ETA_PG = FINV12(ZPG) 
!                                                                       
      OUTVAR(1) = INPVAR(1) 
      OUTVAR(2) = INPVAR(2) 
      OUTVAR(3) = ETA_PG 
      OUTVAR(4) = ETA_NG 
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         ALOADMX.FOR                                          
!    MODULE:       LOADMX                                               
!    TYPE:         LOADMX                                               
!                                                                       
!    PURPOSE:      LOAD THE LOOK-UP TABLE FOR THE MAXWELL CONSTRUCTION  
!                                                                       
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         7/16/90                                              
!                                                                       
!    CALL LINE:    CALL LOADMX                                          
!                                                                       
!    INPUTS:       N/A                                                  
!                                                                       
!    OUTPUTS       N/A                                                  
!                                                                       
!    SUBROUTINE CALLS: EOS_M4A                                          
!                                                                       
!    INCLUDE FILES:  EOS_M4A.INC, MAXWEL.INC                            
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE LOADMX() 
!                                                                       
      IMPLICIT NONE 
!                                                                       
!                                                                       
      INCLUDE 'eos_m4a.inc' 
      INCLUDE 'maxwel.inc' 
!                                                                       
      INTEGER NTMP, NYE, NYE2, NUM_BP 
      INTEGER LUN1, LUN2, KK, KMIN 
      PARAMETER(LUN1=44,LUN2=45) 
!                                                                       
!                                                                       
      INTEGER FNML1, FNML2 
      CHARACTER*60 FNAME1, FNAME2 
!                                                                       
      DOUBLE PRECISION N_SM, SYMM_M, COMP_M, BINDEM, SYM_SM, SIG_SM 
      DOUBLE PRECISION N_SB, SYMM_B, COMP_B, BINDEB, SYM_SB, SIG_SB 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!******below commented out by MBD, always read in max180.atb            
!******and bd180.atb                                                    
!      CALL GETFNM(FNAME1,FNML1,'Enter ASCII Maxwell fname:',26)        
!                                                                       
!      CALL GETFNM(FNAME2,FNML2,'Enter ASCII boundary fname:',27)       
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!        Read the file Maxwell construction data file                   
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
      OPEN(UNIT=LUN1,FILE='max180.atb',STATUS='OLD') 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN1,*) N_SM, SYMM_M 
      READ(LUN1,*) COMP_M,BINDEM 
      READ(LUN1,*) SYM_SM, SIG_SM 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN1,*) NTMP,NYE 
      READ(LUN1,*) T_LOW,T_HI 
      READ(LUN1,*) Y_LOW,Y_HI 
!                                                                       
!                                                                       
!                                                                       
      IF((NTMP.NE.NUMTMP).OR.(NYE.NE.NUMYE)) THEN 
        WRITE(*,*) 'LOADMX:  MXWL TABLE IS INCOMPATIBLE WITH ARRAYS' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      DO 101 J=1,NUMYE,1 
        DO 100 I=1,NUMTMP,3 
          KMIN = MIN0(I+2,NUMTMP) 
          READ(LUN1,*) (BRYLOW(KK,J),KK=I,KMIN,1) 
  100   CONTINUE 
  101 END DO 
!                                                                       
!                                                                       
      DO 103 J=1,NUMYE,1 
        DO 102 I=1,NUMTMP,3 
          KMIN = MIN0(I+2,NUMTMP) 
          READ(LUN1,*) (BRYHI(KK,J),KK=I,KMIN,1) 
  102   CONTINUE 
  103 END DO 
!                                                                       
!                                                                       
!                                                                       
      DO 104 I=1,NUMYE,3 
        KMIN = MIN0(I+2,NUMYE) 
        READ(LUN1,*) (T_H(KK),KK=I,KMIN,1) 
  104 END DO 
!                                                                       
!                                                                       
      DO 105 I=1,NUMYE,3 
        KMIN = MIN0(I+2,NUMYE) 
        READ(LUN1,*) (D_H(KK),KK=I,KMIN,1) 
  105 END DO 
!                                                                       
      READ(LUN1,*) YCUT 
!                                                                       
!                                                                       
      CLOSE(UNIT=LUN1,STATUS='KEEP') 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  MAXWELL CON. TABLE IS INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!        Read the file Boundary data file                               
!-----------------------------------------------------------------------
!                                                                       
!                                                                       
!                                                                       
      OPEN(UNIT=LUN2,FILE='bd180.atb',STATUS='OLD') 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN2,*) N_SB,SYMM_B 
      READ(LUN2,*) COMP_B,BINDEB 
      READ(LUN2,*) SYM_SB,SIG_SB 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
      READ(LUN2,*) NUM_BP,NYE2 
      READ(LUN2,*) LNL,LNH,LNC 
      READ(LUN2,*) Y_LOW2,Y_HI2 
!                                                                       
!                                                                       
      IF((NBPNTS.NE.NUM_BP).OR.(NYE2.NE.NUMYE)) THEN 
        WRITE(*,*) 'LOADMX:  BNDY TABLE IS INCOMPATIBLE WITH ARRAYS' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(LNL-LNLOW).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  LOWER END OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      IF(ABS(LNH-LNHI).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  UPPER END OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      IF(ABS(LNC-LNCUT).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  MID CUT OF PHASE BNDY IS INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(Y_LOW-Y_LOW2).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  LOWER YE LIMITS ARE INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
      IF(ABS(Y_HI-Y_HI2).GT.1.0D-10) THEN 
        WRITE(*,*) 'LOADMX:  UPPER YE LIMITS ARE INCONSIST.' 
        STOP 
      ENDIF 
!                                                                       
!                                                                       
      DO 201 J=1,NUMYE,1 
        DO 200 I=1,NBPNTS,3 
          KMIN = MIN0(I+2,NBPNTS) 
          READ(LUN2,*) (LBOUND(KK,J),KK=I,KMIN,1) 
  200   CONTINUE 
  201 END DO 
!                                                                       
!                                                                       
      DO 203 J=1,NUMYE,1 
        DO 202 I=1,NBPNTS,3 
          KMIN = MIN0(I+2,NBPNTS) 
          READ(LUN2,*) (UBOUND(KK,J),KK=I,KMIN,1) 
  202   CONTINUE 
  203 END DO 
!                                                                       
!                                                                       
!                                                                       
      CLOSE(UNIT=LUN2,STATUS='KEEP') 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  BOUNDARY TABLE IS INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
!                                                                       
!-----------------------------------------------------------------------
!                  All arrays are now loaded so return                  
!-----------------------------------------------------------------------
!                                                                       
      N_S = N_SM 
      NSUBS = N_SM 
      SYMM = SYMM_M 
      COMP = COMP_M 
      BIND_E = BINDEM 
      SYM_S = SYM_SM 
      SIG_S = SIG_SM 
!                                                                       
      SKYRMC=(.3*((HBAR*C)**2)/MASSN)*(1.5*N_S*(PI**2))**OVR23 
      DD = (COMP+2.0*SKYRMC)/(3.0*SKYRMC+9.0*BIND_E) 
      BB = (SKYRMC*(2.0**OVR23-1.0)-SYMM)/N_S 
      AA = (OVR23*SKYRMC-DD*(SKYRMC+BIND_E))/(N_S*(DD-1.0))-BB 
      CC = (COMP+2.0*SKYRMC)/(9.0*DD*(DD-1.0)*N_S**DD) 
!                                                                       
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX:  SKYRME PARAMETERS FOR THIS RUN ARE:>>' 
      WRITE(*,*) 'ABCD: ',AA,BB,CC,DD 
      WRITE(*,*) ' Satur. density, symmetry engy, & compression mod.:' 
      WRITE(*,*) N_SM, SYMM_M, COMP_M 
      WRITE(*,*) N_SB, SYMM_B, COMP_B 
      WRITE(*,*) ' Binding engy, surf. symm. engy, & surface tension:' 
      WRITE(*,*) BINDEM,SYM_SM,SIG_SM 
      WRITE(*,*) BINDEB,SYM_SB,SIG_SB 
!                                                                       
      WRITE(*,*) 
!                                                                       
!                                                                       
      CALL INITFERM() 
!                                                                       
      WRITE(*,*) 
      WRITE(*,*) '<<LOADMX: FERMI INTEGRAL TABLES ARE INITIALIZED>>' 
      WRITE(*,*) 
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       GETFNM.FOR                                           
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         8/5/89                                               
!                                                                       
!    PURPOSE:      OBTAINS A FILE NAME FROM USER                        
!                                                                       
!    CALL LINE:    CALL GETFNM(FNAME,FNML,PROMPT,PROMTL)                
!                                                                       
!    INPUTS:       PROMPT = STRING TO POMPT USER WITH (C*60)            
!                  PROMTL = LENGTH OF PROMPT STRING (I)                 
!                                                                       
!    OUTPUTS:      FNAME = FILE NAME (C*60)                             
!                  FNML = FILE NAME LENGTH (I)                          
!***********************************************************************
!                                                                       
      SUBROUTINE GETFNM(FNAME,FNML,PROMPT,PROMTL) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER FNML, PROMTL 
      CHARACTER*60 FNAME, PROMPT 
!                                                                       
!                       Local variables                                 
!                                                                       
      INTEGER STDOUT, STDIN, I 
      DATA STDOUT/6/, STDIN/5/ 
!                                                                       
!                       Prompt user for file name                       
!                                                                       
      WRITE(STDOUT,'(T2,A,$)') PROMPT(1:PROMTL) 
      READ(STDIN,'(A)') FNAME 
!                                                                       
!                        Figure out input file name length              
      DO 10 I=1,20,1 
        IF(FNAME(I:I).EQ.' ') GOTO 20 
   10 END DO 
!                                                                       
   20 CONTINUE 
      FNML = I-1 
!                                                                       
!                                                                       
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         DMATRIX.FOR                                          
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATINV                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Inverts a N by N Matrix                              
!                                                                       
!    CALL LINE:    CALL MATINV(A,AINV,N)                                
!                                                                       
!    INPUTS:       A = Array to be inverted  (D)                        
!                  N = dimesion of arrays (I)                           
!                                                                       
!    OUTPUTS:      AINV = Inverse of A (D)                              
!                                                                       
!    CALLS :       Numerical recipes routines LUDCMP, LUBKSB            
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATINV(A,AINV,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N,N), AINV(N,N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER NPHYS, I, J 
      PARAMETER(NPHYS=10) 
      DOUBLE PRECISION TEMP(NPHYS,NPHYS), Y(NPHYS,NPHYS), D 
      INTEGER INDEX(NPHYS) 
!                                                                       
!                 Make a copy of the array, and initialize              
!                 the indentity matrix                                  
      DO 20 J=1,N,1 
        DO 10 I=1,N,1 
          Y(I,J) = 0.0 
          TEMP(I,J) = A(I,J) 
   10   CONTINUE 
        Y(J,J) = 1.0 
   20 END DO 
!                                                                       
!                                                                       
!                 LU decompose the matrix                               
      CALL LUDCMP(TEMP,N,NPHYS,INDEX,D) 
!                                                                       
!                                                                       
!                 Back substitute to get inverse                        
      DO 30 J=1,N,1 
        CALL LUBKSB(TEMP,N,NPHYS,INDEX,Y(1,J)) 
   30 END DO 
!                                                                       
!                                                                       
!                 Copy temporary array into the inverse array           
      DO 50 J=1,N,1 
        DO 40 I=1,N,1 
          AINV(I,J) = Y(I,J) 
   40   CONTINUE 
   50 END DO 
!                                                                       
!                                                                       
  999 RETURN 
      END                                           
!***********************************************************************
!                                                                       
!    MODULE:       MATADD                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Adds two N by M Matrices                             
!                                                                       
!    CALL LINE:    CALL MATINV(A,B,C,N,M)                               
!                                                                       
!    INPUTS:       A,B= Arrays to be added  (D)                         
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      C = Array containing A+B (D)                         
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATADD(A,B,C,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N,M 
      DOUBLE PRECISION A(N,M), B(N,M), C(N,M) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
          C(I,J) = A(I,J)+B(I,J) 
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATMUL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiplies two Matrices (LxM)x(MxN)                  
!                                                                       
!    CALL LINE:    CALL MATMUL(A,B,C,L,M,N)                             
!                                                                       
!    INPUTS:       A,B= Arrays to be added  (D)                         
!                  L,M,N = Dimensions of arrays (I)                     
!                                                                       
!    OUTPUTS:      C = Array containing A x B (D)                       
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATMUL(A,B,C,L,M,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER L, M, N 
      DOUBLE PRECISION A(L,M), B(M,N), C(L,N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J, K 
      DOUBLE PRECISION SUM 
!                                                                       
!                 Loop over all elements of the array                   
      DO 30 I=1,L,1 
        DO 20 J=1,N,1 
!                                                                       
!                 Initialize SUM for a new element                      
          SUM = 0.0 
!                 Calculate (i,j)th element                             
          DO 10 K=1,M,1 
            SUM = SUM+A(I,K)*B(K,J) 
   10     CONTINUE 
          C(I,J) = SUM 
!                                                                       
   20   CONTINUE 
   30 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!***********************************************************************
!                                                                       
!    MODULE:       MV_MUL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiplies a Matrix times a vector (NxN)x(N)         
!                                                                       
!    CALL LINE:    CALL MV_MUL(A,V,RV,N)                                
!                                                                       
!    INPUTS:       A = Array to be multiplied  (D)                      
!                  V = Vector to be multiplied (D)                      
!                  N = Dimensions of arrays & vector (I)                
!                                                                       
!    OUTPUTS:      RV = resultant vector (D)                            
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MV_MUL(A,V,RV,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N,N), V(N), RV(N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
      DOUBLE PRECISION SUM 
!                                                                       
!                 Loop over all elements of the array                   
      DO 20 I=1,N,1 
!                                                                       
!                 Initialize SUM for a new element                      
        SUM = 0.0 
!                 Calculate (i)th element                               
        DO 10 J=1,N,1 
          SUM = SUM+A(I,J)*V(J) 
   10   CONTINUE 
        RV(I) = SUM 
!                                                                       
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATSCL                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Multiply a N by M Matrix by a scalar                 
!                                                                       
!    CALL LINE:    CALL MATSCL(A,SCALAR,B,N,M)                          
!                                                                       
!    INPUTS:       A = Array to be scaled  (D)                          
!                  SCALAR = Constant to multiply matrix by (D)          
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      B = Array containing SCALAR x A (D)                  
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATSCL(A,SCALAR,B,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N, M 
      DOUBLE PRECISION A(N,M), B(N,M), SCALAR 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
!                 Loop over all elements of the array                   
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
!                                                                       
          B(I,J) = SCALAR*A(I,J) 
!                                                                       
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       MATCOP                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/3/90                                               
!                                                                       
!    PURPOSE:      Copy one N by M Matrix into another                  
!                                                                       
!    CALL LINE:    CALL MATCOP(A,B,N,M)                                 
!                                                                       
!    INPUTS:       A = Array to be copied  (D)                          
!                  N = Number of rows in arrays (I)                     
!                  M = Number of columns in arrays (I)                  
!                                                                       
!    OUTPUTS:      C = Array to be copied into (D)                      
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE MATCOP(A,B,N,M) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N,M 
      DOUBLE PRECISION A(N,M), B(N,M) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I, J 
!                                                                       
      DO 20 J=1,M,1 
        DO 10 I=1,N,1 
          B(I,J) = A(I,J) 
   10   CONTINUE 
   20 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       VECCOP                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         4/10/90                                              
!                                                                       
!    PURPOSE:      Copy one vector on length N into another             
!                                                                       
!    CALL LINE:    CALL VECCOP(A,B,N)                                   
!                                                                       
!    INPUTS:       A = Array to be copied  (D)                          
!                  N = Number of rows in arrays (I)                     
!                                                                       
!    OUTPUTS:      B = Array to be copied into (D)                      
!                                                                       
!    CALLS :       None                                                 
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE VECCOP(A,B,N) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      INTEGER N 
      DOUBLE PRECISION A(N), B(N) 
!                                                                       
!                                                                       
!                 Local variables                                       
!                                                                       
      INTEGER I 
!                                                                       
      DO 10 I=1,N,1 
        B(I) = A(I) 
   10 END DO 
!                                                                       
  999 RETURN 
!                                                                       
      END                                           
!                                                                       
!                                                                       
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       LUDCMP                                               
!    TYPE:         SUBROUTINE                                           
!                                                                       
!***********************************************************************
!                                                                       
      SUBROUTINE LUDCMP(A,N,NP,INDX,D) 
!                                                                       
!               Make the default double precision                       
      implicit real*8(a-h,o-z) 
!                                                                       
!                                                                       
      PARAMETER (NMAX=100,TINY=1.0E-20) 
      DIMENSION A(NP,NP),INDX(N),VV(NMAX) 
      D=1. 
      DO 12 I=1,N 
        AAMAX=0. 
        DO 11 J=1,N 
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J)) 
   11   CONTINUE 
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.' 
        VV(I)=1./AAMAX 
   12 END DO 
      DO 19 J=1,N 
        IF (J.GT.1) THEN 
          DO 14 I=1,J-1 
            SUM=A(I,J) 
            IF (I.GT.1)THEN 
              DO 13 K=1,I-1 
                SUM=SUM-A(I,K)*A(K,J) 
   13         CONTINUE 
              A(I,J)=SUM 
            ENDIF 
   14     CONTINUE 
        ENDIF 
        AAMAX=0. 
        DO 16 I=J,N 
          SUM=A(I,J) 
          IF (J.GT.1)THEN 
            DO 15 K=1,J-1 
              SUM=SUM-A(I,K)*A(K,J) 
   15       CONTINUE 
            A(I,J)=SUM 
          ENDIF 
          DUM=VV(I)*ABS(SUM) 
          IF (DUM.GE.AAMAX) THEN 
            IMAX=I 
            AAMAX=DUM 
          ENDIF 
   16   CONTINUE 
        IF (J.NE.IMAX)THEN 
          DO 17 K=1,N 
            DUM=A(IMAX,K) 
            A(IMAX,K)=A(J,K) 
            A(J,K)=DUM 
   17     CONTINUE 
          D=-D 
          VV(IMAX)=VV(J) 
        ENDIF 
        INDX(J)=IMAX 
        IF(J.NE.N)THEN 
          IF(A(J,J).EQ.0.)A(J,J)=TINY 
          DUM=1./A(J,J) 
          DO 18 I=J+1,N 
            A(I,J)=A(I,J)*DUM 
   18     CONTINUE 
        ENDIF 
   19 END DO 
      IF(A(N,N).EQ.0.)A(N,N)=TINY 
      RETURN 
      END                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       LUBKSB                                               
!    TYPE:         SUBROUTINE                                           
!                                                                       
!***********************************************************************
!                                                                       
!                                                                       
      SUBROUTINE LUBKSB(A,N,NP,INDX,B) 
!                                                                       
!               Make the default double precision                       
      implicit real*8(a-h,o-z) 
!                                                                       
!                                                                       
      DIMENSION A(NP,NP),INDX(N),B(N) 
      II=0 
      DO 12 I=1,N 
        LL=INDX(I) 
        SUM=B(LL) 
        B(LL)=B(I) 
        IF (II.NE.0)THEN 
          DO 11 J=II,I-1 
            SUM=SUM-A(I,J)*B(J) 
   11     CONTINUE 
        ELSE IF (SUM.NE.0.) THEN 
          II=I 
        ENDIF 
        B(I)=SUM 
   12 END DO 
      DO 14 I=N,1,-1 
        SUM=B(I) 
        IF(I.LT.N)THEN 
          DO 13 J=I+1,N 
            SUM=SUM-A(I,J)*B(J) 
   13     CONTINUE 
        ENDIF 
        B(I)=SUM/A(I,I) 
   14 END DO 
      RETURN 
      END                                           
                                                                        
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    FILE:         EL_EOS.FOR                                           
!                                                                       
!***********************************************************************
!                                                                       
!    MODULE:       EL_EOS                                               
!    TYPE:         SUBROUTINE                                           
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!                                                                       
!    DATE:         2/12/91                                              
!                                                                       
!                  BITNET:  SWESTY@SUNYSBNP or                          
!                  INTERNET: FSWESTY@SBAST1.SUNYSB.EDU or               
!                            fswesty@sbast3.sunysb.edu                  
!                                                                       
!    PURPOSE:      The elctron and photon equation of state             
!                                                                       
!                                                                       
!    CALL LINE:    CALL EL_EOS(T,YE,BRYDNS)                             
!                                                                       
!    INPUTS:       T = TEMPERATURE                                      
!                  YE = ELECTRON FRACTION                               
!                  BRYDNS = BARYON NUMBER DENSITY                       
!                                                                       
!    OUTPUTS:      NONE                                                 
!                                                                       
!                                                                       
!                                                                       
!    INCLUDE FILES:  EL_EOS.INC                                         
!                                                                       
!                                                                       
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
      SUBROUTINE EL_EOS(T,YE,BRYDNS) 
!                                                                       
      IMPLICIT NONE 
!                                                                       
      DOUBLE PRECISION T, YE, BRYDNS 
!                                                                       
      INCLUDE 'el_eos.inc' 
!                                                                       
!                                                                       
!                                                                       
!                           Plancks constant & speed of light           
      DOUBLE PRECISION HBAR, C 
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23) 
!                                                                       
!                           Pi and 1/3                                  
      DOUBLE PRECISION PI, PI2, OVR3, MOVR3, OVR23 
      PARAMETER(PI=3.1415927, PI2=9.8696044) 
      PARAMETER(OVR3=0.33333333, MOVR3=-0.33333333, OVR23=0.66666667) 
!                                                                       
!*****added by MBD                                                      
      double precision me, qcoef 
      me = 5.11d-1 
!                                                                       
!                    Leptons                                            
!                                                                       
!                    Electron number density                            
      NSUBE = BRYDNS*YE 
!                                                                       
!                    Coefficants for chemical potential                 
!                    and thermodynamics quantities                      
      QSUBE = 1.0/( 3.0*(PI**2)*((HBAR*C)**3) ) 
!                                                                       
      ACOEF = 0.5*NSUBE/QSUBE 
!                                                                       
!*****added by MBD new qcoef to include me term                         
      qcoef=((PI*T)**2)/3. - me*me/2. 
      BCOEF= (ACOEF**2+qcoef**3)**0.5 
!      BCOEF = (ACOEF**2+((PI**6)*T**6)/27.0)**0.5                      
!                                                                       
!      DBDT = (PI**6)*(T**5)/(9.0*BCOEF)                                
      DBDT = PI2*T*((PI2*T*T/3. -me*me/2.)**2)/BCOEF 
!                                                                       
      CCOEF = (ACOEF+BCOEF)**OVR3 
!                                                                       
!                                                                       
!                    Electron chemical potential                        
!      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF                            
      MUSUBE = CCOEF-qcoef/CCOEF 
!                                                                       
!*****added by MBD                                                      
!                                                                       
!      write(*,101) qcoef, ACOEF, BCOEF, CCOEF, MUSUBE                  
!                                                                       
!*****now do it again with me=0                                         
!                                                                       
!      qcoef=((PI*T)**2)/3.                                             
!      BCOEF= (ACOEF**2+qcoef**3)**0.5                                  
!      CCOEF = (ACOEF+BCOEF)**OVR3                                      
!      MUSUBE = CCOEF-OVR3*((PI*T)**2)/CCOEF                            
!      write(*,101) qcoef, ACOEF, BCOEF, CCOEF, MUSUBE                  
!                                                                       
!101   format(1x,5(1pd11.4))                                            
!                                                                       
!                                                                       
!                                                                       
!                    Electron pressure for rel. case                    
!*****again me term added                                               
      EPRESS = 0.25*QSUBE*(MUSUBE**4+2.0*(PI*T*MUSUBE)**2+              &
     & 7.0*((PI*T)**4)/15.0 -3.0*(MUSUBE*me)**2 - 0.5*(PI*T*me)**2)     
!                                                                       
!                                                                       
!                    Electron internal energy per baryon                
!*****again me term added                                               
      EU = 0.75*QSUBE*(MUSUBE**4+2.0*(PI*MUSUBE*T)**2+                  &
     & 7.0*((PI*T)**4)/15.0 -(MUSUBE*me)**2 - 0.5*(PI*T*me)**2)         &
     & /BRYDNS                                                          
!                                                                       
!                                                                       
!                    Electron free energy per baryon                    
      FSUBE = ((MUSUBE*NSUBE)-EPRESS)/BRYDNS 
!                                                                       
!                    Electron entropy per baryon                        
      ES = QSUBE*(((PI*MUSUBE)**2)*T+7.0*(PI**4)*(T**3)/                &
     & 15.0 - 0.5*T*me*me)/BRYDNS                                       
!                                                                       
!                    Photons                                            
!                                                                       
!                    Photon pressure                                    
      PPRESS = (PI**2)*(T**4)/(45.0*((HBAR*C)**3)) 
!                    Photon entropy per baryon                          
      PS = 4.0*PPRESS/(T*BRYDNS) 
!                                                                       
!                    Photon internal energy per baryon                  
      PU = 3.0*PPRESS/BRYDNS 
!                                                                       
!                    Photon free energy per baryon                      
      PF = PU-T*PS 
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                    Derivatives of chem. potential w.r.t. T,           
!                    BRYDNS, YE                                         
!                                                                       
!      DEMUDT = DBDT/(3.0*CCOEF**2)-OVR23*(PI**2)*T/CCOEF+              
!     1         DBDT*((PI*T)**2)/(9.0*CCOEF**4)                         
!                                                                       
!      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2)          
!                                                                       
      DEMUDN = (YE*PI2*(HBAR*C)**3)/(MUSUBE**2+OVR3*PI2*T**2            &
     &         - me*me/2.)                                              
!                                                                       
      DEMUDT = -2.*DEMUDN*MUSUBE*T/(3.*YE*(HBAR*C)**3) 
!                                                                       
      DEMUDY = BRYDNS*DEMUDN/YE 
!                                                                       
!                                                                       
!                    Derivatives of pressure w.r.t. BRYDNS,YE,T         
!                                                                       
      DEPDN = BRYDNS*YE*DEMUDN 
!                                                                       
      DEPDY = BRYDNS*DEPDN/YE 
!                                                                       
      DEPDT = BRYDNS*(ES+YE*DEMUDT) 
!                                                                       
!                                                                       
!                    Derivatives of entropy w.r.t. T,BRYDNS,YE          
!                                                                       
      DESDT = ES/T+OVR23*(7.0*PI2*(T**2)/15.0+MUSUBE*T*DEMUDT)/         &
     &        (BRYDNS*(HBAR*C)**3)                                      
!                                                                       
      DESDN = -1.0*DEPDT/(BRYDNS**2) 
!                                                                       
      DESDY = 2.0*T*QSUBE*PI2*MUSUBE*DEMUDY/BRYDNS 
!                                                                       
!                                                                       
!                    Derivatives of internal energy w.r.t.              
!                    T,BRYDNS,YE                                        
      DEUDT = T*DESDT 
!                                                                       
      DEUDN = (YE*(MUSUBE-T*DEMUDT)-EU)/BRYDNS 
!                                                                       
      DEUDY = 3.0*QSUBE*((MUSUBE**3)+PI2*(T**2)*MUSUBE)*                &
     &        DEMUDY/BRYDNS                                             
!                                                                       
!                                                                       
!                               Photons                                 
!                                                                       
!                    Derivatives of photon pressure                     
      DPPDN = 0.0 
      DPPDT = BRYDNS*PS 
      DPPDY = 0.0 
!                                                                       
!                    Derivatives of photon entropy                      
      DPSDN = -PS/BRYDNS 
      DPSDT = 3.0*PS/T 
      DPSDY = 0.0 
!                                                                       
!                    Derivatives of internal energy                     
      DPUDN = -0.75*T*PS/BRYDNS 
      DPUDT = 3.0*PS 
      DPUDY = 0.0 
!                                                                       
!                                                                       
!                                                                       
!23456789012345678901234567890123456789012345678901234567890123456789012
!                                                                       
!                                                                       
  999 RETURN 
!                                                                       
!                                                                       
      END                                           
                                                                        
!      Program to compute spline fits to fermi integrals                
!c  Must provide data file 14                                           
      subroutine initferm 
      implicit real*8(a-h,o-z) 
      parameter (n=201) 
      dimension f32(n),f12(n),fm12(n),eta(n),fr(n) 
      dimension f32a(n),f12a(n),fra(n),fia(n) 
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
!                                                                       
      open(14,file='fermi.atb',status='old') 
!                                                                       
      do 10 i=1,n 
       read(14,*)eta(i),f32(i),f12(i),fm12(i) 
   10  fr(i)=f12(i)/fm12(i) 
!                                                                       
      close(14,status='keep') 
!                                                                       
       call spline(eta,f12,n,f12a) 
!       write(*,1)f12a                                                  
       call spline(eta,f32,n,f32a) 
!       write(*,1)f32a                                                  
       call spline(eta,fr,n,fra) 
!       write(*,1)fra                                                   
       call spline(f12,eta,n,fia) 
!       write(*,1)fia                                                   
    1  format(1p8e10.3) 
       return 
      END                                           
      SUBROUTINE SPLINE(X,Y,N,Y2) 
!  Computes spline coefficients; Y(X) is input data; Y2 is output.      
      implicit real*8(a-h,o-z) 
        DIMENSION X(N),Y(N),Y2(N),U(500) 
      Y2(1)=0. 
      U(1)=0. 
      DO 11 I=2,N-1 
      SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1)) 
      P=SIG*Y2(I-1)+2. 
      Y2(I)=(SIG-1.)/P 
   11    U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))            &
     & /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P                    
      Y2(N)=0. 
      DO 12 K=N-1,1,-1 
   12    Y2(K)=Y2(K)*Y2(K+1)+U(K) 
      RETURN 
      END                                           
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,KLO,KHI) 
!     Computes spline fit of Y(X); YA(XA) is input data, Y2A are spline 
!  coefficents, klo and khi are running indices which bracket X.        
      implicit real*8(a-h,o-z) 
      DIMENSION XA(N),YA(N),Y2A(N) 
!c  Determine the bracketing indices                                    
      IF(KHI-KLO.GT.1)GOTO 1 
      IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
      KHI=KHI-1 
      KLO=KLO-1 
      IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
      KHI=KHI+2 
      KLO=KLO+2 
      IF(XA(KHI).GT.X.AND.XA(KLO).LE.X) GOTO 2 
      KLO=1 
      KHI=N 
    1 IF(KHI-KLO.EQ.1) GOTO 2 
      K=(KHI+KLO)/2 
      IF(XA(K).GT.X)THEN 
      KHI=K 
      ELSE 
      KLO=K 
      ENDIF 
      GOTO 1 
    2 H=XA(KHI)-XA(KLO) 
      IF(H.EQ.0.) PAUSE 'BAD XA INPUT. ' 
!c  Compute spline fit.                                                 
      A=(XA(KHI)-X)/H 
      B=(X-XA(KLO))/H 
      Y=A*YA(KLO)+B*YA(KHI)+                                            &
     & ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*H**2/6.                    
!      write(*,5)klo,khi,x,xa(klo),xa(khi),ya(klo),ya(khi)              
!     > ,y2a(klo),y2a(khi),y                                            
    5   format(2i3,1p8e9.2) 
      RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       F_1_2                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    F_1_2(Y)      (1/2th Fermi Integral)                 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       1/2th Fermi Integral (DOUBLE PRECISION)              
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION F_1_2(Y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n), &
     & fia(n)                                                           
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.17677669&
     &,6.41500299D-02,.4D0,1.32934039D0,.33333333333d0,1,n/             
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,f12,f12a,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=A(6)*y*SQRT(y)*th*(5+(A(1)+(3*A(2)+7*X2*A(3))*X2)*x2) 
      GO TO 100 
   20 F0=DEXP(y) 
      F1=A(7)*th*F0*(2-(4*A(4)-(6*A(5)-.25*F0)*F0)*F0) 
  100 F_1_2=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       F_3_2                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    F_3_2(Y)      (3/2th Fermi Integral)                 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       3/2th Fermi Integral (DOUBLE PRECISION)              
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION F_3_2(y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION A(7),eta(n),f32(n),f32a(n),f12(n),fr(n),f12a(n),fra(n), &
     & fia(n)                                                           
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.176776695D0&
     &,6.41500299D-02,.4D0,1.32934039D0,1,n/                            
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,f32,f32a,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=A(6)*SQRT(y)*(1./X2+A(1)-(A(2)+X2*A(3))*X2) 
      GO TO 100 
   20 F0=DEXP(y) 
      F1=A(7)*F0*(1.-(A(4)-(A(5)-.03125*F0)*F0)*F0) 
  100 F_3_2=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       FINV12                                               
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         11/29/89                                             
!                                                                       
!    CALL LINE:    FINV12(Y)      (Inverse of the 1/2th Fermi Integral) 
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       Inverse of Fermi Integral (DOUBLE PRECISION)         
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION FINV12(y) 
      IMPLICIT REAL*8(A-H,O-Z) 
      parameter(n=201) 
      DIMENSION AI(8),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),        &
     & fra(n),f12a(n)                                                   
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA AI,klo,khi/-.822467032D0,-1.21761363D0,-9.16138616D0,        &
     &.398942281D0,.0732748216D0,-1.310707D0,1.12837917D0,              &
     &8.2810645D-3,1,n/                                                 
      if(y .gt. 109.695) goto 10 
      if(y .lt. 4.0234e-5) goto 20 
      call splint(f12,eta,fia,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=(1.5*y)**(.666666667) 
      X4=1./(X2*X2) 
      F1=X2*(1.+(AI(1)+(AI(2)+AI(3)*X4)*X4)*X4) 
      GO TO 100 
   20 F1=LOG(AI(7)*MAX(y,1.D-20)*(1.+(AI(4)+(AI(5)+AI(8)*y)*y)*y)) 
  100 finv12=F1 
  999 RETURN 
      END                                           
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       FHALF                                                
!    TYPE:         DOUBLE PRECISION FUNCTION                            
!    AUTHOR:       F. DOUGLAS SWESTY, Dpt of Physics SUNY @ Stony Brook 
!    DATE:         3/10/90                                              
!                                                                       
!    CALL LINE:    FHALF(Y)                                             
!                                                                       
!    INPUTS:       Y (DOUBLE PRECISION)   (Argument)                    
!                                                                       
!    RETURN:       Ratio of 1/2th Fermi Integral to the -1/2th Fermi    
!                  Integral (DOUBLE PRECISION)                          
!                                                                       
!***********************************************************************
      DOUBLE PRECISION FUNCTION FHALF(y) 
      IMPLICIT real*8(a-h,o-z) 
      parameter(n=201) 
      DIMENSION A(7),f12(n),fia(n),eta(n),f32(n),fr(n),f32a(n),         &
     & fra(n),f12a(n)                                                   
      common /spl/eta,f32,f12,fr,f32a,f12a,fra,fia 
      DATA A,th,klo,khi/6.16850274D0,1.77568655D0,6.92965606D0,.17677669&
     &,6.41500299D-02,.4D0,1.32934039D0,.3333333333333d0,1,n/           
      IF(y .gt. 30.) goto 10 
      if(y .lt. -10.) goto 20 
      call splint(eta,fr,fra,n,y,f1,klo,khi) 
      GO TO 100 
   10 X2=y**(-2) 
      F1=y*th*(1.+(.2*A(1)+(.6*A(2)+1.4*X2*A(3))*X2)*x2)                &
     & /(1.-(.2*th*a(1)-(a(2)-4.2*x2*a(3))*x2)*x2)                      
      GO TO 100 
   20 F0=EXP(y) 
      F1=(1.-(2*a(4)-(3*a(5)-.125*f0)*f0)*f0)/                          &
     & (2.-(8*a(4)-(18*a(5)-f0)*f0)*f0)                                 
  100 FHALF=F1 
  999 RETURN 
      END                   