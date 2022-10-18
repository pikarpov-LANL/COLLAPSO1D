!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       MAXWEL.INC                                           
!    TYPE:         INCLUDE FILE                                         
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         6/13/90                                              
!                                                                       
!                                                                       
!    CALL LINE:    INCLUDE 'MAXWEL.INC/LIST'                            
!                                                                       
!***********************************************************************
!                                                                       
!                            Number of T points in each boundary        
      INTEGER NUMTMP 
      PARAMETER (NUMTMP=201) 
!c      PARAMETER (NUMTMP=21)                                           
!                                                                       
!                            Number of Ye points in each boundary       
      INTEGER NUMYE 
!c      PARAMETER(NUMYE=49)                                             
!c      PARAMETER(NUMYE=34)                                             
      PARAMETER(NUMYE=49) 
!                                                                       
!                            Number of points in alpha-nuclei           
!                            boundary                                   
      INTEGER NBPNTS, NUMLOW, NUMHI 
      PARAMETER(NBPNTS=101,NUMLOW=51,NUMHI=51) 
!                                                                       
!                                                                       
!                            Baryon densities                           
      DOUBLE PRECISION BRYLOW(NUMTMP,NUMYE), BRYHI(NUMTMP,NUMYE) 
      DOUBLE PRECISION LOWDNS, HIDNS, DNS_1, DNS_2 
      DOUBLE PRECISION DNL_DT, DNH_DT, DNL_DY, DNH_DY 
!                                                                       
!                                                                       
      DOUBLE PRECISION LNLOW, LNHI, LNCUT, LNMINS, LNPLUS 
      DOUBLE PRECISION LOGBRY, LOGBCH, DLTLN1, DLTLN2, LNFRAC 
      PARAMETER(LNLOW=-8.92d0,LNHI=-0.92d0,LNCUT=-2.92d0) 
!                                                                       
!                                                                       
      DOUBLE PRECISION YLOW,YHI,Y_CUT 
!ccc 180 MeV parms      PARAMETER(YLOW=0.03,YHI=0.46,Y_CUT=0.175)       
!ccc 220 MeV parms      PARAMETER(YLOW=0.03,YHI=0.46,Y_CUT=0.135)       
!ccc 375 MeV parms      PARAMETER(YLOW=0.03,YHI=0.46,Y_CUT=0.155)       
      PARAMETER(YLOW=0.03d0,YHI=0.51d0,Y_CUT=0.175d0) 
!                                                                       
!                                                                       
      DOUBLE PRECISION TCHK_B, TCHK_N, T_MXWL, D_MXWL 
!                                                                       
!                                                                       
!                            Highest temperatures where Coexistence     
!                            of bulk and nuclear phases occurs          
      DOUBLE PRECISION T_H(NUMYE), D_H(NUMYE) 
!                                                                       
!                                                                       
!                            Arrays containing number of boundary       
!                            points                                     
      DOUBLE PRECISION LBOUND(NBPNTS,NUMYE), UBOUND(NBPNTS,NUMYE) 
!                                                                       
!                            Minimum density code will work at          
!                                                                       
      DOUBLE PRECISION MINDNS 
      PARAMETER(MINDNS=1.0D-10) 
!                                                                       
!                                                                       
!                                                                       
      DOUBLE PRECISION NEWDNS 
      DOUBLE PRECISION T_LOW, T_HI, TFRAC 
      DOUBLE PRECISION Y_LOW, Y_LOW2, Y_HI, Y_HI2, YFRAC, YCH, YCUT 
      DOUBLE PRECISION LNL, LNH, LNC 
!                                                                       
      DOUBLE PRECISION YMINUS, YPLUS, YINTRP, TMINUS, TPLUS, TINTRP 
      DOUBLE PRECISION DELT_Y, DELT_T 
      INTEGER I_MXWL, J_MXWL 
      INTEGER I_BD, I_BNDY, J_BD, J_BNDY 
!                                                                       
!                            Total chemical potential                   
      DOUBLE PRECISION MUTILD, MUTLOW, MUTHI 
!                                                                       
!                            Total pressure                             
      DOUBLE PRECISION PRTILD, PRLOW, PRHI 
!                                                                       
!                                                                       
!                            Muhat                                      
      DOUBLE PRECISION MUHLOW, MUHHI 
!                                                                       
!                            Electron chemical potential                
      DOUBLE PRECISION MUELOW, MUEHI 
!                                                                       
!                            Total entropy per baryon                   
      DOUBLE PRECISION S_LOW, S_HI 
!                                                                       
!                            Total free energy density                  
      DOUBLE PRECISION F_LOW, F_HI 
!                                                                       
!                                                                       
!                            Phase fraction                             
      DOUBLE PRECISION PHASEF 
!                                                                       
!                                                                       
      COMMON /MAXCOM/ T_LOW,T_HI,Y_LOW,Y_HI,           &
     &       BRYLOW,BRYHI,LBOUND,UBOUND,T_H,D_H,       &
     &       LOWDNS,HIDNS,TCHK_N,TCHK_B,T_MXWL,D_MXWL, &
     &       DNL_DT,DNH_DT,DNL_DY,DNH_DY