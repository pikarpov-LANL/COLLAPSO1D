!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!                                                                       
!    MODULE:       EOS_M4A.INC                                          
!    TYPE:         INCLUDE FILE                                         
!    AUTHOR:       F. DOUGLAS SWESTY                                    
!    DATE:         6/27/90                                              
!                                                                       
!                                                                       
!    CALL LINE:    INCLUDE 'EOS_M4A.INC/LIST'                           
!                                                                       
!***********************************************************************
!                                                                       
!                   INPUT & OUTPUT VARIBLES                             
!                                                                       
!                            Electron fraction                          
      DOUBLE PRECISION YE 
!                            Baryon density                             
      DOUBLE PRECISION BRYDNS 
!                            Previous value of the nuclear proton       
!                            fraction (must be supplied)                
      DOUBLE PRECISION XPREV 
!                                                                       
!                            Previous value of the proton number        
!                            density (must be supplied)                 
      DOUBLE PRECISION P_PREV 
!                                                                       
!                            Temperature, density, neutron chem. pot.,  
!                            and proton chem. pot.                      
      DOUBLE PRECISION INPVAR(4) 
!                                                                       
!                            Flag to tell code what input is being      
!                            supplied (i.e. what INPVAR is)             
!                            1 --> Temperature                          
!                            2 --> Internal energy (NOT IMPLEMENTED)    
!                            3 --> Entropy (NOT IMPLEMENTED)            
      INTEGER IFLAG 
!                                                                       
!                            Flag returned by code to tell user what    
!                            scheme was used                            
!                            1 --> "No Nuclei" scheme                   
!                            2 --> full scheme                          
!                            3 --> bulk scheme (above saturation dnsty) 
      INTEGER EOSFLG 
!                                                                       
      INTEGER FORFLG, NGFLAG, NEWFLG, SSFLAG, FFLAG, ADFLAG 
      INTEGER NF_FLG, RSFLAG, DBFLAG 
!c      DATA NF_FLG/1/                                                  
!c      DATA RSFLAG/0/                                                  
!                                                                       
!***********************************************************************
!                        Local variables                                
!***********************************************************************
!                                                                       
!                            Flag used internally to this routine to    
!                            record whether or not another scheme       
!                            has already been tried                     
      INTEGER SWTFLG 
!                                                                       
!                                                                       
!                            Accuracy in X                              
      DOUBLE PRECISION XACC 
      PARAMETER (XACC=1.0D-8) 
!                                                                       
!                            Maximum number of X iterations             
      INTEGER MAXIT, MAXIT2 
      PARAMETER (MAXIT=200,MAXIT2=30) 
!                                                                       
!                            Update flag                                
      INTEGER UPFLAG 
      PARAMETER(UPFLAG=3) 
!                                                                       
!                            External Interaction flag                  
!                            1 = on, 0 = off                            
      DOUBLE PRECISION EIFLAG 
      PARAMETER(EIFLAG=1.0) 
!                                                                       
!                            Alpha particle enable flag                 
!                            1= alphas, 0 = no alphas                   
      INTEGER ALFLAG 
      PARAMETER(ALFLAG=1) 
!                                                                       
!                            Translational energy scaling factor        
!                            1 = on, 0 = off                            
      DOUBLE PRECISION TRSCAL 
      PARAMETER(TRSCAL=1.0) 
!                                                                       
!                            Surface and Coulomb scaling factor         
!                            1.0 = on, 0 = off                          
      DOUBLE PRECISION CSSCAL 
      PARAMETER(CSSCAL=1.0) 
!                                                                       
!                                                                       
!                                                                       
!                            Cutoff fraction of nuclei (If XH is below  
!                            this the "no nuclei" scheme is used)       
      DOUBLE PRECISION HEAVCT 
      PARAMETER(HEAVCT=1.0D-6) 
!                                                                       
!                            Iteration loop variable                    
      INTEGER I, J, L, IT_NUM 
!                                                                       
!                            Combination of Coulomb and Surface         
!                            coefficients                               
      DOUBLE PRECISION ALPHA, BETA 
      PARAMETER (ALPHA=224.4, BETA=1.3927) 
!                                                                       
!                            Surface tension (MeV/fm**2)                
      DOUBLE PRECISION SIG_0 
      PARAMETER(SIG_0=46.5) 
!                            Binding energy of alpha particles (MeV)    
      DOUBLE PRECISION BALPHA 
      PARAMETER (BALPHA=28.3) 
      DOUBLE PRECISION V_ALFA 
      PARAMETER (V_ALFA=24.0) 
!                                                                       
!                            Nuclear symmetry energy (MeV)              
      DOUBLE PRECISION SV 
      PARAMETER (SV=29.3) 
!                                                                       
!                            Nuclear compression modulus (MeV)          
      DOUBLE PRECISION COMPRS, K 
      PARAMETER (COMPRS=370.0) 
!                                                                       
!                            Nuclear level density (per MeV)            
      DOUBLE PRECISION LEVDNS, AV 
      PARAMETER (LEVDNS=0.0666667) 
!                                                                       
!                            Approximate size of nuclei                 
      DOUBLE PRECISION AZERO 
      PARAMETER (AZERO=60.0) 
!                                                                       
!                                                                       
!                           Plancks constant & speed of light           
      DOUBLE PRECISION HBAR, C 
      PARAMETER (HBAR=6.58217317D-22,C=2.997924581D23) 
!                                                                       
!                           Pi and square root of Pi                    
      DOUBLE PRECISION PI, ROOTPI 
      PARAMETER(PI=3.1415927,ROOTPI=1.7724539) 
!                                                                       
!                           Square root of two                          
      DOUBLE PRECISION ROOT2 
      PARAMETER (ROOT2=1.4142136) 
!                                                                       
!                           Parameters containing powers of 1/3         
!                                                                       
      DOUBLE PRECISION OVR3, OVR23, OVR43, OVR53, OVR29 
      DOUBLE PRECISION MOVR3, M2OVR3, M4OVR3, M5OVR3 
      DOUBLE PRECISION OVR49 
      PARAMETER(OVR49=0.44444444) 
      PARAMETER (OVR3=0.33333333, OVR23=0.66666666, OVR43=1.33333333) 
      PARAMETER (OVR53=1.66666667, OVR29=0.2222222222) 
      PARAMETER (MOVR3=-0.33333333, M2OVR3=-0.66666666) 
      PARAMETER (M4OVR3=-1.33333333, M5OVR3=-1.666666667) 
!                                                                       
!                          Ratio of baryon density to saturation density
      DOUBLE PRECISION Y 
!                                                                       
      DOUBLE PRECISION TFINAL 
!                                                                       
!                          Quantum concentration & Fermi integral       
!                          coefficent                                   
      DOUBLE PRECISION NQ, MQ, LQ, KQ, UQ 
!                                                                       
!                          Estimate of X (used in "no nuclei" scheme)   
      DOUBLE PRECISION XGUESS 
!                                                                       
!                          Cutoff used to determine how to calc the     
!                          number density of heavy nuclei (used in      
!                          the "no nuclei" scheme)                      
      DOUBLE PRECISION MHCUT, MHOVT 
      PARAMETER(MHCUT=30.0) 
!                                                                       
!                          Chem. pot., mass, # density, & quantum       
!                          concentration of heavy nuclei (used in       
!                          the "no nuclei" scheme)                      
      DOUBLE PRECISION MUHEAV, NUCQ, MASSHV, NHEAVY 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                          Bulk nuclear pressures                       
      DOUBLE PRECISION PSUBI, POVRNI 
!                                                                       
!                          Neutron & proton rest masses                 
      DOUBLE PRECISION MASSN, MASSP 
      PARAMETER(MASSN=939.5731,MASSP=938.2796) 
!                                                                       
!                          Saturation density, nuclear density, the     
!                          difference of the two and their ratio        
      DOUBLE PRECISION NSUBS, NSUBI, N_I, DNI, NRATIO, N_IOLD 
!CC      PARAMETER (NSUBS=0.155)                                        
      DOUBLE PRECISION NIOLD 
!                                                                       
!                                                                       
!                          Skyrme interaction coefficants               
      DOUBLE PRECISION AA, BB, CC, DD 
      DOUBLE PRECISION COMP, SYMM, BIND_E, N_S, SKYRMC 
      DOUBLE PRECISION SYM_S, SIG_S 
      DOUBLE PRECISION Q, R_0, ZETA_0 
!c      PARAMETER(AA=-281.0,BB=-111.6,CC=966.0,DD=2.0)                  
!                                                                       
!                                                                       
!                                                                       
!                          Coulomb liquid correction functions          
!                          & derivatives                                
      DOUBLE PRECISION W, DWDX, DWDT, DWDTDX, DWDXDX, DWDTDT 
      DOUBLE PRECISION TZERO, DTZDX, DTZDXX 
!                                                                       
!                                                                       
!                          Coulomb functions and derivatives            
      DOUBLE PRECISION DU, DMU, DUP, DMUP, DUPP, DMUPP 
      DOUBLE PRECISION DUX, DMUX, DUPX, DMUPX 
      DOUBLE PRECISION DUT, DMUT, DUPT, DMUPT, DUPPT 
      DOUBLE PRECISION DUXT, DMUXT 
      DOUBLE PRECISION DUXX, DMUXX, DUTT, DMUTT 
!                                                                       
      DOUBLE PRECISION SCRDU, SCRDUP, SCRDPP, SCRD, SCRD2 
      DOUBLE PRECISION SCRDUT, SCRDUX, SCRDXX, SCRDT, SCRDX, SCRDXT 
      DOUBLE PRECISION SCRDTT 
      DOUBLE PRECISION SCRD2T, SCRD2X, SCRDPT, SCRDPX 
      DOUBLE PRECISION UTILDE, DUTIL, DUTILP, OMEG_R 
!                                                                       
!                                                                       
!                          Finite size energy coefficent                
      DOUBLE PRECISION ZETA, DZDT, DZDX, DZDNI 
!                                                                       
!                          Equilibrium equation finite size terms       
      DOUBLE PRECISION BN, BP, BSUBP 
!                                                                       
!                          Baryon densities                             
      DOUBLE PRECISION NNOUT, ALFDNS, NUCDNS, NOUT 
!                                                                       
!                          Baryon fractions                             
      DOUBLE PRECISION XALFA, XNUT, XPROT, XH, XHCHK 
      DOUBLE PRECISION XALFA2, XNUT2, XPROT2 
!                                                                       
!                          Fraction of space occupied by nuclei         
      DOUBLE PRECISION U_NUC, U_N, COMPU, RY 
!                          Fraction of space available to remaining     
!                          baryons  (1-U_NUC)                           
      DOUBLE PRECISION EXALFA, EXCLU 
!                                                                       
!                                                                       
!                          Surface tension variables                    
      DOUBLE PRECISION COMPX, SIGMA, AHBN, OVRX4, DSIGDX, DSGRDX, SIGSGP 
      DOUBLE PRECISION SIGSG2 
!                                                                       
!                                                                       
!                          Critical temperature variables               
      DOUBLE PRECISION TSUBC, TSC1, DTCDX, DTCDXX, TSC_12 
!                                                                       
!                          Surface & translational temperature variables
      DOUBLE PRECISION DHDX, DHDXX, HX, DHDTDX, CHOVA 
      DOUBLE PRECISION H, HPRIM, HPPRIM 
      DOUBLE PRECISION CAPH, CAPHP 
!                                                                       
!                          Neutron-Proton mass difference               
      DOUBLE PRECISION DELTAM 
      PARAMETER(DELTAM=-1.2935) 
!cccccc      PARAMETER(DELTAM=0.0)                                      
!                                                                       
!                          Chemical potentials                          
      DOUBLE PRECISION MUHAT, MUN, MUNOVT, MUALFA, MUPROT, MUN_I 
      DOUBLE PRECISION MUP_I, MUN_O, MUP_O, MUN_Z, MUP_Z 
      DOUBLE PRECISION ETA_NO, ETA_PO, ETA_NI, ETA_PI, ETAMAX 
      DOUBLE PRECISION ETP_HI, ETP_LO, ETN_HI, ETN_LO 
      DOUBLE PRECISION DELTMN, DELTMP, DELTAP 
!                                                                       
      DOUBLE PRECISION DMNDPI, DMPDPI, DPIDNI 
!                                                                       
      DOUBLE PRECISION OMEGA 
!                                                                       
      DOUBLE PRECISION DETA_P, DETA_N, NETA_N, NETA_P 
      DOUBLE PRECISION CHNG_N, CHNG_P 
!                                                                       
!                                                                       
!                                                                       
!                          N-N interaction variables                    
      DOUBLE PRECISION VNOUT, VPOUT, VNI, VPI 
      DOUBLE PRECISION ZNO, ZPO, ZNI, ZPI 
!                                                                       
!                          Nucleon kinetic energy densities             
      DOUBLE PRECISION TAU_PO, TAU_NO, TAU_NI, TAU_PI 
      DOUBLE PRECISION F32_PO, F32_NO, F32_NI, F32_PI 
      DOUBLE PRECISION FRAT_N, FRAT_P 
!                                                                       
!                          Newton-Raphson equation variables            
      DOUBLE PRECISION G, DGDX, DGDPRT, GOLD 
!                                                                       
!                          Atomic weight & number                       
      DOUBLE PRECISION A, Z 
!                                                                       
!                          Nuclear radius & volume                      
      DOUBLE PRECISION RSUBN, VSUBN 
!                                                                       
!                          Surface, Coulomb, translational, and         
!                          bulk free energies (per baryon)              
      DOUBLE PRECISION FSUBS, FSUBC, FSUBSC, FSUBI, FTRANS, MUSUBT 
      DOUBLE PRECISION DMUTDT 
!                                                                       
      DOUBLE PRECISION FTR, DFTRDT, DFTRDU, DFTRDX, DFTRDN 
      DOUBLE PRECISION F_SC, DFSCDT, DFSCDU, DFSCDX, DFSCDN 
!                                                                       
      DOUBLE PRECISION E_TR, DETRDT, DETRDU, DETRDX, DETRDN 
      DOUBLE PRECISION E_SC, DESCDT, DESCDU, DESCDX, DESCDN 
!                                                                       
      DOUBLE PRECISION S_TR, DSTRDT, DSTRDU, DSTRDX, DSTRDN 
      DOUBLE PRECISION S_SC, DSSCDT, DSSCDU, DSSCDX, DSSCDN 
!                                                                       
      DOUBLE PRECISION S_NUC, S_OUT, S_ALFA 
      DOUBLE PRECISION E_NUC, E_OUT, E_ALFA 
!                                                                       
      DOUBLE PRECISION S_DENS, E_DENS 
!                                                                       
      DOUBLE PRECISION DETPDX, DETNDX, DETPDN, DETNDN 
      DOUBLE PRECISION DTPIDX, DTNIDX, DTPDNI, DTNDNI 
!                                                                       
      DOUBLE PRECISION DSIDT, DSIDX, DSIDN, DEIDT, DEIDX, DEIDN 
!                                                                       
      DOUBLE PRECISION DSODT, DSODEP, DSODEN, DEODT, DEODEP, DEODEN 
!                                                                       
      DOUBLE PRECISION DSADT, DSADEP, DSADEN, DEADT, DEADEP, DEADEN 
!                                                                       
      DOUBLE PRECISION DNPDEP, DNPDEN, DNNDEP, DNNDEN 
      DOUBLE PRECISION DTPDEP, DTPDEN, DTNDEP, DTNDEN 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                            BARYON THERMODYNAMIC VARIABLES             
!                                                                       
!                          Pressures                                    
      DOUBLE PRECISION BPRESS, BPROUT, BPRNUC, BPRALF 
!                                                                       
!                          Entropies (per baryon)                       
      DOUBLE PRECISION BSOUT, BSNUC, BS, BSALFA 
      DOUBLE PRECISION DSFSDU, DSFSDX, DSFSDN, DSFSDT 
!                                                                       
!                                                                       
!                                                                       
      DOUBLE PRECISION DBSODT, DBSNDT, DBSADT, SFS_DT 
!                                                                       
!                          Internal energies (per baryon)               
      DOUBLE PRECISION BUOUT, BUNUC, BU, BUALFA 
!                                                                       
!                          Helmholtz free energies (per baryon)         
      DOUBLE PRECISION BFOUT, BFALFA, BFNUC, BFTOT 
!                                                                       
!                          Adiabatic index                              
      DOUBLE PRECISION GAM_S 
!                                                                       
!                            TOTAL ENERGY, ENTROPY, AND PRESSURE        
!                                                                       
      DOUBLE PRECISION UTOT, STOT, PTOT, FTOT 
      DOUBLE PRECISION FBARY, PBARY, MUBARY, MU_MAT 
!                                                                       
!                            X VARIABLES                                
!                                                                       
!                          Internal proton fraction (in nuclei)         
      DOUBLE PRECISION X 
!                          Outside proton fraction (nucleon vapor)      
      DOUBLE PRECISION XOUT 
!                                                                       
      DOUBLE PRECISION XNEW, XOLD, DX 
!                          Limits on X                                  
      DOUBLE PRECISION X_MIN 
!C      DOUBLE PRECISION XL, XHI                                        
!C      PARAMETER(XL=0.01,XHI=0.52)                                     
!                                                                       
!                                                                       
!                            PROTON VARIABLES                           
      DOUBLE PRECISION NPOUT, PRTNEW, PRTOLD, DPRT 
!                                                                       
!                                                                       
!                            TEMPERATURE VARIABLES                      
!                                                                       
      DOUBLE PRECISION T 
!                                                                       
!                            Diagonal line in the temp-density plane    
      DOUBLE PRECISION TCHECK 
!                                                                       
!                            Temporary storage variables                
      DOUBLE PRECISION TMP1, TMP1P, TMP1PP 
      DOUBLE PRECISION TMP2, TMP2P, TMP2PP, TMP2PT, TMP2PX 
      DOUBLE PRECISION TMP2T, TMP2X, TMP2XX, TMP2XT, TMP2TT 
      DOUBLE PRECISION TMP3, TMP3P, TMP3PP, TMP3PT, TMP3PX 
      DOUBLE PRECISION TMP3T, TMP3X, TMP3XX, TMP3XT, TMP3TT 
      DOUBLE PRECISION TMP4, TMP5 
!                                                                       
!                                                                       
!                                                                       
      DOUBLE PRECISION NSUBIN, NETAP, NETAN 
!                                                                       
      DOUBLE PRECISION NSIACC, PRTACC, NUTACC 
      PARAMETER(NSIACC=1.0D-8,PRTACC=1.0D-9,NUTACC=1.0D-9) 
!ccc      PARAMETER(NSIACC=1.0D-9,PRTACC=1.0D-10,NUTACC=1.0D-10)        
!                                                                       
      DOUBLE PRECISION DETERM 
!                                                                       
      DOUBLE PRECISION DNSUBI, DETAP, DETAN 
!                                                                       
      DOUBLE PRECISION A1, A2, A3 
!                                                                       
      DOUBLE PRECISION A1H, A1L, DA1DN, DA1ETP, DA1ETN 
      DOUBLE PRECISION A2H, A2L, DA2DN, DA2ETP, DA2ETN 
      DOUBLE PRECISION A3H, A3L, DA3DN, DA3ETP, DA3ETN 
!                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!                                                                       
      DOUBLE PRECISION GPI,GNI,GPO,GNO 
      DOUBLE PRECISION DVPIDP, DVPIDN, DVNIDP, DVNIDN 
      DOUBLE PRECISION DVPODP, DVPODN, DVNODP, DVNODN 
      DOUBLE PRECISION MSSCON 
      DOUBLE PRECISION DTPIDP, DTPIDN, DTNIDP, DTNIDN 
      DOUBLE PRECISION DTPODP, DTPODN, DTNODP, DTNODN 
      DOUBLE PRECISION DMPIDP, DMPIDN, DMNIDP, DMNIDN 
      DOUBLE PRECISION DMPODP, DMPODN, DMNODP, DMNODN 
      DOUBLE PRECISION DPIDP, DPIDN, DPODP, DPODN, DPADP, DPADN 
      DOUBLE PRECISION N1, N2, DUDPO, DUDNO, DUDNI 
      DOUBLE PRECISION DXDPO, DXDNO, DXDNI 
      DOUBLE PRECISION DB1DNI, DB1DX, DB1DU 
      DOUBLE PRECISION DB2DNI, DB2DX, DB2DU 
      DOUBLE PRECISION DB3DNI, DB3DX, DB3DU 
      DOUBLE PRECISION DA1ID1, DA1ID2, DA1ID3, DA1OD1, DA1OD2, DA1OD3 
      DOUBLE PRECISION DA2ID1, DA2ID2, DA2ID3, DA2OD1, DA2OD2, DA2OD3 
      DOUBLE PRECISION DA3ID1, DA3ID2, DA3ID3, DA3OD1, DA3OD2, DA3OD3 
      DOUBLE PRECISION DA1D1, DA1D2, DA1D3, DB1D1, DB1D2, DB1D3 
      DOUBLE PRECISION DA2D1, DA2D2, DA2D3, DB2D1, DB2D2, DB2D3 
      DOUBLE PRECISION DA3D1, DA3D2, DA3D3, DB3D1, DB3D2, DB3D3 
      DOUBLE PRECISION DNDETN, DPDETP 
!                                                                       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
!                                                                       
!                                                                       
!                                                                       
!/////////////////////////////////////////////////////////////////////  
!                   Derivative variables                                
!/////////////////////////////////////////////////////////////////////  
!                                                                       
!                                                                       
!                   Temperature derivatives                             
!                                                                       
      DOUBLE PRECISION DMPODT, DMNODT, DMPIDT, DMNIDT 
      DOUBLE PRECISION DNADPO, DNADNO 
!                                                                       
      DOUBLE PRECISION DTPIDT, DTNIDT, DTPODT, DTNODT 
      DOUBLE PRECISION DPODT, DPIDT 
!                                                                       
      DOUBLE PRECISION DNPODT, DNNODT 
!                                                                       
      DOUBLE PRECISION DV_DT 
!                                                                       
      DOUBLE PRECISION DNADT, DPADT, DMUADT 
!                                                                       
!                                                                       
      DOUBLE PRECISION DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT 
!                                                                       
      DOUBLE PRECISION DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN 
!                                                                       
      DOUBLE PRECISION DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY 
!                                                                       
      DOUBLE PRECISION DNA_DT, DNA_DY, DNA_DN 
!                                                                       
      DOUBLE PRECISION DB1DT, DB2DT, DB3DT 
!                                                                       
!                                                                       
!                                                                       
!                   Exterior particle density derivatives               
!                                                                       
!                                                                       
      DOUBLE PRECISION DV_DPO, DV_DNO 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                   ETA derivatives                                     
!                                                                       
!                                                                       
      DOUBLE PRECISION DMPDEP, DMPDEN, DMNDEP, DMNDEN 
!                                                                       
      DOUBLE PRECISION DV_DEP, DV_DEN 
!                                                                       
      DOUBLE PRECISION DNADEP, DNADEN, DPADEP, DPADEN 
!                                                                       
      DOUBLE PRECISION DPODEP, DPODEN 
!                                                                       
      DOUBLE PRECISION DNPODN, DNPIDN, DNNIDN 
!                                                                       
      DOUBLE PRECISION DMADEP, DMADEN 
!                                                                       
!                                                                       
!                                                                       
      DOUBLE PRECISION DFDOM(5,5), DFDOMI(5,5) 
      DOUBLE PRECISION DFDL_1(5), DFDL_2(5), DFDL_3(5) 
!                                                                       
      DOUBLE PRECISION DG1DO1, DG1DO2, DG2DO1, DG2DO2, DET_GT 
      DOUBLE PRECISION DG1DL1, DG1DL2, DG1DL3 
      DOUBLE PRECISION DG2DL1, DG2DL2, DG2DL3 
!                                                                       
      DOUBLE PRECISION A_1X, A_1U, A_1N 
      DOUBLE PRECISION A_2X, A_2U, A_2N 
      DOUBLE PRECISION A_3X, A_3U, A_3N 
      DOUBLE PRECISION DET_AL 
!                                                                       
      DOUBLE PRECISION AT_11, AT_12, AT_13 
      DOUBLE PRECISION AT_21, AT_22, AT_23 
      DOUBLE PRECISION AT_31, AT_32, AT_33 
!                                                                       
      DOUBLE PRECISION B_1P, B_2P, B_3P 
      DOUBLE PRECISION B_1N, B_2N, B_3N 
!                                                                       
!                                                                       
      DOUBLE PRECISION DNIDPO, DNIDNO 
!                                                                       
!                                                                       
      DOUBLE PRECISION C_P, C_N, CP_P, CP_N 
      DOUBLE PRECISION BB_P, BB_N, BP_P, BP_N 
!                                                                       
      DOUBLE PRECISION DMP_DN, DMN_DN 
!                                                                       
      DOUBLE PRECISION DMP_DY, DMN_DY 
!                                                                       
      DOUBLE PRECISION DFDNDN, DFDYDY, DFDNDY, DFDNDT, DFDYDT 
!                                                                       
      DOUBLE PRECISION DFDTDT 
!                                                                       
      DOUBLE PRECISION DBPDT, DBPDN, DBPDY 
!                                                                       
      DOUBLE PRECISION DBSDT, DBSDN, DBSDY 
!                                                                       
      DOUBLE PRECISION DBUDT, DBUDN, DBUDY 
!                                                                       
      DOUBLE PRECISION DBMUDT, DBMUDN, DBMUDY 
!                                                                       
      DOUBLE PRECISION DMUDT, DMUDN, DMUDY, DSDT, DSDN, DSDY 
      DOUBLE PRECISION DPDT, DPDN, DPDY, DUDT, DUDN, DUDY 
!                                                                       
!                                                                       
      DOUBLE PRECISION NPI, NNI, DNPIDT, DNNIDT 
!                                                                       
      DOUBLE PRECISION DXHDT, DXODT, DXADT 
!                                                                       
      DOUBLE PRECISION VP_IN, VN_IN 
!                                                                       
!                                                                       
      DOUBLE PRECISION DUODT, DUIDT, DUADT 
!                                                                       
      DOUBLE PRECISION DBFDT, DBFDY, DBFDN 
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!/////////////////////////////////////////////////////////////////////  
!                                                                       
!                                                                       
!                                                                       
!                                                                       
!                   This common block contains all of the thermodynamic 
!                   and compositional variables that the EOS code       
!                   calculates & it should be included anywhere         
!                   these variables need to be referenced               
!                                                                       
!                                                                       
      COMMON /EOSVAR/ MUN,MUN_I,MUPROT,MUALFA,MUHAT,                    &
     &    BU, BS, BPRESS, UTOT, STOT, PTOT,                             &
     &    XPROT, XNUT, XALFA, XH, X, A,                                 &
     &    BUOUT, BUNUC, BPROUT, BPRNUC, BPRALF, BSOUT,BSNUC,BFOUT,      &
     &    BFNUC, BFALFA, FTOT, BFTOT,U_NUC, NSUBI, DNI, RSUBN,          &
     &    VSUBN, FBARY, PBARY, MUBARY, MU_MAT,                          &
     &    XALFA2, XNUT2, XPROT2, IT_NUM                                 
!                                                                       
!                                                                       
!                    This common contains the variables that are used   
!                    in solving the chemical and pressure equilibrium   
!                    equations                                          
!                                                                       
      COMMON /EQLEQN/ NNOUT, NPOUT, NOUT, VNOUT, VPOUT, F32_NO, F32_PO, &
     & MUN_O, MUP_O, ALFDNS, COMPX, SIGMA, SIGSGP, DU, DMU, DUP, DMUP,  &
     & SCRDU, SCRDUP, SCRD, ZETA, FSUBS, FSUBC, MUSUBT, FTRANS, NRATIO, &
     & VNI, VPI, ETA_NI, ETA_PI, MUP_I, F32_NI, F32_PI, BN, BP, BSUBP,  &
     & H, HPRIM, PSUBI, TSUBC, DHDX, DHDXX, SCRD2, SCRDPP, TFINAL       
!                                                                       
!                                                                       
!                   This common block contains the derivatives of       
!                   the thermodynamic potentials (both baryon and       
!                   the totals)                                         
!                                                                       
      COMMON /DERIVS/ DBMUDT, DBMUDN, DBMUDY, DMUDT, DMUDN, DMUDY,      &
     &                DBPDT, DBPDN, DBPDY, DPDT, DPDN, DPDY,            &
     &                DBSDT, DBSDN, DBSDY, DSDT, DSDN, DSDY,            &
     &                DBUDT, DBUDN, DBUDY, DUDT, DUDN, DUDY, GAM_S,     &
     &                DU_DT, DX_DT, DNI_DT, DEP_DT, DEN_DT,             &
     &                DU_DY, DX_DY, DNI_DY, DEP_DY, DEN_DY,             &
     &                DU_DN, DX_DN, DNI_DN, DEP_DN, DEN_DN              
!                                                                       
!                                                                       
!                   This common block contains the variables that       
!                   specify the particular interaction used             
      COMMON /SKYRME/ AA, BB, CC, DD, COMP, SYMM, SKYRMC, BIND_E,       &
     &                NSUBS, SYM_S, SIG_S                               
!                                                                       
!                                                                       
!                   This common block contains control flags            
      COMMON /FLAGS/ ADFLAG, NF_FLG, RSFLAG, DBFLAG 
!                                                                       
!                                                                       
!                   This common block contains the quality of solution  
!                   variables that determine how "zeroed" the equilbrium
!                   equations are                                       
      COMMON /QCHECK/ A1,A2,A3 
!                                                                       
      double precision ncomp 
      COMMON /DTEST/ DFDOM, DFDL_1,DFDL_2,DFDL_3,DNA_DT,DNA_DN,DNA_DY, &
     & ncomp                                        
