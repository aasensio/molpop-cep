	  ']  �   k820309    k          13.1        �T�Q                                                                                                           
       escape_cep.f90 ESCAPE_CEP                                                    
                                                          
                        �                                 
                @ @                                                 
                &                   &                                                      @ @                                   
                @ @                                                 
                &                                                    @ @                                                 
                &                                                    @ @                                                 
                &                                                    @ @                              	                   
                &                                                    @ @                              
                   
                &                                                                                                    #         @                                                     #SPLIN1%SIZE    #X    #Y    #YP1    #YPN    #Y2                                                     SIZE           
                                                   
               &                                                     
                                                    
 !             &                                                     
                                      
                
                                      
               
                                                    
 "    p          H r      7
S
O p        j            j                                      H r      7
S
O p        j            j                                                             @                                                   
                &                   &                                                    @                                                   
                &                                                                                            
                 
                 PERT�!	@        3.141592654D0                                                                                                                      
                &                                           %         @                                                  
       #EXPINT%ABS    #EXPINT%EXP    #EXPINT%LOG    #N    #X                                                     ABS                                                  EXP                                                  LOG                                                                                                       
                                                                    
                &                                                                                                    #         @                                                      #SPLINE%SIZE !   #XA "   #YA #   #Y2A $   #X %   #Y &                                               !     SIZE           
                                "                   
 %             &                                                     
                                 #                   
 &             &                                                                                     $                   
 (              &                                                     
                                %                   
 '             &                                                     
                                &                   
 $              &                                                                                        '                                                       (                                   &                   &                                                                                     )                   
                &                   &                                                                                     *                   
                &                   &                                                                                       +     
                   
                  {����56                   @ @                               ,                                                         -                     @                                .                   
                &                                                                                     /                   
                &                   &                                                    @                                0                   
                &                                                    @                                1                   
                &                                                                                        2                                                        3     
                  @                                4     
                  @                                5     
                                                 6                   
                &                                                                                     7                   
                &                                                                                       8     
                                                   9     
                    @                              :     P                  @                                ;     
                @                                <                   
                &                   &                                                    @                                =                   
                &                   &                                                    @                                >                   
                &                   &                                                                                     ?                   
                &                                                                                        @                     @                                A                   
                &                                                    @                                B                   
                &                                                    @                                C                   
                &                   &                                                    @                                D                   
                &                                                    @                                E                   
                &                                                    @                                F                   
                &                                                    @                                G                   
                &                   &                                                    @                                H                   
                &                                           #         @                                  I                  #NL -   #COLLIS%TRIM J   #COLLIS%ADJUSTL K   #IP L   #C M                                               J     TRIM                                             K     ADJUSTL           
                                  L                    
                                M                    
       p        5 r -   p          5 r -     5 r -       5 r -     5 r -                                                                  N                                                        O     
                   
                  mx���E                 @ @                              P                   
                &                                           #         @                                  Q                   #LUDCMP%SIZE R   #LUDCMP%MAXVAL S   #LUDCMP%DABS T   #A U   #INDX V   #D W                                               R     SIZE                                             S     MAXVAL                                             T     DABS           
                               U                   
               &                   &                                                     
                                 V                                  &                                                     
                                W     
       #         @                                  X                   #LUBKSB%SIZE Y   #A Z   #INDX [   #B \                                               Y     SIZE           
                                Z                   
 	             &                   &                                                     
                                  [                                 &                                                     
                                \                   
 
              &                                                                                       ]     
                  @                                 ^                       @ @                               _                                                         `                       @                                a     
                  @                                 b                       @                                 c                       @                                 d                       @ @                               e                       @                                 f                                                         g            #         @                                  h                                                                  i                                                      j                   
                &                                           #         @                                  k                   #REGRID%TRIM l   #REGRID%ADJUSTL m   #REGRID%ALLOCATED n   #RESET o                                               l     TRIM                                             m     ADJUSTL                                             n     ALLOCATED                                            o            #         @                                  p                   #LLSLV%DABS q   #LLSLV%DSQRT r   #S s   #X t   #N u   #NR v                                               q     DABS                                             r     DSQRT          
                               s                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                   
                               t                    
     p          5 O p            5 O p                                    
                                 u                     
                                 v           #         @                                   w                    #PRECALCULATE_BETA%LOG10 x   #PRECALCULATE_BETA%ALLOCATED y                                              x     LOG10                                            y     ALLOCATED %         @                               z                   
       #BETA%EXP {   #BETA%ABS |   #BETA%LOG }   #BETA%DSQRT ~   #BETA%SQRT    #TAU_IN �                                              {     EXP                                            |     ABS                                            }     LOG                                            ~     DSQRT                                                 SQRT                                           �     
       %         @                               �                   
       #ALPHAP%EXP �   #ALPHAP%SQRT �   #TAU_IN �                                              �     EXP                                            �     SQRT                                           �     
       %         @                               �                   
       #BETA2%EXP �   #BETA2%ABS �   #BETA2%LOG �   #BETA2%DSQRT �   #BETA2%SQRT �   #BETA2%LOG10 �   #TAU_IN �                                              �     EXP                                            �     ABS                                            �     LOG                                            �     DSQRT                                            �     SQRT                                            �     LOG10           D @                              �     
       %         @                               �                   
       #ALPHAP2%EXP �   #ALPHAP2%ABS �   #ALPHAP2%LOG �   #ALPHAP2%DSQRT �   #ALPHAP2%SQRT �   #ALPHAP2%LOG10 �   #TAU_IN �                                              �     EXP                                            �     ABS                                            �     LOG                                            �     DSQRT                                            �     SQRT                                            �     LOG10           D @                              �     
       #         @                                   �                   #ADAPT_COL_DENSITY_STRATEGY%ADJUSTL �   #ADAPT_COL_DENSITY_STRATEGY%TRIM �   #ADAPT_COL_DENSITY_STRATEGY%SUM �   #ADAPT_COL_DENSITY_STRATEGY%MINVAL �   #ADAPT_COL_DENSITY_STRATEGY%INT �   #ADAPT_COL_DENSITY_STRATEGY%MAXVAL �   #ADAPT_COL_DENSITY_STRATEGY%ABS �   #ADAPT_COL_DENSITY_STRATEGY%LOG10 �   #X �                                              �     ADJUSTL                                            �     TRIM                                            �     SUM                                            �     MINVAL                                            �     INT                                            �     MAXVAL                                            �     ABS                                            �     LOG10                                           �                   
 	              &                                           #         @                                  �                    #X �                                             �                   
 
              &                                           #         @                                  �                   #RT1D_CEP%DABS �   #NP �   #IT �                                              �     DABS                                            �                                                       �            #         @                                  �                   #FUNCV%MATMUL �   #X �   #N �   #FVEC �                                              �     MATMUL          D @                              �                    
     p          5 � p        r �       5 � p        r �                                                                �                     D                                �                    
     p          5 � p        r �       5 � p        r �                     #         @                                  �                   #FUNCV_ANALYTIC%MATMUL �   #FUNCV_ANALYTIC%SUM �   #X �   #N �   #FVEC �   #FJAC �                                              �     MATMUL                                            �     SUM          D @                              �                    
     p          5 � p        r �       5 � p        r �                                                                �                     D                                �                    
     p          5 � p        r �       5 � p        r �                              D                                �                    
       p        5 � p        r �   p          5 � p        r �     5 � p        r �       5 � p        r �     5 � p        r �                     #         @                                  �                   #FDJAC%DABS �   #N �   #X �   #FVEC �   #DF �                                              �     DABS           D @                               �                     D @                              �                    
     p          5 � p        r �       5 � p        r �                                                              �                    
     p          5 � p        r �       5 � p        r �                              D                                �                    
       p        5 � p        r �   p          5 � p        r �     5 � p        r �       5 � p        r �     5 � p        r �                     #         @                                  �                    #X �   #N �   #FVEC �   #FJAC �   #DERIVATIVES �            D @                              �                    
 !    p          5 � p        r �       5 � p        r �                               D @                               �                     D @                              �                    
 "    p          5 � p        r �       5 � p        r �                               D @                              �                   
 #              &                   &                                                                                      �            #         @                                  �                   #MNEWT%DABS �   #MNEWT%MINVAL �   #MNEWT%MAXVAL �   #MNEWT%ABS �   #NTRIAL �   #X �   #N �   #TOLX �   #TOLF �   #DERIVATIVES �   #ERROR �                                              �     DABS                                            �     MINVAL                                            �     MAXVAL                                            �     ABS                                            �                     D @                              �                    
 %    p          5 � p        r �       5 � p        r �                               D @                               �                                                      �     
                                                 �     
                 D @                               �                      D                                 �            #         @                                   �                    #N �   #X �   #FVEC �   #FJAC �   #LDFJAC �   #IFLAG �             D @                               �                     D @                              �                    
 *    p          5 � p        r �       5 � p        r �                              D                                �                    
 +    p          5 � p        r �       5 � p        r �                              D                                �                    
 ,      p        5 � p        r �   p          5 � p        r �     5 � p        r �       5 � p        r �     5 � p        r �                                                                �                                                       �            #         @                                   �                    #NTRIAL �   #X �   #N �   #TOLX �   #TOLF �   #DERIVATIVES �   #ERROR �                                              �                     D                                �                    
 /    p          5 � p        r �       5 � p        r �                                                                �                                                      �     
                                                 �     
                                                  �                      D                                 �            #         @                                   �                   #CEP_SOLVER%MAXVAL �   #CEP_SOLVER%ABS �   #ERROR �                                              �     MAXVAL                                            �     ABS           D @                               �            #         @                                  �                     %         @                                �                          #NG%DABS �   #Y �   #M �   #N �   #YY �                                              �     DABS          
D                                �                    
 9    p          5 � p        r �       5 � p        r �                               
                                  �                     
  @                               �                  B  
D                                �                    
 :      p        5 � p        r �   p          5 � p        r �     1     5 � p        r �     1                      �   "      fn#fn    �   @   J   CONSTANTS_CEP      @   J   FUNCTIONS_CEP    B  @   J   MATHS_CEP #   �  �       SCRATCH+GLOBAL_CEP %   &  @       PRECISION+GLOBAL_CEP $   f  �       TAU_BETA+GLOBAL_CEP )   �  �       BETA_FUNCTION+GLOBAL_CEP '   ~  �       BETA_SPLINE+GLOBAL_CEP +   
  �       ALPHAP_FUNCTION+GLOBAL_CEP )   �  �       ALPHAP_SPLINE+GLOBAL_CEP #   "  @       VERBOSE+GLOBAL_CEP !   b  �       SPLIN1+MATHS_CEP &   �  =      SPLIN1%SIZE+MATHS_CEP #      �   a   SPLIN1%X+MATHS_CEP #   �  �   a   SPLIN1%Y+MATHS_CEP %   8  @   a   SPLIN1%YP1+MATHS_CEP %   x  @   a   SPLIN1%YPN+MATHS_CEP $   �  T  a   SPLIN1%Y2+MATHS_CEP    	  �       TAU+GLOBAL_CEP    �	  �       B+GLOBAL_CEP !   <
  }       PI+CONSTANTS_CEP (   �
  @       N_QUADR_BETA+GLOBAL_CEP     �
  �       W_E3+GLOBAL_CEP !   �  �       EXPINT+MATHS_CEP %     <      EXPINT%ABS+MATHS_CEP %   O  <      EXPINT%EXP+MATHS_CEP %   �  <      EXPINT%LOG+MATHS_CEP #   �  @   a   EXPINT%N+MATHS_CEP #     @   a   EXPINT%X+MATHS_CEP     G  �       X_E3+GLOBAL_CEP 1   �  @       ESCAPE_PROB_ALGORITHM+GLOBAL_CEP !     �       SPLINE+MATHS_CEP &   �  =      SPLINE%SIZE+MATHS_CEP $   �  �   a   SPLINE%XA+MATHS_CEP $   \  �   a   SPLINE%YA+MATHS_CEP %   �  �   a   SPLINE%Y2A+MATHS_CEP #   t  �   a   SPLINE%X+MATHS_CEP #      �   a   SPLINE%Y+MATHS_CEP    �  @       NR+GLOBAL_CEP !   �  �       ITRAN+GLOBAL_CEP !   p  �       DTRAN+GLOBAL_CEP "     �       DLEVEL+GLOBAL_CEP %   �  p       TWOHC2+CONSTANTS_CEP    (  @       NZ+GLOBAL_CEP    h  @       NL+GLOBAL_CEP     �  �       CHIL+GLOBAL_CEP $   4  �       DOPPLERW+GLOBAL_CEP    �  �       SL+GLOBAL_CEP    d  �       DZ+GLOBAL_CEP +   �  @       KTHICK_STRATEGY+GLOBAL_CEP )   0  @       TAU_THRESHOLD+GLOBAL_CEP ,   p  @       FACTOR_ABUNDANCE+GLOBAL_CEP '   �  @       COL_DENSITY+GLOBAL_CEP %   �  �       ABUNDANCE+GLOBAL_CEP    |  �       NH+GLOBAL_CEP !     @       ABUND+GLOBAL_CEP ,   H  @       HYDROGEN_DENSITY+GLOBAL_CEP 4   �  @       FILE_PHYSICAL_CONDITIONS+GLOBAL_CEP )   �  @       COL_THRESHOLD+GLOBAL_CEP !     �       DSLDN+GLOBAL_CEP #   �  �       DCHILDN+GLOBAL_CEP "   P  �       DTAUDN+GLOBAL_CEP     �  �       POPL+GLOBAL_CEP *   �  @       OPTICALLY_THIN+GLOBAL_CEP '   �  �       LSTAR_TOTAL+GLOBAL_CEP &   L  �       JBAR_TOTAL+GLOBAL_CEP )   �  �       DJBAR_TOTALDN+GLOBAL_CEP &   |  �       FLUX_TOTAL+GLOBAL_CEP !     �       LSTAR+GLOBAL_CEP     �  �       JBAR+GLOBAL_CEP #       �       DJBARDN+GLOBAL_CEP     �   �       FLUX+GLOBAL_CEP %   P!  �       COLLIS+FUNCTIONS_CEP *   �!  =      COLLIS%TRIM+FUNCTIONS_CEP -   "  @      COLLIS%ADJUSTL+FUNCTIONS_CEP (   Q"  @   a   COLLIS%IP+FUNCTIONS_CEP '   �"  �   a   COLLIS%C+FUNCTIONS_CEP     e#  @       NACT+GLOBAL_CEP #   �#  p       PI4H+CONSTANTS_CEP    $  �       POP+GLOBAL_CEP !   �$  �       LUDCMP+MATHS_CEP &   6%  =      LUDCMP%SIZE+MATHS_CEP (   s%  ?      LUDCMP%MAXVAL+MATHS_CEP &   �%  =      LUDCMP%DABS+MATHS_CEP #   �%  �   a   LUDCMP%A+MATHS_CEP &   �&  �   a   LUDCMP%INDX+MATHS_CEP #   '  @   a   LUDCMP%D+MATHS_CEP !   _'  q       LUBKSB+MATHS_CEP &   �'  =      LUBKSB%SIZE+MATHS_CEP #   (  �   a   LUBKSB%A+MATHS_CEP &   �(  �   a   LUBKSB%INDX+MATHS_CEP #   =)  �   a   LUBKSB%B+MATHS_CEP )   �)  @       CEP_PRECISION+GLOBAL_CEP 8   	*  @       CEP_CONVERGENCE_MAXIMUM_ITER+GLOBAL_CEP /   I*  @       NEWTON_MAXIMUM_ITER+GLOBAL_CEP (   �*  @       WHICH_SCHEME+GLOBAL_CEP '   �*  @       RELAT_ERROR+GLOBAL_CEP     	+  @       ITER+GLOBAL_CEP "   I+  @       ITRACC+GLOBAL_CEP "   �+  @       IACCEL+GLOBAL_CEP     �+  @       NORD+GLOBAL_CEP $   	,  @       NINTRACC+GLOBAL_CEP #   I,  @       N_ITERS+GLOBAL_CEP 2   �,  H       CORRECT_POPULATIONS+FUNCTIONS_CEP -   �,  @       RELAT_ERROR_LEVEL+GLOBAL_CEP /   -  �       POP_PREVIOUS_REGRID+GLOBAL_CEP %   �-  �       REGRID+FUNCTIONS_CEP *   +.  =      REGRID%TRIM+FUNCTIONS_CEP -   h.  @      REGRID%ADJUSTL+FUNCTIONS_CEP /   �.  B      REGRID%ALLOCATED+FUNCTIONS_CEP +   �.  @   a   REGRID%RESET+FUNCTIONS_CEP     */  �       LLSLV+MATHS_CEP %   �/  =      LLSLV%DABS+MATHS_CEP &   �/  >      LLSLV%DSQRT+MATHS_CEP "   +0  �   a   LLSLV%S+MATHS_CEP "   '1  �   a   LLSLV%X+MATHS_CEP "   �1  @   a   LLSLV%N+MATHS_CEP #   2  @   a   LLSLV%NR+MATHS_CEP "   K2  �       PRECALCULATE_BETA (   �2  >      PRECALCULATE_BETA%LOG10 ,   3  B      PRECALCULATE_BETA%ALLOCATED    Q3  �       BETA    �3  <      BETA%EXP    24  <      BETA%ABS    n4  <      BETA%LOG    �4  >      BETA%DSQRT    �4  =      BETA%SQRT    %5  @   a   BETA%TAU_IN    e5  }       ALPHAP    �5  <      ALPHAP%EXP    6  =      ALPHAP%SQRT    [6  @   a   ALPHAP%TAU_IN    �6  �       BETA2    V7  <      BETA2%EXP    �7  <      BETA2%ABS    �7  <      BETA2%LOG    
8  >      BETA2%DSQRT    H8  =      BETA2%SQRT    �8  >      BETA2%LOG10    �8  @   a   BETA2%TAU_IN    9  �       ALPHAP2    �9  <      ALPHAP2%EXP    :  <      ALPHAP2%ABS    B:  <      ALPHAP2%LOG    ~:  >      ALPHAP2%DSQRT    �:  =      ALPHAP2%SQRT    �:  >      ALPHAP2%LOG10    7;  @   a   ALPHAP2%TAU_IN +   w;  |      ADAPT_COL_DENSITY_STRATEGY 3   �<  @      ADAPT_COL_DENSITY_STRATEGY%ADJUSTL 0   3=  =      ADAPT_COL_DENSITY_STRATEGY%TRIM /   p=  <      ADAPT_COL_DENSITY_STRATEGY%SUM 2   �=  ?      ADAPT_COL_DENSITY_STRATEGY%MINVAL /   �=  <      ADAPT_COL_DENSITY_STRATEGY%INT 2   '>  ?      ADAPT_COL_DENSITY_STRATEGY%MAXVAL /   f>  <      ADAPT_COL_DENSITY_STRATEGY%ABS 1   �>  >      ADAPT_COL_DENSITY_STRATEGY%LOG10 -   �>  �   a   ADAPT_COL_DENSITY_STRATEGY%X #   l?  O       CALCJBAR_LSTAR_CEP %   �?  �   a   CALCJBAR_LSTAR_CEP%X    G@  k       RT1D_CEP    �@  =      RT1D_CEP%DABS    �@  @   a   RT1D_CEP%NP    /A  @   a   RT1D_CEP%IT    oA  r       FUNCV    �A  ?      FUNCV%MATMUL     B  �   a   FUNCV%X    �B  @   a   FUNCV%N    C  �   a   FUNCV%FVEC    �C  �       FUNCV_ANALYTIC &   eD  ?      FUNCV_ANALYTIC%MATMUL #   �D  <      FUNCV_ANALYTIC%SUM !   �D  �   a   FUNCV_ANALYTIC%X !   �E  @   a   FUNCV_ANALYTIC%N $   �E  �   a   FUNCV_ANALYTIC%FVEC $   �F  $  a   FUNCV_ANALYTIC%FJAC    �G  x       FDJAC    $H  =      FDJAC%DABS    aH  @   a   FDJAC%N    �H  �   a   FDJAC%X    UI  �   a   FDJAC%FVEC    	J  $  a   FDJAC%DF    -K  {       USRFUN    �K  �   a   USRFUN%X    \L  @   a   USRFUN%N    �L  �   a   USRFUN%FVEC    PM  �   a   USRFUN%FJAC #   �M  @   a   USRFUN%DERIVATIVES    4N  �       MNEWT    	O  =      MNEWT%DABS    FO  ?      MNEWT%MINVAL    �O  ?      MNEWT%MAXVAL    �O  <      MNEWT%ABS     P  @   a   MNEWT%NTRIAL    @P  �   a   MNEWT%X    �P  @   a   MNEWT%N    4Q  @   a   MNEWT%TOLX    tQ  @   a   MNEWT%TOLF "   �Q  @   a   MNEWT%DERIVATIVES    �Q  @   a   MNEWT%ERROR    4R  �       FCN_NAG    �R  @   a   FCN_NAG%N    �R  �   a   FCN_NAG%X    �S  �   a   FCN_NAG%FVEC    ]T  $  a   FCN_NAG%FJAC    �U  @   a   FCN_NAG%LDFJAC    �U  @   a   FCN_NAG%IFLAG    V  �       MNEWT_NAG !   �V  @   a   MNEWT_NAG%NTRIAL    �V  �   a   MNEWT_NAG%X    �W  @   a   MNEWT_NAG%N    �W  @   a   MNEWT_NAG%TOLX    X  @   a   MNEWT_NAG%TOLF &   GX  @   a   MNEWT_NAG%DERIVATIVES     �X  @   a   MNEWT_NAG%ERROR    �X  ~       CEP_SOLVER "   EY  ?      CEP_SOLVER%MAXVAL    �Y  <      CEP_SOLVER%ABS !   �Y  @   a   CEP_SOLVER%ERROR     Z  H       ACCELERATION    HZ  z       NG    �Z  =      NG%DABS    �Z  �   a   NG%Y    �[  @   a   NG%M    �[  @   a   NG%N    3\  �   a   NG%YY 