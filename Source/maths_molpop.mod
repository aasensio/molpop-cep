	  ud  "  k820309    q          13.1        ��[R                                                                                                           
       maths_molpop.f90 MATHS_MOLPOP                                                    
                                                                                                                 
                @ @                                                 
                &                   &                                                    @ @                                                 
                &                   &                                                                                                        
                &                   &                                                                                                        
                &                   &                                                                                                                                              	                   
                &                   &                                                    @                                
                   
                &                   &                                                                                                             @ @                                                 
                &                   &                                                      @                                     
                                                        
                  @ @                                   
                                                                    
                &                   &                                                      @                                     
                                                                    
                &                   &                                                      @                                                      @                                                                    &                                                    @                                                                    &                                                      @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                     
                  @                                      
                                                   !     
                 
                 ������@        2.725D0         @ @                              "                   
                &                   &                                                    @                                #                   
                &                   &                                                    @                                $                   
                &                   &                                                    @                                %                   
                &                   &                                                    @                                &                   
                &                   &                                                                                     '                   
                &                   &                                           #         @                                  (                    #IUNIT )                                              )            %         @                                 *                          #ERROR_MESSAGE%ADJUSTL +   #ERROR_MESSAGE%TRIM ,   #OPT -   #STR .                                              +     ADJUSTL                                            ,     TRIM            @                              -     �                                  @                              .     �                       #         @                                   /                    #X 0                                             0                    
     p          5 r        5 r                      #         @                                  1                   #BETA_SLAB%DABS 2   #BETA_SLAB%DLOG 3   #BETA_SLAB%DEXP 4   #BETA_SLAB%DSQRT 5   #TAUIN 6   #BETA 7   #DBDX 8                                              2     DABS                                            3     DLOG                                            4     DEXP                                            5     DSQRT                                            6     
                 D                                 7     
                 D                                 8     
       #         @                                  9                   #BETA_LVG%DABS :   #BETA_LVG%DLOG ;   #BETA_LVG%DEXP <   #BETA_LVG%DSQRT =   #TAUIN >   #BETA ?   #DBDX @                                              :     DABS                                            ;     DLOG                                            <     DEXP                                            =     DSQRT                                            >     
                 D                                 ?     
                 D                                 @     
       #         @                                   A                   #MATIN1%DABS B   #A C   #DIM1 D   #N1 E   #DIM2 F   #N2 G   #INDEX H   #NERROR I   #DETERM J                                              B     DABS           D                                 C                   
     p          p            p                                                                     D                                                       E                                                       F                                                       G                      D                                 H                        p          p            p                                    D                                 I                      D                                 J     
       #         @                                   K                    #ARRAY L   #N M   #INDEX N   #JNDEX O   #NBIG P                                              L                   
               &                   &                                                                                      M                      D                                 N                                  &                                                     D                                 O                                  &                                                     D                                 P            #         @                                   Q                    #V R   #N S   #INDEX T   #NBIG U                                             R                    
     p          5 � p        r S       5 � p        r S                                                                S                     D                                 T                         p          5 � p        r U       5 � p        r U                                                                U            %         @                                V                    
       #V W   #N X                                             W                    
 	    p          5 � p        r X       5 � p        r X                                                                X            %         @                                Y                           #LINE Z             D @                               Z                     1 %         @                               [                          #INMIN%LEN \   #STRING ]                                              \     LEN            @                               ]                     1 %         @                               ^                          #INMAX%LEN _   #STRING `                                              _     LEN            @                               `                     1 #         @                                   a                    #S1 b   #S2 c   #S d             D @                               b                     1           D @                               c                     1           D                                 d                     1 #         @                                   e                    #S1 f   #S2 g   #S3 h   #S i             D @                               f                     1           D @                               g                     1           D @                               h                     1           D                                 i                     1 #         @                                   j                    #N k   #X l   #Y m   #Y2 n                                              k                                                       l     d              
 
    p          p d           p d                                                                    m     d              
     p          p d           p d                                   D                                 n     d              
     p          p d           p d                         #         @                                   o                    #N p   #XA q   #YA r   #Y2A s   #X t   #Y u                                              p                                                       q     d              
     p          p d           p d                                                                    r     d              
     p          p d           p d                                                                    s     d              
     p          p d           p d                                                                    t     
                 D                                 u     
       %         @                               v                   
       #BESSI0%DABS w   #BESSI0%DEXP x   #BESSI0%DSQRT y   #X z                                              w     DABS                                            x     DEXP                                            y     DSQRT            @                               z     
       %         @                               {                   
       #BESSK0%DLOG |   #BESSK0%DEXP }   #BESSK0%DSQRT ~   #X                                               |     DLOG                                            }     DEXP                                            ~     DSQRT           D @                                    
       %         @                               �                   
       #BESSI1%DABS �   #BESSI1%DEXP �   #BESSI1%DSQRT �   #X �                                              �     DABS                                            �     DEXP                                            �     DSQRT            @                               �     
       %         @                               �                   
       #BESSK1%DLOG �   #BESSK1%DEXP �   #BESSK1%DSQRT �   #X �                                              �     DLOG                                            �     DEXP                                            �     DSQRT           D @                               �     
       %         @                                �                    
       #N �   #X �                                              �                      D @                               �     
       %         @                                �                           #STRING �             D @                               �                     1 #         @                                   �                    #NN �   #STRING �                                              �                      D                                 �                     1 #         @                                   �                   #FILE_NAME%CHAR �   #STR1 �   #STR2 �   #STR3 �                                              �     CHAR           D @                               �                     1           D @                               �                     1           D                                 �                     1 $         @                                �                           #CARD �   #I �                                                      �     �                                                                  �            %         @                                �                           #CR �                                              �                            %         @                                �                           #CR �                                              �                            %         @                                �                           #CR �                                              �                            %         @                                �                           #CR �                                              �                            %         @                                �                           #CR �                                              �                            %         @                                �                           #CR �                                              �                            %         @                                �                          #IVAL%ICHAR �   #CR �                                              �     ICHAR            @                               �                            #         @                                   �                   #RDINPS%INDEX �   #RDINPS%ADJUSTL �   #RDINPS%TRIM �   #RDINPS%LEN �   #EQUAL �   #IUNIT �   #STR �                                              �     INDEX                                            �     ADJUSTL                                            �     TRIM                                            �     LEN                                            �                      D                                 �                      D                                 �                     1 #         @                                   �                   #RDINPS2%ICHAR �   #RDINPS2%CHAR �   #RDINPS2%INDEX �   #RDINPS2%ADJUSTL �   #RDINPS2%TRIM �   #RDINPS2%LEN �   #EQUAL �   #IUNIT �   #STR �   #LENGTH �   #UCASE �                                              �     ICHAR                                            �     CHAR                                            �     INDEX                                            �     ADJUSTL                                            �     TRIM                                            �     LEN                                            �                      D                                 �                      D                                 �                     1           D                                 �                                                       �            %         @                                �                   
       #RDINP%INDEX �   #RDINP%LEN �   #EQUAL �   #IUNIT �   #OUTUNIT �                                              �     INDEX                                            �     LEN           
                                  �                     D                                 �                                                       �            %         @                               �                    
       #CARD �   #FIRST �   #LAST �   #DECIMAL �   #TERM �                                              �                     1           D                                 �                                                       �                                                       �                      D @                               �                            %         @                                �                   
       #HERMITE%DEXP �   #IT �   #DTAU �   #DELTA �   #I2 �   #J2 �   #M �   #ITER �                                              �     DEXP                                            �                                                       �     
                                                  �     
                                                  �     2                   p          p 2           p 2                                                                    �     2                   p          p 2           p 2                                                                    �                                                       �            %         @                               �                   
       #EXPINT%DABS �   #EXPINT%DLOG �   #EXPINT%DEXP �   #N �   #X �                                              �     DABS                                            �     DLOG                                            �     DEXP                                            �                      D  @                               �     
       #         @                                   �                     %         @                                �                   
       #PART_FUNC%DEXP �   #TT �   #ACON �   #JMAX �                                              �     DEXP                                            �     
                                                  �     
                                                  �            %         @                                �                   
       #THREEJ%SQRT �   #THREEJ%MIN �   #THREEJ%MAX �   #THREEJ%ABS �   #THREEJ%INT �   #THREEJ%DSQRT �   #J1 �   #J2 �   #J �                                              �     SQRT                                            �     MIN                                            �     MAX                                            �     ABS                                            �     INT                                            �     DSQRT           D @                               �                      D @                               �                      D @                               �            %         @                               �                   
       #BINOM%DBLE �   #N �   #R �                                              �     DBLE            @                               �                                                       �            #         @                                   �                    #CONST%DSQRT �                                              �     DSQRT %         @                                �                   
       #PRLOG%DLOG10 �   #X �                                              �     DLOG10            @                               �     
       %         @                               �                   
       #PLEXP%DABS �   #PLEXP%DEXP �   #X �                                              �     DABS                                            �     DEXP            @                               �     
       %         @                               �                   
       #INV_PLEXP%DLOG �   #P �                                              �     DLOG           
                                  �     
      %         @                                �                   
       #TBR_TX%DABS �   #TBR_TX%DEXP �   #TL �   #TX �   #TAUL �                                              �     DABS                                            �     DEXP           
                                  �     
                
                                  �     
                
                                  �     
      %         @                                �                   
       #TBR_I%DABS �   #TBR_I%DEXP    #NU   #I   #TAUL                                              �     DABS                                                 DEXP           
                                      
                
                                      
                
                                      
      #         @                                                     #N   #N1   #N2   #X   #Y 	  #INTEGRAL 
                                                                                                                                                                                                                                 
     p          5 � p        r       5 � p        r                                                               	                   
     p          5 � p        r       5 � p        r                               D                                 
    
       #         @                                                    #INTERPOLATEEXTERNALFILE%SIZE   #INTERPOLATEEXTERNALFILE%ADJUSTL   #INTERPOLATEEXTERNALFILE%TRIM   #FILENAME   #WL   #OUT   #NORM   #ERROR                                                  SIZE                                                ADJUSTL                                                TRIM            @                                                   1           @                                                
               &                   &                                                     D                                                  
               &                   &                                                                                                           D                                            #         @                                                     #RAD_FILE%ADJUSTL   #RAD_FILE%TRIM   #FNAME   #JBOL   #STR   #L   #ERROR                                                  ADJUSTL                                                TRIM            @                                                   1                                               
                                                                      1                                                                 D @                                          #         @                                                     #DUST_RAD%DEXP   #TD   #TAU_D   #STR    #L !                                                 DEXP                                                
                                                      
                                                                       1                                            !              �   &      fn#fn    �   @   j   GLOBAL_MOLPOP       @       N+GLOBAL_MOLPOP     F  @       R+GLOBAL_MOLPOP "   �  �       TAU+GLOBAL_MOLPOP "   *  �       ESC+GLOBAL_MOLPOP #   �  �       TAUX+GLOBAL_MOLPOP "   r  �       GAP+GLOBAL_MOLPOP -     @       DUSTABSORPTION+GLOBAL_MOLPOP $   V  �       QDUST+GLOBAL_MOLPOP !   �  �       XD+GLOBAL_MOLPOP $   �  @       KBETA+GLOBAL_MOLPOP %   �  �       DBDTAU+GLOBAL_MOLPOP %   �  @       ROOTPI+GLOBAL_MOLPOP "   �  @       EPS+GLOBAL_MOLPOP !     @       PI+GLOBAL_MOLPOP     B  �       A+GLOBAL_MOLPOP !   �  @       BK+GLOBAL_MOLPOP #   &  �       FREQ+GLOBAL_MOLPOP "   �  @       NUM+GLOBAL_MOLPOP $   
	  �       UPLIN+GLOBAL_MOLPOP %   �	  �       LOWLIN+GLOBAL_MOLPOP !   "
  @       CL+GLOBAL_MOLPOP "   b
  @       HPL+GLOBAL_MOLPOP "   �
  @       XME+GLOBAL_MOLPOP "   �
  @       XMP+GLOBAL_MOLPOP %   "  @       SOLARL+GLOBAL_MOLPOP %   b  @       SOLARM+GLOBAL_MOLPOP %   �  @       SOLARR+GLOBAL_MOLPOP !   �  @       PC+GLOBAL_MOLPOP $   "  @       TWOPI+GLOBAL_MOLPOP %   b  @       FOURPI+GLOBAL_MOLPOP $   �  @       EITPI+GLOBAL_MOLPOP #   �  w       TCMB+GLOBAL_MOLPOP !   Y  �       WL+GLOBAL_MOLPOP "   �  �       RAD+GLOBAL_MOLPOP '   �  �       RAD_TAU0+GLOBAL_MOLPOP '   E  �       RAD_TAUT+GLOBAL_MOLPOP +   �  �       RAD_INTERNAL+GLOBAL_MOLPOP "   �  �       TIJ+GLOBAL_MOLPOP    1  S       PASS_HEADER "   �  @   a   PASS_HEADER%IUNIT    �  �       ERROR_MESSAGE &   Y  @      ERROR_MESSAGE%ADJUSTL #   �  =      ERROR_MESSAGE%TRIM "   �  P   a   ERROR_MESSAGE%OPT "   &  P   a   ERROR_MESSAGE%STR    v  O       OPTDEP    �  �   a   OPTDEP%X    Y  �       BETA_SLAB      =      BETA_SLAB%DABS    N  =      BETA_SLAB%DLOG    �  =      BETA_SLAB%DEXP     �  >      BETA_SLAB%DSQRT       @   a   BETA_SLAB%TAUIN    F  @   a   BETA_SLAB%BETA    �  @   a   BETA_SLAB%DBDX    �  �       BETA_LVG    z  =      BETA_LVG%DABS    �  =      BETA_LVG%DLOG    �  =      BETA_LVG%DEXP    1  >      BETA_LVG%DSQRT    o  @   a   BETA_LVG%TAUIN    �  @   a   BETA_LVG%BETA    �  @   a   BETA_LVG%DBDX    /  �       MATIN1    �  =      MATIN1%DABS      �   a   MATIN1%A    �  @   a   MATIN1%DIM1    �  @   a   MATIN1%N1    '  @   a   MATIN1%DIM2    g  @   a   MATIN1%N2    �  �   a   MATIN1%INDEX    ;  @   a   MATIN1%NERROR    {  @   a   MATIN1%DETERM    �  z       ORDERA    5  �   a   ORDERA%ARRAY    �  @   a   ORDERA%N      �   a   ORDERA%INDEX    �  �   a   ORDERA%JNDEX    1  @   a   ORDERA%NBIG    q  k       ORDERV    �  �   a   ORDERV%V    �   @   a   ORDERV%N    �   �   a   ORDERV%INDEX    �!  @   a   ORDERV%NBIG    �!  ^       VSUM    ""  �   a   VSUM%V    �"  @   a   VSUM%N    #  Z       EMPTY    p#  L   a   EMPTY%LINE    �#  k       INMIN    '$  <      INMIN%LEN    c$  L   a   INMIN%STRING    �$  k       INMAX    %  <      INMAX%LEN    V%  L   a   INMAX%STRING    �%  _       ATTACH2    &  L   a   ATTACH2%S1    M&  L   a   ATTACH2%S2    �&  L   a   ATTACH2%S    �&  g       ATTACH3    L'  L   a   ATTACH3%S1    �'  L   a   ATTACH3%S2    �'  L   a   ATTACH3%S3    0(  L   a   ATTACH3%S    |(  e       SPLINE    �(  @   a   SPLINE%N    !)  �   a   SPLINE%X    �)  �   a   SPLINE%Y    I*  �   a   SPLINE%Y2    �*  v       SPLINT    S+  @   a   SPLINT%N    �+  �   a   SPLINT%XA    ',  �   a   SPLINT%YA    �,  �   a   SPLINT%Y2A    O-  @   a   SPLINT%X    �-  @   a   SPLINT%Y    �-  �       BESSI0    Z.  =      BESSI0%DABS    �.  =      BESSI0%DEXP    �.  >      BESSI0%DSQRT    /  @   a   BESSI0%X    R/  �       BESSK0    �/  =      BESSK0%DLOG    0  =      BESSK0%DEXP    W0  >      BESSK0%DSQRT    �0  @   a   BESSK0%X    �0  �       BESSI1    `1  =      BESSI1%DABS    �1  =      BESSI1%DEXP    �1  >      BESSI1%DSQRT    2  @   a   BESSI1%X    X2  �       BESSK1    �2  =      BESSK1%DLOG     3  =      BESSK1%DEXP    ]3  >      BESSK1%DSQRT    �3  @   a   BESSK1%X    �3  ^       BESSK    94  @   a   BESSK%N    y4  @   a   BESSK%X    �4  \       FSIZE    5  L   a   FSIZE%STRING    a5  \       CLEAR_STRING     �5  @   a   CLEAR_STRING%NN $   �5  L   a   CLEAR_STRING%STRING    I6  z       FILE_NAME    �6  =      FILE_NAME%CHAR     7  L   a   FILE_NAME%STR1    L7  L   a   FILE_NAME%STR2    �7  L   a   FILE_NAME%STR3    �7  i       CHR    M8  P   a   CHR%CARD    �8  @   a   CHR%I    �8  X       BAR    59  P   a   BAR%CR    �9  X       DIGIT    �9  P   a   DIGIT%CR    -:  X       MINUS    �:  P   a   MINUS%CR    �:  X       SGN    -;  P   a   SGN%CR    };  X       DOT    �;  P   a   DOT%CR    %<  X       EE    }<  P   a   EE%CR    �<  h       IVAL    5=  >      IVAL%ICHAR    s=  P   a   IVAL%CR    �=  �       RDINPS    q>  >      RDINPS%INDEX    �>  @      RDINPS%ADJUSTL    �>  =      RDINPS%TRIM    ,?  <      RDINPS%LEN    h?  @   a   RDINPS%EQUAL    �?  @   a   RDINPS%IUNIT    �?  L   a   RDINPS%STR    4@  �       RDINPS2    "A  >      RDINPS2%ICHAR    `A  =      RDINPS2%CHAR    �A  >      RDINPS2%INDEX     �A  @      RDINPS2%ADJUSTL    B  =      RDINPS2%TRIM    XB  <      RDINPS2%LEN    �B  @   a   RDINPS2%EQUAL    �B  @   a   RDINPS2%IUNIT    C  L   a   RDINPS2%STR    `C  @   a   RDINPS2%LENGTH    �C  @   a   RDINPS2%UCASE    �C  �       RDINP    sD  >      RDINP%INDEX    �D  <      RDINP%LEN    �D  @   a   RDINP%EQUAL    -E  @   a   RDINP%IUNIT    mE  @   a   RDINP%OUTUNIT    �E  �       VAL_DP    3F  L   a   VAL_DP%CARD    F  @   a   VAL_DP%FIRST    �F  @   a   VAL_DP%LAST    �F  @   a   VAL_DP%DECIMAL    ?G  P   a   VAL_DP%TERM    �G  �       HERMITE    /H  =      HERMITE%DEXP    lH  @   a   HERMITE%IT    �H  @   a   HERMITE%DTAU    �H  @   a   HERMITE%DELTA    ,I  �   a   HERMITE%I2    �I  �   a   HERMITE%J2    TJ  @   a   HERMITE%M    �J  @   a   HERMITE%ITER    �J  �       EXPINT    eK  =      EXPINT%DABS    �K  =      EXPINT%DLOG    �K  =      EXPINT%DEXP    L  @   a   EXPINT%N    \L  @   a   EXPINT%X    �L  H       NUMLIN    �L  �       PART_FUNC    dM  =      PART_FUNC%DEXP    �M  @   a   PART_FUNC%TT    �M  @   a   PART_FUNC%ACON    !N  @   a   PART_FUNC%JMAX    aN  �       THREEJ    +O  =      THREEJ%SQRT    hO  <      THREEJ%MIN    �O  <      THREEJ%MAX    �O  <      THREEJ%ABS    P  <      THREEJ%INT    XP  >      THREEJ%DSQRT    �P  @   a   THREEJ%J1    �P  @   a   THREEJ%J2    Q  @   a   THREEJ%J    VQ  n       BINOM    �Q  =      BINOM%DBLE    R  @   a   BINOM%N    AR  @   a   BINOM%R    �R  Y       CONST    �R  >      CONST%DSQRT    S  i       PRLOG    �S  ?      PRLOG%DLOG10    �S  @   a   PRLOG%X     T  w       PLEXP    wT  =      PLEXP%DABS    �T  =      PLEXP%DEXP    �T  @   a   PLEXP%X    1U  k       INV_PLEXP    �U  =      INV_PLEXP%DLOG    �U  @   a   INV_PLEXP%P    V  �       TBR_TX    �V  =      TBR_TX%DABS    �V  =      TBR_TX%DEXP    W  @   a   TBR_TX%TL    _W  @   a   TBR_TX%TX    �W  @   a   TBR_TX%TAUL    �W  �       TBR_I    hX  =      TBR_I%DABS    �X  =      TBR_I%DEXP    �X  @   a   TBR_I%NU    "Y  @   a   TBR_I%I    bY  @   a   TBR_I%TAUL    �Y  {       SIMPSON    Z  @   a   SIMPSON%N    ]Z  @   a   SIMPSON%N1    �Z  @   a   SIMPSON%N2    �Z  �   a   SIMPSON%X    �[  �   a   SIMPSON%Y !   E\  @   a   SIMPSON%INTEGRAL (   �\  �       INTERPOLATEEXTERNALFILE -   j]  =      INTERPOLATEEXTERNALFILE%SIZE 0   �]  @      INTERPOLATEEXTERNALFILE%ADJUSTL -   �]  =      INTERPOLATEEXTERNALFILE%TRIM 1   $^  L   a   INTERPOLATEEXTERNALFILE%FILENAME +   p^  �   a   INTERPOLATEEXTERNALFILE%WL ,   _  �   a   INTERPOLATEEXTERNALFILE%OUT -   �_  @   a   INTERPOLATEEXTERNALFILE%NORM .   �_  @   a   INTERPOLATEEXTERNALFILE%ERROR    8`  �       RAD_FILE !   �`  @      RAD_FILE%ADJUSTL    a  =      RAD_FILE%TRIM    Va  L   a   RAD_FILE%FNAME    �a  @   a   RAD_FILE%JBOL    �a  L   a   RAD_FILE%STR    .b  @   a   RAD_FILE%L    nb  @   a   RAD_FILE%ERROR    �b  ~       DUST_RAD    ,c  =      DUST_RAD%DEXP    ic  @   a   DUST_RAD%TD    �c  @   a   DUST_RAD%TAU_D    �c  L   a   DUST_RAD%STR    5d  @   a   DUST_RAD%L 