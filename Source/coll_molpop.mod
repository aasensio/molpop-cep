	  �  U   k820309    q          13.1        UR                                                                                                           
       coll_molpop.f90 COLL_MOLPOP                  �                                 
                                                          
                  @ @                                           %         @                                                  
       #RDINP%LEN    #RDINP%INDEX    #EQUAL    #IUNIT    #OUTUNIT 	                                                    LEN                                                  INDEX           
                                                                                                                                               	            %         @                                
                          #ERROR_MESSAGE%TRIM    #ERROR_MESSAGE%ADJUSTL    #OPT    #STR                                                     TRIM                                                  ADJUSTL                                                �                                                                      �                                  @ @                                   
              
      p          p 
           p 
                                                                                               @                                   �                                                                #         @                                                     #RDINPS2%LEN    #RDINPS2%TRIM    #RDINPS2%ADJUSTL    #RDINPS2%INDEX    #RDINPS2%CHAR    #RDINPS2%ICHAR    #EQUAL    #IUNIT    #STR    #LENGTH    #UCASE                                                     LEN                                                  TRIM                                                  ADJUSTL                                                  INDEX                                                  CHAR                                                  ICHAR                                                                                                                                                                               1                                                                                                               #         @                                                     #RDINPS%LEN     #RDINPS%TRIM !   #RDINPS%ADJUSTL "   #RDINPS%INDEX #   #EQUAL $   #IUNIT %   #STR &                                                     LEN                                             !     TRIM                                             "     ADJUSTL                                             #     INDEX                                            $                                                       %                                                       &                     1 +           @ @                               '     
       �             p          p 
           p 
                                            @ @                              (     
              
      p          p 
           p 
                                    @                                 )     
                    p          p 
           p 
                         %         @                               *                           #STRING +                                              +                     1            @ @                               ,     
                    p          p 
           p 
                         #         @                                  -                    #S1 .   #S2 /   #S 0                                              .                     1                                            /                     1                                            0                     1              @                               1     �       %         @                               2                    
       #V 3   #N 4                                             3                    
 	    p          5 O p            5 O p                                                                     4                                                         5                     @                                6                   
                &                   &                                                                                       7     
                                                  8                                   &                                           #         @                                  9                    #IUNIT :                                              :                                                        ;     
       %         @                               <                          #INMIN%LEN =   #STRING >                                               =     LEN                                            >                     1 %         @                               ?                          #INMAX%LEN @   #STRING A                                               @     LEN                                            A                     1 #         @                                   B                   #READ_COLLISIONS%INT C   #READ_COLLISIONS%ADJUSTL D   #READ_COLLISIONS%TRIM E   #IEQUAL F                                              C     INT                                            D     ADJUSTL                                            E     TRIM           D @                               F            #         @                                   G                     #         @                                  H                   #LOADCIJ_NORM_TABLE%DSQRT I   #FILENAME J   #FRACTION K   #I_EXTRAP L   #XSEC M                                              I     DSQRT           D @                               J     @                                                                  K     
                                                  L                                                       M     
       #         @                                  N                    #FILENAME O   #WEIGHT P   #XSEC Q             D @                               O     @                                                                  P     
                                                  Q     
       #         @                                  R                    #WEIGHT S   #XSEC T                                              S     
                                                  T     
          �   $      fn#fn    �   @   J   GLOBAL_MOLPOP      @   J   MATHS_MOLPOP $   D  @       N_COL+GLOBAL_MOLPOP #   �  �       RDINP+MATHS_MOLPOP '     <      RDINP%LEN+MATHS_MOLPOP )   S  >      RDINP%INDEX+MATHS_MOLPOP )   �  @   a   RDINP%EQUAL+MATHS_MOLPOP )   �  @   a   RDINP%IUNIT+MATHS_MOLPOP +     @   a   RDINP%OUTUNIT+MATHS_MOLPOP +   Q  �       ERROR_MESSAGE+MATHS_MOLPOP 0   �  =      ERROR_MESSAGE%TRIM+MATHS_MOLPOP 3   #  @      ERROR_MESSAGE%ADJUSTL+MATHS_MOLPOP /   c  P   a   ERROR_MESSAGE%OPT+MATHS_MOLPOP /   �  P   a   ERROR_MESSAGE%STR+MATHS_MOLPOP %     �       FR_COL+GLOBAL_MOLPOP $   �  @       KBETA+GLOBAL_MOLPOP 7   �  @       FILE_PHYSICAL_CONDITIONS+GLOBAL_MOLPOP $     @       I_SUM+GLOBAL_MOLPOP %   W  �       RDINPS2+MATHS_MOLPOP )   E  <      RDINPS2%LEN+MATHS_MOLPOP *   �  =      RDINPS2%TRIM+MATHS_MOLPOP -   �  @      RDINPS2%ADJUSTL+MATHS_MOLPOP +   �  >      RDINPS2%INDEX+MATHS_MOLPOP *   <  =      RDINPS2%CHAR+MATHS_MOLPOP +   y  >      RDINPS2%ICHAR+MATHS_MOLPOP +   �  @   a   RDINPS2%EQUAL+MATHS_MOLPOP +   �  @   a   RDINPS2%IUNIT+MATHS_MOLPOP )   7	  L   a   RDINPS2%STR+MATHS_MOLPOP ,   �	  @   a   RDINPS2%LENGTH+MATHS_MOLPOP +   �	  @   a   RDINPS2%UCASE+MATHS_MOLPOP $   
  �       RDINPS+MATHS_MOLPOP (   �
  <      RDINPS%LEN+MATHS_MOLPOP )   �
  =      RDINPS%TRIM+MATHS_MOLPOP ,   *  @      RDINPS%ADJUSTL+MATHS_MOLPOP *   j  >      RDINPS%INDEX+MATHS_MOLPOP *   �  @   a   RDINPS%EQUAL+MATHS_MOLPOP *   �  @   a   RDINPS%IUNIT+MATHS_MOLPOP (   (  L   a   RDINPS%STR+MATHS_MOLPOP %   t  �       FN_KIJ+GLOBAL_MOLPOP $     �       GXSEC+GLOBAL_MOLPOP (   �  �       I_COL_SRC+GLOBAL_MOLPOP #   8  \       FSIZE+MATHS_MOLPOP *   �  L   a   FSIZE%STRING+MATHS_MOLPOP *   �  �       I_COL_EXSUB+GLOBAL_MOLPOP %   t  _       ATTACH2+MATHS_MOLPOP (   �  L   a   ATTACH2%S1+MATHS_MOLPOP (     L   a   ATTACH2%S2+MATHS_MOLPOP '   k  L   a   ATTACH2%S+MATHS_MOLPOP ,   �  @       PATH_DATABASE+GLOBAL_MOLPOP "   �  ^       VSUM+MATHS_MOLPOP $   U  �   a   VSUM%V+MATHS_MOLPOP $   �  @   a   VSUM%N+MATHS_MOLPOP     9  @       N+GLOBAL_MOLPOP     y  �       C+GLOBAL_MOLPOP !     @       VT+GLOBAL_MOLPOP     ]  �       G+GLOBAL_MOLPOP )   �  S       PASS_HEADER+MATHS_MOLPOP /   <  @   a   PASS_HEADER%IUNIT+MATHS_MOLPOP     |  @       T+GLOBAL_MOLPOP #   �  k       INMIN+MATHS_MOLPOP '   '  <      INMIN%LEN+MATHS_MOLPOP *   c  L   a   INMIN%STRING+MATHS_MOLPOP #   �  k       INMAX+MATHS_MOLPOP '     <      INMAX%LEN+MATHS_MOLPOP *   V  L   a   INMAX%STRING+MATHS_MOLPOP     �  �       READ_COLLISIONS $   F  <      READ_COLLISIONS%INT (   �  @      READ_COLLISIONS%ADJUSTL %   �  =      READ_COLLISIONS%TRIM '   �  @   a   READ_COLLISIONS%IEQUAL    ?  H       LOADCIJ #   �  �       LOADCIJ_NORM_TABLE )   !  >      LOADCIJ_NORM_TABLE%DSQRT ,   _  P   a   LOADCIJ_NORM_TABLE%FILENAME ,   �  @   a   LOADCIJ_NORM_TABLE%FRACTION ,   �  @   a   LOADCIJ_NORM_TABLE%I_EXTRAP (   /  @   a   LOADCIJ_NORM_TABLE%XSEC $   o  l       LOADCIJ_SPEC_TABLE1 -   �  P   a   LOADCIJ_SPEC_TABLE1%FILENAME +   +  @   a   LOADCIJ_SPEC_TABLE1%WEIGHT )   k  @   a   LOADCIJ_SPEC_TABLE1%XSEC '   �  ^       HARD_SPHERE_COLLISIONS .   	  @   a   HARD_SPHERE_COLLISIONS%WEIGHT ,   I  @   a   HARD_SPHERE_COLLISIONS%XSEC 