	  �  H   k820309    k          13.1        ���Q                                                                                                           
       sol_cep.f90 SOLUTION_CEP                                                    
                @       �                                 
                                                          
                                                          
                        �                                 
                                                                           @                                                        @                                     
       #         @                                  	                    #INIT%TRIM 
   #INIT%ADJUSTL    #INIT%ALLOCATED                                                
     TRIM                                                  ADJUSTL                                                  ALLOCATED #         @                                                      #PRECALCULATE_BETA%ALLOCATED    #PRECALCULATE_BETA%LOG10                                                     ALLOCATED                                                  LOG10              @                                   P                                                                                                                                    @                                             #         @                                                      #X                                                                 
 
              &                                                    @ @                                                 
                &                                           #         @                                                                                                                #         @                                                     #ADAPT_COL_DENSITY_STRATEGY%LOG10    #ADAPT_COL_DENSITY_STRATEGY%ABS    #ADAPT_COL_DENSITY_STRATEGY%MAXVAL    #ADAPT_COL_DENSITY_STRATEGY%INT    #ADAPT_COL_DENSITY_STRATEGY%MINVAL    #ADAPT_COL_DENSITY_STRATEGY%SUM    #ADAPT_COL_DENSITY_STRATEGY%TRIM     #ADAPT_COL_DENSITY_STRATEGY%ADJUSTL !   #X "                                                    LOG10                                                  ABS                                                  MAXVAL                                                  INT                                                  MINVAL                                                  SUM                                                   TRIM                                             !     ADJUSTL                                           "                   
 	              &                                                                                     #                   
                &                                                                                     $                   
                &                                                                                     %                   
                &                                           #         @                                  &                   #CEP_SOLVER%ABS '   #CEP_SOLVER%MAXVAL (   #ERROR )                                               '     ABS                                             (     MAXVAL                                            )                       @                                 *                                                         +            #         @                                  ,                     #         @                                  -                   #WRITE_INTERMEDIATE_RESULTS%EXP .   #WRITE_INTERMEDIATE_RESULTS%SQRT /   #WRITE_INTERMEDIATE_RESULTS%SUM 0   #WRITE_INTERMEDIATE_RESULTS%LOG 1   #WRITE_INTERMEDIATE_RESULTS%MAXVAL 2   #X 3   #XLTE 4                                               .     EXP                                             /     SQRT                                            0     SUM                                             1     LOG                                             2     MAXVAL                                           3                   
               &                                                                                     4                   
               &                                           #         @                                  5                   #REGRID%TRIM 6   #REGRID%ADJUSTL 7   #REGRID%ALLOCATED 8   #RESET 9                                               6     TRIM                                             7     ADJUSTL                                             8     ALLOCATED                                            9                                                        :     
                @                                ;                   
                &                                           #         @                                   <                   #LLSLV%DABS =   #LLSLV%DSQRT >   #S ?   #X @   #N A   #NR B                                               =     DABS                                             >     DSQRT          
                               ?                    
       p        5 O p        p          5 O p          5 O p            5 O p          5 O p                                   
                               @                    
     p          5 O p            5 O p                                    
                                 A                     
                                 B           #         @                                   C                   #SOLVE_WITH_CEP%SUM D   #SOLVE_WITH_CEP%ADJUSTL E   #SOLVE_WITH_CEP%TRIM F   #ALGORITHM G                                              D     SUM                                            E     ADJUSTL                                            F     TRIM                                            G               �   !      fn#fn    �   @   J   CONSTANTS_CEP      @   J   GLOBAL_CEP    A  @   J   FUNCTIONS_CEP    �  @   J   ESCAPE_CEP    �  @   J   IO_CEP #     @       VERBOSE+GLOBAL_CEP (   A  @       WHICH_SCHEME+GLOBAL_CEP ,   �  @       FACTOR_ABUNDANCE+GLOBAL_CEP #   �  }       INIT+FUNCTIONS_CEP (   >  =      INIT%TRIM+FUNCTIONS_CEP +   {  @      INIT%ADJUSTL+FUNCTIONS_CEP -   �  B      INIT%ALLOCATED+FUNCTIONS_CEP -   �  �       PRECALCULATE_BETA+ESCAPE_CEP 7   �  B      PRECALCULATE_BETA%ALLOCATED+ESCAPE_CEP 3   �  >      PRECALCULATE_BETA%LOG10+ESCAPE_CEP '     @       OUTPUT_FILE+GLOBAL_CEP &   C  @       START_MODE+GLOBAL_CEP .   �  @       READ_PREVIOUS_FLAG+GLOBAL_CEP *   �  @       OPTICALLY_THIN+GLOBAL_CEP .     O       CALCJBAR_LSTAR_CEP+ESCAPE_CEP 0   R  �   a   CALCJBAR_LSTAR_CEP%X+ESCAPE_CEP    �  �       POP+GLOBAL_CEP 2   j  H       CORRECT_POPULATIONS+FUNCTIONS_CEP +   �  @       KTHICK_STRATEGY+GLOBAL_CEP 6   �  |      ADAPT_COL_DENSITY_STRATEGY+ESCAPE_CEP <   n	  >      ADAPT_COL_DENSITY_STRATEGY%LOG10+ESCAPE_CEP :   �	  <      ADAPT_COL_DENSITY_STRATEGY%ABS+ESCAPE_CEP =   �	  ?      ADAPT_COL_DENSITY_STRATEGY%MAXVAL+ESCAPE_CEP :   '
  <      ADAPT_COL_DENSITY_STRATEGY%INT+ESCAPE_CEP =   c
  ?      ADAPT_COL_DENSITY_STRATEGY%MINVAL+ESCAPE_CEP :   �
  <      ADAPT_COL_DENSITY_STRATEGY%SUM+ESCAPE_CEP ;   �
  =      ADAPT_COL_DENSITY_STRATEGY%TRIM+ESCAPE_CEP >     @      ADAPT_COL_DENSITY_STRATEGY%ADJUSTL+ESCAPE_CEP 8   [  �   a   ADAPT_COL_DENSITY_STRATEGY%X+ESCAPE_CEP    �  �       DZ+GLOBAL_CEP %   s  �       ABUNDANCE+GLOBAL_CEP    �  �       NH+GLOBAL_CEP &   �  ~       CEP_SOLVER+ESCAPE_CEP *   	  <      CEP_SOLVER%ABS+ESCAPE_CEP -   E  ?      CEP_SOLVER%MAXVAL+ESCAPE_CEP ,   �  @   a   CEP_SOLVER%ERROR+ESCAPE_CEP '   �  @       GLOBAL_ITER+GLOBAL_CEP 0     @       PRINTINGS_PER_DECADE+GLOBAL_CEP %   D  H       WRITE_RESULTS+IO_CEP 2   �        WRITE_INTERMEDIATE_RESULTS+IO_CEP 6   �  <      WRITE_INTERMEDIATE_RESULTS%EXP+IO_CEP 7   �  =      WRITE_INTERMEDIATE_RESULTS%SQRT+IO_CEP 6     <      WRITE_INTERMEDIATE_RESULTS%SUM+IO_CEP 6   R  <      WRITE_INTERMEDIATE_RESULTS%LOG+IO_CEP 9   �  ?      WRITE_INTERMEDIATE_RESULTS%MAXVAL+IO_CEP 4   �  �   a   WRITE_INTERMEDIATE_RESULTS%X+IO_CEP 7   Y  �   a   WRITE_INTERMEDIATE_RESULTS%XLTE+IO_CEP %   �  �       REGRID+FUNCTIONS_CEP *   s  =      REGRID%TRIM+FUNCTIONS_CEP -   �  @      REGRID%ADJUSTL+FUNCTIONS_CEP /   �  B      REGRID%ALLOCATED+FUNCTIONS_CEP +   2  @   a   REGRID%RESET+FUNCTIONS_CEP )   r  @       COL_THRESHOLD+GLOBAL_CEP "   �  �       POPOLD+GLOBAL_CEP     >  �       LLSLV+MATHS_CEP %   �  =      LLSLV%DABS+MATHS_CEP &     >      LLSLV%DSQRT+MATHS_CEP "   ?  �   a   LLSLV%S+MATHS_CEP "   ;  �   a   LLSLV%X+MATHS_CEP "   �  @   a   LLSLV%N+MATHS_CEP #     @   a   LLSLV%NR+MATHS_CEP    _  �       SOLVE_WITH_CEP #     <      SOLVE_WITH_CEP%SUM '   ?  @      SOLVE_WITH_CEP%ADJUSTL $     =      SOLVE_WITH_CEP%TRIM )   �  @   a   SOLVE_WITH_CEP%ALGORITHM 