module io_molpop
use global_molpop
use maths_molpop, only: Pass_Header, error_message, inmin, inmax, attach2, ordera, optdep
use cep_molpop_interface
implicit none
contains
		
	SUBROUTINE INPUT(error)
!     enters model information from input file
   use global_molpop
	use maths_molpop
	use coll_molpop
      
   integer i,j, npr, iunit, L, L2, L3
	integer k, npoints_local
	integer :: j_max, vib_max
	integer unit,unit2,mol_maxlev
	double precision aux,xm4, w,Tbb,tau_d, T_d, Lbol,dist, temp
	double precision mol_const, temp1, temp2, temp3, temp4
	character*128 str, Method, Option, OPT,header
	character*168 fn_dusty
	logical iequal, error, UCASE,NoCase, stat
	integer LMeth, LOpt, nradiat, ComingFrom, n_modified
	character*72 TableName, TableInfo
	data UCASE/.TRUE./, NoCase/.FALSE./
			
!      start reading with negative unit to clean residual info in rdinps
!      UCASE and NoCase signal whether to convert input strings
!      to UCASE to trap possible keyboard entry problems    
   	iequal = .true.
      iunit = -15
                  
!         read type of escape probability approach
      call rdinps2(iequal,15,str,L,UCASE)

!       IF(str(1:L) .eq. 'LVG') THEN
!       	kbeta = 0
!       ELSE IF(str(1:L) .eq. 'LVGPP') THEN
! 			kbeta = -1
!       ELSE IF(str(1:L) .eq. 'SLAB') THEN
!         	kbeta = 2
! 		ELSE IF(str(1:L) .eq. 'CEP-NEWTON') THEN
!         	kbeta = 3			
! 		ELSE IF(str(1:L) .eq. 'CEP-ALI') THEN
!         	kbeta = 4			
! 		ELSE IF(str(1:L) .eq. 'CEP-NAG') THEN
!         	kbeta = 5					
! 		ELSE
! 	   	OPT = 'escape probability'
! 	   	error = error_message(opt,str)
!       	return
!       END IF
      
      IF(str(1:L) .eq. 'LVG') THEN
      	kbeta = 0
      ELSE IF(str(1:L) .eq. 'LVGPP') THEN
			kbeta = -1
      ELSE IF(str(1:L) .eq. 'SLAB') THEN
        	kbeta = 2
		ELSE IF(str(1:L) .eq. 'CEP') THEN
        	kbeta = 3		
		ELSE
	   	OPT = 'escape probability'
	   	error = error_message(opt,str)
      	return
      END IF
            
!     Decide whether to output PLOTTING file
!       call rdinps2(iequal,15,str,L,UCASE)      
!       if(str(1:L) .eq. 'NO') then
!       	i_sum = 0
!       else if(str(1:L) .eq. 'YES') then
!         	i_sum = 1
!         	open(17, file=fn_sum, status='unknown')
!       else
!         	OPT = 'output PLOT'        
!         	error = error_message(opt,str)
!         	return
!       end if
      
      i_sum = 0
      
!     read what molecule is used
      call rdinps2(iequal,iunit,STR,L,Nocase)
! Input the root directory for the database
      call rdinps2(iequal,15,path_database,L2,Nocase)

      
!
! ********* ME: From now (Oct 27, 2005) on, all data files are in subdirectory DataBase 
!
!     Read molecular information from file: mol_list.dat
! 		inquire(file=trim(adjustl(path_database))//'/mol_list.dat',exist=stat)
!       if (.not.stat) then
!       	print *, 'Error when opening file ', trim(adjustl(path_database))//'/mol_list.dat'
!       	stop
!       endif
      
      mol_name = STR(1:L)
      
!       open(11, file=trim(adjustl(path_database))//'/mol_list.dat', status='old')      
! !     First read Header lines
!       call Pass_Header(11)      
!       do while(.not. (mol_name(1:size(mol_name)).eq.STR(1:L)))
!          k=0
!          read(11, *,iostat=k) mol_name, N_max, mol_mass
!          
!          if (k .eq. -1) then
! !          No match:
!          	OPT = 'molecule name'
!            	error = error_message(opt,str)
!         		return
!          end if
! !           read(11, *) N_max, mol_mass, mol_const, j_max
!       end do
! 5    	close(11)

!
1   	call clear_string(len(s_mol),s_mol)
      s_mol(1:L) = str(1:L)
!
! Now it's safe to attach the directory to the molecule name 
!
      call attach2(trim(adjustl(path_database))//'/',mol_name,str)
      mol_name = str(1: size(str))
      
!     the number of levels to be used
      n = rdinp(iequal,15)
                        
! Allocate all the arrays that depend on the number of included energy levels		
		allocate(tau(n,n))
		allocate(esc(n,n))
		allocate(dbdtau(n,n))
		allocate(a(n,n))
		allocate(tij(n,n))
		allocate(taux(n,n))
		allocate(c(n,n))
		allocate(rad(n,n))
		allocate(we(n))
		allocate(gap(n,n))
		allocate(ems(n,n))
		allocate(boltz(n))
		allocate(rad_internal(n,n))
		allocate(rad_tau0(n,n))
		allocate(rad_tauT(n,n))
		allocate(freq(n,n))
		allocate(fr(n))
		allocate(wl(n,n))
		allocate(ti(n))
		allocate(g(n))
		allocate(pop(n))
		allocate(coolev(n))
		allocate(xp(n))
		allocate(ledet(n))
		allocate(imaser(n))
		allocate(jmaser(n))
		
!     Is there a correction for finite number of rotation levels in ground vib stat?
		vib_max = rdinp(iequal,15)
		if (vib_max .gt. 0) then
	   	j_max     = rdinp(iequal,15)
	   	mol_const = rdinp(iequal,15)
		endif
	           
      
! Read collisional partners 
      call read_collisions(iequal)

      
! Read the filename with the spatial variation of the physical conditions
      call rdinps2(iequal,15,file_physical_conditions,L2,Nocase)

! If "none", then use the standard inputs given by nH2 and T
      if (trim(adjustl(file_physical_conditions)) /= 'none') then
			write(16,"(6x,'Using file with physical conditions : ', a)")&
				trim(adjustl(file_physical_conditions))
		else
			T = rdinp(iequal,15)			
	   	nh2 = rdinp(iequal,15)
	   	xmol = rdinp(iequal,15)
	   	v   = rdinp(iequal,15)
      	eps = rdinp(iequal,15)
		endif
		      		
		if (trim(adjustl(file_physical_conditions)) /= 'none' .and. kbeta < 3) then			
			WRITE (16,*)' *** Cannot use varying physical conditions without CEP'
	      OPT = 'varying physical conditions without CEP'
	      error = error_message(opt,trim(adjustl(file_physical_conditions)))
        	return
		endif
      
				   										   
	   if (int(T) > 10000) then
      	write(6, "(' *** T = ',F10.2,' is a bit much!')") T
	      goto 1
      end if
      
 !     Load molecular data:		
      call data(error)
      if (error) return
      
!     Test whether the desired number of levels is larger than the number of tabulated levels
      if (n .gt. N_max) then
      	WRITE (16,"(/,' *** Cannot specify N_levels = ', I3,/,&
                &' *** Maximum allowed for ',a,' is ',I3)") n, s_mol(1:size(s_mol)), N_max
			OPT = 'number of levels of'
			str = s_mol(1:size(s_mol))
			error = error_message(opt,str)
			return
      endif

      
!       temp  = rdinp(iequal,15)
      
      nmol = xmol*nh2
      xm4  = 1.0e4*(xmol/2.0)
      unit2 = 16
      if(i_sum .eq. 1) unit2 = 17
      DO unit = 16,unit2
         write(unit,"(6x,'Molecule --- ',a)")s_mol(1:L)
         write(unit,"(6x,'Number of levels = ',I2)")n
         if (trim(adjustl(file_physical_conditions)) == 'none') then
         	write(unit,'(6x,a,f6.1,a)') 'Tgas = ',T, ' K'
         	write(unit,'(6x,a,1pe9.2,a)') 'n_H2 = ',nh2, ' cm-3'
         	write(unit,'(6x,3a,1pe9.2)') 'n_', s_mol(1:size(s_mol)),'/n_H2 = ',xmol
         	write(unit,'(6x,3a,1pe9.2,a)') 'n_',s_mol(1:size(s_mol)),' = ', nmol, ' cm-3'
         endif
      end do
            
!   read overlap flag (overlaptest)
      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'ON') then
         overlaptest =.true.
         write(16,*)'     LINE OVERLAP USED'
         if(i_sum .eq. 1) write(17,*)'     LINE OVERLAP USED'
      else
      	overlaptest =.false.
         write(16,*)'     NO LINE OVERLAP USED' 
        	if(i_sum .eq. 1) write(17,*)'     NO LINE OVERLAP USED'
      end if
      
! Test that line overlap is not used with CEP
      if (overlaptest .and. kbeta >= 3) then
      	OPT = 'overlap and CEP'
        	error = error_message(opt,str)
        	return
      endif
      
      
      vt = 1.0d-5*dsqrt(2.0*bk*T/(mol_mass*xmp))
      
		if (trim(adjustl(file_physical_conditions)) == 'none') then
			write(16,"(6x,'Thermal speed = ',f5.1,' km/s')") vt
		endif

		
		IF(kbeta == 0 .or. kbeta == -1) THEN      	
         IF (VT.GT.V) THEN
		  		WRITE (16,"(/,' *** Cannot use LVG approximation with V =',F5.1,&
		  			&' km/s and thermal speed of',F5.1,' km/s')") v, vt
	      	OPT = 'escape probability input for'
	      	error = error_message(opt,str)
        		return
	   	END IF
         IF (overlaptest) THEN
		   	WRITE (16,*)' *** Cannot use LVG with Line Overlap'
	      	OPT = 'LVG with Line Overlap'
	      	error = error_message(opt,str)
        		return
	   	END IF
	   	
         WRITE(16,"(6x,'LVG escape probability.  V = ',F5.1,' km/s,  dlogV/dlogr= ', F7.2)") V, eps
        	if (i_sum .eq. 1) then
        		WRITE(17,"(6x,'LVG escape probability.  V = ',F5.1,  &
        			&' km/s,  dlogV/dlogr= ', F7.2)") V, eps 
        	endif
      ELSE IF(kbeta == 2) THEN        	
        	WRITE (16,"(6x,'Slab escape probability ')")
        	if (i_sum .eq. 1) WRITE (17,"(6x,'Slab escape probability ')")
        	IF (V.GT.VT) THEN
	     		WRITE (16,"(6X,'Using Entered Doppler Width = ',F6.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Entered Doppler Width = ',F6.1,' km/s')") V
	   	ELSE
           	V = VT
	     		WRITE (16,"(6X,'Using Thermal Doppler Width = ',F6.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	   	END IF
	   	write(16,"(6x,'FWHM = ',f5.1,' km/s')")V*1.665
	   	if (i_sum .eq. 1)write(17,"(6x,'Full Width half-max = ',f5.1,' km/s')")V*1.665
		ELSE IF(kbeta == 3) THEN
        	kbeta = 3			
        	WRITE (16,"(6x,'Coupled escape probability with Newton')")
        	if (i_sum .eq. 1) WRITE (17,"(6x,'Coupled escape probability with Newton')")
        	IF (V.GT.VT) THEN
	     		WRITE (16,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	   	ELSE
           	V = VT
	     		WRITE (16,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	   	END IF
	   	write(16,"(6x,'FWHM = ',f5.1,' km/s')")V*1.665
	   	if (i_sum .eq. 1)write(17,"(6x,'Full Width half-max = ',f5.1,' km/s')")V*1.665
		ELSE IF(kbeta == 4) THEN        		
        	WRITE (16,"(6x,'Coupled escape probability with ALI')")
        	if (i_sum .eq. 1) WRITE (17,"(6x,'Coupled escape probability with ALI')")
        	IF (V.GT.VT) THEN
	     		WRITE (16,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	   	ELSE
           	V = VT
	     		WRITE (16,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	   	END IF
	   	write(16,"(6x,'FWHM = ',f5.1,' km/s')")V*1.665
	   	if (i_sum .eq. 1)write(17,"(6x,'Full Width half-max = ',f5.1,' km/s')")V*1.665
		ELSE IF(kbeta == 5) THEN        		
        	WRITE (16,"(6x,'Coupled escape probability with NAG-NEWTON')")
        	if (i_sum .eq. 1) WRITE (17,"(6x,'Coupled escape probability with ALI')")
        	IF (V.GT.VT) THEN
	     		WRITE (16,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Entered Doppler Width = ',F5.1,' km/s')") V
	   	ELSE
	         V = VT
	     		WRITE (16,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	  			if (i_sum .eq. 1)WRITE (17,"(6X,'Using Thermal Doppler Width = ',F5.1,' km/s')") V
	   	END IF
	   	write(16,"(6x,'FWHM = ',f5.1,' km/s')")V*1.665
	   	if (i_sum .eq. 1)write(17,"(6x,'Full Width half-max = ',f5.1,' km/s')")V*1.665
		ELSE
	   	OPT = 'escape probability'
	   	error = error_message(opt,str)
      	return
      END IF
      

! Read the filename with the spatial variation of the physical conditions
      auxiliary_functions = 'KROLIK-MCKEE'
            	  				      		      
!     convert velocities to cm/sec for the internal working
      V  = V*1.D5
      VT = dsqrt(2.0*bk*T/(mol_mass*xmp))
!     Decide whether to include maser effects on populations
      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'NO') then
      	sat = 0
        	write(16,"(6x,'NO SATURATION EFFECTS')") 
        	if(i_sum .eq. 1)write(17,"(6x,'NO SATURATION EFFECTS')")
      else if(str(1:L) .eq. 'YES') then
        	sat = 1
        	write(16,"(6x,'INCLUDES MASER SATURATION')") 
        	if(i_sum .eq. 1)write(17,"(6x,'INCLUDES MASER SATURATION')")
      else
        	OPT = 'saturation parameter'
        	error = error_message(opt,str)
        	return
      end if
      
      if (kbeta >= 3 .and. sat == 1) then
      	sat = 0
      	print *, 'Maser saturation has been turned off in CEP'
      endif
      
!     read the maser aspect ratio
!       aspect  = rdinp(iequal,15)
      aspect = 1.d0
      if(aspect .gt.  1.0) then      
!     Find beaming angle from eq 3.26 Elitzur et al.,ApJ,367,333,1991
      	omega = 1.0/(11.*aspect**2)
        	write(16,'(6x,a,f5.0)')'Filamentary Masers; aspect ratio = ', aspect
      	if(i_sum .eq. 1)write(17,'(6x,a,f5.0)')'Filamentary Masers; Aspect ratio = ', aspect
      else if(aspect .eq. 1.) then
!     Use beaming angle of 4*pi
        	omega = 1.0
        	write(16,"(6x,'Aspect Ratio = 1')") 
        	if(i_sum .eq. 1)write(17,"(6x,'Aspect ratio = 1')")
      else
         OPT = 'Aspect less than one'
         error = error_message(opt,str)
        	return
      end if
      
! Test whether we want dust absorption
		dustAbsorption = rdinp(iequal,15)
                        
!     and load collision rates
		if (trim(adjustl(file_physical_conditions)) == 'none') then
! If everything is constant
      	call loadcij
      else
      
! If varying physical conditions
! First read the number of zones of the file with the physical conditions
			inquire(file=trim(adjustl(file_physical_conditions)),exist=stat)
      	if (.not.stat) then
      		print *, 'Error when opening file ',&
      			trim(adjustl(file_physical_conditions))
      		stop
      	endif
      	open(unit=45,file=trim(adjustl(file_physical_conditions)),&
      		action='read',status='old')
			call Pass_Header(45)
			n_zones_slab = rdinp(iequal,45)
			read(45,*)
			read(45,*)		
			read(45,*)
			if (allocated(collis_all)) deallocate(collis_all)
      	allocate(collis_all(n,n,n_zones_slab))
      	do i = 1, n_zones_slab

! Read the T and the relative abundance of each collider
				read(45,*) temp1, temp2, T, temp3, temp4, (fr_col(j),j=1,n_col)

! In case hard sphere collision rates are used, change the thermal velocity				
				vt = dsqrt(2.0*bk*T/(mol_mass*xmp))
				
				call loadcij
				collis_all(:,:,i) = c(1:n,1:n)
			enddo
			close(45)	
      	
      endif
      
!     Radiation Field;
		rad_tau0 = 0.d0
		rad_tauT = 0.d0
		rad_internal = 0.d0
		
      DO unit = 16,unit2
      	write(unit,'(/6x,a)') 'External Radiation field includes:'
      	write(unit,'(8x,a)') '- 3K Cosmic Background'
      end do
      do i = 2, n
        	do j = 1, i-1
          	rad(i,j)  = plexp(tij(i,j)/2.725d0)

! Include CMB radiation in both sides of the slab
          	rad_tau0(i,j)  = plexp(tij(i,j)/2.725d0)
          	rad_tauT(i,j)  = plexp(tij(i,j)/2.725d0)
        	end do
      end do
      
!     optionally add diluted black bodies with temperatures
!     Tbb and dilution coefficients W:
!	
		W = 1.
		do while (W.gt.0.)
      	W   = rdinp(iequal,15)
      	
         if(W .gt. 0) then
         	Tbb = rdinp(iequal,15)
         	call rdinps2(iequal,15,str,L,UCASE)
         	
		 		write(16,"(8x,'- Black-body with T =',f7.1,'K  Dilution factor = ',1pe9.2)") Tbb, W 
       		if (i_sum .eq. 1) then
       			write(17,"(8x,'- Black-body with T =',f7.1,'K  Dilution factor = ',1pe9.2)") Tbb, W
       		endif
       	
         	do i = 2, n
         		do j = 1, i-1
            		rad(i,j)  = rad(i,j) + W*plexp(tij(i,j)/Tbb)
            		
! Coming from the left
               	if (str(1:L) == 'LEFT') then
               		rad_tau0(i,j) = rad_tau0(i,j) + W*plexp(tij(i,j)/Tbb)
               	endif
               	
! Coming from the right
               	if (str(1:L) == 'RIGHT') then
	               	rad_tauT(i,j) = rad_tauT(i,j) + W*plexp(tij(i,j)/Tbb)
               	endif
               	
! Coming from both sides
               	if (str(1:L) == 'BOTH') then
               		rad_tau0(i,j) = rad_tau0(i,j) + W*plexp(tij(i,j)/Tbb)
	               	rad_tauT(i,j) = rad_tauT(i,j) + W*plexp(tij(i,j)/Tbb)
               	endif
               	
! Coming from both sides
               	if (str(1:L) == 'INTERNAL') then
               		rad_internal(i,j) = rad_internal(i,j) + W*plexp(tij(i,j)/Tbb)
               	endif
            	end do
         	end do
         end if         
      end do

!     optionally add dust radiation
		tau_d = 1.		
		do while (tau_d.gt.0.)
      	tau_d = rdinp(iequal,15)      	
         if(tau_d .gt. 0) then
         	T_d   = rdinp(iequal,15)
         	call rdinps2(iequal,15,str,L,UCASE)
           	write(16,"(8x,'- Dust with Td ='f7.1,' K and Tau(v) = ',1pe9.2)") T_d, tau_d
        		if (i_sum .eq. 1) then
        			write(17,"(8x,'- Dust with Td ='f7.1, ' K and Tau(v) = ',1pe9.2)") T_d, tau_d 
        		endif
           	call dust_rad(T_d,tau_d,str,L)
	   	end if 
		end do 
		
!     Radiation from DUSTY output file
      call rdinps2(iequal,15,fn_DUSTY,L,Nocase)
      if (fn_DUSTY(1:L) /= 'none') then
			Lbol = rdinp(iequal,iunit)
			dist = rdinp(iequal,iunit)
			call rdinps2(iequal,15,str,L2,UCASE)		  		
	  		call rad_file(fn_DUSTY(1:L),Lbol,dist,error,str,L2)
        	if (error) return
 		end if
 		

! Method of solution. If single-zone, we don't care because we always use Newton
! If CEP, we can use NEWTON or ALI
		call rdinps2(iequal,15,str,L,UCASE)
		if (kbeta == 3) then
			IF(str(1:L) .eq. 'NEWTON') THEN
         	kbeta = 3			
			ELSE IF(str(1:L) .eq. 'ALI') THEN
         	kbeta = 4
         endif
		endif

! ACC    - ACCURACY REQUIRED IN THE SOLUTION
      ACC     = rdinp(iequal,iunit)
      itmax   = rdinp(iequal,iunit)

!   Solution strategy: 
!   "increasing" -- start from optically thin solution for R = 0
!   by solving the linear rate equations. Then find R such that all
!   optical depths are smaller than the input TAUM. Solve for that based
!   on the linear solution. Increase R until it exceeds the input RM or the
!   column exceeds the input COLM, whichever comes first
!   "decreasing" -- start from thermal equilibrium level populations.
!   Then find R such that all 
!   optical depths are larger than the input TAUM. Solve for that based
!   on the thermal populations. Decrease R until it is smaller then the 
!   input RM or the column is less then the input COLM, whichever comes first.
!   "fixed" -- solve the fixed given problem and stop

! Precision in the CEP grid convergence	   
      call rdinps2(iequal,15,str,L,UCASE)            
      TAUM  = rdinp(iequal,iunit)
      NPR     = rdinp(iequal,iunit)
      NRMAX   = rdinp(iequal,iunit)
      NMAX    = rdinp(iequal,iunit)
      RM    = rdinp(iequal,iunit)
      COLM  = rdinp(iequal,iunit)
      if(str(1:L) .eq. 'INCREASING') then
      	KTHICK = 0
      	write(16,*) 'Using INCREASING strategy'
      else if(str(1:L) .eq. 'DECREASING') then
         KTHICK = 1
         write(16,*) 'Using DECREASING strategy'
      else if(str(1:L) .eq. 'FIXED') then      	
      	KTHICK = 2
! Since no physical dimensions are given in the case of constant physical conditions
! in needs to be done in the increasing/decreasing mode
      	if (trim(adjustl(file_physical_conditions)) == 'none') then
      		OPT = 'FIXED mode. It needs a file with physical conditions'
      		error = error_message(opt,str)
        		return
      	endif
      else
	    	OPT = 'solution strategy'
	    	error = error_message(opt,str)
        	return
      end if
      
      cep_precision = rdinp(iequal,15)

! Output value of mu
		mu_output = rdinp(iequal,15)

!
! NPR    - # OF STEPS PER DECADE FOR PRINTING
! NR     - # OF STEPS PER DECADE FOR INCREASING R; DEFINES THE STEP SIZE AND
!          STARTS EQUAL TO NPR
! NMAX   - MAXIMUM # OF STEPS TO REACH RM
! NRMAX  - MAXIMUM # OF STEPS PER DECADE; EQUIVALENT TO MINIMUM STEP
! ACC    - ACCURACY REQUIRED IN THE SOLUTION
!
!
      NR      = NPR
      NR0     = NPR
      STEP    = 10.**(1./NR)
		PRSTEP  = INT(STEP)
      IF (KTHICK.EQ.1) STEP = 1./STEP

!
!             print control
!
!     Printing molecular data
!
      do i=1,6
        	call rdinps2(iequal,15,str,L,UCASE)
        	if(str(1:L) .eq. 'ON') then
         	ipr_lev(i)=.true.
        	else
          	ipr_lev(i)=.false.
        	end if
      end do

      do i=1,6
      	call rdinps2(iequal,15,str,L,UCASE)       
        	if(str(1:L) .eq. 'ON') then
         	ipr_tran(i)=.true.
        	else
         	ipr_tran(i)=.false.
        	end if
      end do

      if(ipr_lev(1) .or. ipr_tran(1)) call print_mol_data

!     stop without running after printing molecular data
!
      call rdinps2(iequal,15,str,L,UCASE)       
      if(str(1:L) .eq. 'ON') then
!        That's it, so:
	   	write(16,"(/6x,'Done with printing molecular data')")
         error = .true.
         return
      end if

!     print messages
      call rdinps2(iequal,15,str,L,UCASE)       
      if (str(1:L) .eq. 'OFF') then
        	newtpr = 0
!      else if(str(1:L) .eq. 'TERSE') then
!        newtpr = 1
      else if(str(1:L) .eq. 'ON') then
        	newtpr = 2
      else
        	OPT = 'Newton messages'
        	error = error_message(opt,str)
        	return
      end if

      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'OFF') then
        	ksolpr = 0
      else if(str(1:L) .eq. 'ON') then
        	ksolpr = 1
      else
        	OPT = 'step-size messages'
        	error = error_message(opt,str)
        	return
      end if

      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'OFF') then
        	kfirst = 0
      else if(str(1:L) .eq. 'ON') then
        	kfirst = 2
      else
        	OPT = 'initial guess'
        	error = error_message(opt,str)
        	return
      end if


      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'OFF') then
        	kprt = 0
      else if(str(1:L) .eq. 'ON') then
        	kprt = 1
      else
        	OPT = 'information on each step'
        	error = error_message(opt,str)
        	return
      end if

      call rdinps2(iequal,15,str,L,UCASE)
      if(str(1:L) .eq. 'OFF') then
        	kprt = kprt + 0
      else if(str(1:L) .eq. 'ON') then
        	kprt = kprt + 2
      else
        	OPT = 'population printing'
        	error = error_message(opt,str)
        	return
      end if

      nbig    = rdinp(iequal,15)
      		
      allocate(final(12,nmax))

!     input parameters for selected transitions to output
!
      n_tr = rdinp(iequal,15)
                  
		if (n_tr .gt. 0) then
					
			allocate(itr(n_tr))
      	allocate(jtr(n_tr))
      	allocate(in_tr(n_tr))
      	allocate(f_tr(n_tr))
      	allocate(fin_tr(n_tr,12,nmax))
      	
      	do i=1, n_tr
         	itr(i) = rdinp(iequal,15)
            jtr(i) = rdinp(iequal,15)
         end do
		end if

! ANDRES: if n_tr=-1, output for all the lines	
		if (n_tr .eq. -1) then
			nradiat = 0
			do i = 1, n
				do j = 1, i-1
					if (a(i,j) .ne. 0.d0) then
						nradiat = nradiat + 1
					endif
				enddo
			enddo
			n_tr = nradiat
						
			allocate(itr(n_tr))
      	allocate(jtr(n_tr))
      	allocate(in_tr(n_tr))
      	allocate(f_tr(n_tr))
      	allocate(fin_tr(n_tr,12,1000))
      	
			nradiat = 0
			do i = 1, n
				do j = 1, i-1
					if (a(i,j) .ne. 0.d0) then
						nradiat = nradiat + 1
						itr(nradiat) = i
						jtr(nradiat) = j
					endif
				enddo
			enddo
		endif
					
      write(16,'(78(''-'')/)')

!     Finish up: calculate the partition function and renormalize density
!
      if(vib_max .gt. 0) then
       	write(16,*)'Correct for Finite Number of Rotational Levels'
       	write(16,*)'Using Method of Lockett & Elitzur,ApJ,399,704'
	 		write(16,123)nmol
123  		format(' Original Molecular Density = ',1pe9.2)
        	nmol=nmol*part_func(t,mol_const,j_max)
       	write(16,124)nmol
124  		format(' Recalculated Molecular Density = ',1pe9.2)
      end if
      
!    
!     SCALE RATES WITH WEIGHT FACTORS AND DEFINE THE REST OF CONSTANTS
!
		AUX    = (NMOL/V)*CL**3/EITPI
		
!  Change slab optical depth to line center (11 March 05)

! Write intermediate file in case a CEP calculations is being done		
		if (kbeta > 2) then
			call generate_intermediate_data
		endif
	
		if (kbeta .eq. 2 .or. kbeta .eq. 0 .or. kbeta .eq. -1) aux = aux/rootpi 
      DO I = 2,N
        	DO J = 1,I-1
        		A(I,J)      =  A(I,J)*WE(I) 
        		C(I,J)      =  C(I,J)*WE(I)*NH2
        		TAUX(I,J)   =  AUX*A(I,J)/FREQ(I,J)**3
        		EMS(I,J)    =  HPL*FREQ(I,J)*A(I,J)
        	END DO
      END DO

      return

	end SUBROUTINE INPUT
		
						
	subroutine print_mol_data
!     outputs molecular properties
   integer i,j,len1,len2,len3,len4,len5,len6
   character*80 str1,str2,str3,str4,str5,str6

      if(ipr_lev(1)) then
      	if(.not. ipr_lev(2) .and. .not. ipr_lev(3) .and..not. ipr_lev(4) .and. &
      		.not. ipr_lev(5) .and..not. ipr_lev(6)) ipr_lev(6)=.true.
        	write(16,'(31(''-''),'' Energy levels '',32(''-'')/)')
        	str1='   i'
        	len1=len('   i')
        	if(ipr_lev(2)) then
          	str2='     g'
          	len2=len('     g')
        	else
          	str2=' '
          	len2=len(' ')
        	end if
        	if(ipr_lev(3)) then
          	str3='    cm^{-1}  '
          	len3=len('    cm^{-1}  ')
        	else
          	str3=' '
          	len3=len(' ')
        	end if
        	if(ipr_lev(4)) then
          	str4='      GHz    '
          	len4=len('      GHz    ')
        	else
          	str4=' '
          	len4=len(' ')
        	end if
        	if(ipr_lev(5)) then
          	str5='       K     '
          	len5=len('       K     ')
        	else
          	str5=' '
          	len5=len(' ')
        	end if
        	if(ipr_lev(6)) then
          	str6='   quantum numbers'
          	len6=len('   quantum numbers')
        	else
          	str6=' '
          	len6=len(' ')
        	end if

        	write(16,'(a/)')str1(1:len1) // str2(1:len2) // str3(1:len3) // str4(1:len4) //&
        		str5(1:len5) // str6(1:len6)
        
        	do i=1,n
          	write(str1,'(i4)') i
          	len1=4
          	if(ipr_lev(2)) then
            	write(str2,'(i6)') g(i)
            	len2=6
          	else
            	str2=' '
            	len2=len(' ')
          	end if
          	if(ipr_lev(3)) then
            	write(str3,'(1pe13.5)') fr(i)
            	len3=13
          	else
	            str3=' '
            	len3=len(' ')
          	end if
          	if(ipr_lev(4)) then
	            write(str4,'(1pe13.5)') 1.0e-9*cl*fr(i)
            	len4=13
          	else
	            str4=' '
            	len4=len(' ')
          	end if
          	if(ipr_lev(5)) then
	            write(str5,'(1pe13.5)') ti(i)
            	len5=13
          	else
	            str5=' '
            	len5=len(' ')
          	end if
          	if(ipr_lev(6)) then
	            write(str6,'(3x,a)') ledet(i)(inmin(ledet(i)):inmax(ledet(i)))
            	len6=3+len(ledet(i)(inmin(ledet(i)):inmax(ledet(i))))
          	else
	            str6=' '
            	len6=len(' ')
          	end if
          	write(16,'(a)') str1(1:len1) // str2(1:len2) // str3(1:len3) // str4(1:len4) //&
          		str5(1:len5) // str6(1:len6)
        	end do
        	write(16,'(78(''-'')/)')
     	end if
!
      if(ipr_tran(1)) then
      	if(.not. ipr_tran(2) .and. .not. ipr_tran(3) .and. .not. ipr_tran(4) .and. &
      		.not. ipr_tran(5) .and..not. ipr_tran(6)) ipr_tran(5)=.true.
        	write(16,'(33(''-''),'' Transitions '',32(''-'')/)')
        	str1='  i -> j'
        	len1=len('  i -> j')
        	if(ipr_tran(2)) then
          	str2='     micron   '
          	len2=len('     micron   ')
        	else
          	str2=' '
          	len2=len(' ')
        	end if
        	if(ipr_tran(3)) then
          	str3='       GHz    '
          	len3=len('       GHz    ')
        	else
          	str3=' '
          	len3=len(' ')
        	end if
        	if(ipr_tran(4)) then
          	str4='        K     '
          	len4=len('        K     ')
        	else
          	str4=' '
          	len4=len(' ')
        	end if
        	if(ipr_tran(5)) then
          	str5='   Aij, s^{-1}'
          	len5=len('   Aij, s^{-1}')
        	else
          	str5=' '
          	len5=len(' ')
        	end if
        	if(ipr_tran(6)) then
          	str6='   Cij, s^{-1}'
          	len6=len('   Cij, s^{-1}')
        	else
          	str6=' '
          	len6=len(' ')
        	end if

        	write(16,'(a/)')str1(1:len1) // str2(1:len2) // str3(1:len3) // str4(1:len4) //&
        		str5(1:len5) // str6(1:len6)

        	do i = 1, n
          	do j = 1, i-1
!            if(ipr_tran(5) .and. (a(i,j) .ne. 0) .or.
!     *          ipr_tran(6) .and. (c(i,j) .ne. 0)) then
            	if(.true.) then
              		write(str1,'(i3,i5)') i,j
              		len1=8
              		if(ipr_tran(2)) then
                		write(str2,'(1pe14.5)') wl(i,j)
                		len2=14
              		else
                		str2=' '
                		len2=len(' ')
              		end if
              		if(ipr_tran(3)) then
                		write(str3,'(1pe14.5)') 1.0e-9*freq(i,j)
                		len3=14
              		else
                		str3=' '
                		len3=len(' ')
              		end if
              		if(ipr_tran(4)) then
                		write(str4,'(1pe14.5)') tij(i,j)
                		len4=14
              		else
                		str4=' '
                		len4=len(' ')
              		end if
              		if(ipr_tran(5)) then
                		write(str5,'(1pe14.5)') a(i,j)
                		len5=14
              		else
                		str5=' '
                		len5=len(' ')
              		end if
              		if(ipr_tran(6)) then
                		write(str6,'(1pe14.5)') c(i,j)*nh2
                		len6=14
              		else
                		str6=' '
                		len6=len(' ')
              		end if
              		write(16,'(a)') str1(1:len1) // str2(1:len2) // str3(1:len3) // &
              			str4(1:len4) // str5(1:len5) // str6(1:len6)
            	end if
          	end do
        end do
        write(16,'(78(''-'')/)')
		end if
      return
	end subroutine print_mol_data
		
		
	subroutine data(error)
!     Get molecular data: level properties and Einstein A-coefficients
   character*80 line
	character*64 fn_lev,fn_aij
	integer i,j,in,k,l,ii,i1,j1,i2,j2
	double precision aux, temp
	logical error, stat

!      include 'common__data.inc'
!     molecular data
!
!     g     - statistical weights;
!     a     - Einstein coefficients;
!     c     - collision rates;
!     we    - weight factors
!
!             a and c are scaled by we; note -
!                 n(i) a(i,j) = nmol x(i) (we(i) a(i,j))
!
!     fr    -  state energy (in cm**-1)
!     ti    -  state energy (in K)
!     freq  -  transition frequencies in Hz
!     wl    -  wavelength in microns
!     tij   -  transition frequencies in K
!     boltz -  Boltzmann factor for state
!     gap   -  Boltzmann factor for transition
!     taux  -  multipliers for optical depths: tau = r*taux*delta(x)
!     ems   -  multipliers for line emissivities
!     rad   -  external radiation intensity in the lines
!
!     inprt - printing index:
!
!             0 - no print;
!             1 - print level assignments;
!             2 - print also einstein a and collision rates;
!            <0 - terminate execution after printing

		do i = 1,n
      	g(i)  = 0
        	fr(i) = 0.
        	ti(i) = 0.
        	do j = 1,n
          	freq(i,j) = 0.
          	wl(i,j)   = 0.
          	tij(i,j)  = 0.
          	a(i,j)    = 0.
          	c(i,j)    = 0.
          	ems(i,j)  = 0.
          	taux(i,j) = 0.
          	rad(i,j)  = 0.
          	rad_tau0(i,j)  = 0.
          	rad_tauT(i,j)  = 0.
        	end do
      end do
!
!         read atomic data (energy levels and Einstein coefficients)
!
      call attach2(mol_name, '.molecule', fn_lev)
      inquire(file=fn_lev,exist=stat)
      if (.not.stat) then
      	print *, 'Error when opening file ', fn_lev
      	stop
      endif
		open(4, err=800, file=fn_lev, status='unknown')
		call Pass_Header(4)
		read(4,*) N_max, mol_mass
		call Pass_Header(4)
		do i=1,n
			read(4,*,end=5) in,g(i),fr(i),ledet(i)
		end do

      do i=1,n
        	if(i .gt. 1) then
          	do j = 1,i-1
	            freq(i,j) = cl*(fr(i) - fr(j))
            	wl(i,j)=1.0e4*cl/freq(i,j)
          	end do
        	end if
      end do

      !     Get past Header lines
		read(4,*)
      call Pass_Header(4)
      do while(.true.)
      	read(4,*,end=5)  i,j,temp
      	
! Only read the transitions between levels included in the model
      	if (i <= n .and. j <= n) then
      		a(i,j) = temp
      	endif
      end do
5     close(4)		
!
!                    define weight and Boltzmann factors:
!
      aux    = hpl/bk
      do i = 1,n
        	ti(i)    = aux*cl*fr(i)
        	BOLTZ(I) = dexp(-TI(I)/T)
        	WE(I)    = G(I)*BOLTZ(I)
        	do j = 1,i
          	tij(i,j) =  aux*freq(i,j)
          	GAP(I,J) =  DEXP(TIJ(I,J)/T)
        	end do
      end do

      return

800 	write(6,'(6x,3a)') 'File ',fn_lev(inmin(fn_lev):inmax(fn_lev)),' is missing! Stop.'
      error = .true.
      return

810 	write(6,'(6x,3a)') 'File ',fn_aij(inmin(fn_aij):inmax(fn_aij)),' is missing! Stop.'
      error = .true.
      return

	end subroutine data
		
		
   subroutine finish(ier)
!     Output summary of the run   
   integer ier,i,j,k,unit,unit2,n1,n2,dn
   
   	if(ier .eq. -1) write(16,"(/6x,'Terminated. Step size was ',1pe10.3/)") step
      if(ier .eq. 0) write(16,"(/6x,'Terminated. Number of steps reached',I5/)") Nmax
      if(ier .eq. 1) write(16,"(/6x,'Normal completion. Dimension is ',1pe10.2,' cm'/)") r
      if(ier .eq. 2) write(16,"(/6x,'Normal completion. H2 column is',1pe10.2,' cm-2')") Hcol

      IF (KTHICK.EQ.0) THEN
	   	n1 = 1
	   	n2 = nprint
	   	dn = 1
      ELSE
!     print in revesrse order when KTHICK = 1
	   	n1 = nprint
	   	n2 = 1
	   	dn = -1
      END IF
      
!     Attach original Summary to Output file
		unit = 16
      write(unit,"(/T35,'SUMMARY')")
      write(unit,FMT='(/6x,A,8x,A,3x,A,5x,A,3X,A,/6x,A,10x,A,8x,A,5x,A)') &
      	&'R','H2 column','mol column','emission',&
	   	&'xi','cm','cm-2','cm-2','erg/s/mol'
      do i = n1, n2, dn
      	write(unit,'(6(1pe12.3))') (final(k,i), k = 1,5)
      end do
      
      if(n_tr .eq. 0) return
      	write (unit,"(/20x,'Parameters for selected transitions')")
         do j = 1, n_tr
         	if (1.0e-4*wl(itr(j),jtr(j)) .ge. 1) then
            	write(unit,'(/5x,f8.3,'' cm ('',f7.3,'' GHz) transition between levels '',&
            		&i2,'' and '', i2)') &
						1.0e-4*wl(itr(j),jtr(j)),1.e-9*freq(itr(j),jtr(j)),itr(j), jtr(j)
            else
            	write(unit,'(/5x,f8.2,'' mic ('',f7.3,'' GHz) transition between levels '',&
            		&i2,'' and '', i2)') &
						wl(itr(j),jtr(j)),1.e-9*freq(itr(j),jtr(j)),itr(j), jtr(j)
            end if
            write(unit,"(5x,'Upper Level: ',a)")ledet(itr(j))
            write(unit,"(5x,'Lower Level: ',a)")ledet(jtr(j))
 	       	write(unit,"(/5x,'R',7x,'Tex',6x,'tau',5x,'Flux',5x,'eta',7x,'p1',7x,'p2',5x,&
 	       		&'Gamma1',3x,'Gamma2',4x,'Gamma'/,5x,'cm',7x,'K',16x,'Jy',13x,'cm-3s-1',&
 	       		&2x,'cm-3s-1',3(5x,'s-1 '))") 
            do i = n1, n2, dn
            	write(unit,'(11(1pe9.2))') (fin_tr(j,k,i), k = 1,10)
            end do
	    	end do
	    	
!      Now Output Plotting Summary
        	if(i_sum .eq. 1) then 
        		do j = 1, n_tr
            	write(17,'(60(''*''))')
             	if (1.0e-4*wl(itr(j),jtr(j)) .ge. 1) then
               	write(17,'(''*'',f8.3,'' cm ('',f7.3,'' GHz) transition between levels '',i2,&
               		&'' and '', i2)') &
							1.0e-4*wl(itr(j),jtr(j)),1.e-9*freq(itr(j),jtr(j)),itr(j), jtr(j)
             	else
               	write(17,'(''*'',f8.2,'' mic ('',f7.3,'' GHz) transition between levels '',i2,&
               		&'' and '', i2)') &
							wl(itr(j),jtr(j)),1.e-9*freq(itr(j),jtr(j)),itr(j), jtr(j)
             	end if
             	write(17,"('*',5x,'Upper Level',a33)")ledet(itr(j))
             	write(17,"('*',5x,'Lower Level',a33)")ledet(jtr(j))    
        			write(17,"(7x,'R',8x,'nmol*R/V',8x,'Tex',8x,'Tau')")
       			do i = n1, n2, dn
        				write(17,'(4(1pe12.2))')fin_tr(j,1,i),fin_tr(j,1,i)*nmol/V*1.e5,fin_tr(j,2,i),&
        					fin_tr(j,3,i)
       			end do      
       		end do 
       	end if    
		return
	end subroutine finish
		
		
	subroutine output(x,cool,kpr)
!              kpr  = 0 - print only summary
!                     1 - print each step; info on the nbig strongest
!                         cooling lines
!                     2 - also level populations
!
	integer kpr,i,j,k,kool
	integer, allocatable :: index(:),jndex(:)
	double precision x(n),bright,cool(:,:)
	double precision arg,rprev,test,tex,eta

		kool(arg) = 100.0*arg/tcool + 0.5
				
		allocate(index(nbig+1))
		allocate(jndex(nbig+1))		
		      
!     First check whether R changed enough for printing:
      if(nprint .gt. 1) then
        	test = r/rprev
        	IF (KTHICK.EQ.1) TEST = 1./TEST        	
!         	if(test .lt. prstep) return
      end if
      
!     OK, this is a printing step. Prepare output in LINES
!     and proceed with printing if detailed printout is requested by kpr
      rprev = r
      call lines(x,cool)
      if(kpr .eq. 0) return

!     kpr > 0;  more detailed printing
      write(16,FMT='(/2X,A,1pe9.2,A,5X,A,1pe9.2,A,5X,A,1pe9.2,A,/2X,A,1pe9.2,A)') &
      	&'R = ',R,'cm', 'H2 column = ',HCOL,'cm-2','mol column = ',MCOL,&
      	&' cm-2','Total emission = ',TCOOL,' erg/s/mol '
			
!     print populations, if required:
      if (kpr .ge. 2) call printx(x)
      if (kpr .eq. 1 .or. kpr .eq. 3) then
        	if(nbig .gt. 0) then
!         find the main cooling lines and print:
         	call ordera(cool,n,index,jndex,nbig)
          	write(16,'(2x,''The top '',i2,'' emitting lines are:''/)') nbig
          	write(16,"(8x,'i',4x,'j',2x,'wavelength',4x,'tau',4x,'beta',7x,'Tex',6x,&
          		&'emission',3x,'%')")
          	write(16,'(18x,''micron'',25x,''K'',7x,''erg/s/mol'')')
          	do k = 1, nbig
            	i = index(k)
            	j = jndex(k)
            	tex  = tij(i,j)/dlog(pop(j)/pop(i))
            	if(kool(cool(i,j)) .gt. 0) then
            		write(16,'(6x,i3,2x,i3,1pe11.3,4(1pe10.2),i4)')i,j,wl(i,j),tau(i,j),&
							esc(i,j),tex,cool(i,j),kool(cool(i,j))
					endif
          	end do
        	end if
!       nmaser - # of inverted transitions;  
!                imaser and jmaser are the (i,j) of these lines.
        	if (nmaser .gt. 0) then
!         print basic maser output for inverted transitions:
         	write(16,'(/2x,''Inverted lines:''/)')
          	write(16,'(8x,''i'',4x,''j'',2x,''wavelength'',4x,''tau'',5x,''Tex'',8x,''eta'')')
          	write(16,"(18x,'micron',15x,'K',/)")
          	do k = 1, nmaser
            	i = imaser(k)
            	j = jmaser(k)
            	TEX  = TIJ(I,J)/DLOG(POP(J)/POP(I))
            	eta = (pop(i) - pop(j))/(pop(i) + pop(j))
!            IF (KBETA.EQ.0) TAU(I,J) = TAU(I,J)/EPS
            	WRITE (16,"(6X,I3,2X,I3,1PE11.3,6(1PE10.2))")I,J,WL(I,J),TAU(I,J),TEX,eta
          	end do
        	else
          	write(16,'(/2x,''No inversions'')')
        	end if
      end if
      write(16,'(78(''-'')/)')
      
      deallocate(index)
		deallocate(jndex)
		
      return
	end subroutine output


	subroutine printx(x)
!     print populations:
	integer kool,i
	double precision x(:),tex,xni,bi

		write(16,'(2x,a/,6x,a,4x,a,9x,a,7x,a,5x,a,5x,a,5x,a/,26x,a,22x,a,9x,a)')&
			&'Populations:','level','x(i)',&
			&'n(i)','n(i)/n*','T_ex(i,1)','emission','%','cm^{-3}','K','erg/s/mol'
      do i=1,n
        	if(i .eq. 1) then
          	tex = 0.0
        	else
          	tex  = tij(i,1)/dlog(pop(1)/pop(i))
        	end if
        	xni = nmol*x(i)*we(i)
        	bi  = x(i)*sumw
        	if(tcool .ne. 0.0) then
          	kool = 100.0*coolev(i)/tcool +0.5
        	else
          	kool=0
        	end if
        	write(16,'(i9,5(1pe13.3),i5)') i,x(i),xni,bi,tex,coolev(i),kool
      end do
      write(16,'()')
      return
	end subroutine printx

 
   subroutine lines(x,cool)
!     Prepare all line output
!     Calculate line emissivities, optical depths and cooling.
!     cool(i,j) is the cooling of the (i,j) transition; coolev(i) the
!     overall cooling of the i-th level;
!      
!  	use sol_molpop, only: optdep
   integer m(2),i,j,k,i_m
   double precision x(n),cool(:,:),taumaser,pop1,pop2,eta
   double precision Gamma_av,Gamma(2),p1,p2,q,aux,flux,Tex

      call optdep(x)
      nmaser  = 0
      tcool   = 0.0
!
      do i = 1, n
        	coolev(i) = 0.0
        	do j = 1, i
          	cool(i,j) = 0.0
          	if(a(i,j) .ne. 0.0) then
            	cool(i,j) = ems(i,j)*esc(i,j)*x(i)
            	if(tau(i,j) .lt. 0.0) then
!             this is an inverted transition:
              		nmaser = nmaser + 1
              		imaser(nmaser) = i
              		jmaser(nmaser) = j
!             forget the emission from maser transitions
	        			cool(i,j) = 0.
            	end if
            	coolev(i) = coolev(i) + cool(i,j)
            	tcool = tcool + cool(i,j)
          	end if
        	end do
      end do
		if (R.eq.0.0) return

!     for final printing
      nprint = nprint + 1
      hcol = nh2*r
      mcol = nmol*r
      AUX = 1.D23*MCOL*CL/V
!     scaling parameter xi, as defined in EHM (xmol is normalized
!     with H-nuclei number)
      xi  = (1.0e-18*(2*nh2)**2)*(xmol/2.0e-4)*(r/1.0e13)*(1.0e5/v)
      final( 1,nprint) = r
      final( 2,nprint) = hcol
      final( 3,nprint) = mcol
      final( 4,nprint) = tcool
      final( 5,nprint) = xi
      if(n_tr .eq. 0) return
 
!  Printing elements for every selected transition:
!     2 - Excitation temperature
!     3 - tau
!     4 - flux density in Jy, obtained from COOL(i,j) which is in erg/s/mol
!  For masers only:
!     5 - inversion efficiency
!     6 - pump rate of lower level; obtained from loss through p_i = n_i*Gamma_i
!     7 - pump rate of upper level; obtained from loss through p_i = n_i*Gamma_i
!     8 - loss rate of lower level
!     9 - loss rate of upper level
!    10 - mean loss rate from average with statistical weights
!    11 - q: mean pump rate coefficient per sub-level
!         Follow definition by CFM, eq A11 of May 29, 02, Appendix

      do k = 1, n_tr
      	m(2)= itr(k)
         m(1)= jtr(k)
         Tex = TIJ(m(2),m(1))/DLOG(POP(m(1))/POP(m(2)))
         fin_tr(k, 1,nprint) = r
         fin_tr(k, 2,nprint) = Tex
         fin_tr(k, 3,nprint) = tau(m(2),m(1))
         fin_tr(k, 4,nprint) = aux*cool(m(2),m(1))/freq(m(2),m(1))
         If (tau(m(2),m(1)) .lt. 0.0) then
            eta = (pop(m(2)) - pop(m(1)))/(pop(m(2)) + pop(m(1)))
            call loss(m,Gamma)
            Gamma_av = (Gamma(1)*g(m(1)) + Gamma(2)*g(m(2)))/(g(m(1))+g(m(2)))
            q  = (pop(m(2)) + pop(m(1)))*Gamma_av/(4.*NH2)
            p2 = nmol*x(m(2))*we(m(2))*Gamma(2)
            p1 = nmol*x(m(1))*we(m(1))*Gamma(1)
            fin_tr(k, 5,nprint) = eta
            fin_tr(k, 6,nprint) = p1
            fin_tr(k, 7,nprint) = p2
            fin_tr(k, 8,nprint) = Gamma(1)
            fin_tr(k, 9,nprint) = Gamma(2)
            fin_tr(k,10,nprint) = Gamma_av
            fin_tr(k,11,nprint) = q
	   	else
	      	do i = 5,11
               fin_tr(k,i,nprint) = 0.0
	      	end do
         end if
      end do

		return
	end subroutine lines



   subroutine loss(m,Gamma)
!     calculate the loss rate of each maser level using eq. 7.1.1 of
!	Astronomical Masers   
   integer m(2),i,j,l
   double precision Gamma(2)
!      include 'common__loss.inc'

		DO i = 1, 2
      	l = m(i)
         Gamma(i) = 0.
         do j = 1, m(1) - 1
         	Gamma(i) = Gamma(i) + (A(l,j)*ESC(l,j)*(1 + RAD(l,j)) + C(l,j))/WE(l)
         end do   
         do j = m(2) + 1, N
            Gamma(i) = Gamma(i) + (A(j,l)*ESC(j,l)*RAD(j,l) + C(j,l)/GAP(j,l))*G(j)/(G(l)*WE(j))
         end do   
      End DO
      return
   end subroutine loss



   double precision function bright(flux,i,j)
!     calculates brightness temperature
!     enter with flux density in ERG/CM**2/SEC/HZ   
   integer i,j
   double precision intensity,flux,arg,lambda

!     DIVIDE BY SOLID ANGLE OMEGA TO GET THE INTENSITY; FORGET FILAMENTS:
!     IF (TAU(I,J).LT.0.) FLUX = FLUX/OMEGA
		intensity =	FLUX/fourpi
		lambda = 1.d4*wl(i,j)
      BRIGHT = intensity*lambda**2/(2.*BK)
      IF (BRIGHT.GE.TIJ(I,J) .OR. BRIGHT.EQ.0.) RETURN
!     Can't use the Rayleigh Jeans limit, so use full expression:
      ARG    = 2.*HPL*CL/(INTENSITY*LAMBDA**3)
		BRIGHT = TIJ(I,J)/DLOG(1. + ARG)
      RETURN
   END function bright
end module io_molpop
