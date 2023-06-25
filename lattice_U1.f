      program U1
      
C==========================================================================
C Preamble with useful definitions
C==========================================================================
      
      include "par.f"                           ! file with all parameters
      
      character(len=30) file_name, formatstring ! for output file
      character :: cr                           ! for percentage loading
      cr = char(13)    
      
	call cpu_time(start)
	call ranstart
	
	formatstring = "(A13,i2,A4)"
      write(file_name, formatstring) "raw_data/dati",NX,".dat"
      
	open(2, file=file_name, status='unknown')	! file for results
      
      write(2,*) measures  ! useful information for the analysis
      write(2,*) beta
      write(2,*) NX, NT

C==========================================================================
C Start simulation
C==========================================================================
	
	call init()   
      call boundary()
      
      print*, "Loop for tremalizzation..."
      do k = 1, i_term
          call metropolis()      
      enddo
	
	print*, "Start of measurament..."
	do i = 1, measures
	    write(*, '(A1)', advance='no')  cr !I write cr needed to clean the shell
	    
          do j = 1, i_dec
              call metropolis()
              call over_relaxation()
              call over_relaxation()
              call over_relaxation()
              call over_relaxation()
          enddo
          
          call charge(q)
          write(2,*) q
          
          !percentage loading
          perc = i/float(measures)*100
          write(*,'(f8.1, a2)', advance='NO') perc, " %"
          flush(6)
      enddo

C==========================================================================
C End of simulation
C==========================================================================	
	
	call ranfinish
	
	write(*,*)
	call cpu_time(finish)
	time = finish-start
	
	i_hr = time/3600          ! hour
	imin = time/60            ! minutes
	secd = mod(time, 60.0)    ! seconds
	
	write(*,'(i2, a3, i2, a3, f8.1, a3)') 
     & i_hr," h ", imin, " m ", secd, " s "
     
	
	end program
	
c=========================================================================
C Initial condition : Random configuration
c=========================================================================

      subroutine init()

      include "par.f"         ! file with all parameters
      
      do j = 1, NX            ! loop over lattice
            do i = 1, NT
              lv(i,j) = 1.0 - 2.*ran2() 
              lh(i,j) = 1.0 - 2.*ran2()
          enddo
      enddo
        
      return
      end

c=============================================================================
C Periodic boundary
c=============================================================================

      subroutine boundary()

      include "par.f"             ! file with all parameters
      
      integer, dimension(2) :: d
      d = (/ NX, NT /)	          ! array of spatial and temporal extents
      
      do j = 1, 2                 ! loop over number of dimensions
          do i = 1, d(j)          ! loop over lattice
          
              nl(i, j) = i + 1    ! nearest neighbors
              ns(i, j) = i - 1
              
          enddo
          nl(d(j), j) = 1	    ! fix boundary with
          ns(1, j) = d(j)	    ! periodic conditions
      enddo
      
      return
      end

c============================================================================
C METROPOLIS
c============================================================================

      subroutine metropolis()
       
      include "par.f"    ! file with all parameters
      
      complex :: c
      c = (0.0, 1.)      ! immaginary units
      
      pi = 4 * atan(1.0) ! pi greco
      
      do i_t = 1, NT     ! loop over lattice
          do i_s = 1, NX
            
C****************************** VERTICAL LINK  *******************************

              c_l = lv(i_t, i_s)                     ! value of current link
              t_l = c_l + (ran2()*2. - 1.)*delta     ! new trial value of link
              
              f1 = - lh(i_t, i_s)
     &             - lv(i_t, nl(i_s, 1)) + lh(nl(i_t, 2), i_s)      ! first staple

              f2 = lh(i_t, ns(i_s, 1)) -
     &             lv(i_t, ns(i_s, 1)) - lh(nl(i_t, 2), ns(i_s, 1)) ! second staple
              
              F_c = real(exp(c*(f1 + c_l)) + exp(c*(f2 + c_l)))   ! current action
              F_p = real(exp(c*(f1 + t_l)) + exp(c*(f2 + t_l)))   ! trial action
              r = beta*(F_p - F_c)
              prob = log(ran2())
                
              if (prob < r) then                       ! acceptance test
                    
                  t_l = mod(t_l, 2.0*pi)               ! the value must be
                  if (t_l < (-pi)) t_l = t_l + 2.0*pi  ! between -pi and pi
                  if (t_l >   pi)  t_l = t_l - 2.0*pi
                    
                  lv(i_t, i_s) = t_l                   ! update link

              endif                       
            
C****************************** HORIZZONTAL LINK *******************************

              c_l = lh(i_t, i_s)                     ! value of current link
              t_l = c_l + (ran2()*2. - 1.)*delta     ! new trial value of link
                
              f1 = - lv(i_t,i_s)
     &             - lh(nl(i_t, 2), i_s) + lv(i_t, nl(i_s, 1))    ! frist staple

              f2 = lv(ns(i_t, 2), i_s) - lh(ns(i_t, 2), i_s) -    ! second staple
     &             lv(ns(i_t, 2), nl(i_s, 1))
                
              F_c = real(exp(c*(f1 + c_l)) + exp(c*(f2 + c_l)))   ! current action
              F_p = real(exp(c*(f1 + t_l)) + exp(c*(f2 + t_l)))   ! trial action
              r = beta*(F_p - F_c)
              prob = log(ran2())
                
              if (prob < r) then                       ! acceptance test
                    
                  t_l = mod(t_l, 2.0*pi)               ! the value must be
                  if (t_l < (-pi)) t_l = t_l + 2.0*pi  ! between -pi and pi
                  if (t_l >   pi)  t_l = t_l - 2.0*pi
                    
                  lh(i_t, i_s) = t_l                   ! update link

              endif 
          enddo
      enddo
        
      return
      end

c============================================================================
C OVER RELAXATION
c============================================================================

      subroutine over_relaxation()
       
      include "par.f"    ! file with all parameters
      
      complex :: c, ff
      c = (0.0, 1.)      ! immaginary units
      
      pi = 4 * atan(1.0) ! pi greco
      
      do i_t = 1, NT     ! loop over lattice
          do i_s = 1, NX
            
C****************************** VERTICAL LINK  *******************************

              c_l = lv(i_t, i_s)                   ! value of current link

              f1 = - lh(i_t, i_s)
     &             - lv(i_t, nl(i_s, 1)) + lh(nl(i_t, 2), i_s)       ! first staple

              f2 = lh(i_t, ns(i_s, 1)) -
     &             lv(i_t, ns(i_s, 1)) - lh(nl(i_t, 2), ns(i_s, 1))  ! second staple
                
              ff = exp(-c*f1) + exp(-c*f2)
              xf = imag(log(ff))
              
		  xf = mod(xf, 2.0*pi)                 ! the value must be
              if (xf < (-pi)) xf = xf + 2.0*pi     ! between -pi and pi
              if (xf >   pi)  xf = xf - 2.0*pi 
              t_l = 2*xf - c_l
              
              t_l = mod(t_l, 2.0*pi)               ! the value must be
              if (t_l < (-pi)) t_l = t_l + 2.0*pi  ! between -pi and pi
              if (t_l >   pi)  t_l = t_l - 2.0*pi
              
              lv(i_t, i_s) = t_l                   ! update link
            
C****************************** HORIZZONTAL LINK *******************************

              c_l = lh(i_t, i_s)                   ! value of current link

              f1 = - lv(i_t,i_s)
     &             - lh(nl(i_t, 2), i_s) + lv(i_t, nl(i_s, 1))   ! first staple

              f2 = lv(ns(i_t, 2), i_s) - lh(ns(i_t, 2), i_s) -   ! second staple
     &             lv(ns(i_t, 2), nl(i_s, 1))
                
              ff = exp(-c*f1) + exp(-c*f2)
              xf = imag(log(ff))
              
		  xf = mod(xf, 2.0*pi)                 ! the value must be
              if (xf < (-pi)) xf = xf + 2.0*pi     ! between -pi and pi
              if (xf >   pi)  xf = xf - 2.0*pi 
              t_l = 2*xf - c_l
              
              t_l = mod(t_l, 2.0*pi)               ! the value must be
              if (t_l < (-pi)) t_l = t_l + 2.0*pi  ! between -pi and pi
              if (t_l >   pi)  t_l = t_l - 2.0*pi
                            
              lh(i_t, i_s) = t_l                   ! update link
              
          enddo
      enddo
        
      return
      end

c============================================================================
C Charge measurement
c============================================================================

      subroutine charge(q)
       
      include "par.f"    ! file with all parameters
      
      pi = 4 * atan(1.0) ! pi greco
      q = 0.0            ! charge
      plaq = 0.0         ! plaquette
      
      do i_t = 1, NT     ! loop over lattice
          do i_s = 1, NX
          
              square = lh(i_t, i_s) + lv(i_t, nl(i_s, 1)) -
     &                 lh(nl(i_t, 2), i_s) - lv(i_t, i_s)
              
              square = mod(square, 2.0*pi)                   ! the value must be
              if (square < (-pi)) square = square + 2.0*pi   ! between -pi and pi
              if (square >   pi ) square = square - 2.0*pi
              
              plaq = plaq + square
              
          enddo
      enddo

      q = int(plaq/(2.0*pi) + 0.1)
      !q = int(plaq/(2.0*pi))
      return
      end

c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================

