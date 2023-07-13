	program analisi
	
C======================================================================	
C Program for analysis and calculation of relevant quantities
C======================================================================
	
	real*8, dimension(:), allocatable :: Q     ! array for markov chain
	real*8:: dq, dx                            ! error on charge an susceptibility
	real*8:: aq, ax                            ! average charge an susceptibility
	integer :: R                               ! number of resampling 
	integer(16) :: Db                         ! Size of block
	integer, dimension(5) :: cp                ! axuliar array
      character(len=30) file_name, formatstring  ! for loop over files
	
	call cpu_time(start)
	call ranstart
	
	open(unit=2, file='plot_data/datiplot.dat',status='unknown')
	
	R  = 100				! number of resampling
	Db = 1000                     ! size of block
	cp = [36, 48, 56, 60, 64]     ! size of lattice
	
	write(2,*) "#N, beta, charge^2, susc, err_charge^2, err_susc"
	
	do j = 1, 5
	
	    formatstring = "(A13,i2,A4)"
          write(file_name, formatstring) "raw_data/dati",cp(j),".dat"
          open(unit=0, file=file_name, status="old")
            
	    read(0, *) N			! I read the first values ​​that are 
	    read(0, *) beta		! necessary for the analysis
	    read(0, *) NX, NT
	      
	    nvol = NX*NT 
	      
	    allocate(Q(N))
	    
	    do i = 1, N
	        read(0, *) Q(i)      ! read data
	    enddo
	     
	    ! computation of phisical quantities
	    aq = sum(Q)/float(N*nvol)     ! it should be zero
	          
	    !ax = (sum(Q**2)/float(N) - (sum(Q)/float(N))**2)/float(nvol)
	    ax = sum(Q**2)/float(N*nvol)
	    
	    !computation of errors     
          call bootstrap(N, Q, dq, 0, R, nvol, Db)
          call bootstrap(N, Q, dx, 1, R, nvol, Db)
          
          close(0)
          deallocate(Q)

          write(2,*) cp(j), beta, aq, ax, dq, dx	! save on file
          write(*,*) cp(j), beta, aq, ax, dq, dx      ! write on standad output

	enddo
	
	call ranfinish
	
	call cpu_time(finish)
	print '("tempo di esecuzione= ", f8.4," secondi.")', finish-start
	
	end program analisi
	
C=============================================================================
C Bootstrap
C=============================================================================

	subroutine bootstrap(N, x, dx, ics, R, nvol, Db)
C=============================================================================
C     Subroutine for calculating errors using binned bootstrap
C     
C     Parameters
C     N : int
C         size of data, length(x)
C     x : one dimensional array of  length(x) = N
C         data, markov chain
C     dx : float
C         variable to which the error will be written
C     ics : int
C         flag if 0 computere error on mean, if 1 compute error on variance
C     R : int
C         number of resampling
C     nvol : float
C         volume of lattice
C     Db : int
C         seize of the blocks
C=============================================================================
	real*8, dimension(:), allocatable :: z, a   ! axuliar array
	real*8, dimension(N) :: x                   ! initial array
	integer(16) :: nb, Db, i, j, l              ! parameter and indices
	real*8 :: media_x, dx                       ! final results
	integer :: g, R, ics                        ! others, pamaters
	
	allocate(z(N), a(R))
	
	! The calculation of the error is done through the binned bootstrap since
      ! it is not possible to know a priori which is the best decorellation time
      ! to insert in the simulation, and in order not to waste machine time,
      ! a possible correlation between the data must be taken into account
	
	nb = N/Db		      ! number of blocks
	
	do l = 1, R			! loop of resampling
	
	    do i = 1, nb
		  j = int(ran2()*N +1) ! I choose site at random
		  do g = 1, Db		 	
			z((i-1)*Db+g) = x(mod(j+g-2,N)+1) ! block's resampling	
		  enddo	
	    enddo
		
	    ! calculation of the mean or the variance of the resamplings
	    
	    if (ics==0) then
	        a(l) = sum(z)/float(N*nvol)
	    endif
	    
	    if (ics==1) then 
	        a(l) = sum(z**2)/float(N) !- (sum(z)/float(N))**2
	        a(l) = a(l)/float(nvol)
	    endif
	    
	enddo
	
	media_x  = sum(a)/float(R)         ! mean
	dx = 0
	do i=1, R
	    dx = dx + (a(i) - media_x )**2 ! standard deviation
	enddo
	
	dx  = sqrt(dx/float(R - 1))		
	! I take the error on the sample because we made R
      ! resamples so I only divide by R-1 and not R(R-1)
	
	return
	end
	
c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*4 ran2,am,eps,rnmx
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
