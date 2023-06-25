C================================================================
C File with simulations' parameters     
C================================================================
      
      parameter(NX       = 64   )         ! Spatial size of lattice
      parameter(NT       = NX   )         ! Impose square lattice
      parameter(measures = 10000)         ! number of measures
      parameter(i_dec    = 35   )         ! updating between measures
      parameter(i_term   = 10000)         ! initial steps for termalizzation
      parameter(beta     = 6.4  )         ! parameter of the action, "temperature"
      parameter(delta    = 0.5  )         ! metropolis step's size
      
      ! Some value of N = NX = NT and beta with N**2/beta = 640 
      ! We choose a set of N and after compute beta
      ! We also report other values ​​used in the simulation 
      
      ! N     = 36,    48,  56,  60,    64
      ! beta  = 2.025, 3.6, 4.9, 5.625, 6.4
      ! i_dec = 20,    20,  20,  30,    35 
      
      real lv, lh ! lattice
      
      common/lattice/lv(NT, NX),lh(NT, NX)  ! lattice common for all routine link vertical and orizzontal
      common/neighbors/nl(NX, 2),ns(NX, 2)  ! matrix of nearest neighbors
      
C============================================================================
C This file makes it easier to change the parameters of the simulation
C by doing it once here rather than several times within the simulation code.
C Performing a cycle on the size lattice immediately would mean a very
C long simulation with the risk of not being able to verify that everything
C is working properly and therefore a potential waste of time 
C============================================================================
