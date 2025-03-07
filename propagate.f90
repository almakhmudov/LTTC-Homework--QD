!-----------------!
 module parameters
!-----------------!

   real(8), parameter :: length=5.12d0                                      ! length of the box (in Bohr)
   real(8), parameter :: mass = 1822.88839d0                                ! = 1 atomic mass unit (=1 g/mol)
   real(8), parameter :: pi=3.141592653589793d0
   real(8), parameter :: au2kcalmol=627.509d0
   real(8), parameter :: fs2au=41.341373336561d0
   real(8) :: angfreq, barrier
   character(10) :: potentialtype

   interface
      subroutine initpsi(npoints, dx, alpha, x0, coeff, psi0, mass, angfreq, pi)
         integer, intent(in) :: npoints
         real(8), intent(in) :: dx, alpha, x0, mass, angfreq, pi
         real(8), intent(in) :: coeff(:)
         complex(8), intent(out) :: psi0(npoints)
      end subroutine initpsi
   end interface

 end module parameters
!---------------------!
   
!-----------------!
 program propagate
!-----------------!
   use parameters
   implicit none

   integer :: npoints, ntime, snapshot, i, nmax
   real(8) :: alpha, dt, t, dx, x0
   real(8), allocatable :: pot(:), kin(:), psisquare(:)
   complex(8), allocatable :: psi(:), psi0(:), exppot(:), expkin(:)
   real(8), allocatable :: coeff(:)

   open(unit=10, file='wavepacket')
   read(10,*) npoints                                                       ! Number of lattice points
   read(10,*) x0                                                            ! Initial position
   read(10,*) alpha                                                         ! Governs the initial width of the wave packet
   read(10,*) dt                                                            ! Propagation time step
   read(10,*) ntime                                                         ! Number of propagation steps
   read(10,*) snapshot                                                      ! Snapshot frequency
   read(10,*) nmax                                                          ! Used for allocating memory for coefficients
   allocate(coeff(nmax))                                                    ! Allocate memory for coefficients
   read(10,*) coeff(:)                                                      ! Read coefficients from file
   close(10)

   open(unit=11, file='potential')
   read(11,*) potentialtype                                                 ! Harmonic or double well potential
   read(11,*) angfreq                                                       ! Angular frequency for harmonic potential
   read(11,*) barrier                                                       ! Height of barrier in double well potential (in kcal/mol)
   close(11)

   dt = dt * fs2au                                                          ! Convert femtoseconds to atomic units
   angfreq = angfreq / fs2au                                                ! Convert femtoseconds to atomic units

   allocate(psi(npoints), psi0(npoints))                                    ! Allocate memory for wavepacket psi
   allocate(pot(npoints), exppot(npoints))                                  ! Allocate memory for potential and its exponential
   allocate(kin(npoints), expkin(npoints))                                  ! Allocate memory for kinetic and its exponential
   allocate(psisquare(npoints))                                             ! Allocate memory for psi squared

   dx = length / dble(npoints)

   call initpsi(npoints, dx, alpha, x0, coeff, psi0, mass, angfreq, pi)     ! Obtain initial wavepacket psi0
   call fourier(0, npoints, psi0)                                           ! Initialize the FFT
   call operators(npoints, dx, dt, pot, kin, exppot, expkin)                ! Calculate the kinetic and potential operators

   psi = psi0                                                               ! Set the wavepacket psi at t=0 equal psi0
   do i = 0, ntime                                                          ! Start propagation
      t = i * dt
      if (i > 0) then
         psi = psi * exppot                                                 ! Multiply psi with exp(-i*dt*potential)
         call fourier(1, npoints, psi)                                      ! Forward FFT to momentum space
         psi = psi * expkin                                                 ! Multiply psi with exp(-i*dt*kinetic operator)
         call fourier(-1, npoints, psi)                                     ! Backward FFT to position space
      endif
      if (mod(i, snapshot) == 0) then                                       ! Take a snapshot if i/snapshot remainder is 0
         call initgraph(i / snapshot, t)                                    ! Initialize graph
         psisquare = (abs(psi))**2
         call graphpot(dx, npoints)                                         ! Plot the potential
         call graphpsi(dx, npoints, psisquare)                              ! Plot |psi|^2
      endif
   end do                                                                   ! End propagation

   deallocate(psi, psi0, coeff)
   deallocate(pot, exppot)
   deallocate(kin, expkin)
   deallocate(psisquare)

 end program propagate
!---------------------!
   
!-------------------------------------------------------------------------!
 subroutine initpsi(npoints, dx, alpha, x0, coeff, psi0, mass, angfreq, pi)
!-------------------------------------------------------------------------!
   implicit none

   integer :: i, j, npoints, n, nmax
   real(8) :: alpha, x, x0, dx, norm, fact, mass, angfreq, pi, hermite
   real(8) :: coeff(:)
   complex(8) :: psi0(npoints)

   norm = 0.0d0                                                             ! Initialise the norm
   nmax = size(coeff) - 1                                                   ! The max order of the Hermite polynomial
   do i = -npoints/2 + 1, npoints/2                                         ! Loop over all lattice points to get the initial wavepacket
      x = dble(i) * dx
      if (i > 0) then
         j = i
      else     
         j = i + npoints
      endif
      psi0(j) = 0.0d0                                                       ! Initialise psi0
      do n = 0, nmax                                                        ! Loop through all eigenstates
         call factorial(n, fact)                                            ! Calculate the factorial of n
         psi0(j) = psi0(j) + &                                              ! Superposition of all eigenstates
                   coeff(n+1) * &
                   (1.0d0 / sqrt(2.0d0**n * fact)) * &
                   (mass * angfreq / pi)**0.25d0 * &
                   exp(-mass * angfreq * alpha * (x-x0)**2.0d0 / 2.0d0) * &
                   hermite(n, (x-x0) * sqrt(mass * angfreq)) 
      end do
      norm = norm + abs(psi0(j))**2.0d0 * dx                                ! Calculate the norm
   end do

   norm = 1.0d0 / sqrt(norm)                                                ! Calculate the normalisation factor
   do i = -npoints/2 + 1, npoints/2
      x = dble(i) * dx
      if (i > 0) then
         j = i
      else     
         j = i + npoints
      endif
      psi0(j) = norm * psi0(j)                                              ! Normalise the wavepacket
   end do 

 end subroutine initpsi
!----------------------!

!------------------------------!
 subroutine factorial(n, result)
!------------------------------!
   implicit none
   integer, intent(in) :: n
   real(8), intent(out) :: result
   integer :: i

   result = 1.0d0
   do i = 1, n
      result = result * i
   end do

 end subroutine factorial
!------------------------!

!--------------------------------!
 function hermite(n, y) result(Hn)
!--------------------------------!
   implicit none
   
   integer, intent(in) :: n
   real(8), intent(in) :: y
   real(8) :: Hn, H0, H1, H_prev
   integer :: i

   if (n == 0) then                                                         ! Hermite polynomial of order 0
      Hn = 1.0d0
   elseif (n == 1) then                                                     ! Hermite polynomial of order 1
      Hn = 2.0d0 * y
   else                                                                     ! Hermite polynomial of higher order
      H0 = 1.0d0
      H1 = 2.0d0 * y
      do i = 2, n                                                           ! The solution is generalised up to the order n
         H_prev = H1
         H1 = 2.0d0 * y * H1 - 2.0d0 * (i - 1) * H0                         ! Calculated using the recursion relation
         H0 = H_prev
      end do
      Hn = H1
   endif

 end function hermite
!--------------------!

!--------------------------------------------------------!
 subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!--------------------------------------------------------!
   use parameters
   implicit none

   integer :: i,j,npoints
   real(8) :: x,p,b,dt,dx,dp,pot(npoints),kin(npoints)
   complex(8) :: exppot(npoints),expkin(npoints)

   dp=2.d0*pi/length
   do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      p=dble(i-1)*dp
      if (i>0) then
         j=i
      else
         j=i+npoints
      endif
      if (potentialtype=='harmonic') then
         pot(j)=0.5d0*mass*angfreq**2*x**2
      elseif (potentialtype=='doublewell') then
         pot(j)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
      endif
      kin(j)=0.5d0*p**2/mass
      exppot(j)=exp(-dt*(0,1)*pot(j))
      expkin(j)=exp(-dt*(0,1)*kin(j))
   end do

 end subroutine operators
!------------------------!
   
!----------------------------------!
 subroutine fourier(dir,npoints,psi)
!----------------------------------!
   implicit none

   integer :: i,npoints,dir
   real(8) :: nr
   complex(8) :: psi(npoints)
   real(8), allocatable, save :: wsave(:)       

   if (dir==1) then
      call dcfftf(npoints,psi,wsave)
      nr=1.d0/dble(npoints)
      do i=1,npoints
         psi(i)=psi(i)*nr
      end do
   elseif (dir==-1) then
      call dcfftb(npoints,psi,wsave)
   elseif (dir==0) then
      if (allocated(wsave)) deallocate(wsave)
      allocate(wsave(4*npoints+20))
      call dcffti(npoints,wsave)
   endif

 end subroutine fourier
!----------------------!