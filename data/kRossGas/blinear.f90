program main
  implicit none
  real(8) :: Tgas,nHtot,x,y,kRoss
  Tgas = 1500.0
  nHtot = 1.E+14
  x=LOG10(Tgas)
  y=LOG10(nHtot)
  call GAS_OPACITY(x,y,kRoss)
  print'(" n<H>,T,kRoss=",0pF9.2,1pE11.4,1pE11.4)',Tgas,nHtot,kRoss
  call GAS_OPACITY(x,y,kRoss)
  print'(" n<H>,T,kRoss=",0pF9.2,1pE11.4,1pE11.4)',Tgas,nHtot,kRoss
end program main

!-------------------------------------------------------------------------
subroutine GAS_OPACITY(x,y,kRoss)
!-------------------------------------------------------------------------
! ***  Calculates the Rosseland mean gas opacity for solar abundances  ***  
! ***  (c) Paola Marigo 2023, private communication.                   ***
! ***  https://ui.adsabs.harvard.edu/abs/2023arXiv231014588M/abstract  ***
! ***  explains the basics, the paper has gas+dust opacities.          ***
! ***  Here, this is only the gas opacity without dust                 ***
!-----------------------------------------------------------------------  
  implicit none
  real(8),intent(in)  :: x      ! LOG10(T)    [K]
  real(8),intent(in)  :: y      ! LOG10(n<H>) [cm-3]
  real(8),intent(out) :: kRoss  ! Rosseland opacity [cm2/g]
  integer, parameter  :: num_temperatures = 231
  integer, parameter  :: num_opacities = 141
  integer, parameter  :: num_hydrogen = 141
  integer :: size, i, j, jt, jnh,iread
  real(8) :: y1, y2, y3,y4, t, u, opac
  real(8) :: start, stop, step
  real(8),save,dimension(num_temperatures) :: temperature_values
  real(8),save,dimension(num_temperatures,num_opacities) :: opacity_values
  real(8),save,dimension(num_hydrogen) :: hydrogen_values
  logical,save :: firstCall=.true.
  
  if (firstCall) then
    !--- Define the range and step of log10(nH) ---
    start = 6.0d0
    stop = 20.0d0
    step = 0.1d0
    size = int((stop - start) / step) + 1
    !print*, size  
    !--- Populate the vector using a DO loop ---
    do i = 1, size
      hydrogen_values(i) = start + (i - 1) * step
      !print*,hydrogen_values(i)
    enddo
    call READ_OPACITY_TABLE(temperature_values, opacity_values)
    firstCall = .false.
  endif
  
  call LOCATE(temperature_values, num_temperatures, x, jt)
  jt = max(1,jt)
  jt = min(num_temperatures-1,jt)
  call LOCATE(hydrogen_values, num_hydrogen, y, jnh)
  jnh = max(1,jnh)
  jnh = min(num_hydrogen-1,jnh)
  
  ! Perform bilinear interpolation  
  t=(x-temperature_values(jt))/(temperature_values(jt+1)-temperature_values(jt))
  u=(y-hydrogen_values(jnh))/(hydrogen_values(jnh+1)-hydrogen_values(jnh))
  
  y1=opacity_values(jt,jnh)
  y2=opacity_values(jt+1,jnh)
  y3=opacity_values(jt+1,jnh+1)
  y4=opacity_values(jt,jnh+1)
  
  opac = (1.d0 - t)*(1.d0 - u)*y1 + t*(1.d0 - u)*y2 + t*u*y3 + (1.d0 - t)*u*y4
  kRoss = 10.d0**opac
  
  ! Display the result
  ! print *, "Interpolated opacity value at (", x, ",", y, ") =", opac

  contains

  subroutine READ_OPACITY_TABLE(temperature_values, opacity_values)
  implicit none
  integer, parameter :: num_temperatures = 231
  integer, parameter :: num_opacities = 141
  real(8),intent(out),dimension(num_temperatures) :: temperature_values
  real(8),intent(out),dimension(num_temperatures, num_opacities) :: opacity_values
  character(len=50) :: filename
  character(len=1) :: s
  integer :: i, j
  logical :: file_exists
  ! Specify the file name
  filename = "aesopus2.0_gas_AGSS09_Z0.0134_X0.7374.tab"
  ! Check if the file exists
  inquire(file=filename, exist=file_exists)
  if (.not.file_exists) then
    print *, "Error: File '", filename, "' not found."
    stop
  end if
  print*,"reading Rosseland opacity table from "//trim(filename)
  open(unit=10, file=filename, status='old', action='read')
  do i=1,330
    read(10,*) s
  enddo
  ! Read temperature values and opacity data from the file
  do i=1,num_temperatures
    read(10,*) temperature_values(i),(opacity_values(i,j),j=1,num_opacities)
  end do
  close(unit=10)
  end subroutine read_opacity_table
  
  subroutine LOCATE(xx, n, x, j)
  integer,intent(in) :: n
  real*8,intent(in) :: x, xx(n)
  integer,intent(out) :: j
  integer :: jl, jm, ju
  jl = 0
  ju = n + 1
10 continue
  if (ju - jl > 1) then
    jm = (ju + jl) / 2
    if ((xx(n) > xx(1)).eqv.(x > xx(jm))) then
      jl = jm
    else
      ju = jm
    endif
    goto 10
  endif
  j = jl
  end subroutine LOCATE

end subroutine GAS_OPACITY
