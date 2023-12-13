
!-------------------------------------------------------------------------
real*8 function GAS_OPACITY(Tgas,nH)
!-------------------------------------------------------------------------
! ***  Calculates the Rosseland mean gas opacity for solar abundances  ***  
! ***  (c) Paola Marigo 2023, private communication, in [g/cm2].       ***
! ***  https://ui.adsabs.harvard.edu/abs/2023arXiv231014588M/abstract  ***
! ***  explains the basics, the paper has gas+dust opacities.          ***
! ***  Here, this is only the gas opacity without dust.                ***
!-------------------------------------------------------------------------  
  implicit none
  real(8),intent(in)  :: Tgas         ! T[K]
  real(8),intent(in)  :: nH           ! n<H>[cm-3]
  integer, parameter  :: num_temperatures = 231
  integer, parameter  :: num_opacities = 141
  integer, parameter  :: num_hydrogen = 141
  integer :: size, i, j, jt, jnh,iread
  real(8) :: x, y, y1, y2, y3,y4, t, u, opac
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
    !--- Populate the vector using a DO loop ---
    do i = 1, size
      hydrogen_values(i) = start + (i - 1) * step
    enddo
    call READ_OPACITY_TABLE(temperature_values, opacity_values)
    firstCall = .false.
  endif

  x = LOG10(Tgas)
  y = LOG10(nH)
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
  GAS_OPACITY = 10.d0**opac
  return
  
  contains

  subroutine READ_OPACITY_TABLE(temperature_values, opacity_values)
  implicit none
  integer, parameter :: num_temperatures = 231
  integer, parameter :: num_opacities = 141
  real(8),intent(out),dimension(num_temperatures) :: temperature_values
  real(8),intent(out),dimension(num_temperatures,num_opacities) :: opacity_values
  character(len=100) :: filename
  character(len=1) :: s
  integer :: i, j
  logical :: file_exists
  
  ! Specify the file name
  filename = "data/kRossGas/aesopus2.0_gas_AGSS09_Z0.0134_X0.7374.tab"
  ! Check if the file exists
  inquire(file=filename, exist=file_exists)
  if (.not.file_exists) then
    print *, "Error: File '", filename, "' not found."
    stop
  end if
  print'(A100)',"reading Rosseland opacity table "//trim(filename)
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
  implicit none  
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

end function GAS_OPACITY
