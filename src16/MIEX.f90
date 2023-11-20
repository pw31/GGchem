! ******************************************************************
! MIEX: MIE SCATTERING CODE FOR LARGE GRAINS
!  _____________________________________________________
!  Contact information:   swolf@mpia.de (Sebastian Wolf)
! ==================================================================

module datatype
  implicit none
  integer, parameter, public ::  r1=selected_real_kind(1)  ! real*4
  integer, parameter, public ::  r2=selected_real_kind(9)  ! real*8 (double precision)
end module datatype


! ====================================================================================================
! Collection of subroutines
! ====================================================================================================
module mie_routines
  private :: aa2
  public  :: shexqnn2
contains
  ! ==================================================================================================
  ! Subroutine for calculations of the ratio of derivative to the function for Bessel functions
  ! of half order with complex argument: J'(n)/J(n). The calculations are given by the recursive
  ! expression ``from top to bottom'' beginning from n=num.
  ! *  a=1/x (a=2*pi*a(particle radius)/lambda - size parameter).
  ! *  ri - complex refractive index.
  ! *  ru-array of results.
  ! - this routine is based on the routine 'aa' published by
  !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
  !                           Astrophysics and Space Physics Review 12,  1 (2002)
  ! ==================================================================================================
  subroutine aa2( a, ri, num, ru )
    use datatype

    implicit none

    ! variables for data exchange.....................................................................
    real(kind=r2), intent(in)                   :: a
    complex(kind=r2), intent(in)                :: ri
    integer, intent(in)                         :: num
    complex(kind=r2), dimension(:), intent(out) :: ru

    ! local variables.................................................................................
    integer :: i, i1, j, num1
    complex(kind=r2) :: s, s1
    !-------------------------------------------------------------------------------------------------
    ! initialisierung: not necessary (+ slowes the code down remarkably)
    ! ru(:) = (0.0, 0.0)

    s       = a / ri
    ru(num) = real(num+1,kind=r2) * s
    num1    = num - 1
    do j=1, num1
       i     = num - j
       i1    = i + 1
       s1    = i1 * s
       ru(i) = s1 - 1.0_r2 / (ru(i1) + s1)
    end do
  end subroutine aa2


  !===================================================================================================
  ! shexqnn2
  ! --------
  ! - for a given size parameter 'x' and (complex) refractive index 'ri' the following quantities
  !   are determined:
  !   * Qext     - extinction effiency
  !   * Qsca     - scattering effiency
  !   * Qabs     - absorption effiency
  !   * Qbk      - backscattering effiency
  !   * Qpr      - radiation pressure effiency
  !   * albedo   - Albedo
  !   * g        - g scattering assymetry factor
  !   * SA1, SA2 - scattering amplitude function
  ! - further input parameters
  !   * doSA = .true.  ->  calculation of the scattering amplitudes
  !   * nang ... half number of scattering angles theta in the intervall 0...PI/2
  !              (equidistantly distributed)
  ! - this routine is based on the routine 'shexqnn' published by
  !       N.V.Voshchinnikov: "Optics of Cosmic Dust",
  !                           Astrophysics and Space Physics Review 12,  1 (2002)
  !===================================================================================================
  subroutine shexqnn2( ri, x, Qext, Qsca, Qabs, Qbk, Qpr, albedo, g, ier, SA1, SA2, doSA, nang )
    use datatype

    implicit none
    
    ! variables for data exchange.....................................................................
    complex(kind=r2), intent(in)                :: ri
    real(kind=r2), intent(in)                   :: x
    real(kind=r2), intent(out)                  :: Qext, Qsca, Qabs, Qbk, Qpr, albedo, g
    integer, intent(out)                        :: ier
    complex(kind=r2), dimension(:), intent(out) :: SA1, SA2
    logical, intent(in)                         :: doSA
    integer, intent(in)                         :: nang

    ! local variables.................................................................................
    integer       :: iterm, nterms, num, iu0, iu1, iu2, iang2, iang
    real(kind=r2) :: r_iterm, factor, eps, pi, ax, besJ0, besJ1, besJ2, besY0, besY1, besY2, b, an, &
         y, ass, w1, qq, fac, an2, P, T, Si, Co, z, xmin
    complex(kind=r2) :: ra0, rb0, ra1, rb1, r, ss, s1, s2, s3, s, rr

    real(kind=r2),dimension(0:1)                   :: fact
    real(kind=r2),dimension(:),allocatable,save    :: mu, fpi, fpi0, fpi1, ftau
    complex(kind=r2),dimension(:),allocatable,save :: ru
    integer                                        :: istat
    logical,external :: is_nan
    !-------------------------------------------------------------------------------------------------
    ! Maximum number of terms to be considered
    nterms = 20000000
    ! this works for x up to 1.d9, but needs a hell more of memory!!
    !nterms= 550000000

    ! Accuracy to be achieved
    eps    = 1.0e-15_r2

    ! Minimum size parameter
    xmin   = 1.0e-6_r2

    !-------------------------------------------------------------------------------------------------
    ! initialization
    if (.not.allocated(ru)) then
      allocate( ru(1:nterms), mu(1:nang), fpi(1:nang), &
              & fpi0(1:nang), fpi1(1:nang), ftau(1:nang),stat=istat)
      if (istat.ne.0) then
        write(*,*) "nterms: ",nterms,"nang: ",nang
        stop "*** Not enough memory for MIEX in shexqnn2."
      endif
    endif
    ier     = 0
    Qext    = 0.0_r2
    Qsca    = 0.0_r2
    Qabs    = 0.0_r2
    Qbk     = 0.0_r2
    Qpr     = 0.0_r2
    albedo  = 0.0_r2
    g       = 0.0_r2
    fact(0) = 1.0_r2
    fact(1) = 1.0e+250_r2
    factor  = 1.0e+250_r2

    ! null argument
    if (x <= xmin) then
       ier = 1
       print *, "<!> Error in subroutine shexqnn2:"
       print *, "    - Mie scattering limit exceeded:"
       print *, "      current size parameter: ", x
    else
       pi = 4.0_r2 * atan(1.0_r2) ! PI = 3.14...
       ax = 1.0_r2 / x
       b  = 2.0_r2 * ax**2
       ss = (0.0_r2, 0.0_r2)
       s3 = (0.0_r2,-1.0_r2)
       an = 3.0_r2
       
       ! define the number for subroutine aa2 [Loskutov (1971)]
       y   = sqrt( RI * conjg(ri) )  *  x
       num = 1.25 * y + 15.5
       
       if      ( y<1.0_r2 ) then
          num = 7.5 * y + 9.0
       else if ( (y>100.0_r2) .and. (y<50000.0_r2) ) then
          num = 1.0625 * y + 28.5
       else if ( y>=50000.0_r2 ) then
          num=1.005*y+50.5
       end if

       if(num > nterms) then
          ier = 2
          print *, "<!> Error in subroutine shexqnn2:"
          print *, "    - Maximum number of terms  : ", nterms
          print *, "    - Number of terms required : ", num
          print *, "    ** Solution: Increase default value of the variable 'nterm' **"
       else
          ! logarithmic derivative to Bessel function (complex argument)
          call aa2(ax,ri,num,ru)

          ! ------------------------------------------------------------------------------------------
          ! FIRST TERM
          ! ------------------------------------------------------------------------------------------
          ! initialize term counter
          iterm = 1
          
          ! Bessel functions
          ass = sqrt( pi / 2.0_r2 * ax )
          w1  = 2.0_r2/pi * ax
          Si  = sin(x)/x
          Co  = cos(x)/x
          
          ! n=0
          besJ0 =  Si / ass
          besY0 = -Co / ass
          iu0   = 0
          
          ! n=1
          besJ1 = ( Si * ax - Co) / ass
          besY1 = (-Co * ax - Si) / ass
          iu1   = 0
          iu2   = 0
          
          ! Mie coefficients
          s   = ru(1) / ri + ax
          s1  = s * besJ1 - besJ0
          s2  = s * besY1 - besY0
          ra0 = s1 / (s1 - s3 * s2)   ! coefficient a_1
          
          s   = ru(1) * ri + ax
          s1  = s * besJ1 - besJ0
          s2  = s * besY1 - besY0
          rb0 = s1 / (s1 - s3 * s2)   ! coefficient b_1
          
          ! efficiency factors
          r    = -1.5_r2 * (ra0-rb0)
          Qext = an * (ra0 + rb0)
          Qsca = an * (ra0 * conjg(ra0)  +  rb0 * conjg(rb0))
          
          ! scattering amplitude functions
          if (doSA) then
             do iang=1, nang
                mu(iang) = cos( (real(iang,kind=r2)-1.0_r2) * (pi/2.0_r2)/real(nang-1,kind=r2) )
             end do

             fpi0(:) = 0.0_r2
             fpi1(:) = 1.0_r2
             SA1(:)  = cmplx( 0.0_r2, 0.0_r2 )
             SA2(:)  = cmplx( 0.0_r2, 0.0_r2 )
             
             r_iterm = real(iterm,kind=r2)  ! double precision             
             fac     = (2.0*r_iterm + 1.0_r2) / (r_iterm * (r_iterm+1.0_r2))
             
             do iang=1, nang
                iang2      = 2 * nang - iang
                
                fpi(iang)  = fpi1(iang)
                ftau(iang) = r_iterm * mu(iang) * fpi(iang)  -  (r_iterm+1.0) * fpi0(iang)
                
                P          = (-1.0)**(iterm-1)
                SA1(iang)  = SA1(iang)   +   fac * (ra0*fpi(iang)  + rb0*ftau(iang))      
                
                T          = (-1.0)**iterm
                SA2(iang)  = SA2(iang)   +   fac * (ra0*ftau(iang) + rb0*fpi(iang) )
                
                if  ( iang /= iang2 )  then
                   SA1(iang2) = SA1(iang2)   +   fac * (ra0*fpi( iang)*P + rb0*ftau(iang)*T)
                   SA2(iang2) = SA2(iang2)   +   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                end if
             end do
             
             iterm   = iterm + 1
             r_iterm = real(iterm, kind=r2)
             
             do iang=1, nang
                fpi1(iang) = ((2.0*r_iterm-1.0) / (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                fpi1(iang) = fpi1(iang)   -   r_iterm * fpi0(iang)/(r_iterm-1.0)
                fpi0(iang) = fpi(iang)
             end do
          else
             ! start value for the next terms
             iterm = 2
          end if

          ! ------------------------------------------------------------------------------------------
          ! 2., 3., ... num
          ! ------------------------------------------------------------------------------------------
          z = -1.0_r2

          do        
             an  = an + 2.0_r2
             an2 = an - 2.0_r2
             
             ! Bessel functions
             if(iu1 == iu0) then
                besY2 = an2 * ax * besY1 - besY0
             else
                besY2 = an2 * ax * besY1 - besY0 / factor
             end if
             if(dabs(besY2) > 1.0e+300_r2) then
                besY2 = besY2 / factor
                iu2   = iu1 + 1
             end if
             besJ2 = (w1 + besY2 * besJ1) / besY1
             
             ! Mie coefficients
             r_iterm = real(iterm,kind=r2)
             
             s   = ru(iterm) / ri + r_iterm * ax
             if(iu1>1) then
                ier=1
                return
             endif
             if(iu2>1) then
                ier=1
                return
             endif  
             s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1) ! Subscript #1 of the array FACT has value 2 which is greater than the upper bound of 1
             s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
             ra1 = s1 / (s1 - s3 * s2)                        ! coefficient a_n, (n=iterm)
             
             s   = ru(iterm) * ri + r_iterm * ax
             s1  = s * besJ2 / fact(iu2) - besJ1 / fact(iu1)
             s2  = s * besY2 * fact(iu2) - besY1 * fact(iu1)
             rb1 = s1 / (s1 - s3 * s2)                        ! coefficient b_n, (n=iterm)
             
             ! efficiency factors
             z  = -z
             rr = z * (r_iterm + 0.5_r2) * (ra1 - rb1)
             r  = r + rr
             ss = ss + (r_iterm - 1.0_r2) * (r_iterm + 1.0_r2) / r_iterm * (ra0 * conjg(ra1)  &
                  + rb0 * conjg(rb1)) &
                  + an2 / r_iterm / (r_iterm - 1.0_r2) * (ra0 * conjg(rb0))
             qq   = an * (ra1 + rb1)
             Qext = Qext + qq
             Qsca = Qsca + an * (ra1 * conjg(ra1) + rb1 * conjg(rb1))
             if(IS_NAN(Qext)) then
                 ier=1
                 return
             endif    
             ! leaving-the-loop criterion
             if ( dabs(qq / qext) < eps ) then
                exit
             end if
             
             ! Bessel functions
             besJ0 = besJ1
             besJ1 = besJ2
             besY0 = besY1
             besY1 = besY2
             iu0   = iu1
             iu1   = iu2
             ra0   = ra1
             rb0   = rb1
             
             ! scattering amplitude functions
             if (doSA) then
                r_iterm = real(iterm,kind=r2)       
                fac      = (2.0 * r_iterm+1.0) / (r_iterm * (r_iterm+1.0))
                
                do iang=1, nang
                   iang2      = 2 * nang - iang
                   
                   fpi(iang)  = fpi1(iang)
                   ftau(iang) = r_iterm * mu(iang) * fpi(iang)  -  (r_iterm+1.0) * fpi0(iang)
                   
                   P          = (-1.0)**(iterm-1)
                   SA1(iang)  = SA1(iang)   +   fac * (ra0*fpi(iang) + rb0*ftau(iang))      
                   
                   T          = (-1.0)**iterm
                   SA2(iang)  = SA2(iang)   +   fac * (ra0*ftau(iang) + rb0*fpi(iang))
                   
                   if  ( iang /= iang2 ) then
                      SA1(iang2) = SA1(iang2)   +   fac * (ra0*fpi(iang)*P  + rb0*ftau(iang)*T)
                      SA2(iang2) = SA2(iang2)   +   fac * (ra0*ftau(iang)*T + rb0*fpi( iang)*P)
                   end if
                end do
                
                iterm   = iterm + 1
                r_iterm = real(iterm,kind=r2)
                
                do iang=1, nang
                   fpi1(iang) = ((2.0*r_iterm-1.0) / (r_iterm-1.0))   *   mu(iang)  *  fpi(iang)
                   fpi1(iang) = fpi1(iang)   -   r_iterm * fpi0(iang)/(r_iterm-1.0)
                   fpi0(iang) = fpi(iang)
                end do
             else
                iterm = iterm + 1
             endif

             if ( iterm==num ) then
                exit
             else
                cycle
             end if
          end do
          
          ! efficiency factors (final calculations)
          Qext   = b * Qext
          Qsca   = b * Qsca
          Qbk    = 2.0_r2 * b * r * conjg(r)
          Qpr    = Qext - 2.0_r2 * b * ss
          Qabs   = Qext - Qsca
          albedo = Qsca / Qext
          g      = (Qext - Qpr) / Qsca
       end if
    end if
    !deallocate( ru, mu, fpi, fpi0, fpi1, ftau )
    ier=0
    return
  end subroutine shexqnn2
end module mie_routines
