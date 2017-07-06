************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use PARAMETERS,ONLY: elements,abund_pick,model_dim,model_pconst,
     >                     model_struc,model_eqcond,Npoints,
     >                     Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax
      use CHEMISTRY,ONLY: NewChemIt,NewBackIt,dispol_file
      use DUST_DATA,ONLY: bar
      implicit none
      integer :: iarg,iline,i
      character(len=200) :: ParamFile,line

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      dispol_file(1) = 'dispol_new.dat'
      dispol_file(2) = ''
      dispol_file(3) = ''
      dispol_file(4) = ''
      elements     = 'H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li el'
      abund_pick   = 3
      model_eqcond = .false.
      model_dim    = 1
      model_pconst = .true.
      model_struc  = .false.
      Npoints      = 200
      Tfast        = 750.d0
      Tmin         = 100.d0
      Tmax         = 6000.d0
      pmin         = 1.d0*bar
      pmax         = 1.d0*bar
      nHmin        = 4.d+19
      nHmax        = 4.d+19
      NewChemIt    = .true.
      NewBackIt    = 5

      !-------------------------------------------
      ! ***  change parameters via input file  ***
      !-------------------------------------------
      iarg = iargc()
      if (iarg==0) then
        print*,"using default parameters"
        return
      endif  
      call getarg(1,ParamFile)
      open(unit=1,file=ParamFile,status='old')
      iline = 0
      do 
        read(1,'(A200)',end=100) line
        if (line(1:1).eq.'#') cycle           ! ignore comment lines
        if (len(trim(line))==0) cycle         ! ignore comment lines
        iline = iline+1
        print*,trim(line)
        if (iline.eq.1) then
          elements = ' '//trim(line)//' '     ! selection of element 
        else if (index(line,"! abund_pick")>0) then   
          read(line,*) abund_pick
        else if (index(line,"! model_eqcond")>0) then   
          read(line,*) model_eqcond
        else if (index(line,"! model_dim")>0) then   
          read(line,*) model_dim
        else if (index(line,"! model_pconst")>0) then   
          read(line,*) model_pconst
        else if (index(line,"! model_struc")>0) then   
          read(line,*) model_struc
        else if (index(line,"! Tmax")>0) then   
          read(line,*) Tmax
        else if (index(line,"! Tmin")>0) then   
          read(line,*) Tmin
        else if (index(line,"! Tfast")>0) then   
          read(line,*) Tfast
        else if (index(line,"! pmax")>0) then   
          read(line,*) pmax
          pmax = pmax*bar
        else if (index(line,"! pmin")>0) then   
          read(line,*) pmin
          pmin = pmin*bar
        else if (index(line,"! nHmax")>0) then 
          read(line,*) nHmax
        else if (index(line,"! nHmin")>0) then 
          read(line,*) nHmin
        else if (index(line,"! Npoints")>0) then 
          read(line,*) Npoints
        else if (index(line,"! NewChemIt")>0) then 
          read(line,*) NewChemIt
        else if (index(line,"! NewBackIt")>0) then 
          read(line,*) NewBackIt
        else if (index(line,"! dispol_file2")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) dispol_file(2)
        else if (index(line,"! dispol_file3")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) dispol_file(3)
        else if (index(line,"! dispol_file4")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) dispol_file(4)
        else if (index(line,"! dispol_file")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) dispol_file(1)
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      end
