************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use PARAMETERS,ONLY: elements,abund_pick,model_dim,model_pconst,
     >                     model_struc,model_eqcond,Npoints,
     >                     Tmin,Tmax,pmin,pmax,nHmin,nHmax
      use DUST_DATA,ONLY: bar
      implicit none
      integer :: iarg,iline
      character(len=200) :: ParamFile,line

      iarg = iargc()
      if (iarg==0) then
        print*,"*** syntax error."
        print*,"expected syntax is, for example, ./ggchem default.in"
        stop
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
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      end
