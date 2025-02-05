************************************************************************
      subroutine READ_PARAMETER
************************************************************************
      use PARAMETERS,ONLY: elements,abund_pick,model_dim,model_pconst,
     >     model_struc,model_eqcond,Npoints,useDatabase,verbose,
     >     Tfast,Tmin,Tmax,pmin,pmax,nHmin,nHmax,pick_mfrac,use_SiO,
     >     abund_file,struc_file,remove_condensates,phyllosilicates,
     >     metal_sulphates,output_dispol, 
     >     initchem_info,auto_atmos,stop_after_init,Mpl,Rpl,gamma,
     >     adapt_cond,adapt_file,method_eqcond,Nseq,Tseq,Tmin_atmos,
     >     disk_model 
      use CHEMISTRY,ONLY: NewBackIt,NewFullIt,NewBackFac,NewPreMethod,
     >                    NewFastLevel,Natmax,dispol_file
      use DUST_DATA,ONLY: DustChem_file,bar,MEarth,REarth
      implicit none
      integer :: iarg,iline,i,dispol_set
      character(len=200) :: ParamFile,line

      !-------------------------
      ! ***  default values  ***
      !-------------------------
      dispol_file(1) = 'dispol_BarklemCollet.dat'
      dispol_file(2) = 'dispol_StockKitzmann_withoutTsuji.dat'
      dispol_file(3) = 'dispol_WoitkeRefit.dat'
      dispol_file(4) = ''
      DustChem_file  = 'DustChem.dat'
      elements       = ' H He C N O Na Mg Si Fe Al Ca Ti S Cl K Li el '
      pick_mfrac         = .false.
      model_eqcond       = .false.
      remove_condensates = .false.
      phyllosilicates    = .true.
      metal_sulphates    = .true.
      use_SiO            = .true.
      auto_atmos         = .false.
      adapt_cond         = .false.
      stop_after_init    = .false.
      model_pconst       = .true.
      initchem_info      = .true.
      disk_model         = .true.
      output_dispol      = .false.
      abund_pick   = 3
      model_dim    = 0
      model_struc  = 0
      verbose      = 0
      Npoints      = 100
      Tfast        = 1000.d0
      Tmin         = 100.d0
      Tmax         = 1000.d0
      pmin         = 1.d0*bar
      pmax         = 1.d0*bar
      nHmin        = 4.d+19
      nHmax        = 4.d+19
      Mpl          = MEarth
      Rpl          = REarth
      gamma        = 7.0/5.0
      Tmin_atmos   = 0.0
      method_eqcond= 2
      UseDataBase  = .true.
      NewFullIt    = .true.
      NewBackIt    = 5
      NewBackFac   = 1.E+2
      NewFastLevel = 1
      NewPreMethod = 3
      Natmax       = 16
      Nseq         = 1

      !-------------------------------------------
      ! ***  change parameters via input file  ***
      !-------------------------------------------
      iarg = iargc()
      if (iarg<=0) then
        print*,"using default parameters"
        return
      endif  
      call getarg(1,ParamFile)
      open(unit=1,file=ParamFile,status='old')
      iline = 0
      dispol_set = 0
      do 
        read(1,'(A200)',end=100) line
        if (line(1:1).eq.'#') cycle           ! ignore comment lines
        if (len(trim(line))==0) cycle         ! ignore comment lines
        iline = iline+1
        print*,trim(line)
        if (iline.eq.1) then
          elements = ' '//trim(line)//' '     ! selection of elements
        else if (index(line,"! abund_pick")>0) then   
          read(line,*) abund_pick
          if (abund_pick==0) read(1,'(A200)',end=100) abund_file
        else if (index(line,"! pick_mfrac")>0) then   
          read(line,*) pick_mfrac
        else if (index(line,"! model_eqcond")>0) then   
          read(line,*) model_eqcond
        else if (index(line,"! remove_condensates")>0) then   
          read(line,*) remove_condensates
        else if (index(line,"! phyllosilicates")>0) then   
          read(line,*) phyllosilicates
        else if (index(line,"! metal_sulphates")>0) then   
          read(line,*) metal_sulphates
        else if (index(line,"! use_SiO")>0) then   
          read(line,*) use_SiO
        else if (index(line,"! disk_model")>0) then   
          read(line,*) disk_model
        else if (index(line,"! output_dispol")>0) then   
          read(line,*) output_dispol
        else if (index(line,"! model_dim")>0) then   
          read(line,*) model_dim
        else if (index(line,"! model_pconst")>0) then   
          read(line,*) model_pconst
        else if (index(line,"! model_struc")>0) then   
          read(line,*) model_struc
          if (model_struc>0) read(1,'(A200)') struc_file
        else if (index(line,"! Tmax")>0) then   
          read(line,*) Tmax
        else if (index(line,"! Tmin_atmos")>0) then 
          read(line,*) Tmin_atmos
        else if (index(line,"! Tmin")>0) then   
          read(line,*) Tmin
        else if (index(line,"! Tfast")>0) then   
          read(line,*) Tfast
        else if (index(line,"! verbose")>0) then   
          read(line,*) verbose
        else if (index(line,"! initchem_info")>0) then   
          read(line,*) initchem_info
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
        else if (index(line,"! auto_atmos")>0) then 
          read(line,*) auto_atmos
        else if (index(line,"! adapt_cond")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) adapt_file
          adapt_cond = .true.
          if (adapt_file=='.false.') adapt_cond = .false.
        else if (index(line,"! stop_after_init")>0) then 
          read(line,*) stop_after_init
        else if (index(line,"! Mplanet")>0) then 
          read(line,*) Mpl
          Mpl = Mpl*MEarth
        else if (index(line,"! Rplanet")>0) then 
          read(line,*) Rpl
          Rpl = Rpl*REarth
        else if (index(line,"! gamma")>0) then 
          read(line,*) gamma
        else if (index(line,"! Npoints")>0) then 
          read(line,*) Npoints
        else if (index(line,"! NewBackIt")>0) then 
          read(line,*) NewBackIt
        else if (index(line,"! NewBackFac")>0) then 
          read(line,*) NewBackFac
        else if (index(line,"! NewFullIt")>0) then 
          read(line,*) NewFullIt
        else if (index(line,"! NewFastLevel")>0) then 
          read(line,*) NewFastLevel
        else if (index(line,"! NewPreMethod")>0) then 
          read(line,*) NewPreMethod
        else if (index(line,"! Natmax")>0) then 
          read(line,*) Natmax
        else if (index(line,"! Nseq")>0) then 
          read(line,*) Nseq
          read(1,*) Tseq(1:Nseq)
        else if (index(line,"! useDatabase")>0) then 
          read(line,*) useDatabase
        else if (index(line,"! method_eqcond")>0) then 
          read(line,*) method_eqcond
        else if (index(line,"! dispol_file2")>0) then 
          i = index(line,"!")
          dispol_file(2) = line(1:i-1)
          dispol_set = 2
        else if (index(line,"! dispol_file3")>0) then 
          i = index(line,"!")
          dispol_file(3) = line(1:i-1)
          dispol_set = 3
        else if (index(line,"! dispol_file4")>0) then 
          i = index(line,"!")
          dispol_file(4) = line(1:i-1)
          dispol_set = 4
        else if (index(line,"! dispol_file")>0) then 
          i = index(line,"!")
          dispol_file(1) = line(1:i-1)
          dispol_set = 1
        else if (index(line,"! DustChem_file")>0) then 
          i = index(line,"!")
          read(line(1:i-1),*) DustChem_file
        else
          print*,"*** syntax error in "//trim(ParamFile)//":"
          print*,trim(line)
          stop
        endif  
      enddo  
 100  continue
      if (dispol_set>0.and.dispol_set<4) dispol_file(4)=""
      if (dispol_set>0.and.dispol_set<3) dispol_file(3)=""
      if (dispol_set>0.and.dispol_set<2) dispol_file(2)=""
      end
