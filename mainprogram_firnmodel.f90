!------------------------------------------------------
!- MAIN PROGRAM OF THE FIRN DENSIFICATION MODEL
!--------------
!- to compile: pgf90 -I/include mainprogram.f90 openNetCDF.f90 ini_model.f90 &
!-		time_loop.f90 output.f90 subprogram.f90 -o firnmodel.x &
!-		-L/lib -lnetcdf
!--------------
!- Cleaned up to become IMAU-FDM v1.0 (SL: 12-2014)
!- Adapted for ECMWF use (01-2012)
!- Made by: Stefan Ligtenberg (02-2010) 
!- based on the firn model of Michiel Helsen (2004-2007)
!------------------------------------------------------

    program mainprogram
    
    IMPLICIT NONE
    
    integer :: lon,lat,numTimes,numLats,numLons,numSteps,numPoints,numPointsSU,dtSnow
    integer :: ImpExp,kkinit,dtmodel,dtmodelImp,dtmodelExp,dtobs
    integer :: numberofrepeat,writeinprof,writeinspeed,writeindetail
    integer :: proflayers,detlayers
    integer :: NoR,BeginT,startasice,IceShelf
    integer :: i,j,a,b,nyears,nyearsSU
    integer,dimension(2) :: Points 
    integer, parameter :: kk=20000
    
    character*255 :: username,settingsfile,domain,fname_p1,ini_fname
    
    double precision :: dzmax,initdepth,th,maxpore,rho0max,rhoi,R,pi,detthick
    double precision :: tsav,acav,ffav,sumMelt
    
    double precision, dimension(2) :: Grid
    double precision, dimension(:), allocatable :: SnowMelt,PreTot,PreSol,PreLiq
    double precision, dimension(:), allocatable :: Sublim,SnowDrif,TempSurf,FF10m
    double precision, dimension(:,:), allocatable :: AveTsurf,AveAcc,AveWind,AveMelt
    double precision, dimension(:,:), allocatable :: ISM,LSM,Latitude,Longitude
      
    
    write (*,*) "------------------------------------"
    write (*,*) "----- FIRN DENSIFICATION MODEL -----"
    write (*,*) "------------------------------------"
    
    ! Get command line arguments:
    ! 1: ECMWF username (e.g. nmg)
    ! 2: Simulation number, corresponding to the line number in the IN-file.
    ! 3: Domain name, similar to the RACMO forcing (e.g. ANT27)
    ! 4: General part of output filename
    ! 5: Optional, name of the initialization file
    call getarg(1, username)
    call getarg(2, settingsfile)       
    call getarg(3, domain)
    call getarg(4, fname_p1)
    call getarg(5, ini_fname)
    
    print *, "TEST"
    
    ! Read in the model settings, input settings and constants
    call model_settings(dtSnow,nyears,nyearsSU,dtmodelImp,dtmodelExp,ImpExp,dtobs, &
        kkinit,NoR,startasice,beginT,writeinprof,writeinspeed,writeindetail, &
        proflayers,detlayers,detthick,dzmax,initdepth,th,Grid,settingsfile, &
        username,domain)
    call input_settings(numLons,numLats,numTimes,domain,username)
    call constants(maxpore,rho0max,rhoi,R,pi)

    ! Adjust matrices to correct dimensions according to the netCDF files	
    allocate(SnowMelt(numTimes))
    allocate(PreTot(numTimes))
    allocate(PreSol(numTimes))
    allocate(PreLiq(numTimes))
    allocate(Sublim(numTimes))
    allocate(TempSurf(numTimes))
    allocate(SnowDrif(numTimes))
    allocate(FF10m(numTimes))
    
    allocate(AveTsurf(numLons,numLats))
    allocate(LSM(numLons,numLats))
    allocate(ISM(numLons,numLats))
    allocate(Latitude(numLons,numLats))
    allocate(Longitude(numLons,numLats))
    allocate(AveAcc(numLons,numLats))
    allocate(AveWind(numLons,numLats))
    allocate(AveMelt(numLons,numLats))
    
    ! Calculate averages, T-amplitude and Melt-or-no-melt for all points
    call GetNetCDFAverage(AveTsurf,AveAcc,AveWind,AveMelt,LSM,numLats,numLons, &
        numTimes,nyears,Latitude,Longitude,ISM,username,domain)
    print *, "Read all averaged values"
    print *, " "
    
    
    ! Find corresponding grid point for given Lat and Lon	
    call FindGrid(Points,Grid,Latitude,Longitude,LSM,numLons,numLats)

    lon = Points(1)
    lat = Points(2)
    numberofrepeat = NoR
    IceShelf = ISM(lon,lat)

    write (*,*) "------ Point number: ", a 
    write (*,*) " Run for Lon: ",Grid(1)," and Lat: ",Grid(2)
    write (*,*) " Gridpoint: ",Points(1),", ",Points(2)
    write (*,*) " Number of spin up times: ",numberofrepeat
    write (*,*) " Grounded (0) or Floating (1) ice: ",IceShelf
    write (*,*) " Implicit (1) or Explicit (2) scheme: ",ImpExp
    write (*,*) "------------------------------------"
    write (*,*) " "
    
!	! Check if gridpoint has a positive annual accumulation rate
!	if (AveAcc(lon,lat).lt.0.) goto 101

    ! Get variables from the NetCDF files
    call GetVarNetCDF(SnowMelt,PreTot,PreSol,PreLiq,Sublim,SnowDrif,TempSurf, &
        FF10m,numTimes,lon,lat,username,domain,dtobs)
    print *, "Got all variables from the NetCDF files"
    print *, "---------------------------------------"

!	! Check if snowmelt is present: implicit or explicit temperature calculation
!	call CheckForMelt(SnowMelt,sumMelt,numTimes,nyears,dtmodel,dtmodelImp, &
!			dtmodelExp,ImpExp)
    
    if (ImpExp .eq. 2) then
        dtmodel = dtmodelExp
    else
        dtmodel = dtmodelImp
    endif
    
    numSteps = dtobs/dtmodel
    numPoints = numTimes*numSteps
    numPointsSU = numPoints*(REAL(nyearsSU)*365.25*24.*3600.)/(REAL(numTimes*dtobs))
    if (numPointsSU.gt.numPoints) numPointsSU = numPoints
    
    ! Read averages for the current grid point		
    call read_averages(AveTsurf,AveAcc,AveWind,AveMelt,tsav,acav,ffav, &
        numLons,numLats,lon,lat,nyears)
    if (acav.lt.0) acav = 0.1
    
    ! Call subprogram for spin-up and the time-integration			
    call subprogram(beginT,dtmodel,dtobs,ImpExp,numTimes,numSteps, &
        numPoints,numPointsSU,dtSnow,nyears,nyearsSU,kk,kkinit,initdepth,numberofrepeat,th,R, &
        pi,dzmax,rhoi,maxpore,writeinprof,writeinspeed,writeindetail, &
        proflayers,detlayers,detthick,acav,tsav,TempSurf,PreSol,PreLiq, &
        Sublim,SnowMelt,SnowDrif,FF10m,IceShelf,settingsfile,fname_p1,username, &
        domain,ini_fname)
    
101 continue

    end program

