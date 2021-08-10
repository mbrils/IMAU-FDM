mainprogram_firnmodel.f90                                                                           000644  000765  000024  00000012156 14104232652 016454  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         !------------------------------------------------------
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
		

	! Get variables from the NetCDF files
	call GetVarNetCDF(SnowMelt,PreTot,PreSol,PreLiq,Sublim,SnowDrif,TempSurf, &
			FF10m,numTimes,lon,lat,username,domain,dtobs)
	print *, "Got all variables from the NetCDF files"
	print *, "---------------------------------------"

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
		numPoints,numPointsSU,dtSnow,nyears,kk,kkinit,initdepth,numberofrepeat,th,R, &
		pi,dzmax,rhoi,maxpore,writeinprof,writeinspeed,writeindetail, &
		proflayers,detlayers,detthick,acav,tsav,TempSurf,PreSol,PreLiq, &
		Sublim,SnowMelt,SnowDrif,FF10m,IceShelf,settingsfile,fname_p1,username, &
		domain,ini_fname)
		
101	continue	

	
	end program
                                                                                                                                                                                                                                                                                                                                                                                                                  openNetCDF.f90                                                                                      000644  000765  000024  00000064742 14104232747 014043  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         !------------------------------------
!  SUBROUTINES TO OPEN NETCDF FILES
!------------------------------------
	
	subroutine GetNetCDFAverage(AveTsurf,AveAcc,AveWind,AveMelt,LSM,numLats, &
		numLons,numTimes,nyears,Latitude,Longitude,ISM,username,domain)

	use netcdf, only : nf90_open,nf90_inq_varid,nf90_close,nf90_get_var,nf90_noerr


	implicit none
	
	integer :: status,ncid(50),ID(50),numLats,numLons,numTimes,i,j,nyears
	double precision, dimension(numLons,numLats) :: AveTsurf,AveAcc,AveWind,AveSubl, &
		AveSnowDrif,AveMelt,LSM,ISM,Latitude,Longitude
		
	character*255 :: add,pad,username,domain

        
	if (domain .eq. "ANT27") then    
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/ANT27_averages/"	
	 add = "_ANT27_79-16_ave.nc"
	elseif (domain .eq. "XPEN055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/XPEN055_averages/"	
	 add = "_XPEN055_79-16_ave.nc"
	elseif (domain .eq. "FGRN11") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/data/input/era_files/averages/"	
	 add = "_FGRN11_60-79_ave.nc"
	elseif (domain .eq. "FGRN055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/data/input/era055_files/averages/"	
	 add = "_FGRN055_60-80_ave.nc"
	elseif (domain .eq. "PAT055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/PAT055_averages/"	
	 add = "_PAT055_79-12_ave.nc"
	elseif (domain .eq. "XDML055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/XDML055_averages/"	
	 add = "_XDML055_79-15_ave.nc"
	elseif (domain .eq. "ASE055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/ASE055_averages/"	
	 add = "_ASE055_79-15_ave.nc"
	elseif (domain .eq. "DMIS055") then
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/DMIS055_averages/"	
	 add = "_DMIS055_79-17_ave.nc"
	else
	 call handle_err(41,'no valid domain') 
 	endif
	
	if (domain .eq. "ANT27") then

	 status = nf90_open(trim(pad)//"../lsm_ANT27.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"LSM",ID(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"Lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"Lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! Open Ice Shelf Mask

	 status = nf90_open(trim(pad)//"../ism_ANT27.nc",0,ncid(6))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open6')
	 status = nf90_inq_varid(ncid(6),"ISM",ID(6))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	 status  = nf90_get_var(ncid(6),ID(6),ISM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_ism')


	 status = nf90_close(ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_close1')
	 status = nf90_close(ncid(6))
	 if(status /= nf90_noerr) call handle_err(status,'nf_close6')

	elseif (domain .eq. "XPEN055") then

	 status = nf90_open(trim(pad)//"../Height_latlon_XPEN055.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"mask2d",ID(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! Open Ice Shelf Mask

	 status = nf90_inq_varid(ncid(1),"iceshelves",ID(6))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	 status  = nf90_get_var(ncid(1),ID(6),ISM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_ism')

	 status = nf90_close(ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_close1')

	elseif (domain .eq. "FGRN11") then

	 status = nf90_open(trim(pad)//"../mask/FGRN11_Masks_wholedomain.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"icemask",ID(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid_lsm')	
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! No Ice Shelves in Greenland
	 ISM(:,:) = 0

	elseif (domain .eq. "FGRN055") then

	 status = nf90_open(trim(pad)//"../mask/FGRN055_Masks.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"Icemask_GR",ID(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid_lsm')	
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! No Ice Shelves in Greenland
	 ISM(:,:) = 0

	elseif (domain .eq. "PAT055") then

	 status = nf90_open(trim(pad)//"../lsm_PAT055.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"mask",ID(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! No Ice Shelves in Patagonia
	 ISM(:,:) = 0

	elseif (domain .eq. "XDML055") then

	 status = nf90_open(trim(pad)//"../lsm_XDML055.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"ism",ID(1))   !ice sheet mask....
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	! No Ice Shelves file yet
	 ISM(:,:) = 0

	elseif (domain .eq. "ASE055") then

	 status = nf90_open(trim(pad)//"../Masks_ASE055.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"LSM",ID(1))   !ice sheet mask....
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"ISM",ID(2))   !ice shelf mask....
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	 status  = nf90_get_var(ncid(1),ID(2),ISM,start=(/1,1/), &
		count=(/numLons,numLats/))

	elseif (domain .eq. "DMIS055") then

	 status = nf90_open(trim(pad)//"../Masks_DMIS055.nc",0,ncid(1))
	 if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	 status = nf90_inq_varid(ncid(1),"LSM",ID(1))   !ice sheet mask....
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')	
	 status = nf90_inq_varid(ncid(1),"ISM",ID(2))   !ice shelf mask....
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	 status = nf90_inq_varid(ncid(1),"lat",ID(11))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid11')
	 status = nf90_inq_varid(ncid(1),"lon",ID(12))
	 if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid12')
	 
	 status  = nf90_get_var(ncid(1),ID(1),LSM,start=(/1,1/), &
		count=(/numLons,numLats/))
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lsm')
	 status  = nf90_get_var(ncid(1),ID(11),Latitude)		
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lat')
	 status  = nf90_get_var(ncid(1),ID(12),Longitude)	
	 if(status /= nf90_noerr) call handle_err(status,'nf_get_var_lon')

	 status  = nf90_get_var(ncid(1),ID(2),ISM,start=(/1,1/), &
		count=(/numLons,numLats/))

	else

	 call handle_err(42,'no valid domain')

	endif 

	print *, trim(pad)//"snowmelt"//trim(add)

	status = nf90_open(trim(pad)//"snowmelt"//trim(add),0,ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	status = nf90_open(trim(pad)//"precip"//trim(add),0,ncid(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_open2')
	status = nf90_open(trim(pad)//"ff10m"//trim(add),0,ncid(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_open3')
	status = nf90_open(trim(pad)//"tskin"//trim(add),0,ncid(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_open4')
	status = nf90_open(trim(pad)//"evap"//trim(add),0,ncid(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_open5')
	status = nf90_open(trim(pad)//"sndiv"//trim(add),0,ncid(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_open7')

	status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	status = nf90_inq_varid(ncid(2),"precip",ID(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	status = nf90_inq_varid(ncid(3),"ff10m",ID(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	status = nf90_inq_varid(ncid(4),"tskin",ID(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	status = nf90_inq_varid(ncid(5),"evap",ID(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')
	status = nf90_inq_varid(ncid(7),"sndiv",ID(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid')

	status  = nf90_get_var(ncid(1),ID(1),AveMelt,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var_avemelt')
	status  = nf90_get_var(ncid(2),ID(2),AveAcc,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var_aveacc')
	status  = nf90_get_var(ncid(3),ID(3),AveWind,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_varavewind')
	status  = nf90_get_var(ncid(4),ID(4),AveTsurf,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_varavetsurf')
	status  = nf90_get_var(ncid(5),ID(5),AveSubl,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_varavesubl')
	status  = nf90_get_var(ncid(7),ID(7),AveSnowDrif,start=(/1,1,1,1/), &
		count=(/numLons,numLats,1,1/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_varavesndiv')

	! Close all netCDF files	

	status = nf90_close(ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_close1')
	status = nf90_close(ncid(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_close2')
	status = nf90_close(ncid(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_close3')
	status = nf90_close(ncid(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_close4')
	status = nf90_close(ncid(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_close5')
	status = nf90_close(ncid(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_close7')

        ! Convert units from [mm w.e./s] to [mm w.e./yr]
        do i=1,numLons
         do j=1,numLats
            AveAcc(i,j) = (AveAcc(i,j)+AveSubl(i,j)-AveSnowDrif(i,j)) &
                * (365.*24.*3600.)
            AveMelt(i,j) = AveMelt(i,j) * (365.*24.*3600.)
          end do
        end do
	
	end subroutine


!----------------------
	subroutine GetVarNetCDF(SnowMelt,PreTot,PreSol,PreLiq,Sublim,SnowDrif,TempSurf, &
		FF10m,numTimes,lon,lat,username,domain,dtobs)

	use netcdf, only : nf90_open,nf90_close,nf90_inq_varid,nf90_get_var,nf90_noerr

	implicit none
	
	integer :: status,k,ncid(50),numTimes,lat,lon,ID(50),dtobs
	double precision :: wegPsol,wegPliq,wegPtot,wegMelt
	double precision, dimension(numTimes) :: SnowMelt,PreTot,PreSol,PreLiq,Sublim,TempSurf, &
		SnowDrif,FF10m

	integer :: latfile,lonfile,fnumb_i
	character*255 :: add,pad,fnumb,username,domain
	
	if (domain .eq. "ANT27") then    
	 latfile = mod(lat,15)
	 if (latfile.eq.0) latfile = 15
	 fnumb_i = floor((real(lat)-0.001)/15.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 lonfile = lon

	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/ANT27_files/"	
	 add = "_ANT27_79-16_p"//trim(fnumb)//".nc"
	 
	elseif (domain .eq. "XPEN055") then
	 lonfile = mod(lon,5)
	 if (lonfile.eq.0) lonfile = 5
	 fnumb_i = floor((real(lon)-0.001)/5.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 latfile = lat
	
	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/XPEN055_files/"	
	 add = "_XPEN055_79-16_p"//trim(fnumb)//".nc"
	
	elseif (domain .eq. "FGRN11") then
	 lonfile = mod(lon,5)
	 if (lonfile.eq.0) lonfile = 5
	 fnumb_i = floor((real(lon)-0.001)/5.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 latfile = lat

	 pad = "/scratch/ms/nl/"//trim(username)//"/data/input/era_files/files/"	
	 add = "_FGRN11_60-16_p"//trim(fnumb)//".nc"

	elseif (domain .eq. "FGRN055") then
	 lonfile = mod(lon,6)
	 if (lonfile.eq.0) lonfile = 6
	 fnumb_i = floor((real(lon)-0.001)/6.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 latfile = lat

	 pad = "/scratch/ms/nl/"//trim(username)//"/data/input/era055_files/files/"	
	 add = "_FGRN055_57-20_p"//trim(fnumb)//".nc"
	
	elseif (domain .eq. "PAT055") then
	 latfile = mod(lat,4)
	 if (latfile.eq.0) latfile = 4
	 fnumb_i = floor((real(lat)-0.001)/4.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 lonfile = lon

	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/PAT055_files/"	
	 add = "_PAT055_79-12_p"//trim(fnumb)//".nc"

	elseif (domain .eq. "XDML055") then
	 latfile = mod(lat,5)
	 if (latfile.eq.0) latfile = 5
	 fnumb_i = floor((real(lat)-0.001)/5.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 lonfile = lon

	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/XDML055_files/"	
	 add = "_XDML055_79-15_p"//trim(fnumb)//".nc"

	elseif (domain .eq. "ASE055") then
	 latfile = mod(lat,5)
	 if (latfile.eq.0) latfile = 5
	 fnumb_i = floor((real(lat)-0.001)/5.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 lonfile = lon

	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/ASE055_files/"	
	 add = "_ASE055_79-15_p"//trim(fnumb)//".nc"

	elseif (domain .eq. "DMIS055") then
	 lonfile = mod(lon,6)
	 if (lonfile.eq.0) lonfile = 6
	 fnumb_i = floor((real(lon)-0.001)/6.)+1
	 if (fnumb_i.le.9) write(fnumb,'(I1)') fnumb_i
	 if (fnumb_i.ge.10) write(fnumb,'(I2)') fnumb_i

	 latfile = lat

	 pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/DMIS055_files/"	
	 add = "_DMIS055_79-17_p"//trim(fnumb)//".nc"
	
	else
	 call handle_err(43,'no valid domain') 
 	endif
	
	print *, numTimes
	print *, lat,latfile,lon,lonfile,fnumb,numTimes
	print *, pad,add

	! Open the snowmelt netCDF file
	status = nf90_open(trim(pad)//"snowmelt"//trim(add),0,ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_open1')
	status = nf90_open(trim(pad)//"precip"//trim(add),0,ncid(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_open2')
	status = nf90_open(trim(pad)//"snowfall"//trim(add),0,ncid(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_open3')
	status = nf90_open(trim(pad)//"evap"//trim(add),0,ncid(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_open4')
	status = nf90_open(trim(pad)//"tskin"//trim(add),0,ncid(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_open5')
	status = nf90_open(trim(pad)//"sndiv"//trim(add),0,ncid(6))
	if(status /= nf90_noerr) call handle_err(status,'nf_open6')
	status = nf90_open(trim(pad)//"ff10m"//trim(add),0,ncid(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_open7')
	
	!Get ID of the variables in the NetCDF files
	status = nf90_inq_varid(ncid(1),"snowmelt",ID(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid1')
	status = nf90_inq_varid(ncid(2),"precip",ID(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid2')
	status = nf90_inq_varid(ncid(3),"snowfall",ID(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid3')
	status = nf90_inq_varid(ncid(4),"evap",ID(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid4')
	status = nf90_inq_varid(ncid(5),"tskin",ID(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid5')
	status = nf90_inq_varid(ncid(6),"sndiv",ID(6))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid6')
	status = nf90_inq_varid(ncid(7),"ff10m",ID(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid7')

	! Get all variables from the netCDF files
	status  = nf90_get_var(ncid(1),ID(1),SnowMelt,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var1')
	print *, "Read snowmelt..."
	status  = nf90_get_var(ncid(2),ID(2),PreTot,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var2')
	print *, "Read precipitation..."
	status  = nf90_get_var(ncid(3),ID(3),PreSol,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var3')
	print *, "Read snowfall..."
	status  = nf90_get_var(ncid(4),ID(4),Sublim,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var4')
	print *, "Read sublimation..."
	status  = nf90_get_var(ncid(5),ID(5),TempSurf,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var5')
	print *, "Read skin temperature..."
	status  = nf90_get_var(ncid(6),ID(6),SnowDrif,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var6')
	print *, "Read snow drift..."
	status  = nf90_get_var(ncid(7),ID(7),FF10m,start=(/lonfile,latfile,1,1/), &
		count=(/1,1,1,numTimes/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var7')
	print *, "Read wind speed..."
	
	! Close all netCDF files	
	status = nf90_close(ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_close1')
	status = nf90_close(ncid(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_close2')
	status = nf90_close(ncid(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_close3')
	status = nf90_close(ncid(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_close4')
	status = nf90_close(ncid(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_close5')
	status = nf90_close(ncid(6))
	if(status /= nf90_noerr) call handle_err(status,'nf_close6')
	status = nf90_close(ncid(7))
	if(status /= nf90_noerr) call handle_err(status,'nf_close7')

	wegPsol = 0.
	wegPliq = 0.
	wegPtot = 0.
	wegMelt = 0.

        do k=1,numTimes
          SnowMelt(k) = SnowMelt(k) * dtobs
          PreTot(k) = PreTot(k) * dtobs
          PreSol(k) = PreSol(k) * dtobs
          Sublim(k) = Sublim(k) * dtobs
          SnowDrif(k) = SnowDrif(k) * dtobs
          PreLiq(k) = PreLiq(k) * dtobs
        end do

	do k=1,numTimes
	 if (SnowMelt(k) .lt.1.e-04) then
	  wegMelt = wegMelt + SnowMelt(k)
	  SnowMelt(k) = 0.	 
	 endif
	
         if (PreSol(k) .lt. 1.e-04) then
      	  wegPsol = wegPsol + PreSol(k)
      	  PreSol(k) = 0.
         endif
	 
         if (PreTot(k) .lt. 1.e-04) then
      	  wegPtot = wegPtot + PreTot(k)
      	  PreTot(k) = 0.
         endif
	 
         if (TempSurf(k) .gt. 267.) then
          PreLiq(k) = PreTot(k)-PreSol(k)
          if (PreLiq(k) .lt. 1.e-04) then
       	   wegPliq = wegPliq + PreLiq(k)
      	   PreLiq(k) = 0.
           PreSol(k) = PreTot(k)
          endif
         else 
          PreSol(k) = PreTot(k)
         endif
        end do
	
	print *, "wegPsol: ",wegPsol
	print *, "wegPliq: ",wegPliq
	print *, "wegPtot: ",wegPtot
	print *, "wegMelt: ",wegMelt



	end subroutine
	
	
!---------------------------------
	subroutine CheckForMelt(SnowMelt,sumMelt,numTimes,nyears,dtmodel,dtmodelImp, &
		dtmodelExp,ImpExp)

	implicit none
	
	integer :: k,numTimes,nyears,dtmodel,dtmodelImp,dtmodelExp,ImpExp
	double precision :: sumMelt,sumMeltY
	double precision, dimension(numTimes) :: SnowMelt


  	sumMelt = 0.
	do k = 1, numTimes
         if (SnowMelt(k) .lt. 1e-04) SnowMelt(k) = 0.
         sumMelt = sumMelt + SnowMelt(k)
        end do
	
	sumMeltY = sumMelt/nyears	
	
	if (sumMeltY .lt. 0.5) then
	  ImpExp = 1 !=implicit
	  dtmodel = dtmodelImp
	  print *, "No melt, implicit scheme!    Melt sum = ",sumMelt," mm yr-1"
	  SnowMelt = 0.
	else
	  ImpExp = 0 !=explicit
	  dtmodel = dtmodelExp
	  print *, "Melt, explicit scheme!      Melt sum = ",sumMelt," mm yr-1"
	endif
		
	end subroutine


!---------------------------------
	subroutine ini_fromfile(kk,kUL,Rho,M,T,Depth,Mlwc,DZ,username,domain, &
		ini_fname,settingsfile)

	use netcdf, only : nf90_open,nf90_close,nf90_inq_varid,nf90_inq_dimid, &
		nf90_inquire_dimension,nf90_get_var,nf90_noerr
	
	implicit none
	
	integer :: kk,kUL
	integer :: k,status,ncid(50),ID(50),LayerID
	
	double precision, dimension(kk) :: Rho,M,T,Depth,Mlwc,DZ 
	
	character*255 :: charac_kUL,fname,pad,username,domain,ini_fname,settingsfile
	
	pad = "/scratch/ms/nl/"//trim(username)//"/FM_Data/INPUT/ini_files/"//trim(domain)//"/"	
	fname = trim(ini_fname)//"_ini_"//trim(settingsfile)//".nc"
	
	print *, trim(pad)//trim(fname)

	! Open the snowmelt netCDF file
	status = nf90_open(trim(pad)//trim(fname),0,ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_open1')

	!Get ID of the variables in the NetCDF files
	status = nf90_inq_varid(ncid(1),"dens",ID(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid1')
	status = nf90_inq_varid(ncid(1),"temp",ID(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid2')
	status = nf90_inq_varid(ncid(1),"mass",ID(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid3')
	status = nf90_inq_varid(ncid(1),"depth",ID(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid4')
	status = nf90_inq_varid(ncid(1),"lwc",ID(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_varid5')

	! Get dimension of the array
	status = nf90_inq_dimid(ncid(1),"layer",LayerID)
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_dimid1')
	status = nf90_inquire_dimension(ncid(1),LayerID,len=kUL)
	if(status /= nf90_noerr) call handle_err(status,'nf_inq_dim1')
	
	! Get all variables from the netCDF files
	status  = nf90_get_var(ncid(1),ID(1),Rho(1:kUL),start=(/1/),count=(/kUL/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var1')
	status  = nf90_get_var(ncid(1),ID(2),T(1:kUL),start=(/1/),count=(/kUL/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var2')
	status  = nf90_get_var(ncid(1),ID(3),M(1:kUL),start=(/1/),count=(/kUL/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var3')
	status  = nf90_get_var(ncid(1),ID(4),Depth(1:kUL),start=(/1/),count=(/kUL/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var4')
	status  = nf90_get_var(ncid(1),ID(5),Mlwc(1:kUL),start=(/1/),count=(/kUL/))
	if(status /= nf90_noerr) call handle_err(status,'nf_get_var5')
	
	! Close all netCDF files	
	status = nf90_close(ncid(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_close1')

	DZ(kUL) = Depth(kUL) * 2.	
	do k = kUL-1,1,-1
	  DZ(k) = (Depth(k) - Depth(k+1) - 0.5*DZ(k+1)) * 2.
	end do
	
	end subroutine
	
	
!---------------------------------	
	subroutine handle_err(stat,msg)
	!--------------------------------------------------------------------------
	! error stop for netCDF    
	!--------------------------------------------------------------------------
	use netcdf, only : nf90_noerr, nf90_strerror
	
	implicit none

	integer,intent(in) :: stat       
	character(len=*),intent(in) :: msg

	if(stat /= nf90_noerr) then
  	  print*, 'netCDF error (',msg,'): ', nf90_strerror(stat)
	  stop
	endif

	end subroutine handle_err

                              output.f90                                                                                          000644  000765  000024  00000075524 14104233010 013436  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         !------------------------------------
!  SUBROUTINES TO OUTPUT THE RESULTS TO NETCDF
!------------------------------------
!----------------------------
!----------------------------
	subroutine write_initial(kk,kUL,Rho,M,T,Depth,Mlwc,Year,settingsfile, &
		fname_p1,username)
	
	use netcdf, only: nf90_create,nf90_def_dim,nf90_def_var,nf90_real,nf90_noerr, &
		nf90_enddef,nf90_put_var,nf90_close
	
	implicit none
	
	integer :: status,ncid,ID,varID(10),kk,kkinit,kUL
	integer :: nlayers,skUL,ekUL
	double precision, dimension(kk) ::  M,T,Depth,Mlwc
	double precision, dimension(kk) :: Rho,Year
	character*255 :: pad,settingsfile,fname_p1,username
	
	pad = "/scratch/ms/nl/"//trim(username)//"/data/output/era055/"

	! CREATE NETCDF FILES
	
	status = nf90_create(trim(pad)//trim(fname_p1)//"_ini_"//trim(settingsfile)// &
		".nc",0,ncid)
	if(status /= nf90_noerr) call handle_err(status,'nf_create')


	! DEFINE DIMENSIONS
	status = nf90_def_dim(ncid,"layer",kUL,ID)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim1')
	
	
	! DEFINE VARIABLES
	status = nf90_def_var(ncid,"dens",nf90_real,ID,varID(1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
	status = nf90_def_var(ncid,"temp",nf90_real,ID,varID(2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var2')
	status = nf90_def_var(ncid,"mass",nf90_real,ID,varID(3))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var3')
	status = nf90_def_var(ncid,"depth",nf90_real,ID,varID(4))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var4')
	status = nf90_def_var(ncid,"lwc",nf90_real,ID,varID(5))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var5')
	status = nf90_def_var(ncid,"year",nf90_real,ID,varID(6))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var6')	
	
	! END OF DEFINING FILES
	status = nf90_enddef(ncid)
	if(status /= nf90_noerr) call handle_err(status,'nf_enddef')
	
	! PUT VARIABLES
	status = nf90_put_var(ncid,varID(1),Rho(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var1')	
	status = nf90_put_var(ncid,varID(2),T(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var2')
	status = nf90_put_var(ncid,varID(3),M(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var3')
	status = nf90_put_var(ncid,varID(4),Depth(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var4')
	status = nf90_put_var(ncid,varID(5),Mlwc(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var5')
	
!	print *, Year
	
	status = nf90_put_var(ncid,varID(6),Year(1:kUL))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var6')
	
	! CLOSE NETCDF FILE
	status = nf90_close(ncid)
	if(status /= nf90_noerr) call handle_err(status,'nf_close')
	
	end subroutine

!---------------------------------
    subroutine save_out_restart(time,kk,kUL,Rho,DenRho,M,T,Depth,Mlwc,Year,Refreeze,DZ,settingsfile, &
        fname_p1,username)
    
    use netcdf, only: nf90_create,nf90_def_dim,nf90_def_var,nf90_real,nf90_int,nf90_noerr, &
        nf90_enddef,nf90_put_var,nf90_close
    
    implicit none
    
    integer :: status,ncid,dimID(2),varID(11),kk,kkinit,kUL
    integer :: nlayers,time
    double precision, dimension(kk) ::  M,T,Depth,Mlwc,DenRho
    double precision, dimension(kk) :: Rho,Year,Refreeze,DZ
    character*255 :: pad,settingsfile,fname_p1,username
    
    pad = "/scratch/ms/nl/"//trim(username)//"/data/output/restart/"
    
    ! CREATE NETCDF FILES
    
    status = nf90_create(trim(pad)//trim(fname_p1)//"_restart_"//trim(settingsfile)// &
        ".nc",0,ncid)
    if(status /= nf90_noerr) call handle_err(status,'nf_create')
    
    
    ! DEFINE DIMENSIONS
    status = nf90_def_dim(ncid,"layer",kk,dimID(1))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_dim1')
    status = nf90_def_dim(ncid,"constant",1,dimID(2))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_dim2')
    
    
    ! DEFINE VARIABLES
    status = nf90_def_var(ncid,"dens",nf90_real,dimID(1),varID(1))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
    status = nf90_def_var(ncid,"temp",nf90_real,dimID(1),varID(2))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var2')
    status = nf90_def_var(ncid,"mass",nf90_real,dimID(1),varID(3))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var3')
    status = nf90_def_var(ncid,"depth",nf90_real,dimID(1),varID(4))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var4')
    status = nf90_def_var(ncid,"lwc",nf90_real,dimID(1),varID(5))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var5')
    status = nf90_def_var(ncid,"year",nf90_real,dimID(1),varID(6))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var6')
    status = nf90_def_var(ncid,"refreeze",nf90_real,dimID(1),varID(7))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var7')
    status = nf90_def_var(ncid,"dz",nf90_real,dimID(1),varID(8))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var8')
    status = nf90_def_var(ncid,"drho",nf90_real,dimID(1),varID(9))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var9')
    
    status = nf90_def_var(ncid,"kUL",nf90_int,dimID(2),varID(10))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var10')
    status = nf90_def_var(ncid,"time",nf90_int,dimID(2),varID(11))
    if(status /= nf90_noerr) call handle_err(status,'nf_def_var11')
    
    ! END OF DEFINING FILES
    status = nf90_enddef(ncid)
    if(status /= nf90_noerr) call handle_err(status,'nf_enddef')
    
    ! PUT VARIABLES
    status = nf90_put_var(ncid,varID(1),Rho)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var1')
    status = nf90_put_var(ncid,varID(2),T)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var2')
    status = nf90_put_var(ncid,varID(3),M)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var3')
    status = nf90_put_var(ncid,varID(4),Depth)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var4')
    status = nf90_put_var(ncid,varID(5),Mlwc)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var5')
    status = nf90_put_var(ncid,varID(6),Year)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var6')
    status = nf90_put_var(ncid,varID(7),Refreeze)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var7')
    status = nf90_put_var(ncid,varID(8),DZ)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var8')
    status = nf90_put_var(ncid,varID(9),DenRho)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var9')
    
    status = nf90_put_var(ncid,varID(10),kUL)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var10')
    status = nf90_put_var(ncid,varID(11),time)
    if(status /= nf90_noerr) call handle_err(status,'nf_put_var11')
    
    
    ! CLOSE NETCDF FILE
    status = nf90_close(ncid)
    if(status /= nf90_noerr) call handle_err(status,'nf_close')
    
    end subroutine



!-----------------------
	subroutine to_out_1D(time,numOutputSpeed,zs,Totvice,Totvfc,Totvacc,Totvsub, &
		Totvsnd,Totvmelt,Totvbouy,Totv,TotRunoff,FirnAir,TotLwc, &
		Totrefreeze,TotRain,TotSurfmelt,TotSolIn,IceMass,Rho0out, &
		out_1D,outputSpeed)
	
	implicit none
	
	integer :: time,plek,numOutputSpeed,outputSpeed
	
  	double precision :: zs
  	double precision :: Totvice,Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy,Totv
  	double precision :: TotRunoff,FirnAir,TotLwc
  	double precision :: Totrefreeze,TotRain,TotSurfmelt,TotSolIn,IceMass,Rho0out 

  	real,dimension((outputSpeed+50),18) :: out_1D 
	
	plek = time/numOutputSpeed
	! zs, vice, vacc, vfc, vmelt, vbouy, vsub, vsnd, vtotal, Runoff, FirnAir, TotLwc, refreeze, rain, surfmelt, solin, icemass
	out_1D(plek,1) = REAL(zs)
	out_1D(plek,2) = REAL(Totvice)
	out_1D(plek,3) = REAL(Totvacc)
	out_1D(plek,4) = REAL(Totvfc)
	out_1D(plek,5) = REAL(Totvmelt)
	out_1D(plek,6) = REAL(Totvbouy)
	out_1D(plek,7) = REAL(Totvsub)
	out_1D(plek,8) = REAL(Totvsnd)
	out_1D(plek,9) = REAL(Totv)
	out_1D(plek,10) = REAL(TotRunoff)
	out_1D(plek,11) = REAL(FirnAir)
	out_1D(plek,12) = REAL(TotLwc)
	out_1D(plek,13) = REAL(TotRefreeze)
	out_1D(plek,14) = REAL(TotRain)
	out_1D(plek,15) = REAL(TotSurfmelt)
	out_1D(plek,16) = REAL(TotSolIn)
	out_1D(plek,17) = REAL(IceMass)
	out_1D(plek,18) = REAL(Rho0out)
	
    	Totvice = 0.
    	Totvacc = 0.
    	Totvsub = 0.
    	Totvsnd = 0.
    	Totvmelt = 0.
    	Totvfc = 0.
    	Totvbouy = 0.
    	Totv = 0.
    	TotRunoff = 0.
    	TotRefreeze = 0. !####################### modification
    	TotRain = 0. !####################### modification
    	TotSurfmelt = 0. !####################### modification
    	TotSolIn = 0. !####################### modification
	Rho0out = 0.
	
	
	end subroutine
	
	


!-----------------------
	subroutine to_out_2D(kk,kUL,time,dtmodel,numOutputProf,outputProf,proflayers, &
		Rho,T,Mlwc,Depth,DenRho,Year,out_2D_dens,out_2D_temp,out_2D_lwc, &
		out_2D_depth,out_2D_dRho,out_2D_year)

	implicit none
	
	integer :: kk,kUL,plek,skUL,ekUL
	integer :: time,dtmodel,numOutputProf,outputProf,proflayers
	double precision :: factor3
	
	double precision, dimension(kk) :: Rho,T,Mlwc,Depth,DenRho,Year
	real, dimension(outputProf+50,proflayers) :: out_2D_dens,out_2D_temp
	real, dimension(outputProf+50,proflayers) :: out_2D_lwc,out_2D_depth
	real, dimension(outputProf+50,proflayers) :: out_2D_dRho,out_2D_year
		
	plek = time/numOutputProf
        skUL = 1    !kUL
        ekUL = kUL
        if (kUL .gt. proflayers) then 
          skUL = kUL - proflayers + 1
          ekUL = proflayers
        endif
	
	factor3 = 3600. * 24. * 365. / dtmodel / numOutputProf
	DenRho(:) = DenRho(:) * factor3
	
	out_2D_dens(plek,1:ekUL) = REAL(Rho(skUL:kUL))
	out_2D_temp(plek,1:ekUL) = REAL(T(skUL:kUL))
	out_2D_lwc(plek,1:ekUL) = REAL(Mlwc(skUL:kUL))
	out_2D_depth(plek,1:ekUL) = REAL(Depth(skUL:kUL))
	out_2D_dRho(plek,1:ekUL) = REAL(DenRho(skUL:kUL))
	out_2D_year(plek,1:ekUL) = REAL(Year(skUL:kUL))
	
	DenRho(:) = 0.
	
	end subroutine



!-----------------------
	subroutine to_out_2Ddetail(kk,kUL,time,detlayers,detthick,numOutputDetail, &
		outputDetail,Rho,T,Mlwc,Refreeze,DZ,out_2D_det_dens, &
		out_2D_det_temp,out_2D_det_lwc,out_2D_det_refreeze)
	
	implicit none
	
	integer :: kk,kUL,plek,time,detlayers
	integer :: ind_orig,ind_int,numOutputDetail,outputDetail
	
	double precision :: dist,part,detthick,refreeze_ts
	
	double precision, dimension(detlayers) :: IntRho,IntT,IntMlwc,IntRefreeze
	double precision, dimension(kk) :: Rho,T,Mlwc,Refreeze,DZ,DZ_mod
		
	real, dimension(outputDetail+50,detlayers) :: out_2D_det_dens,out_2D_det_temp
	real, dimension(outputDetail+50,detlayers) :: out_2D_det_lwc,out_2D_det_refreeze	
		
	IntRho(:) = 0.
	IntT(:) = 0.
	IntMlwc(:) = 0.
	IntRefreeze(:) = 0.

	DZ_mod = DZ
	dist = 0.
	ind_orig = kUL
	ind_int = 1
	do while(ind_int .le. detlayers)
	  if ((dist + DZ_mod(ind_orig)) .lt. (detthick * ind_int)) then
	    IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
	    IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
	    IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
	    dist = dist + DZ_mod(ind_orig)
	    ind_orig = ind_orig - 1
	  else if ((dist + DZ_mod(ind_orig)) .eq. (detthick * ind_int)) then
	    IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (DZ_mod(ind_orig) / detthick)
	    IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (DZ_mod(ind_orig) / detthick)
	    IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
	    dist = dist + DZ_mod(ind_orig)
	    ind_orig = ind_orig - 1
	    ind_int = ind_int + 1
	  else
	    part = (detthick * ind_int) - dist
	    IntRho(ind_int) = IntRho(ind_int) + Rho(ind_orig) * (part / detthick)
	    IntT(ind_int) = IntT(ind_int) + T(ind_orig) * (part / detthick)
	    IntMlwc(ind_int) = IntMlwc(ind_int) + Mlwc(ind_orig) * (part / DZ(ind_orig))		
	    dist = dist + part
	    DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
	    ind_int = ind_int + 1
	  end if
	end do


	refreeze_ts = SUM(Refreeze)
	if (refreeze_ts .gt. 1e-05) then
	  DZ_mod = DZ
	  dist = 0.
	  ind_orig = kUL
	  ind_int = 1
	  do while(ind_int .le. detlayers)
	    if ((dist + DZ_mod(ind_orig)) .lt. (detthick * ind_int)) then
	      IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
	      dist = dist + DZ_mod(ind_orig)
	      ind_orig = ind_orig - 1
	    else if ((dist + DZ_mod(ind_orig)) .eq. (detthick * ind_int)) then
	      IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (DZ_mod(ind_orig) / DZ(ind_orig))
	      dist = dist + DZ_mod(ind_orig)
	      ind_orig = ind_orig - 1
	      ind_int = ind_int + 1
	    else
	      part = (detthick * ind_int) - dist
	      IntRefreeze(ind_int) = Intrefreeze(ind_int) + Refreeze(ind_orig) * (part / DZ(ind_orig))
	      dist = dist + part
	      DZ_mod(ind_orig) = DZ_mod(ind_orig) - part
   	      ind_int = ind_int + 1
	    end if
	  end do
	end if

	plek = time/numOutputDetail
	
	out_2D_det_dens(plek,:) = REAL(IntRho)
	out_2D_det_temp(plek,:) = REAL(IntT)
	out_2D_det_lwc(plek,:) = REAL(IntMlwc)
	out_2D_det_refreeze(plek,:) = REAL(IntRefreeze)

	Refreeze(:) = 0.

	
	end subroutine
	


!----------------------
	subroutine save_out_1D(outputSpeed,settingsfile,fname_p1,username,out_1D)
	
	use netcdf, only: nf90_create,nf90_noerr,nf90_def_dim,nf90_def_var, &
			  nf90_enddef,nf90_put_var,nf90_real,nf90_int,nf90_unlimited, &
			  nf90_close,nf90_put_att
	
	implicit none

	integer :: outputSpeed
	integer :: status,ncid(50),IDs(50,5),varID(50,20)
	real,dimension((outputSpeed+50),18) :: out_1D
	character*255 :: pad,settingsfile,fname_p1,username
	
	pad = "/scratch/ms/nl/"//trim(username)//"/data/output/era055/"
	
	ncid(:) = 0
	IDs(:,:) = 0
	varID(:,:) = 0

	! CREATE NETCDF FILES
	status = nf90_create(trim(pad)//trim(fname_p1)//"_1D_"//trim(settingsfile)// &
		".nc",0,ncid(32))

	! DEFINE DIMENSIONS
	status = nf90_def_dim(ncid(32),"time",outputSpeed+50,IDs(32,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim1')
	
	! DEFINE VARIABLES
	status = nf90_def_var(ncid(32),"zs",nf90_real,(/IDs(32,1)/),varID(32,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
	status = nf90_def_var(ncid(32),"vice",nf90_real,(/IDs(32,1)/),varID(32,2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var2')
	status = nf90_def_var(ncid(32),"vacc",nf90_real,(/IDs(32,1)/),varID(32,3))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var3')
	status = nf90_def_var(ncid(32),"vfc",nf90_real,(/IDs(32,1)/),varID(32,4))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var4')
	status = nf90_def_var(ncid(32),"vmelt",nf90_real,(/IDs(32,1)/),varID(32,5))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var5')
	status = nf90_def_var(ncid(32),"vbouy",nf90_real,(/IDs(32,1)/),varID(32,6))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var6')
	status = nf90_def_var(ncid(32),"vsub",nf90_real,(/IDs(32,1)/),varID(32,7))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var7')
	status = nf90_def_var(ncid(32),"vsnd",nf90_real,(/IDs(32,1)/),varID(32,8))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var8')
	status = nf90_def_var(ncid(32),"vtotal",nf90_real,(/IDs(32,1)/),varID(32,9))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var9')
	status = nf90_def_var(ncid(32),"Runoff",nf90_real,(/IDs(32,1)/),varID(32,10))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var10')
	status = nf90_def_var(ncid(32),"FirnAir",nf90_real,(/IDs(32,1)/),varID(32,11))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var11')
	status = nf90_def_var(ncid(32),"TotLwc",nf90_real,(/IDs(32,1)/),varID(32,12))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var12')
	status = nf90_def_var(ncid(32),"refreeze",nf90_real,(/IDs(32,1)/),varID(32,13)) 
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var13') 
	status = nf90_def_var(ncid(32),"rain",nf90_real,(/IDs(32,1)/),varID(32,14))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var14')
	status = nf90_def_var(ncid(32),"surfmelt",nf90_real,(/IDs(32,1)/),varID(32,15)) 
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var15') 
	status = nf90_def_var(ncid(32),"solin",nf90_real,(/IDs(32,1)/),varID(32,16))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var16') 
	status = nf90_def_var(ncid(32),"icemass",nf90_real,(/IDs(32,1)/),varID(32,17))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var17') 
	status = nf90_def_var(ncid(32),"Rho0",nf90_real,(/IDs(32,1)/),varID(32,18))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var18')
	
	! DEFINE ATTRIBUTES (unit could also be defined here)
	status = nf90_put_att(ncid(32),varID(32,1),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_1')
	status = nf90_put_att(ncid(32),varID(32,2),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_2')
	status = nf90_put_att(ncid(32),varID(32,3),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_3')
	status = nf90_put_att(ncid(32),varID(32,4),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_4')
	status = nf90_put_att(ncid(32),varID(32,5),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_5')
	status = nf90_put_att(ncid(32),varID(32,6),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_6')
	status = nf90_put_att(ncid(32),varID(32,7),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_7')
	status = nf90_put_att(ncid(32),varID(32,8),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_8')
	status = nf90_put_att(ncid(32),varID(32,9),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_9')
	status = nf90_put_att(ncid(32),varID(32,10),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_10')
	status = nf90_put_att(ncid(32),varID(32,11),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_11')
	status = nf90_put_att(ncid(32),varID(32,12),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_12')
	status = nf90_put_att(ncid(32),varID(32,13),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_13')
	status = nf90_put_att(ncid(32),varID(32,14),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_14')
	status = nf90_put_att(ncid(32),varID(32,15),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_15')
	status = nf90_put_att(ncid(32),varID(32,16),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_16')
	status = nf90_put_att(ncid(32),varID(32,17),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_17')	
	status = nf90_put_att(ncid(32),varID(32,18),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_18')
	
	! END OF DEFINING FILES
	status = nf90_enddef(ncid(32))
	if(status /= nf90_noerr) call handle_err(status,'nf_enddef')

	! SAVE DATA
	status = nf90_put_var(ncid(32),varID(32,1),out_1D(:,1), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed1')	
	status = nf90_put_var(ncid(32),varID(32,2),out_1D(:,2), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var2')
	status = nf90_put_var(ncid(32),varID(32,3),out_1D(:,3), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed3')
	status = nf90_put_var(ncid(32),varID(32,4),out_1D(:,4), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed4')
	status = nf90_put_var(ncid(32),varID(32,5),out_1D(:,5), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed5')
	status = nf90_put_var(ncid(32),varID(32,6),out_1D(:,6), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed6')
	status = nf90_put_var(ncid(32),varID(32,7),out_1D(:,7), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed7')
	status = nf90_put_var(ncid(32),varID(32,8),out_1D(:,8), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed8')
	status = nf90_put_var(ncid(32),varID(32,9),out_1D(:,9), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed9')
	status = nf90_put_var(ncid(32),varID(32,10),out_1D(:,10), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed10')
	status = nf90_put_var(ncid(32),varID(32,11),out_1D(:,11), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed11')
	status = nf90_put_var(ncid(32),varID(32,12),out_1D(:,12), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed12')
	status = nf90_put_var(ncid(32),varID(32,13),out_1D(:,13), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed13')
	status = nf90_put_var(ncid(32),varID(32,14),out_1D(:,14), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed14')
	status = nf90_put_var(ncid(32),varID(32,15),out_1D(:,15), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed15')
	status = nf90_put_var(ncid(32),varID(32,16),out_1D(:,16), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed16')
	status = nf90_put_var(ncid(32),varID(32,17),out_1D(:,17), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed17')
	status = nf90_put_var(ncid(32),varID(32,18),out_1D(:,18), &
	start=(/1/),count=(/(outputSpeed+50)/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_speed18')

	
	! CLOSE NETCDF-FILE
	status = nf90_close(ncid(32))
	if(status /= nf90_noerr) call handle_err(status,'nf_close')
	
	end subroutine	




!----------------------
	subroutine save_out_2D(outputProf,proflayers,out_2D_dens,out_2D_temp, &
		out_2D_lwc,out_2D_depth,out_2D_dRho,out_2D_year,settingsfile,fname_p1,username)
	
	use netcdf, only: nf90_create,nf90_noerr,nf90_def_dim,nf90_def_var, &
			  nf90_enddef,nf90_put_var,nf90_real,nf90_int,nf90_unlimited, &
			  nf90_close,nf90_put_att
	
	implicit none

	integer :: outputProf, proflayers
	integer :: status,ncid(50),IDs(50,5),varID(50,20)
	real,dimension((outputProf+50),proflayers) :: out_2D_dens,out_2D_temp,out_2D_lwc
	real,dimension((outputProf+50),proflayers) :: out_2D_depth,out_2D_dRho,out_2D_year
	character*255 :: pad,settingsfile,fname_p1,username
	
	pad = "/scratch/ms/nl/"//trim(username)//"/data/output/era055/"
	
	ncid(:) = 0
	IDs(:,:) = 0
	varID(:,:) = 0

	! CREATE NETCDF FILES
	status = nf90_create(trim(pad)//trim(fname_p1)//"_2D_"//trim(settingsfile)// &
		".nc",0,ncid(31))

	! DEFINE DIMENSIONS
	status = nf90_def_dim(ncid(31),"time",outputProf+50,IDs(31,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim2')	
	status = nf90_def_dim(ncid(31),"layer",proflayers,IDs(31,2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim3')

	! DEFINE VARIABLES
	status = nf90_def_var(ncid(31),"dens",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
	status = nf90_def_var(ncid(31),"temp",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var2')
	status = nf90_def_var(ncid(31),"year",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,3))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var3')
	status = nf90_def_var(ncid(31),"lwc",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,4))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var4')
	status = nf90_def_var(ncid(31),"depth",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,5))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var5')
	status = nf90_def_var(ncid(31),"dRho",nf90_real,(/IDs(31,1),IDs(31,2)/), &
		varID(31,6))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var6')

	! DEFINE ATTRIBUTES (unit could also be defined here)
	status = nf90_put_att(ncid(31),varID(31,1),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_1')
	status = nf90_put_att(ncid(31),varID(31,2),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_2')
	status = nf90_put_att(ncid(31),varID(31,3),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_3')
	status = nf90_put_att(ncid(31),varID(31,4),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_4')
	status = nf90_put_att(ncid(31),varID(31,5),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_5')
	status = nf90_put_att(ncid(31),varID(31,6),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_6')

	! END OF DEFINING FILES
	status = nf90_enddef(ncid(31))
	if(status /= nf90_noerr) call handle_err(status,'nf_enddef')

	! SAVE DATA
	status = nf90_put_var(ncid(31),varID(31,1),out_2D_dens,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid1')	
	status = nf90_put_var(ncid(31),varID(31,2),out_2D_temp,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid2')
	status = nf90_put_var(ncid(31),varID(31,3),out_2D_year,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid3')
	status = nf90_put_var(ncid(31),varID(31,4),out_2D_lwc,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid4')
	status = nf90_put_var(ncid(31),varID(31,5),out_2D_depth,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid5')
	status = nf90_put_var(ncid(31),varID(31,6),out_2D_dRho,start=(/1,1/), &
		count=(/(outputProf+50),proflayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_grid6')
	
	! CLOSE NETCDF-FILE
	status = nf90_close(ncid(31))
	if(status /= nf90_noerr) call handle_err(status,'nf_close')
	
	end subroutine	




!----------------------
	subroutine save_out_2Ddetail(outputDetail,detlayers,detthick,out_2D_det_dens,out_2D_det_temp, &
		out_2D_det_lwc,out_2D_det_refreeze,settingsfile,fname_p1,username)
	
	use netcdf, only: nf90_create,nf90_noerr,nf90_def_dim,nf90_def_var, &
			  nf90_enddef,nf90_put_var,nf90_real,nf90_int,nf90_unlimited, &
			  nf90_close,nf90_put_att
	
	implicit none

	integer :: outputDetail, detlayers, dd
	integer :: status,ncid(50),IDs(50,5),varID(50,20)
	double precision :: detthick
	double precision, dimension(detlayers) :: DetDepth,DetDZ
	real,dimension((outputDetail+50),detlayers) :: out_2D_det_dens,out_2D_det_temp
	real,dimension((outputDetail+50),detlayers) :: out_2D_det_lwc,out_2D_det_refreeze
	character*255 :: pad,settingsfile,fname_p1,username
	
	pad = "/scratch/ms/nl/"//trim(username)//"/data/output/era055/"

	ncid(:) = 0
	IDs(:,:) = 0
	varID(:,:) = 0

	DetDZ(1) = detthick
	DetDepth(1) = DetDZ(1) / 2.
	do dd = 2, detlayers
	  DetDZ(dd) = detthick
	  DetDepth(dd) = DetDepth(dd-1) + DetDZ(dd)
	end do  

	! CREATE NETCDF FILES
	status = nf90_create(trim(pad)//trim(fname_p1)//"_2Ddetail_"//trim(settingsfile)// &
		".nc",0,ncid(33))	

	! DEFINE DIMENSIONS
	status = nf90_def_dim(ncid(33),"time",outputDetail+50,IDs(33,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim4')	
	status = nf90_def_dim(ncid(33),"layer",detlayers,IDs(33,2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_dim5')

	! DEFINE VARIABLES
	status = nf90_def_var(ncid(33),"dens",nf90_real,(/IDs(33,1),IDs(33,2)/), &
		varID(33,1))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
	status = nf90_def_var(ncid(33),"temp",nf90_real,(/IDs(33,1),IDs(33,2)/), &
		varID(33,2))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var2')
	status = nf90_def_var(ncid(33),"lwc",nf90_real,(/IDs(33,1),IDs(33,2)/), &
		varID(33,3))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var3')
	status = nf90_def_var(ncid(33),"depth",nf90_real,(/IDs(33,2)/),varID(33,4))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var4')
	status = nf90_def_var(ncid(33),"dz",nf90_real,(/IDs(33,2)/),varID(33,5))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var5')
	status = nf90_def_var(ncid(33),"refreeze",nf90_real,(/IDs(33,1),IDs(33,2)/), & 
		varID(33,6))
	if(status /= nf90_noerr) call handle_err(status,'nf_def_var1')
	
	! DEFINE ATTRIBUTES (unit could also be defined here)
	status = nf90_put_att(ncid(33),varID(33,1),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_1')
	status = nf90_put_att(ncid(33),varID(33,2),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_2')
	status = nf90_put_att(ncid(33),varID(33,3),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_3')
	status = nf90_put_att(ncid(33),varID(33,4),"missing_value",9.96921e+36)
	if(status /= nf90_noerr) call handle_err(status,'nf_def_att_miss_val_4')	
	
	! END OF DEFINING FILES
	status = nf90_enddef(ncid(33))
	if(status /= nf90_noerr) call handle_err(status,'nf_enddef')
	
	! PUT IN THE CONSTANT Z-AXIS in 2Ddetail
	status = nf90_put_var(ncid(33),varID(33,4),DetDepth)
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var1_1')
	status = nf90_put_var(ncid(33),varID(33,5),DetDZ)
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var2_1')
	
	! SAVE DATA
	status = nf90_put_var(ncid(33),varID(33,1),out_2D_det_dens, &
		start=(/1,1/),count=(/(outputDetail+50),detlayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_detail1')	
	status = nf90_put_var(ncid(33),varID(33,2),out_2D_det_temp, &
		start=(/1,1/),count=(/(outputDetail+50),detlayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_detail2')
	status = nf90_put_var(ncid(33),varID(33,3),out_2D_det_lwc, &
		start=(/1,1/),count=(/(outputDetail+50),detlayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_detail3')
	status = nf90_put_var(ncid(33),varID(33,6),out_2D_det_refreeze, &
		start=(/1,1/),count=(/(outputDetail+50),detlayers/))
	if(status /= nf90_noerr) call handle_err(status,'nf_put_var_detail6')

	! CLOSE NETCDF-FILE
	status = nf90_close(ncid(33))
	if(status /= nf90_noerr) call handle_err(status,'nf_close')
	
	end subroutine	
                                                                                                                                                                            subprogram.f90                                                                                      000644  000765  000024  00000014274 14104233175 014266  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         !------------------------------------
!  SUBROUTINE WITH TIME LOOP OF FIRN MODEL
!------------------------------------

!----------------------------
	subroutine subprogram(beginT,dtmodel,dtobs,ImpExp,numTimes,numSteps, &
		numPoints,numPointsSU,dtSnow,nyears,kk,kkinit,initdepth,numberofrepeat,th,R, &
		pi,dzmax,rhoi,maxpore,writeinprof,writeinspeed,writeindetail, &
		proflayers,detlayers,detthick,acav,tsav,TempSurf,PreSol,PreLiq, &
		Sublim,SnowMelt,SnowDrif,FF10m,IceShelf,settingsfile,fname_p1,username, &
		domain,ini_fname)

	IMPLICIT NONE
	
	integer :: dtmodel,dtobs,ImpExp,IceShelf,kk,kkinit,kUL,initdepth
	integer :: numTimes,numPoints,numPointsSU,numSteps,dtSnow,nyears,beginT,i,time
	integer :: numberofrepeat,writeinprof,writeinspeed,writeindetail
	integer :: proflayers,detlayers
	double precision :: th,R,pi,dzmax,rho0,rhoi,maxpore,acav,tsav,detthick
	
	integer :: numOutputSpeed,numOutputProf,numOutputDetail
	integer :: outputSpeed,outputProf,outputDetail
	double precision :: zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy
	double precision :: Totvice,Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy
  	double precision :: Mrunoff,TotRunoff,Mrefreeze,Totrefreeze,FirnAir,TotLwc,IceMass
  	double precision :: Mrain,TotRain,Msurfmelt,TotSurfmelt,Msolin,TotSolIn,Rho0out

	character*255 :: settingsfile,fname_p1,username,domain,ini_fname
	
	double precision,dimension(detlayers) :: DetDZ
	double precision,dimension(kk) :: Rho,M,T,Depth,Mlwc,DZ,DenRho,Refreeze,Year
    	double precision,dimension(numTimes) :: TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif,FF10m
    	double precision,dimension(numPoints) :: TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM,Rho0FM    

	real, dimension(:,:), allocatable :: out_1D
	real, dimension(:,:), allocatable :: out_2D_dens,out_2D_temp,out_2D_lwc,out_2D_depth
	real, dimension(:,:), allocatable :: out_2D_dRho,out_2D_year
	real, dimension(:,:), allocatable :: out_2D_det_dens,out_2D_det_temp,out_2D_det_lwc,out_2D_det_refreeze

	call spinup_var(kUL,kkinit,zs,Mrunoff,Mrefreeze,Mrain,Msurfmelt,Msolin)


	Rho(:) = 0.
	M(:) = 0.
	T(:) = 0.
	Depth(:) = 0.
	Mlwc(:) = 0.
	DZ(:) = 0.
	DenRho(:) = 0.
	Refreeze(:) = 0.
	Year(:) = 0.


! Interpolate the RACMO-data to firn model timestep
    	call interpol(TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif,FF10m,TempFM,PsolFM, &
        PliqFM,SublFM,MeltFM,DrifFM,Rho0FM,numTimes,numSteps,numPoints,dtSnow,dtmodel,domain)

! Construct an initial firn layer (T-, rho-, dz-, and M-profile)
	rho0 = Rho0FM(1)
	
	print *, rho0
	
	
	call ini_dens(kk,kUL,initdepth,dzmax,rho0,rhoi,R,acav,tsav, &
            zs,DZ,Depth,Rho,DenRho,Mlwc,Refreeze,domain)
	
	print *, Rho(1:10)
	
	call ini_temp(kk,kUL,beginT,tsav,pi,DZ,M,T,Rho,Depth,Year,rhoi)

	numberofrepeat = numberofrepeat*2 + 5
! Spin up the model to a 'steady state'



    	if (numberofrepeat .ge. 1) then		
	  call spinup(numberofrepeat,numTimes,numPoints,numPointsSU,kk,kUL,dtmodel, &
    		R,rhoi,acav,tsav,th,dzmax,rho0,maxpore,zs,Msurfmelt,Mrain,Msolin, &
		Mrunoff,Mrefreeze,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year,TempFM, &
		PSolFM,PLiqFM,SublFM,MeltFM,DrifFM,Rho0FM,IceShelf,ImpExp,nyears,domain)

! Or read in the initial firn profiles from a file
	else
	  print *, "Load from ini-file"
	  call ini_fromfile(kk,kUL,Rho,M,T,Depth,Mlwc,DZ,username,domain, &
	  	ini_fname,settingsfile)
	  Year(:) = -999.
	  	
	endif


! Calculate the densification	  
	  call densific(kk,kUL,dtmodel,R,rhoi,acav,tsav,Rho,T,domain)

! Re-calculate the Temp-profile (explicit or implicit)		  
	  if (ImpExp .eq. 1) call sub_temp_imp(time,kk,kUL,dtmodel,th,tsav, &
	  	TempFM,numPoints,T,Rho,DZ,Depth,rhoi)
	  if (ImpExp .eq. 2) call sub_temp_exp(time,kk,kUL,dtmodel,numPoints, &
	  	TempFM,T,Rho,DZ,rhoi)

! Re-caluclate DZ/M-values according to new Rho-/T-values
! And calculate all liquid water processes
	  rho0 = Rho0FM(time)
	  
	  if (mod(time,200000) .eq. 0) then
	    print *, time, zs
	  endif  
	  
	  call vertgrid(time,kk,kUL,dtmodel,numPoints,nyears,dzmax,rho0,rhoi,acav, &
   		maxpore,zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy,TempFM,PSolFM,PLiqFM, &
		SublFM,MeltFM,DrifFM,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year, &
		ImpExp,IceShelf,Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze)

! Calculate the firn air content and total liquid water content of the firn column
	  FirnAir = 0.
	  TotLwc = 0.
	  IceMass = 0.
	  do i = 1,kUL
	    if (Rho(i).le.910.) FirnAir = FirnAir + DZ(i)*(rhoi-Rho(i))/(rhoi)
	    TotLwc = TotLwc + Mlwc(i)
	    IceMass = IceMass + M(i)
	  end do
  
	  call SpeedComp(time,numOutputSpeed,dtmodel,zs,dzd,vmelt,vacc,vsub,vsnd, &
  		vfc,vbouy,Totvice,Totvacc,Totvsub,Totvsnd,Totvfc,Totvmelt,Totvbouy, &
		Mrunoff,TotRunoff,Mrefreeze,Totrefreeze,FirnAir,TotLwc,Mrain,TotRain, &
		Msurfmelt,TotSurfmelt,Msolin,TotSolIn,IceMass,rho0,Rho0out,out_1D,outputSpeed)

! Write output to file if needed  
          if (mod(time,numOutputProf) .eq. 0.) then
	    call to_out_2D(kk,kUL,time,dtmodel,numOutputProf,outputProf,proflayers, &
		Rho,T,Mlwc,Depth,DenRho,Year,out_2D_dens,out_2D_temp,out_2D_lwc, &
		out_2D_depth,out_2D_dRho,out_2D_year)
	  endif

	  if (mod(time,numOutputDetail) .eq. 0.) then
	    call to_out_2Ddetail(kk,kUL,time,detlayers,detthick,numOutputDetail, &
		outputDetail,Rho,T,Mlwc,Refreeze,DZ,out_2D_det_dens,out_2D_det_temp, &
		out_2D_det_lwc,out_2D_det_refreeze)
	  endif  
	  
	end do
	
	Year(:) = Year(:) - nyears
	
	print *, "End of Time Loop"
	print *, " "

! Finished time loop

	print *, fname_p1

	call save_out_1D(outputSpeed,settingsfile,fname_p1,username,out_1D)
	call save_out_2D(outputProf,proflayers,out_2D_dens,out_2D_temp, &
		out_2D_lwc,out_2D_depth,out_2D_dRho,out_2D_year,settingsfile,fname_p1,username)
	call save_out_2Ddetail(outputDetail,detlayers,detthick,out_2D_det_dens,out_2D_det_temp, &
		out_2D_det_lwc,out_2D_det_refreeze,settingsfile,fname_p1,username)
    call save_out_restart(time,kk,kUL,Rho,DenRho,M,T,Depth,Mlwc,Year,Refreeze,DZ,settingsfile, &
        fname_p1,username)
	print *, "Written output data to files"
	print *, " "

	
	call ClearAll(numPoints,numTimes,kk,M,T,DZ,Mlwc,Rho,Year, &
		TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif, &
		TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM)

	print *, "Cleared all variables"
	print *, " "
	print *, "___________________________________"

103	continue

	end subroutine
	
                                                                                                                                                                                                                                                                                                                                    time_loop.f90                                                                                       000644  000765  000024  00000062376 14104233356 014103  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         	
!-------------------------------------------------------
! 	SUBROUTINES OF THE FIRN DENSIFICATION MODEL 
!-------------------------------------------------------
    
    subroutine spinup(numberofrepeat,numTimes,numPoints,numPointsSU,kk,kUL,dtmodel, &
    		R,rhoi,acav,tsav,th,dzmax,rho0,maxpore,zs,Msurfmelt,Mrain,Msolin, &
		Mrunoff,Mrefreeze,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year,TempFM, &
		PSolFM,PLiqFM,SublFM,MeltFM,DrifFM,Rho0FM,IceShelf,ImpExp,nyears,domain)
    
    implicit none
    
    integer :: numb,time,numberofrepeat,numPoints,numPointsSU
    integer :: kk,kUL,dtmodel,numTimes,nyears,i
    integer :: IceShelf,ImpExp
    
    double precision :: R,rhoi,acav,tsav,th,dzmax,rho0,maxpore,zs
    double precision :: vacc,vsub,vsnd,vmelt,vfc,vbouy,dzd
    double precision :: Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze
    double precision :: FirnAir,TotLwc,IceMass
    
    double precision, dimension(kk) :: Rho,M,T,Depth,Mlwc,DZ,DenRho,Refreeze,Year
    double precision, dimension(numPoints) :: TempFM,PSolFM,PLiqFM,SublFM
    double precision, dimension(numPoints) :: MeltFM,DrifFM,Rho0FM

    character*255 :: domain


    do numb = 1, numberofrepeat
      zs = 0
      do time = 1, numPointsSU
   
	! Calculate the densification	  
  	call densific(kk,kUL,dtmodel,R,rhoi,acav,tsav,Rho,T,domain)

	! Re-calculate the Temp-profile (explicit or implicit)		  
	if (ImpExp .eq. 1) call sub_temp_imp(time,kk,kUL,dtmodel,th,tsav, &
		TempFM,numPoints,T,Rho,DZ,Depth,rhoi)
	if (ImpExp .eq. 2) call sub_temp_exp(time,kk,kUL,dtmodel,numPoints, &
		TempFM,T,Rho,DZ,rhoi)

	rho0 = Rho0FM(time)
	    
	! Re-caluclate DZ/M-values according to new Rho-/T-values
	call vertgrid(time,kk,kUL,dtmodel,numPoints,nyears,dzmax,rho0,rhoi,acav, &
   		maxpore,zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy,TempFM,PSolFM,PLiqFM, &
		SublFM,MeltFM,DrifFM,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year, &
		ImpExp,IceShelf,Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze)

	Msurfmelt = 0.
	Mrain = 0.
	Msolin = 0.
	Mrunoff = 0.
	Mrefreeze = 0.
	DenRho(:) = 0.
	Refreeze(:) = 0.

      end do
      
      Year(:) = Year(:) - DBLE(nyears)
      
      ! Calculate the firn air content and total liquid water content of the firn column
      FirnAir = 0.
      TotLwc = 0.
      IceMass = 0.
      do i = 1,kUL
        if (Rho(i).le.910.) FirnAir = FirnAir + DZ(i)*(917.-Rho(i))/(917.)
        TotLwc = TotLwc + Mlwc(i)
        IceMass = IceMass + M(i)
      end do
      
      print *, "After spin-up #",numb,Rho(200),T(200),Year(200),kUL,zs,FirnAir,IceMass
    end do
    
    end subroutine
    

!---------------------------    
    subroutine sub_temp_exp(time,kk,kUL,dtmodel,numPoints,TempFM,T,Rho,DZ,rhoi)

    implicit none

    integer :: time,k,kk,kUL,dtmodel,numPoints
    double precision :: rhoint,kice,rhoi,kice_ref,kcal,kf,kair,kair_ref,theta
    
    double precision, dimension(numPoints) :: TempFM
    double precision, dimension(kk) :: ci,ki,G
    double precision, dimension(kk) :: Rho,T,DZ

    do k=1,kUL
      ci(k) = 152.5+7.122*T(k)
      kice = 9.828 * exp(-0.0057*T(k))              ! Paterson et al., 1994
      kice_ref = 9.828 * exp(-0.0057*270.15)              ! Paterson et al., 1994
      kcal = 0.024 - 1.23E-4*Rho(k) + 2.5E-6*Rho(k)**2.
      kf = 2.107 + 0.003618*(Rho(k)-rhoi)           ! Calonne (2019)
      kair = (2.334E-3*T(k)**(3/2))/(164.54 + T(k)) ! Reid (1966)
      kair_ref = (2.334E-3*270.15**(3/2))/(164.54 + 270.15)
      theta = 1./(1.+exp(-0.04*(Rho(k)-450.)))
      ki(k) = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)
    enddo

    G(1) = 0.
    G(kUL) = -ki(kUL)*(TempFM(time)-T(kUL))/DZ(kUL)
    do k=kUL-1,2,-1
      G(k) = -ki(k)*(T(k+1)-T(k))/DZ(k)
    enddo

    T(1) = T(2)
    do k=2,kUL
      rhoint = (Rho(k-1)*DZ(k-1)+Rho(k)*DZ(k))/(DZ(k-1)+DZ(k))
      T(k) = T(k)-dtmodel/(rhoint*ci(k))*(G(k)-G(k-1))/(0.5*DZ(k)+ &
      	0.5*DZ(k-1))
    enddo
    
    end subroutine


!-----------------------------------
    subroutine sub_temp_imp(time,kk,kUL,dtmodel,th,tsav,TempFM,numPoints, &
    	T,Rho,DZ,Depth,rhoi)
        
    implicit none
       
    integer :: k,kUL,kk,dtmodel,numPoints,time
    double precision :: tsav,th,kice,rhoi,kice_ref,kcal,kf,kair,kair_ref,theta
    double precision :: ci(kk),ki(kk),D(kk),l(kk),b(kk)
    double precision :: alpha(kk),beta(kk),gamma(kk),ar(kk),br(kk),cr(kk)
    double precision :: TempFM(numPoints)
    double precision,dimension(kk) :: Rho,T,DZ,Depth

    do k=1,kUL-1
!      if (Depth(k).gt.20) T(k)=tsav
      ci(k) = 152.5+7.122*T(k)

      kice = 9.828 * exp(-0.0057*T(k))              ! Paterson et al., 1994
      kice_ref = 9.828 * exp(-0.0057*270.15)              ! Paterson et al., 1994
      kcal = 0.024 - 1.23E-4*Rho(k) + 2.5E-6*Rho(k)**2.
      kf = 2.107 + 0.003618*(Rho(k)-rhoi)           ! Calonne (2019)
      kair = (2.334E-3*T(k)**(3/2))/(164.54 + T(k)) ! Reid (1966)
      kair_ref = (2.334E-3*270.15**(3/2))/(164.54 + 270.15)
      theta = 1./(1.+exp(-0.04*(Rho(k)-450.)))
      ki(k) = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)

      D(k) = ki(k)/(Rho(k)*ci(k))
      l(k) = 2.*D(k)*dtmodel/((DZ(k))**2.)

      alpha(k) = 1.+th*l(k)
      beta(k) = -th*l(k)/2.
      gamma(k) = -th*l(k)/2.
      br(k) = (1.-th)*l(k)/2.
      ar(k) = 1.-((1.-th)*l(k))
      cr(k) = (1.-th)*l(k)/2.
    enddo

    ci(kUL) = 152.5+7.122*T(kUL) 

    !ki(kUL) = 0.021+2.5*(Rho(kUL)/1000.)**2           ! Anderson (1976)
    kice = 9.828 * exp(-0.0057*T(kUL))              ! Paterson et al., 1994
    kice_ref = 9.828 * exp(-0.0057*270.15)              ! Paterson et al., 1994
    kcal = 0.024 - 1.23E-4*Rho(kUL) + 2.5E-6*Rho(kUL)**2.
    kf = 2.107 + 0.003618*(Rho(kUL)-rhoi)           ! Calonne (2019)
    kair = (2.334E-3*T(kUL)**(3/2))/(164.54 + T(kUL)) ! Reid (1966)
    kair_ref = (2.334E-3*270.15**(3/2))/(164.54 + 270.15)
    theta = 1./(1.+exp(-0.04*(Rho(kUL)-450.)))
    ki(kUL) = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)


    D(kUL) = ki(kUL)/(Rho(kUL)*ci(kUL))
    l(kUL) = 2.*D(kUL)*dtmodel/((DZ(kUL))**2.)

    alpha(kUL) = 1.+th*l(kUL)
    beta(kUL) = -th*l(kUL)/2.
    gamma(kUL) = -th*l(kUL)/2.
    br(kUL) = (1.-th)*l(kUL)/2.
    ar(kUL) = 1.-((1.-th)*l(kUL))
    cr(kUL) = (1.-th)*l(kUL)/2.
    
    !creating the matrix b:
    do k=2,kUL-1
      b(k)=br(k)*T(k-1)+ar(k)*T(k)+cr(k)*T(k+1)
    enddo
    b(kUL)=br(kUL)*T(kUL-1)+ar(kUL)*T(kUL)+cr(kUL)*TempFM(time)
    b(1) = br(1)*T(1)+ar(1)*T(1)+cr(1)*T(2)
                
    call tridiag(alpha,beta,gamma,b,T,kUL,TempFM,numPoints,time)
    T(1)=T(2)

    end subroutine


!----------------------------------
    subroutine tridiag(alpha,beta,gamma,b,T,kUL,TempFM,numPoints,time)
!	this subroutine solves a system of linear equations that is 
!	tridiagonal	

    implicit none

    integer :: kUL,numPoints,time
    double precision :: alpha(kUL),beta(kUL),gamma(kUL),b(kUL),T(kUL)
    double precision :: betadot(kUL),bdot(kUL)
    double precision :: TempFM(numPoints)
    integer :: k

    betadot(1) = 2*beta(1)/alpha(1)
    bdot(1) = b(1)/alpha(1)
    do k=2,kUL      !decomposition and forward substitution
      betadot(k)=beta(k)/(alpha(k)-gamma(k)*betadot(k-1))
      bdot(k)=(b(k)-gamma(k)*bdot(k-1))/(alpha(k)-gamma(k)* &
      	betadot(k-1))
    enddo
    T(kUL)=bdot(kUL)-TempFM(time)*betadot(kUL)
    if (T(kUL).gt.272.65) T(kUL)=272.65
    do k=kUL-1,1,-1     !backsubstitution
      T(k)=bdot(k)-T(k+1)*betadot(k)
      if (T(k).gt.272.65) T(k)=272.65
    enddo

    end subroutine
    

!----------------------------------
    subroutine densific(kk,kUL,dtmodel,R,rhoi,acav,tsav,Rho,T,domain)

    implicit none
    
    integer :: k,kk,kUL,dtmodel
    double precision :: Ec,Eg,g,cons,part1,Krate
    double precision :: rhoi,R,acav,tsav
    
    character*255 :: domain    
    
    double precision, dimension(kk) :: T,Rho
	
    	do k=1,kUL
      	 Ec = 60000.
      	 Eg = 42400.
      	 g = 9.81
	 if (Rho(k) .le. 550.) then
	  if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
	   !cons = 0.3867 + 0.0415*log(acav)      ! Snow density mean annual temp FGRN055, new fit 1957-2020 run, rho=315
	   cons = 0.6569 + 0.0067*log(acav)      ! Snow density mean annual temp FGRN055, new fit 1957-2020 run
	   !cons = 0.6545 + 0.0022*log(acav)      ! Snow density mean annual temp FGRN055
	   !cons = 1.042 - 0.0916*log(acav)     ! Old snow density
       !cons = 1.
	  else
	   cons = 1.435 - 0.151*log(acav)
	  endif
	  if (cons.lt.0.25) cons = 0.25
	  !part1 = exp((-Ec/(R*T(k)))+(Eg/(R*tsav)))
	  part1 = exp((-Ec/(R*T(k)))+(Eg/(R*T(1))))
	  Krate = 0.07*cons*acav*g*part1 
	 else
	  if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
	    !cons = 1.7231 - 0.2009*log(acav)   ! Snow density mean annual temp FGRN055, new fit 1957-2020 run, rho=315
	    cons = 1.7243 - 0.2011*log(acav)   ! Snow density mean annual temp FGRN055, new fit 1957-2020 run
	    !cons = 1.6942 - 0.1941*log(acav)   ! Snow density mean annual temp FGRN055
	    !cons = 1.734 - 0.2093*log(acav)    ! Old snow denisty
        !cons = 1.
	  else
       cons = 2.366 - 0.293*log(acav)
	  endif
	  if (cons.lt.0.25) cons = 0.25
	  !part1 = exp((-Ec/(R*T(k)))+(Eg/(R*tsav)))
	  part1 = exp((-Ec/(R*T(k)))+(Eg/(R*T(1))))
	  Krate = 0.03*cons*acav*g*part1
	 endif
         Rho(k)=Rho(k)+dtmodel/(3600.*24.*365.)*Krate*(rhoi-Rho(k))
         if (Rho(k).gt.rhoi) Rho(k) = rhoi
      	enddo
    
    end subroutine
      

!------------------------------------
    subroutine vertgrid(time,kk,kUL,dtmodel,numPoints,nyears,dzmax,rho0,rhoi,acav, &
   		maxpore,zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy,TempFM,PSolFM,PLiqFM, &
		SublFM,MeltFM,DrifFM,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year, &
		ImpExp,IceShelf,Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze)

    implicit none
    
    integer :: kUL,k,kk,time,dtmodel,numPoints,ImpExp,IceShelf,nyears
    double precision :: dzmax,rho0,rhoi,acav,maxpore
    double precision :: mdiff,macc,mice,mrun
    double precision :: Ts,Psol,Pliq,Su,Sd,Me,Ws,Mmelt
    double precision :: Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze
    double precision :: vmelt,vzd,dzd,zs,vacc,vsub,vsnd,vfc,vbouy,oldDZ

    double precision,dimension(numPoints) :: TempFM,PSolFM,PLiqFM,SublFM
    double precision,dimension(numPoints) :: MeltFM,DrifFM
    double precision,dimension(kk) :: Rho,M,T,DZ,Depth,Mlwc,Refreeze,DenRho,Year

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare climate input for timestep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Ts = TempFM(time)
    Psol = PsolFM(time)
    Pliq = PliqFM(time)
    Su = SublFM(time)
    Sd = DrifFM(time)
    Me = MeltFM(time)
    
    if (ImpExp.eq.2) then
      if (Me .gt. 1e-06) then
	Mmelt = Me + Pliq
	Msurfmelt = Msurfmelt + Me
	Mrain = Mrain + Pliq
      elseif (Ts .gt.273.15) then
	Mmelt = Me + Pliq
	Msurfmelt = Msurfmelt + Me
	Mrain = Mrain + Pliq
      else
        Psol = Psol + Pliq
	Pliq = 0.
        Mmelt = 0.
	Me = 0.
      endif
    else
      Psol = Psol + Pliq
      Mmelt = 0.
      Pliq = 0.
      Me = 0.
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate mass change in kUL and velocity conponents of the surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   vfc = 0.
    do k=1,kUL
      oldDZ = DZ(k)
      DZ(k) = M(k)/Rho(k) 	!update layer thickness for new density
      vfc = vfc - (oldDZ - DZ(k))
      DenRho(k) = DenRho(k) - (oldDZ - DZ(k))   
    enddo

    M(kUL) = M(kUL)+(Psol+Su-Sd) 		!add mass to upper layer
    Msolin = Msolin + (Psol+Su-Sd)

    DZ(kUL) = DZ(kUL)+((Psol)/rho0)	!recalculate height of upper layer
    vacc = (Psol)/rho0

    DZ(kUL) = DZ(kUL)+(Su/Rho(kUL))
    vsub = Su/Rho(kUL)

    if (Sd.gt.0.) then
      DZ(kUL) = DZ(kUL)-(Sd/Rho(kUL))
      vsnd = -1. * (Sd/Rho(kUL))  	
    else
      DZ(kUL) = DZ(kUL)-(Sd/rho0)
      vsnd = -1. * (Sd/rho0)  	
    endif

    Rho(kUL) = M(kUL)/DZ(kUL)		!recalculate density of upper layer

    vzd = acav/(3600.*24.*365.*rhoi)     		! vertical (downward) velocity of lowest firn layer [m/s]
    dzd = vzd*dtmodel                       	! vertical displacement of lowest firn layer in time step

    vmelt = -1 * Me/Rho(kUL)
    if (ImpExp.eq.2) then
      call melt(time,kk,kUL,dtmodel,maxpore,rhoi,Me,Mmelt,T,M,Rho,DZ,Mlwc, &
      	Refreeze,Mrunoff,Mrefreeze)	
      call LWrefreeze(kk,kUL,rhoi,Mrefreeze,T,M,Rho,DZ,Mlwc,Refreeze)
    endif
     
    if (IceShelf.eq.1) then   ! Ice Shelf = on, so vbouy has to be calculated
      macc = Psol + Su + Pliq
      mice = dzd * rhoi
      mrun = Mrunoff
      mdiff = macc - mice - mrun
      vbouy = -1. * (mdiff/1027.)   ! rho-water=1027 kg m-3
    else
      vbouy = 0.
    endif
     
    zs = zs - dzd + vfc + vacc + vsub + vsnd + vmelt + vbouy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if the vertical grid is still valid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (DZ(kUL).gt.dzMAX) call SplitLayers(kk,kUL,dzmax,Rho,M,T,Mlwc,DZ,DenRho,Refreeze, &
    	Year,time,numPoints,nyears)
    if (DZ(kUL).le.dzMAX/3.) call MergeLayers(kk,kUL,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year)
    if (kUL.gt.2500 .and. MINVAL(Rho(1:200)).ge.rhoi-7.) call DeleteLayers(kk,kUL,Rho,M,T,Mlwc,DZ, &
    	DenRho,Refreeze,Year)
    if (kUL.lt.200) call AddLayers(kk,kUL,dzmax,rhoi,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year) 
                 
    Depth(kUL)=DZ(kUL)/2.
    do k=kUL-1,1,-1
      Depth(k) = Depth(k+1) + DZ(k+1)/2. + DZ(k)/2.
    enddo   
     
    return
    end subroutine


!--------------------------------
    subroutine melt(time,kk,kUL,dtmodel,maxpore,rhoi,Me,Mmelt,T,M,Rho, &
    	DZ,Mlwc,Refreeze,Mrunoff,Mrefreeze)

	implicit none

	integer :: time,k,kk,kUL,dtmodel,icecount
	double precision :: Me,Mmelt,Mrunoff,Mrefreeze,rhoi,maxpore
	double precision :: cp,cp0,lh,poro,Mporeavail,Mporespace
	double precision :: MavailCol,MavailMax
	double precision :: Efreeze,Mfreeze,Mavail,Madd,toomuch
        double precision, dimension(kk) :: Rho,T,DZ,M,Mlwc,Refreeze

	lh = 333500.
	cp0 = 152.5+7.122*273.15

	icecount = 0

	M(kUL) = M(kUL) - Me				!Substract melted snow from upper layer
	DZ(kUL) = DZ(kUL) - (Me/Rho(kUL))	!Recalculate the height of the upper layer

	do k = kUL,1,-1
         if (Rho(k) .lt. 917.) then
	 
          cp = 152.5+7.122*T(k)	
	  Efreeze = (273.15-T(k)) * M(k) * cp		!Calculate the energy available for freezing in the layer
	  Mfreeze = Efreeze / lh			!the mass that can be frozen with that energy
	  
	  ! Maximum available capacity for liquid water according to Coleou, 1998
  	  poro = (rhoi-Rho(k))/rhoi
      maxpore = 0.017 + 0.057 * (poro/(1.-poro))
	  MavailCol = maxpore * M(k)

	  ! Maximum available capacity for liquid water for high density firn
	  ! Coleou, 1998 parameterization still has 1.7% of water for the density of ice...
	  MavailMax = (rhoi * DZ(k)) - M(k)

	  Mavail = MIN(MavailCol,MavailMax)	!the available pore space (in kg) in the layer
	  

	  !Check if all the water can be refrozen in this layer
	  !if YES, all the water is refrozen in this layer. T/M/dz recalculated
	  if (Mfreeze.ge.Mmelt.and.Mavail.ge.Mmelt) then
	  
	    Madd = Mmelt	
	    Mmelt = 0.			
	    T(k) = T(k) + ((Madd*lh) / (M(k)*cp0))
	    M(k) = M(k) + Madd
	    Mrefreeze = Mrefreeze + Madd
	    Refreeze(k) = Refreeze(k) + Madd
	    Rho(k) = M(k) / DZ(k)
          
	  ! if NO, check which of the two is the (most) limiting factor.
	  ! If it is the available energy, the refreezable part is refrozen. T .eq. 273.15
	  ! If it is the available pore space, the pore space is filled. But T .ne. 273.15
	  ! The remainder of the meltwater will percolate further
	  else
	    if (Mfreeze .gt. Mavail) then
	      Madd = Mavail
	      M(k) = M(k) + Madd
	      Mrefreeze = Mrefreeze + Madd
	      Refreeze(k) = Refreeze(k) + Madd	      
	      Rho(k) = M(k) / DZ(k)
	      T(k) = T(k) + ((Madd*lh) / (M(k)*cp0))	    
	    
	    else
	      Madd = Mfreeze
	      M(k) = M(k) + Madd
	      Mrefreeze = Mrefreeze + Madd
	      Refreeze(k) = Refreeze(k) + Madd	      
	      Rho(k) = M(k) / DZ(k)
	      T(k) = 273.15	    
	    	    
	    endif
	   
	    Mmelt = Mmelt - Madd	    
	    call LWcontent(k,kk,rhoi,maxpore,Mmelt,M,Rho,DZ,Mlwc) !Possible LWC will remain in the layer as liquid water	  
	  
	   endif	  
	  endif  
	  
	  ! Supersaturation check
	  ! Extra check that densification does not lead to supersaturated layers
	  
	  ! Maximum available capacity for liquid water according to Coleou, 1998
  	  poro = (rhoi-Rho(k))/rhoi
      maxpore = 0.017 + 0.057 * (poro/(1.-poro))
      !maxpore = 0.07 / ( 0.07 + (Rho(k)*rhoi/(1000.*(rhoi-Rho(k)))) )
      !maxpore = 0.02
	  MavailCol = maxpore * M(k)

	  ! Maximum available capacity for liquid water for high density firn
	  ! Coleou, 1998 parameterization still has 1.7% of water for the density of ice...
	  MavailMax = (rhoi * DZ(k)) - M(k)

	  Mavail = MIN(MavailCol,MavailMax)	!the available pore space (in kg) in the layer

	  if (Mavail .lt. Mlwc(k)) then	!Check if densification did not cause too much LWC
	    toomuch = Mlwc(k) - Mavail
	    Mlwc(k) = Mlwc(k) - toomuch
	    Mmelt = Mmelt + toomuch
          endif
          
	  ! End supersaturation check	  
	  
	enddo
	
	if (Mmelt .ne. 0.) then
	  Mrunoff = Mmelt
	  Mmelt = 0.
	endif
	  
46 	continue
      
    return
    end subroutine


!-----------------------------
    subroutine LWcontent(k,kk,rhoi,maxpore,Mmelt,M,Rho,DZ,Mlwc)

    implicit none

    integer :: k,kk
    double precision :: rhoi,Mmelt,maxpore,poro
    double precision :: Mavail,MavailCol,MavailMax,Mporeavail,toomuch
    double precision,dimension(kk) :: M,Rho,DZ,Mlwc
    
    ! Maximum available capacity for liquid water according to Coleou, 1998
    poro = (rhoi-Rho(k))/rhoi
    maxpore = 0.017 + 0.057 * (poro/(1.-poro))
    MavailCol = maxpore * M(k)

    ! Maximum available capacity for liquid water for high density firn
    ! Coleou, 1998 parameterization still has 1.7% of water for the density of ice...
    MavailMax = (rhoi * DZ(k)) - M(k)

    Mavail = MIN(MavailCol,MavailMax)	!the available pore space (in kg) in the layer

    if (Mavail .lt. Mlwc(k)) then	!Check if densification did not cause to much LWC
	  toomuch = Mlwc(k) - Mavail
	  Mlwc(k) = Mlwc(k) - toomuch
	  Mmelt = Mmelt + toomuch
	  goto 47
    endif
    
    Mporeavail = Mavail - Mlwc(k)
    if (Mporeavail .ge. Mmelt) then
      Mlwc(k) = Mlwc(k) + Mmelt	!Enough pore space available: all water stored
	  Mmelt = 0.
	  goto 47
    else
	  Mlwc(k) = Mlwc(k) + Mporeavail	!Not enough pore space available: 2% pore space stored
	  Mmelt = Mmelt - Mporeavail
	  Mporeavail = 0.
    endif

47	continue
    
    return
    end subroutine

!------------------------------
    subroutine LWrefreeze(kk,kUL,rhoi,Mrefreeze,T,M,Rho,DZ,Mlwc,Refreeze)

    implicit none

    integer :: kk,kUL,k
    double precision :: rhoi,cp,Mfreeze,cp0,lh,Mrefreeze
    double precision,dimension(kk) :: Rho,T,M,DZ,Mlwc,Refreeze


	lh = 333500.
	cp0 = 152.5+7.122*273.15

	do k = 1,kUL
	  if (Mlwc(k).gt.0 .and. T(k).ne.273.15) then
        cp = 152.5+7.122*T(k)
	    Mfreeze = ((273.15-T(k)) * M(k) * cp) / lh	!Available energy (mass) for refreezing
        if (Mfreeze .ge. Mlwc(k)) then
	      M(k) = M(k) + Mlwc(k)	!Enough energy: all LWC is refrozen
	      Mrefreeze = Mrefreeze + Mlwc(k)
	      Refreeze(k) = Refreeze(k) + Mlwc(k)
	      Rho(k) = M(k) / DZ(k)
	      T(k) = T(k) + ((Mlwc(k)*lh) / (M(k)*cp0))
	      Mlwc(k) = 0.
	    else 
	      M(k) = M(k) + Mfreeze	!Not enough energy: all energy is used for refreezing. Rest will remain LWC
	      Mrefreeze = Mrefreeze + Mfreeze
	      Refreeze(k) = Refreeze(k) + Mfreeze
	      Rho(k) = M(k) / DZ(k)
	      T(k) = 273.15
	      Mlwc(k) = Mlwc(k) - Mfreeze
	    endif
	  endif
	enddo
      
    return
    end subroutine
    

!--------------------------------
  subroutine SpeedComp(time,numOutputSpeed,dtmodel,zs,dzd,vmelt,vacc,vsub,vsnd, &
  	vfc,vbouy,Totvice,Totvacc,Totvsub,Totvsnd,Totvfc,Totvmelt,Totvbouy, &
	Mrunoff,TotRunoff,Mrefreeze,Totrefreeze,FirnAir,TotLwc,Mrain,TotRain, &
	Msurfmelt,TotSurfmelt,Msolin,TotSolIn,IceMass,rho0,Rho0out,out_1D,outputSpeed)

  implicit none

  integer :: outputSpeed,numOutputSpeed,dtmodel,time
  double precision :: zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy,factor
  double precision :: Totvice,Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy,Totv
  double precision :: Mrunoff,TotRunoff,Mrefreeze,Totrefreeze,FirnAir,TotLwc,IceMass
  double precision :: Mrain,TotRain,Msurfmelt,TotSurfmelt,Msolin,TotSolIn,rho0,Rho0out

  real, dimension(outputSpeed+50,18) :: out_1D	

  Totvice = Totvice + (dzd/dtmodel) 		! m s-1
  Totvacc = Totvacc + (vacc/dtmodel)		! m s-1
  Totvsub = Totvsub + (vsub/dtmodel)
  Totvsnd = Totvsnd + (vsnd/dtmodel)	
  Totvfc = Totvfc + (vfc/dtmodel)			! m s-1
  Totvmelt = Totvmelt + (vmelt/dtmodel)		! m s-1
  Totvbouy = Totvbouy + (vbouy/dtmodel)     ! m s-1
  
  TotRunoff = TotRunoff + Mrunoff
  TotRefreeze = TotRefreeze + Mrefreeze
  TotRain = TotRain + Mrain
  TotSurfMelt = TotSurfMelt + Msurfmelt
  TotSolIn = TotSolIn + Msolin
  Rho0out = Rho0out + rho0
  
  Mrunoff = 0.
  Mrefreeze = 0.
  Mrain = 0.
  Msurfmelt = 0.
  Msolin = 0.
  
  if (mod(time,numOutputSpeed) .eq. 0) then
    factor = (3600.*24.*365.)/numOutputSpeed
    Totvice = -1 * (Totvice * factor)
    Totvacc = Totvacc * factor
    Totvsub = Totvsub * factor
    Totvsnd = Totvsnd * factor
    Totvmelt = Totvmelt * factor
    Totvfc = Totvfc * factor
    Totvbouy = Totvbouy * factor
    
    Totv = Totvacc + Totvice + Totvmelt + Totvfc + Totvbouy
    
    Rho0out = Rho0out / numOutputSpeed

    call to_out_1D(time,numOutputSpeed,zs,Totvice,Totvfc,Totvacc,Totvsub,Totvsnd, &
    	Totvmelt,Totvbouy,Totv,TotRunoff,FirnAir,TotLwc,TotRefreeze,TotRain, &
	TotSurfmelt,TotSolIn,IceMass,Rho0out,out_1D,outputSpeed)

  endif

  end subroutine


!-----------------------------------
  subroutine ClearAll(numPoints,numTimes,kk,M,T,DZ,Mlwc,Rho,Year, &
  	TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif, &
	TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM)
  
  integer :: numPoints,numTimes,kk
  integer :: p,tt,k
  
  double precision, dimension(kk) :: Rho,M,T,DZ,Depth,Mlwc,Year
  double precision, dimension(numTimes) :: TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif
  double precision, dimension(numPoints) :: TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM
  
  do k = 1, kk
   M(k) = 0.
   T(k) = 0.
   DZ(k) = 0.
   Depth(k) = 0.
   Mlwc(k) = 0.
   Rho(k) = 0.
   Year(k) = 0.
  end do
  
  do tt = 1, numTimes
   TempSurf(tt) = 0.
   PreSol(tt) = 0.
   PreLiq(tt) = 0.
   Sublim(tt) = 0.
   SnowMelt(tt) = 0.
   SnowDrif(tt) = 0.
  end do
  
  do p = 1, numPoints
   TempFM(p) = 0.
   PsolFM(p) = 0.
   PliqFM(p) = 0.
   SublFM(p) = 0.
   MeltFM(p) = 0.
   DrifFM(p) = 0.
  end do
    
  end subroutine


!-----------------------------------
  subroutine SplitLayers(kk,kUL,dzmax,Rho,M,T,Mlwc,DZ,DenRho,Refreeze, &
  	Year,time,numPoints,nyears)
  
  integer :: kk,kUL,time,numPoints,nyears
  double precision :: dzmax
  double precision, dimension(kk) :: Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year
  

  T(kUL+1) = T(kUL)
  Rho(kUL+1) = Rho(kUL)
  DZ(kUL+1) = DZ(kUL)/2.
  M(kUL+1) = M(kUL)/2.
  Mlwc(kUL+1) = Mlwc(kUL)/2.
  DenRho(kUL+1) = DenRho(kUL)/2.
  Refreeze(kUL+1) = Refreeze(kUL)/2.
  Year(kUL) = (DBLE(time) * DBLE(nyears)) / DBLE(numPoints)


  DZ(kUL) = DZ(kUL+1)
  M(kUL) = M(kUL+1)
  Mlwc(kUL) = Mlwc(kUL+1)
  DenRho(kUL) = DenRho(kUL+1)
  Refreeze(kUL) = Refreeze(kUL+1)
  
  kUL = kUL + 1
  
    
  end subroutine


!------------------------------------
  subroutine MergeLayers(kk,kUL,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year)
  
  integer :: kk,kUL
  double precision, dimension(kk) :: Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year
  
  
  T(kUL-1) = (T(kUL-1)*M(kUL-1) + T(kUL)*M(kUL)) / (M(kUL-1)+M(kUL))
  M(kUL-1) = M(kUL-1)+M(kUL)
  DZ(kUL-1) = DZ(kUL-1)+DZ(kUL)
  Rho(kUL-1) = M(kUL-1)/DZ(kUL-1)
  DenRho(kUL-1) = DenRho(kUL-1)+DenRho(kUL)
  Mlwc(kUL-1) = Mlwc(kUL-1) + Mlwc(kUL)
  Refreeze(kUL-1) = Refreeze(kUL-1) + Refreeze(kUL)

  T(kUL) = 0.
  M(kUL) = 0.
  DZ(kUL) = 0.
  Rho(kUL) = 0.
  DenRho(kUL) = 0.
  Mlwc(kUL) = 0.
  Refreeze(kUL) = 0.
  Year(kUL) = 0.

  kUL = kUL-1


  end subroutine


!-----------------------------------
  subroutine DeleteLayers(kk,kUL,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year)
  
  integer :: ikk,kk,kUL
  double precision, dimension(kk) :: Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year
  
  do ikk = 1, kUL-200
    T(ikk) = T(ikk+200)
    M(ikk) = M(ikk+200)
    DZ(ikk) = DZ(ikk+200)
    Rho(ikk) = Rho(ikk+200)
    DenRho(ikk) = DenRho(ikk+200)
    Mlwc(ikk) = Mlwc(ikk+200)
    Refreeze(ikk) = Refreeze(ikk+200)
    Year(ikk) = Year(ikk+200)
  end do
      
  kUL = kUL - 200
      
  T(kUL+1:kUL+201) = 0.
  M(kUL+1:kUL+201) = 0.
  DZ(kUL+1:kUL+201) = 0.
  Rho(kUL+1:kUL+201) = 0.
  DenRho(kUL+1:kUL+201) = 0.
  Mlwc(kUL+1:kUL+201) = 0.
  Refreeze(kUL+1:kUL+201) = 0.
  Year(kUL+1:kUL+201) = 0.
    
  

  end subroutine


!------------------------------------
  subroutine AddLayers(kk,kUL,dzmax,rhoi,Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year)
  
  integer :: ikk,kk,kUL
  double precision :: dzmax,rhoi
  double precision, dimension(kk) :: Rho,M,T,Mlwc,DZ,DenRho,Refreeze,Year
  
  
  M(101:100+kUL)    = M(1:kUL)
  T(101:100+kUL)    = T(1:kUL)
  DZ(101:100+kUL)   = DZ(1:kUL)
  Mlwc(101:100+kUL) = Mlwc(1:kUL)
  Rho(101:100+kUL)  = Rho(1:kUL)
  DenRho(101:100+kUL) = DenRho(1:kUL)
  Refreeze(101:100+kUL) = Refreeze(1:kUL)
  Year(101:100+kUL) = Year(1:kUL)
      
  kUL = kUL + 100
      
  ! Properties of the 100 new layers:
  do ikk = 1,100
    DZ(ikk)  = dzmax
    Rho(ikk) = rhoi
    DenRho(ikk) = 0.
    M(ikk)   = Rho(ikk) * DZ(ikk)
    Mlwc(ikk) = 0.
    T(ikk) = T(101)
    Refreeze(ikk) = 0.
    Year(ikk) = -999.
  end do
      
  end subroutine      

                                                                                                                                                                                                                                                                  IMAU-FDM_np.x                                                                                       000644  000765  000024  00001547520 14104224004 013607  0                                                                                                    ustar 00brils001                        staff                           000000  000000                                                                                                                                                                         ELF          >    @     @       �         @ 8  @ ) &       @       @ @     @ @     �      �                           @      @                                          @       @     �     �                   �     �d     �d     @     `                   (�     (�d     (�d                                    @     @     8       8              P�td   ��     ��D     ��D     t      t             Q�td                                                  /lib64/ld-linux-x86-64.so.2          GNU                    SuSESuSE     
    %   <   -   9   7   /       3   #   0      5       !   8      :      )                   2   '   4      1   6      ;         .          (                                                                                    	                     
                                                         "   *             %                 $       +       ,      &                                  �                     %                       �                     �                      l                     �     �dD            0                     X                     n    (�e             K                     z    H�e             %                     A                       �                     z                     �                     �                     �                     g    (�e             E                     P                       "                     �                     /                     :                     �                     �                     �                     ~     @             �                     �                      �                      �                     �                                          T                     T   ��  @             8                     ~                     t                     �                     �                     d                       l                     �                      #                     �                     Y                     �                                          ]                     �                     s                     �                     4                     �     �cD             	                     �      �e     @       e                      libeccodes.so libAtpSigHandler.so.0 _ITM_deregisterTMCloneTable __gmon_start__ _Jv_RegisterClasses _ITM_registerTMCloneTable _init _fini __atpHandlerInstall _ATP_Text_Globals _ATP_Data_Globals librca.so.0 libnetcdff.so.6 _F90_LEN_TRIM_ __cray_sset_HSW nf_def_dim_ nf90_def_var_onedim$netcdf_ nf90_put_att_one_fourbytereal$netcdf_ _DEALLOC nf_open_ nf90_inq_dimid$netcdf_ nf90_def_var_manydims$netcdf_ nf_put_vara_double_ nf90_inq_varid$netcdf_ nf_put_var1_int_ __ALLOCATE nf90_strerror$netcdf_ nf90_def_dim$netcdf_ nf90_inquire_dimension$netcdf_ nf__enddef_ nf_create_ nf_enddef_ nf90_close$netcdf_ nf_get_vara_double_ nf_put_vara_real_ __cray_dset_HSW libeccodes_f90.so __pgas_register_dv _COS_V _COS libpgas-dmapp.so.2 memmove memset libquadmath.so.0 libomp.so.1 libcraymp.so.1 libmodules.so.1 _FWF _STOP2 libfi.so.1 _F90_STRING_COMPARE libcraymath.so.1 __cray2_EXP __cray2_EXP_W _COS_W __cray2_EXP_V __cray2_ALOG libf.so.1 _OPEN _END _CLOSE _FRF libu.so.1 getarg_ __cray_dcopy_HSW libcsup.so.1 libtcmalloc_minimal.so.1 mallopt libstdc++.so.6 libpthread.so.0 libc.so.6 strncmp __libc_start_main libm.so.6 __executable_start _edata __bss_start _end /opt/cray/cce/8.5.8/craylibs/x86-64:/opt/cray/cce/8.5.8/CC/x86-64/lib/x86-64:/opt/cray/gcc-libs:/usr/lib64:/lib64:/opt/cray/dmapp/default/lib64:/opt/cray/mpt/7.5.3/gni/mpich-cray/8.4/lib:/opt/cray/libsci/16.11.1/CRAY/8.3/x86_64/lib:/opt/cray/hdf5/1.10.0.1/CRAY/8.3/lib:/usr/local/apps/netcdf4/4.7.4/CRAY/85/lib:/usr/local/apps/eccodes/2.21.0/CRAY/85/lib/pkgconfig/../../lib:/usr/local/apps/jasper/1.900.1//lib:/opt/cray/rca/1.0.0-2.0502.60530.1.62.ari/lib64:/opt/cray/atp/2.1.0/libApp:/opt/cray/cce/8.5.8/craylibs/x86-64/pkgconfig/../ GLIBC_2.2.5                                                                                                                       &         ui	   �      ��d                   8�d                   ��d        ,           @�d        ,           @e        ,           @e        ,            �e        ,           ��e        ,           ��e        ,            �e        ,           @�e        ,           X�d                   `�d                   h�d                   p�d                   x�d                   ��d        
           ��d                   ��d                   ��d                   ��d                   ��d                   ��d                   ��d                   ��d                   ȫd                   Ыd                   ثd                   �d                   �d                   �d                   ��d                    �d                   �d                    �d        !           �d        "            �d        #           (�d        $           0�d        &           8�d        '           @�d        (           H�d        )           P�d        *           X�d        -           `�d        .           h�d        /           p�d        0           x�d        1           ��d        2           ��d        3           ��d        4           ��d        5           ��d        6           ��d        7           ��d        9           ��d        ;           H���/  �*  �N H����5"�$ �%$�$ @ �%"�$ h    ������%�$ h   ������%�$ h   ������%
�$ h   �����%�$ h   �����%��$ h   �����%�$ h   �����%�$ h   �p����%�$ h   �`����%ڔ$ h	   �P����%Ҕ$ h
   �@����%ʔ$ h   �0����%$ h   � ����%��$ h   �����%��$ h   � ����%��$ h   ������%��$ h   ������%��$ h   ������%��$ h   ������%��$ h   �����%��$ h   �����%z�$ h   �����%r�$ h   �����%j�$ h   �p����%b�$ h   �`����%Z�$ h   �P����%R�$ h   �@����%J�$ h   �0����%B�$ h   � ����%:�$ h   �����%2�$ h   � ����%*�$ h   ������%"�$ h    ������%�$ h!   ������%�$ h"   ������%
�$ h#   �����%�$ h$   �����%��$ h%   �����%�$ h&   �����%�$ h'   �p����%�$ h(   �`����%ړ$ h)   �P����%ғ$ h*   �@����%ʓ$ h+   �0����%$ h,   � ����\$��L$�@�  �T$�Ð1�I��^H��H���PTI�� dD H��dD H���@ �����f�H��H��$ H��t��H���f�     �/�e UH-(�e H��H��w]ø    H��t�]�(�e ���    �(�e UH-(�e H��H��H��H��?H�H��u]ú    H��t�]H�ƿ(�e ���    �=i�%  u_UH��S��d H���d H��H�S�% H��H��H9�s$fD  H��H�5�% ���d H�'�% H9�r��5�����% H��[]��    H�=؍$  t�    H��tU� �d H����]�+��� �#��� UH��H��0�}�H�u�H�U��E�    �E�    �E�    �E�    H�}� ��   �   H�E�H� � ����Ctl��Mt�   �}� u*H�E�H� �   �eD H���v�����u�E�   �E��0�}� u*H�E�H� �   �eD H���F�����u�E�   �E��3�1�}� u*H�E�H� �   �.eD H��������u�E�   �E�� �H�E��}�H�E�H� H���6����}� u*�}� u�    ������L����}� u�    ������7������f.�     f.�     UH��AWAVAUATSH���H���	  I��I�Ơ  �@�d �@�d L���������d ���d L���z����@�d �@�d L���h���Ǆ$p     HǄ$`  �e HǄ$h  �   H��$`  H��H��p  ��   0�����Ǆ$t     HǄ$P  �e HǄ$X  �   H��$P  H��H��t  ��   0������Ǆ$x     HǄ$@  ��e HǄ$H  �   H��$@  H��H��x  ��   0�����Ǆ$|     HǄ$0  �e HǄ$8  �   H��$0  H��H��|  ��   0��D���Ǆ$�     HǄ$   0�e HǄ$(  �   H��$   H��H�ǀ  ��   0��������d ���d L������Z�% �% �% 5�% &�% ��% (�% ��% ��% 3�% ܕ% M�% N�% 7�% �% HǄ$  �e HǄ$  �   H��$  HǄ$   �e HǄ$  �   H��$   HǄ$�  ��e HǄ$�  �   H��$�  H��$�   H��$�   H�D$xHǄ$�   �   HǄ$�   �   HǄ$�   �   H�D$p �e H�D$h��e H�D$`(�e H�D$Xشe H�D$P��e H�D$H��e H�D$@ȷe H�D$8�e H�D$0�e H�D$(�e H�D$ ��e H�D$�e H�D$X�e H�D$8�e H�$ȴe �дe ���e ���e ���e A���e A� �e 0��"_ K�% <�% ]�% HǄ$�  ��e HǄ$�  �   H��$�  HǄ$�  �e HǄ$�  �   L��$�  H�$�   �p�e �h�e ���e A��   0��k ��% 3�% 4�% �% �% �P�e �طe ��e �зe A���e 0��s H��������H���$ H!�H��H�~�$ E1�H���$    Hc��% H��IH�H���$ I��H�}�$    H�      H��$�  I��HǄ$�  h�d HǄ$�      HǄ$�      H��$�  H��$�  1�1�E1�0������h�d 0�����H�S�$ H!�H��H�E�$ H�Z�$    L�%[�$ H�X�$    L��$�  HǄ$�  ذd HǄ$�      HǄ$�      H��$�  H��H���  1�1�E1�0��$����ذd 0�����H���$ H!�H��H�s�$ H���$    L�%��$ H���$    L��$   HǄ$(  ��d HǄ$�      HǄ$�      H��$�  H��H��   1�1�E1�0��������d 0��~���H���$ H!�H��H���$ H���$    L�%��$ H���$    L��$`  HǄ$h  H�d HǄ$�      HǄ$�      H��$�  H��H��`  1�1�E1�0������H�d 0������H���$ H!�H��H��$ H���$    L�%��$ H���$    L��$�  HǄ$�  ��d HǄ$�      HǄ$�      H��$�  H��H�Ǡ  1�1�E1�0��������d 0��j���H�K�$ H!�H��H�=�$ H�R�$    L�%S�$ H�P�$    L��$�  HǄ$�  ��d HǄ$p      HǄ$x      H��$p  H��H���  1�1�E1�0���������d 0������H��$ H!�H��H�ۍ$ H���$    L�%�$ H��$    L��$   HǄ$(   �d HǄ$`      HǄ$h      H��$`  H��H��   1�1�E1�0��r���� �d 0��V���H���$ H!�H��H���$ H�Ɗ$    L�%Ǌ$ H�Ċ$    L��$`  M��HǄ$h  ��d HǄ$P      HǄ$X      H��$P  H��H��`  1�1�E1�0���������d 0������H�r�$ H!�H��H�d�$ H�y�$    Hc�% H��IH�H�l�$ H�i�$    H�f�$    Hc��% H��IH�H�Y�$ I��H�W�$ I��L��$�  HǄ$�  ��d HǄ$@      HǄ$H      H��$@  H��$�  1�1�E1�0��"������d 0�����H�׊$ H!�H��H�Ɋ$ H�ފ$    L�-ߊ$ L��H�ي$    H�֊$    L�=׊$ H�؊$ L��$�  HǄ$�  �d HǄ$0      HǄ$8      H��$0  H��H���  1�1�E1�0��|�����d 0��`���H��$ H!�H��H��$ H��$    L�-�$ L��H��$    H��$    L�=�$ H��$ L��$   HǄ$(  Ȯd HǄ$       HǄ$(      H��$   H��H��   1�1�E1�0�������Ȯd 0�����H�ˈ$ H!�H��H���$ H�҈$    L�-ӈ$ L��H�͈$    H�ʈ$    L�=ˈ$ H�̈$ L��$`  HǄ$h  (�d HǄ$      HǄ$      H��$  H��H��`  1�1�E1�0��0����(�d 0�����H���$ H!�H��H�w�$ H���$    L�-��$ L��H���$    H���$    L�=��$ H���$ L��$�  HǄ$�  ��d HǄ$       HǄ$      H��$   H��H�Ǡ  1�1�E1�0��������d 0��n���H�W�$ H!�H��H�I�$ H�^�$    L�-_�$ L��H�Y�$    H�V�$    L�=W�$ H�X�$ L��$�  HǄ$�   �d HǄ$�      HǄ$�      H��$�  H��H���  1�1�E1�0������� �d 0������H�х$ H!�H��H�Å$ H�؅$    L�-م$ L��H�Ӆ$    H�Ѕ$    L�=х$ H�҅$ L��$   HǄ$(   �d HǄ$�      HǄ$�      H��$�  H��H��   1�1�E1�0��>���� �d 0��"���H#k�$ H��H�`�$ H�u�$    L�-v�$ L��H�p�$    H�m�$    L�=n�$ H�o�$ L��$`  HǄ$h  `�d HǄ$�      HǄ$�      H��$�  H��H��`  1�1�E1�0������`�d 0�����H��$ H�N�$ H�d�$ H���$ H��$ ��% ��% ��% ��% H�2�$ H���$ H���$ H�=��$ H�5�$ H���$ H�7�$ L���$ L��$ L�J�$ L���$ HǄ$�  �e HǄ$�  �   H��$�  HǄ$�  ��e HǄ$�  �   H��$�  H�D$8H�\$0L�\$(L�T$ L�L$H�D$H�   H�D$@�   H�D$��e H�D$��e H�$p�e A�h�e 0���  �@�d �@�d L���������d ���d L��������% ��% H��$ H�n�$ H�Ą$ B�% 3�% H��$ H�E�$ L���$ H�$h�e ���e � �e A�p�e 0���g Lc=P�% D�=ы% HcF�% ���% �΋% �Ћ% H��H��H���$ L�H� �$ ��,D���d�% �@�d �@�d L������H�3�% H��$�  H�,�% H��$�  H���$ H��$ 	  ��(u�$ ��)�$	  ��(T�$ ��)�$ 	  ��(3�$ ��)�$�  ��(�$ ��)�$�  ��(�$ ��)�$�  ��(Б$ ��)�$�  ��(��$ ��)�$�  ��(��$ ��)�$�  H��H�  H��$�  H��H�  H��$	  H��H�Ơ  � �d L�������D��$�  ��$�  H���$ H��$�	  ��(֓$ ��)�$�	  ��(��$ ��)�$�	  ��(��$ ��)�$�	  ��(s�$ ��)�$�	  ��(R�$ ��)�$p	  ��(1�$ ��)�$`	  ��(�$ ��)�$P	  ��(�$ ��)�$@	  H��H�  H��$x	  H��H�  H��$�	  H��H��@	  �@�d L����������d ���d L��������@�d �@�d L�������� �d � �d L���������d ���d L�������@�d �@�d L������H�t�$ H�ڂ$ H���$ H�6�$ H���$ H���$ H�ȃ$ H�F�$ L�% ��% �% o�% H�=�$ H�5q�$ H�"�$ H�Ӂ$ L�4�$ L���$ L�n�$ L��$ HǄ$�  �e HǄ$�  �   H��$�  HǄ$�  ��e HǄ$�  �   H��$�  H�\$0H�D$(L�\$L�$H�D$H�   H�D$@�   H�D$8ȴe H�D$ @�e H�D$H�e H�D$��e 0��?�  ���d ���d L���M����@�d �@�d L���;����5u�% �g�% �=Ȇ% u�Ɖ5N�% �`�% �ș����% ��% ��*���*���Y���% ��W���*�% ����Y�5 ��*���Y���^���,�9�O����% H��}$ H�/}$ H�E~$ H�{}$ �% ��% �% l�% ]�% 6�% '�% x�% H�=�}$ H�5�|$ H��}$ H�$}$ H�D$(��e H�D$ @�e H�D$H�e H�D$h�e H�D$p�e H�$�e A� �e A���e 0���� ����% ��W�������T����4 ��U���V���ׂ% Ǆ$�   N  Ղ% �% ��% H�% ��% ��% ��% ��% Ճ% ��% /�% �% I�% ڇ% ��% L�̀$ L�]$ L�$ L�o�$ L�5 �$ L�=�$ H�*}$ HǄ$�  �e HǄ$�  �   H��$�  HǄ$p  �e HǄ$x  �   H��$p  HǄ$`  �e HǄ$h  �   H��$`  HǄ$P  ��e HǄ$X  �   H��$P  HǄ$@  0�e HǄ$H  �   H��$@  H��$  H��$  H��$   H��$�   H��$�   H��$�   L��$�   L��$�   L��$�   L��$�   L��$�   L��$�   H��H�  H�D$ HǄ$8  �   HǄ$0  �   HǄ$(  �   HǄ$   �   HǄ$  �   HǄ$�   �e HǄ$�    �e HǄ$�   ��e HǄ$�   ��e HǄ$�   ��e HǄ$�   ȷe HǄ$�   �e H�D$x�e H�D$p�e H�D$hP�e H�D$`�e H�D$Xشe H�D$P��e H�D$Hзe H�D$@��e H�D$8`�e H�D$0(�e H�D$(8�e H�D$��e H�D$дe H�D$��e H�$x�e ���e ���e �ȴe � �e A���e A���e 0���� 0�����H��H���[A\A]A^A_]�fD  SH��H���UUH�kH�l$H��AWAVAUATH��@  I��H������M��L��X���I��H������H��H@���H��0���Hǅ8���	   H��H ���H�����Hǅ���	   H��H ���H������Hǅ����	   H��H����H������Hǅ����	   H��H����H������Hǅ����	   H��H����H������Hǅ����	   H��Hp���H��`���Hǅh���   H��H ���H������Hǅ����   H��H����H��p���Hǅx���   H��H ���H������Hǅ����   H��H����H��p���Hǅx���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H{���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H<���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hz���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H8���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hx���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��HX���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H}���H��`���Hǅh���   H��HX���H��@���HǅH���   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��HX���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H}���H��`���Hǅh���   H��HX���H��@���HǅH���   H��H0���H�� ���Hǅ(���   H��H����H������Hǅ����   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H}���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H ���H�����Hǅ���   H��H����H��p���Hǅx���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H|���H��`���Hǅh���   H��HX���H��@���HǅH���   H��H0���H�� ���Hǅ(���   H��H����H������Hǅ����   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H}���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����
   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H ���H�����Hǅ���   H��H����H��p���Hǅx���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hy���H��`���Hǅh���   H��HX���H��@���HǅH���   H��H0���H�� ���Hǅ(���	   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����
   H��Hp���H��`���Hǅh���   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H}���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H:���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����	   H��H����H������Hǅ����	   H��H����H������Hǅ����   H��HX���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H=���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hp���H��`���Hǅh���   H��HP���H��@���HǅH����   H��H����H������H�CHcH��P���1�H��HH�IcU H��H���H��HI�H��Hǅ�����   H��H��H��H���H��H)�H���H������H��H�� ���H��H�{H��eD ��   �   �����L�s@H���w  �fD ��   �   H�{H����H���  �efD ��   �   H�{H����H����  ��fD ��   �   H�{H�s���H����  �gD ��   �   H�{H�R���H����  �dgD ��   �   H�{H�1���H����  ��gD ��   �   H�{H����H����  �hD ��   �   H�{H�����H����  L�����L��`���L��h����q% ��n���f�b% f��l����Q% ��h���H�<% H��`���ǅ����)   H��H`���H������Hǅ����   H������H��H�ǈ����   0���8 �	  L�����L��`���L��h�����   L�������1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L������N��-O���Hǅp����eD Hǅx���   H��p���L��L���T���M��    L��L��������! ��N���f��! f��L�����! ��H���H��! H��@�����! ������f��! f��������({! ��)�����H�        H��s�����y�  ���c������S������C������3������#����������������������������������������������������������������ǅ{���    ��  L�����L��`���L��h�����   L���8���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�    I�� MN�M��LH�M)�LH�H��O���L��L�������N��-O���Hǅ`��� fD Hǅh���    H��`���L��L������M��    L��L���e����m  ��N���f�^  f��L����M  ��H���H�8  H��@����O  �������?  ��������(!  �L  L�����L��`���L��h�����   L���*���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L�������N��-O���HǅP���pfD HǅX���   H��P���L��L������M��    L��L���W����� ��N���f�� f��L����� ��H���H�z H��@����� ��������(o �  L�����L��`���L��h�����   L���(���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�"   I��"MN�M��LH�M)�LH�H��O���L��L�������N��-O���Hǅ@����fD HǅH���"   H��@���L��L������M��    L��L���U���� ��N���f�� f��L����� ��H���H�� H��@����� �������� ��������(� �<  L�����L��`���L��h�����   L������1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L������N��-O���Hǅ0��� gD Hǅ8���   H��0���L��L������M��    L��L���G����O ��N���f�@ f��L����/ ��H���H� H��@����- ��������( �  L�����L��`���L��h�����   L������1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�    I�� MN�M��LH�M)�LH�H��O���L��L������N��-O���Hǅ ����gD Hǅ(���    H�� ���L��L������M��    L��L���E����� ��N���f�� f��L����� ��H���H�x H��@����� ������� ��������(a �,  L�����L��`���L��h�����   L���
���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L������N��-O���Hǅ����gD Hǅ���   H�����L��L���z���M��    L��L���7����� ��N���f�� f��L����� ��H���H�� H��@����� ��������(� ��)�����H��H�ǔ����    ��   ������  L�����L��`���L��h�����   L�������1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�    I�� MN�M��LH�M)�LH�H��O���L��L������N��-O���Hǅ ��� hD Hǅ���    H�� ���L��L���W���M��    L��L������� ��N���f� f��L����� ��H���H�� H��@����� �������� ��������(� ��)�����H�        H��u�����y� ���e������U������E������5������%����������������������������������������������������������������fǅ}���  L�sH�whD ��   �   L������H���T  �QiD ��   �   L���z���H���`  �)jD ��   �   L���Z���H���@  ��jD ��   �   L���:���H���y$  ��kD ��   �   L������H����-  �lD ��   �   L�������H���97  ��lD ��   �   L�������H����@  �>mD ��   �   L������H����K  �� ��N���f�� f��L����� ��H���H�� H��@���ǅ����*   H��H@���H������Hǅ����   H������H��H�ǌ����   �  H��H��@�����   I���������H��0���Lc�1�M��L��HH�H��H���I��I)�I���L���   L���������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L���;���O�4&Hǅ�����hD Hǅ����   H������L��L������M��    L��L�������ǅL���    L������H��0���H������H������H�� H�D$H�D$    H�D$    H�$    H��H��L���H��H��p���1�E1�E1��by  H�� ��\�����tIH�nf_open1H������H��H����H������Hǅ����   H������H��H��\����   0��+, �V ������f�G f������H��H����H������Hǅ����   H������H��H��p���H��H�����   �������\�����L�����tR�� ������H�� H������H��H����H������Hǅ����   H������H��H��\����   0��n+ �� ������f�� f������H��H����H������Hǅ����   H������H��H��p���H��H�¸����   ������\�����t`f�W f������F �����H�1 H�� ���H��H ���H������Hǅ����   H������H��H��\����   0��* �� ��"���f�� f�� ���H��H ���H��p���Hǅx���   H��p���H��H��p���H��H�¼����   �@�����\�����t`f�� f��L����� ��H���H�� H��@���H��H@���H��`���Hǅh���   H��`���H��H��\����   0���) H��P�����`���H��H�����d���ǅh���   ǅl���   H��X���H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CHc 1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������Hǅ���� �d Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��`���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H��p���H��H������I��I��P���E1��u  H����\�����t`f�k f�������Z ������H�E H������H��H����H��P���HǅX���   H��P���H��H��\����   0��' H�C(H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   H�CHc 1�H��HH�H�����Hǅ���   Hǅ���   IcH��HH�H�� ���H��(���H��H�$    H��H��p���H��H�Ƹ���H��H������1�E1�E1���t  H����\�����t`f�8 f��\����' ��X���H� H��P���H��HP���H��@���HǅH���   H��@���H��H��\����   0��S& H�C0H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CHcH��@���1�H��HH�H������Hǅ����   Hǅ����   IcH��8���H��HH�H������H������H��H�$    H��H��p���H��H�Ƽ���H��H��p���1�E1�E1��pt  H����\�����t`f�� f�������� ������H�� H������H��H����H��0���Hǅ8���   H��0���H��H��\����   0��% H��H��@�����   I��������H��(���Lc�1�M��L��HH�H��H���I��I)�I���L���   L���п�����    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L���p���O�4&Hǅ ��� iD Hǅ(���   H�� ���L��L���D���M��    L��L������ǅP���    L�����H��(���H�����H�����H�� H�D$H�D$    H�D$    H�$    H��H��P���H��H����1�E1�E1��s  H�� ��\�����tIH�nf_open6H�� ���H��H ���H�� ���Hǅ���   H�� ���H��H��\����   0��`# � ��"���f�� f�� ���H��H ���H������Hǅ����   H������H��H�Ǆ���H��H�¤����   �������\�����L�����tR�� ��H���H�� H��@���H��H@���H������Hǅ����   H������H��H��\����   0��" H��@�����h���H��8�����l���H�C8H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CHc 1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������Hǅ����@�d Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��h���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H�Ǆ���H��H�Ƥ���H��H��p���H��H������I��I��P���E1���p  H����\�����t`f�� f�������� ������H�q H������H��H����H������Hǅ����   H������H��H��\����   0��j  H��H��p���������\�����tR�+ ������H� H������H��H����H������Hǅ����	   H������H��H��\����	   0��� H��H�Ǆ���耺����\�������H  �� ������H�� H������H��H����H������Hǅ����	   H��������  H��H��@�����   I��誺����H�� ���Lc�1�M��L��HH�H��H���I��I)�I���L���   L���o������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L������O�4&Hǅ����`iD Hǅ����   H������L��L������M��    L��L��蠹��ǅT���    L������H�� ���H������H������H�� H�D$H�D$    H�D$    H�$    H��H��T���H��H��p���1�E1�E1��vo  H�� ��\�����tIH�nf_open1H�� ���H��H ���H������Hǅ����   H������H��H��\����   0��� f� f��$���� �� ���H��H ���H��p���Hǅx���   H��p���H��H��p���H��H�����   蕸����\�����L�����L�stR��
 ��H���H��
 H��@���H��H@���H��`���Hǅh���   H��`���H��H��\����   0��> �n
 ��b���f�_
 f��`���H��H`���H��P���HǅX���   H��P���H��H��p���H��H�¸����   �Է����\�����t`f�
 f�������
 ������H��	 H������H��H����H��@���HǅH���   H��@���H��H��\����   0��z ��	 ������f��	 f������H��H����H��0���Hǅ8���   H��0���H��H��p���H��H�¼����   ������\�����t`f�k	 f�������Z	 ������H�E	 H������H��H����H�� ���Hǅ(���   H�� ���H��H��\����   0�� A�������A�$������H��X���H������Hǅ����@   H�q     H�� ���H�    H�����Hǅ���    Hǅ���    Hǅ ���   Ic1�H��HH�H��(���Hǅ0���   Hǅ8���   Ic$H��HH�H��@���H��H���Hǅp�����d Hǅx���    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H�$    H��H��p���H��H�Ɛ���H��H������H��H��p���I��I������E1��2k  H����\�����t`f�I f��<����8 ��8���H�# H��0���H��H0���H�����Hǅ���   H�����H��H��\����   0�� H�C(H��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   Ic$H��HH�H������H������H��H�$    H��H��p���H��H�Ƹ���H��H��P���1�E1�E1��k  H����\�����t`f� f������� ������H�� H������H��H����H�� ���Hǅ���   H�� ���H��H��\����   0��D H�C0H������Hǅ����@   H�q     H�� ���H�    H�����Hǅ���    Hǅ���    Hǅ ���   Mc61�M��L��HH�H��(���Hǅ0���   Hǅ8���   Mc<$M��L��HH�H��@���H��H���H��H�$    H��H��p���H��H�Ƽ���H��H������1�E1�E1��j  H����\�����t`f�� f��l����� ��h���H�� H��`���H��H`���H������Hǅ����   H������H��H��\����   0��� f�� f������H�y H������H��H����H������Hǅ����
   H������H��H��p���H��H�¤����
   蒱����\�����tR�: ������H�% H������H��H����H������Hǅ����   H������H��H��\����   0��F D������D������H�C8H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   H�CHc 1�H��HH�H�����Hǅ���   Hǅ���   Ic$H��HH�H�� ���H��(���HǅP�����d HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ƥ���H��H������H��H��P���I��I������E1��h  H����\�����t`f�- f������ �����H� H�� ���H��H ���H������Hǅ����   H������H��H��\����   0�� H��H��p���虮����\�������<  �� ��(���H�� H�� ���H��H ���H������Hǅ����	   H������H��H��\����	   0���w� �D<  H��H��@�����   I��誮����#H�����Lc�1�M��L��HH�H��H���I��I)�I���L���   L���o������    I�Lc�E9�MO�M)��    LH�A�#   I��#MN�M��LH�1�M)�LH�L��H��@���L������O�4&Hǅ����@jD Hǅ����#   H������L��L������M��    L��L��蠭��ǅX���    L������H�����H������H������H�� H�D$H�D$    H�D$    H�$    H��H��X���H��H��p���1�E1�E1��g  H�� ��\�����tIH�nf_open1H��@���H��H@���H������Hǅ����   H������H��H��\����   0��� �  ��f���f��� f��d������ ��`���H��H`���H��p���Hǅx���   H��p���H��H��p���H��H�����   艬����\�����L�����tH��(�� ��)�����H��H����H��`���Hǅh���   H��`���H��H��\����   0��@ �\� ������f�M� f������H��H����H��P���HǅX���   H��P���H��H��p���H��H�¸����   �֫����\�����t`f�	� f��������� ������H��� H������H��H����H��@���HǅH���   H��@���H��H��\����   0��| ��� ������f��� f������H��H����H��0���Hǅ8���   H��0���H��H��p���H��H�¼����   ������\�����t`f�]� f������L� �����H�7� H�� ���H��H ���H�� ���Hǅ(���   H�� ���H��H��\����   0�� H�CH�����(���A���,���H��X���H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Hc1�H��HH�H��h���Hǅp���   Hǅx���   IcH��HH�H������H������Hǅ���� �d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��(���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H��p���H��H�Ɛ���H��H��0���H��H������I��I�����E1��Pc  H����\�����t`f�7� f��|����&� ��x���H�� H��p���H��Hp���H�����Hǅ���   H�����H��H��\����   0�� H�C(H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CHc 1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������H��H�$    H��H��p���H��H�Ƹ���H��H����1�E1�E1��-c  H����\�����t`f�� f�������� �����H��� H�����H��H���H�� ���Hǅ���   H�� ���H��H��\����   0��? H�C0H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   H�CLc81�M��L��HH�H��h���Hǅp���   Hǅx���   Mc6M��L��HH�H������H������H��H�$    H��H��p���H��H�Ƽ���H��H��0���1�E1�E1���b  H����\�����t`f��� f��������� ������H��� H������H��H����H������Hǅ����   H������H��H��\����   0���
 L��H��L��H��H����3  M��I��  �  M���j3  1�I��H�s8��  L��H���H��H��1�H���+1���������D ��D@��D`H��H��H���|�H����  @ f.�     ������H�� H��x��[  H��H��@�����   I���Q�����H�����Lc�1�M��L��HH�H��H���I��I)�I���L���   L���������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L��趥��O�4&Hǅ�����jD Hǅ����   H������L��L��芥��M��    L��L���G���ǅ\���    L������H�����H������H������H�� H�D$H�D$    H�D$    H�$    H��H��\���H��H��p���1�E1�E1��=a  H�� ��\�����tIH�nf_open1H������H��H����H������Hǅ����   H������H��H��\����   0�� f�O� f������H�9� H������H��H����H������Hǅ����
   H������H��H��p���H��H�����
   �:�����\�����L�����tH��(�� ��)� ���H��H ���H������Hǅ����   H������H��H��\����   0��� ��� ��"���f��� f�� ���H��H ���H������Hǅ����   H������H��H��p���H��H�¸����   臢����\�����t`f�j� f��L����Y� ��H���H�D� H��@���H��H@���H������Hǅ����   H������H��H��\����   0��- �� ��b���f� � f��`���H��H`���H��p���Hǅx���   H��p���H��H��p���H��H�¼����   �á����\�����t`f��� f��������� ������H��� H������H��H����H��`���Hǅh���   H��`���H��H��\����   0��i H�CH���������A�������H��X���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Hc1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H�� ���H�����Hǅ0���@�d Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H�°���H��H��0���I��I������E1��]  H����\�����t`f��� f��������� ������H�r� H������H��H����H��P���HǅX���   H��P���H��H��\����   0��3 H�C(H�����Hǅ���@   H�q     H�� ���H�    H��(���Hǅ0���    Hǅ8���    Hǅ@���   H�CHc 1�H��HH�H��H���HǅP���   HǅX���   IcH��HH�H��`���H��h���H��H�$    H��H��p���H��H�Ƹ���H��H�����1�E1�E1���\  H����\�����t`f�e� f�������T� ������H�?� H������H��H����H��@���HǅH���   H��@���H��H��\����   0��� H�C0H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CLc81�M��L��HH�H������Hǅ����   Hǅ����   Mc6M��L��HH�H�� ���H�����H��H�$    H��H��p���H��H�Ƽ���H��H�°���1�E1�E1��u\  H����\�����t`f�,� f��,����� ��(���H�� H�� ���H��H ���H��0���Hǅ8���   H��0���H��H��\����   0�� L��H��L��H��H���5*  M��I��  ��  M���*  1�I���   L��H���H��H��1�H���H�s8<1��    f.�     ��������D ��D@��D`H��H��H���|�H��y%�     f.�     ������H�� H��x�Hc�L��H)�H��|Hc���W�H�S8��ʃ�Hc�I9��b)  H�H�K8H��    �O)  H��H��@�����   I��赛����H�����Lc�1�M��L��HH�H��H���I��I)�I���L���   L���z������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L������O�4&Hǅ ����kD Hǅ(���   H�� ���L��L������M��    L��L��諚��ǅ`���    L�����H�����H�����H�����H�� H�D$H�D$    H�D$    H�$    H��H��`���H��H��p���1�E1�E1��Z  H�� ��\�����tIH�nf_open1H��@���H��H@���H�� ���Hǅ���   H�� ���H��H��\����   0��
�  ǅ`���maskH��H`���H������Hǅ����   H������H��H��p���H��H�����   谙����\�����L�����tR�� ������H��� H������H��H����H������Hǅ����   H������H��H��\����   0��]�  ��� ������f��� f������H��H����H������Hǅ����   H������H��H��p���H��H�¸����   ������\�����t`f�^� f�������M� ������H�8� H������H��H����H������Hǅ����   H������H��H��\����   0���  �� ������f��� f������H��H����H������Hǅ����   H������H��H��p���H��H�¼����   �/�����\�����t`f��� f�������� �����H��� H�� ���H��H ���H������Hǅ����   H������H��H��\����   0����  H�CH�����(���A���,���H��X���H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Hc1�H��HH�H��h���Hǅp���   Hǅx���   IcH��HH�H������H������Hǅ������d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��(���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H��p���H��H�Ɛ���H��H��0���H��H������I��I�����E1��mV  H����\�����t`f��� f��|����{� ��x���H�f� H��p���H��Hp���H������Hǅ����   H������H��H��\����   0���  H�C(H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CHc 1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������H��H�$    H��H��p���H��H�Ƹ���H��H����1�E1�E1��JV  H����\�����t`f�Y� f������H� �����H�3� H�����H��H���H������Hǅ����   H������H��H��\����   0��\�  H�C0H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   H�CLc81�M��L��HH�H��h���Hǅp���   Hǅx���   Mc6M��L��HH�H������H������H��H�$    H��H��p���H��H�Ƽ���H��H��0���1�E1�E1���U  H����\�����t`f� � f�������� ������H��� H������H��H����H��p���Hǅx���   H��p���H��H��\����   0���  L��H��L��H��H����   M��I��  �9
  M����   1�I��H�s8��	  L��H���H��H��1�H���H1�f�     f.�     f.�     ��������D ��D@��D`H��H��H���|�H����	  @ f.�     ������H�� H��x��[	  H��H��@�����   I���Q�����H�� ���Lc�1�M��L��HH�H��H���I��I)�I���L���   L���������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L��趒��O�4&Hǅ`��� lD Hǅh���   H��`���L��L��芒��M��    L��L���G���ǅd���    L��P���H�� ���H��X���H��P���H�� H�D$H�D$    H�D$    H�$    H��H��d���H��H��p���1�E1�E1��=T  H�� ��\�����tIH�nf_open1H������H��H����H��@���HǅH���   H��@���H��H��\����   0���  �s� ������f�d� f������H��H����H��0���Hǅ8���   H��0���H��H��p���H��H�����   �<�����\�����L�����tR�� �����H� � H�� ���H��H ���H�� ���Hǅ(���   H�� ���H��H��\����   0����  ��� ��"���f��� f�� ���H��H ���H�����Hǅ���   H�����H��H��p���H��H�¸����   ������\�����t`f�r� f��L����a� ��H���H�L� H��@���H��H@���H�� ���Hǅ���   H�� ���H��H��\����   0��%�  �� ��b���f�� f��`���H��H`���H������Hǅ����   H������H��H��p���H��H�¼����   軎����\�����t`f��� f��������� ������H��� H������H��H����H������Hǅ����   H������H��H��\����   0��a�  H�CH���������A�������H��X���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Hc1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H�� ���H�����Hǅ0�����d Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H�°���H��H��0���I��I������E1���O  H����\�����t`f��� f��������� ������H�z� H������H��H����H������Hǅ����   H������H��H��\����   0��+�  H�C(H�����Hǅ���@   H�q     H�� ���H�    H��(���Hǅ0���    Hǅ8���    Hǅ@���   H�CHc 1�H��HH�H��H���HǅP���   HǅX���   IcH��HH�H��`���H��h���H��H�$    H��H��p���H��H�Ƹ���H��H�����1�E1�E1���O  H����\�����t`f�m� f�������\� ������H�G� H������H��H����H������Hǅ����   H������H��H��\����   0����  H�C0H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   H�CLc81�M��L��HH�H������Hǅ����   Hǅ����   Mc6M��L��HH�H�� ���H�����H��H�$    H��H��p���H��H�Ƽ���H��H�°���1�E1�E1��mO  H����\�����t`f�4� f��,����#� ��(���H�� H�� ���H��H ���H������Hǅ����   H������H��H��\����   0���  L��H��L��H��H���-  M��I��  ��   M���  1�I��H�s8|wL��H���H��H��1�H���81� f.�     ��������D ��D@��D`H��H��H���|�H��y%�     f.�     ������H�� H��x�Hc�L��H)�H��|Hc���W���΃�Hc�I9��f  H�H��    �W  �@�e 0�H�{8L���1����?  H��H��@�����   I��襈����H������Lc�1�M��L��HH�H��H���I��I)�I���L���   L���j������    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L���
���O�4&Hǅ�����lD Hǅ����   H������L��L���ވ��M��    L��L��蛇��ǅh���    L������H������H������H������H�� H�D$H�D$    H�D$    H�$    H��H��h���H��H��p���1�E1�E1��M  H�� ��\�����tIH�nf_open1H��@���H��H@���H������Hǅ����   H������H��H��\����   0����  �X� ��b���f�I� f��`���H��H`���H��p���Hǅx���   H��p���H��H��p���H��H�����   萆����\�����L�����L�stR��� ������H��� H������H��H����H��`���Hǅh���   H��`���H��H��\����   0��9�  ��� ������f��� f������H��H����H��P���HǅX���   H��P���H��H��p���H��H�����   �υ����\�����tR�O� ������H�:� H������H��H����H��@���HǅH���   H��@���H��H��\����   0���  �� ������f��� f������H��H����H��0���Hǅ8���   H��0���H��H��p���H��H�¸����   ������\�����t`f��� f�������� �����H��� H�� ���H��H ���H�� ���Hǅ(���   H�� ���H��H��\����   0���  �Q� ��"���f�B� f�� ���H��H ���H�����Hǅ���   H�����H��H��p���H��H�¼����   �U�����\�����t`f� � f��L������ ��H���H��� H��@���H��H@���H�� ���Hǅ���   H�� ���H��H��\����   0����  A���h���A�$��l���H��X���H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   Ic$H��HH�H������H������Hǅ���� �d Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��h���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H��p���H��H������I��I��P���E1��H  H����\�����t`f��� f��������� ������H��� H������H��H����H������Hǅ����   H������H��H��\����   0����  H�C(H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   Hǅ���   Ic$H��HH�H�� ���H��(���H��H�$    H��H��p���H��H�Ƹ���H��H������1�E1�E1��wH  H����\�����t`f��� f��\������ ��X���H��� H��P���H��HP���H������Hǅ����   H������H��H��\����   0���  H�C0H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Mc61�M��L��HH�H������Hǅ����   Hǅ����   Mc<$M��L��HH�H������H������H��H�$    H��H��p���H��H�Ƽ���H��H��p���1�E1�E1��H  H����\�����t`f�x� f�������g� ������H�R� H������H��H����H������Hǅ����   H������H��H��\����   0��C�  D�����D�����H�C8H�����Hǅ���@   H�q     H�� ���H�    H��(���Hǅ0���    Hǅ8���    Hǅ@���   H�CHc 1�H��HH�H��H���HǅP���   HǅX���   Ic$H��HH�H��`���H��h���Hǅ����@�d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�����H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��p���H��H�Ɣ���H��H�����H��H������I��I������E1���F  H���  H��H��@�����   I���}����H�����Lc�1�M��L��HH�H��H���I��I)�I���L���   L���H}�����    I�Lc�E9�MO�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L��H��@���L����}��O�4&Hǅ����PmD Hǅ����   H������L��L���}��M��    L��L���y|��ǅl���    L������H�����H������H������H�� H�D$H�D$    H�D$    H�$    H��H��l���H��H��p���1�E1�E1��F  H�� ��\�����tIH�nf_open1H��@���H��H@���H������Hǅ����   H������H��H��\����   0����  ��� ��b���f��� f��`���H��H`���H������Hǅ����   H������H��H��p���H��H�����   �n{����\�����L�����L�stR�s� ������H�^� H������H��H����H������Hǅ����   H������H��H��\����   0���  �'� ������f�� f������H��H����H��p���Hǅx���   H��p���H��H��p���H��H�����   �z����\�����tR��� ������H��� H������H��H����H��`���Hǅh���   H��`���H��H��\����   0��a�  ��� ������f�r� f������H��H����H��P���HǅX���   H��P���H��H��p���H��H�¸����   ��y����\�����t`f�*� f������� �����H�� H�� ���H��H ���H��@���HǅH���   H��@���H��H��\����   0���  ��� ��"���f��� f�� ���H��H ���H��0���Hǅ8���   H��0���H��H��p���H��H�¼����   �3y����\�����t`f�~� f��L����m� ��H���H�X� H��@���H��H@���H�� ���Hǅ(���   H�� ���H��H��\����   0����  A���h���A�$��l���H��X���H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   Ic$H��HH�H������H������Hǅ������d Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��h���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H��p���H��H������I��I��P���E1��A  H����\�����t`f�\� f�������K� ������H�6� H������H��H����H�����Hǅ���   H�����H��H��\����   0���  H�C(H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   Hǅ���   Ic$H��HH�H�� ���H��(���H��H�$    H��H��p���H��H�Ƹ���H��H������1�E1�E1��uA  H����\�����t`f�,� f��\����� ��X���H�� H��P���H��HP���H�� ���Hǅ���   H�� ���H��H��\����   0��g�  H�C0H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Mc61�M��L��HH�H������Hǅ����   Hǅ����   Mc<$M��L��HH�H������H������H��H�$    H��H��p���H��H�Ƽ���H��H��p���1�E1�E1��A  H����\�����t`f��� f��������� ������H��� H������H��H����H������Hǅ����   H������H��H��\����   0��!�  D�����D�����H�C8H�����Hǅ���@   H�q     H�� ���H�    H��(���Hǅ0���    Hǅ8���    Hǅ@���   H�CHc 1�H��HH�H��H���HǅP���   HǅX���   Ic$H��HH�H��`���H��h���Hǅ������d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�����H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��p���H��H�Ɣ���H��H�����H��H������I��I������E1���?  H��H��H��@�����   I���fr��A��H��H�ǀ�����   �Or��A�DH������Lc�1�M��L��HH�H��H���I��I)�I���L������L����   L���r�����    H�Lc�E9�MO�M)�LH�   I��M��LO�M��LH�M)�LH�L��H��@���L���r��K�<7H������Hǅ�����mD Hǅ����   H������L���r����   L������L���|q�����    H�Lc�M9�MN�M)�LH�L�����L��L��L���<r��M��    L��L����p��H�R,$ H��������(3,$ ��)�������(,$ ��)�����H������H������H������H������H������H������H������H��H�Ɛ���H��H���������d �q����   I��I��@���L���p��A�ƾ�   H��H�ǀ����p��A�DH������Lc�M��L��A�    IH�H��H���I��I)�I���L������L���   L���Fp����AH�Lc�E9�MO�M)��    LH�I��A�   MN�M��LH�M)�LH�L��H��@���L����p��O�</Hǅ�����mD Hǅ����   H������L��L����p����   H�������o�����    I�Lc�M9�MN�M)�    LH�M�L��H��H�ƀ���L���wp��M�    L��L���4o��ǅp���    H������H������H������H������H������H�� H�D$H�D$    H�D$    H�$    H��H��p���H��H��p���1�E1�E1��c=  H�� ��\�����tLH�nf_open1H��`���H��H`���H������Hǅ����   H������H��H��\����   0���w��  H��H��@�����   I����w�n��A��H��H�ǀ�����   �xn��A�DH������Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���:n����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L����n��K�|% H������Hǅ���� nD Hǅ����   H������L���n����   H�������m�����    I�Lc�M9�MN�M)�    LH�L�����L��H��H�ƀ���L���bn��M�    L��L���m��ǅt���    L��p���H������H��x���H��p���H�� H�D$H�D$    H�D$    H�$    H��H��t���H��H��t���1�E1�E1��u;  H�� ��\�����tLH�nf_open2H������H��H����H��`���Hǅh���   H��`���H��H��\����   0���w�{�  H��H��@�����   I����w�l��A��H��H�ǀ�����   �jl��A�DH������Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���,l����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L����l��K�|% H������HǅP���nD HǅX���   H��P���L���l����   H�������k�����    I�Lc�M9�MN�M)�    LH�L�����L��H��H�ƀ���L���Tl��M�    L��L���k��ǅx���    L��@���H������H��H���H��@���H�� H�D$H�D$    H�D$    H�$    H��H��x���H��H��x���1�E1�E1��9  H�� ��\�����tLH�nf_open3H������H��H����H��0���Hǅ8���   H��0���H��H��\����   0���w�m�  H��H��@�����   I����w�sj��A��H��H�ǀ�����   �\j��A�DH������Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���j����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L����j��K�|% H������Hǅ ���nD Hǅ(���   H�� ���L���j����   H�������i�����    I�Lc�M9�MN�M)�    LH�L�����L��H��H�ƀ���L���Fj��M�    L��L���i��ǅ|���    L�����H������H�����H�����H�� H�D$H�D$    H�D$    H�$    H��H��|���H��H��|���1�E1�E1��7  H�� ��\�����tLH�nf_open4H������H��H����H�� ���Hǅ���   H�� ���H��H��\����   0���w�_�  H��H��@�����   I����w�eh��A��H��H�ǀ�����   �Nh��A�DH������Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���h����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L���h��K�|% H��x���Hǅ����nD Hǅ����   H������L���h����   H�������g�����    I�Lc�M9�MN�M)�    LH�L�x���L��H��H�ƀ���L���8h��M�    L��L����f��ǅ����    L������H������H������H������H�� H�D$H�D$    H�D$    H�$    H��H�ƀ���H��H����1�E1�E1��5  H�� ��\�����tLH�nf_open5H������H��H����H������Hǅ����   H������H��H��\����   0���w�Q�  H��H��@�����   I����w�Wf��A��H��H�ǀ�����   �@f��A�DH������Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���f����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L���f��K�|% H��p���Hǅ����nD Hǅ����   H������L���wf����   H�������ve�����    I�Lc�M9�MN�M)�    LH�L�p���L��H��H�ƀ���L���*f��M�    L��L����d��ǅ����    L������H������H������H������H�� H�D$H�D$    H�D$    H�$    H��H�Ƅ���H��H����1�E1�E1��3  H�� ��\�����tLH�nf_open7H�� ���H��H ���H������Hǅ����   H������H��H��\����   0���w�C�  H�snowmeltH�� ���H��H ���H������Hǅ����   H������H��H��p���H��H�����   ��w��c����\�����L��h���L��`���L�����L�stU��� ��H���H�y� H��@���H��H@���H������Hǅ����   H������H��H��\����   0���w�w�  f�@� f��d����/� ��`���H��H`���H��p���Hǅx���   H��p���H��H��t���H��H�����   ��w�
c����\�����tU�� ������H�չ H������H��H����H��`���Hǅh���   H��`���H��H��\����   0���w��  ��� ��������� ������H��H����H��P���HǅX���   H��P���H��H��x���H��H�����   ��w�Pb����\�����tU�H� ������H�3� H������H��H����H��@���HǅH���   H��@���H��H��\����   0���w��  ��� �������� ������H��H����H��0���Hǅ8���   H��0���H��H��|���H��H�����   ��w�a����\�����tU��� �����H��� H�� ���H��H ���H�� ���Hǅ(���   H�� ���H��H��\����   0���w�G�  ǅ ���evapH��H ���H�����Hǅ���   H�����H��H�ǀ���H��H� ����   ��w��`����\�����tU�
� ��H���H��� H��@���H��H@���H�� ���Hǅ���   H�� ���H��H��\����   0���w��  ��� ��d������ ��`���H��H`���H������Hǅ����   H������H��H�ǈ���H��H�¨����   ��w�0`����\�����tU�h� ������H�S� H������H��H����H������Hǅ����   H������H��H��\����   0���w���  A�������A�E ������ǅ����   ǅ����   L������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   IcU H��HH�H�� ���H�����Hǅ0��� �d Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɛ���H��H�°���H��H��0���I��I������E1���w�M-  H����\�����tYf�H� f�� �����()� ��)�����H��H����H������Hǅ����   H������H��H��\����   0���w��  A��� ���A�E ��$���ǅ(���   ǅ,���   L��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Ic1�H��HH�H��h���Hǅp���   Hǅx���   IcU H��HH�H������H������Hǅ����@�d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�� ���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H��t���H��H�Ɣ���H��H��0���H��H������I��I�����E1���w�2,  H����\�����tW�.� ��������(� ��)�p���H��Hp���H������Hǅ����   H������H��H��\����   0���w�m�  A�������A�E ������ǅ����   ǅ����   H������H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   IcU H��HH�H�� ���H�����Hǅ0�����d Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��x���H��H�Ƙ���H��H�°���H��H��0���I��I������E1���w�+  H����\�����tW�� �� �����(� ��)�����H��H����H������Hǅ����   H������H��H��\����   0���w�-�  A��� ���A�E ��$���ǅ(���   ǅ,���   H������H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Ic1�H��HH�H��h���Hǅp���   Hǅx���   IcU H��HH�H������H������Hǅ������d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�� ���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H��|���H��H�Ɯ���H��H��0���H��H������I��I�����E1���w��)  H����\�����tYf��� f��������(ή ��)�p���H��Hp���H������Hǅ����   H������H��H��\����   0���w��  A�������A�E ������ǅ����   ǅ����   H�� ���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic1�H��HH�H������Hǅ����   Hǅ����   IcU H��HH�H�� ���H�����Hǅ0��� �d Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H�ǀ���H��H�Ơ���H��H�°���H��H��0���I��I������E1���w��(  H����\�����tW�̬ �� �����(�� ��)�����H��H����H������Hǅ����   H������H��H��\����   0���w諹  A��� ���A�E ��$���ǅ(���   ǅ,���   H������H��0���L��L��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Hc1�H��HH�H��h���Hǅp���   Hǅx���   IcU H��HH�H������H������Hǅ����@�d Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�� ���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H�ǈ���H��H�ƨ���H��H��0���H��H������I��I�����E1���w�'  H����\�����tYf��� f��p�����(�� ��)�`���H��H`���H������Hǅ����   H������H��H��\����   0���w�_�  H��H��p�����w��Q����\�����tU�=� ������H�(� H������H��H����H��p���Hǅx���	   H��p���H��H��\����	   0���w��  H��H��t�����w�lQ����\�����tU�ܩ ������H�ǩ H������H��H����H��`���Hǅh���	   H��`���H��H��\����	   0���w�}�  H��H��x�����w��P����\�����tU�{� ������H�f� H������H��H����H��P���HǅX���	   H��P���H��H��\����	   0���w��  H��H��|�����w�P����\�����tU�� ������H�� H������H��H����H��@���HǅH���	   H��@���H��H��\����	   0���w蛵  H��H�ǀ�����w�P����\�����tU��� �����H��� H�����H��H���H��0���Hǅ8���	   H��0���H��H��\����	   0���w�*�  H��H�ǈ�����w�O����\�����tU�X� ��8���H�C� H��0���H��H0���H�� ���Hǅ(���	   H�� ���H��H��\����	   0���w蹴  H�CHc H���S  IcM H���F  H��1�H��H�� �����   H��H���H��H��1�H���r1���}m� D  ��|��X��u\��Y���|��}Y��|��|L ��XL ��u\L ��Y���|L ��}YL ��|L H��@H��H���|�H��y0��|��X��}\��}� ��Y���|��uY��|H��H)�H��|3��x���X���y\���(Ԝ ��Y���x���qY���x�H��H9�~/��{���X���{\���t� ��Y���{���sY���{�H��H���A\A]A^A_H��]H��[��w��    f.�     f.�     �L$ �    H�|$�H�L$�H�|$�0��pN��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0���M���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��M��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���L��H��H���[A\A]A^A_]Ðf.�     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L����L���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��~L��L�A$A$�$�  AH��H���  0�L��L��L���(K��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L����K���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��K��L�A$A$�$�  AH��H���  0�L��L��L���HJ��H��H���[A\A^A_]��     �L$ �    H�|$�H�L$�H�|$�0��pK��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0���J���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��J��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���I��H��H���[A\A]A^A_]�f��L$ �    H�|$�H�L$�H�|$�0��0J��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��I���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��HI��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L����G��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L���H���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��>H��L�A$A$�$�  AH��H���  0�L��L��L����F��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L���G���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��^G��L�A$A$�$�  AH��H���  0�L��L��L���F��H��H���[A\A^A_]��     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��F���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��hF��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L����D��H��H���[A\A]A^A_]�f��L$ �    H�|$�H�L$�H�|$�0��F��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��sE���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��(E��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���C��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L���iD���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��D��L�A$A$�$�  AH��H���  0�L��L��L����B��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L���C���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��>C��L�A$A$�$�  AH��H���  0�L��L��L����A��H��H���[A\A^A_]��     �L$ �    H�|$�H�L$�H�|$�0��C��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��sB���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��(B��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���@��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$���d �   0�L���iA���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��A��L�A$A$�$�  AH��H���  0�L��L��L����?��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L���@���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��>@��L�A$A$�$�  AH��H���  0�L��L��L����>��H��H���[A\A^A_]��     �L$ �    H�|$�H�L$�H�|$�0��@��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$�@�d �   0��s?���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��(?��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���=��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$���d �   0�L���i>���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��>��L�A$A$�$�  AH��H���  0�L��L��L����<��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$���d �   0�L���=���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��>=��L�A$A$�$�  AH��H���  0�L��L��L����;��H��H���[A\A^A_]��     �L$ �    H�|$�H�L$�H�|$�0��=��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �d �   0��s<���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��(<��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���:��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L���i;���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��;��L�A$A$�$�  AH��H���  0�L��L��L����9��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$���d �   0�L���:���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��>:��L�A$A$�$�  AH��H���  0�L��L��L����8��H��H���[A\A^A_]��     �L$ �    H�|$�H�L$�H�|$�0��:��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��s9���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��(9��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���7��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L���i8���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��8��L�A$A$�$�  AH��H���  0�L��L��L����6��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L���7���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��>7��L�A$A$�$�  AH��H���  0�L��L��L����5��H��H���[A\A^A_]��     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��6���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��H6��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L����4��H��H���[A\A]A^A_]�f��L$ �    H�|$�H�L$�H�|$�0���5��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��S5���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��5��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���3��H��H���[A\A]A^A_]�f�UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$� �d �   0�L���I4���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0���3��L�A$A$�$�  AH��H���  0�L��L��L���2��H��H���[A\A^A_]��     UH��AWAVATSH���H��@!  H��I��I��H�C0   H�CH   L�<$�@�d �   0�L���i3���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��3��L�A$A$�$�  AH��H���  0�L��L��L����1��H��H���[A\A^A_]��     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��s2���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��(2��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���0��H��H���[A\A]A^A_]�f��L$ �    H�|$�H�L$�H�|$�0���1���L$ �    H�|$�H�L$�H�|$�0��1���L$ �    H�|$�H�L$�H�|$�0��1���L$ �    H�|$�H�L$�H�|$�0��p1���L$ �    H�|$�H�L$�H�|$�0��P1���L$ �    H�|$�H�L$�H�|$�0��01��UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��0���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0��H0��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L����.��H��H���[A\A]A^A_]Ðf.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �d �   0��s/���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��(/��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���-��H��H���[A\A]A^A_]Ðf.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$�@�d �   0��S.���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�d ��  0��.��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���,��H��H���[A\A]A^A_]Ðf.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��3-���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0���,��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���m+��H��H���[A\A]A^A_]Ðf.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���d �   0��,���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���d ��  0���+��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���M*��H��H���[A\A]A^A_]Ðf.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �d �   0���*���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �d ��  0��*��A�EtI�E �� ��)$A�D$tI�$�� ��)�$�  L�$�$�  AH�$H��H���  0�L��L���-)��H��H���[A\A]A^A_]Ðf.�     SH��H���UUH�kH�l$H��AWAVAUATH��p"  I��I��L��X���L��P���H��@���H��8���H��H����H������Hǅ����	   H��H����H������Hǅ����	   H��Hp���H��`���Hǅh���	   H��HP���H��@���HǅH���	   H��H0���H�� ���Hǅ(���	   H��H���H�� ���Hǅ���	   H��H����H������Hǅ����	   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H+���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��Hl���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H(���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��Hh���H��P���HǅX���   H��HH���H��0���Hǅ8���   H��H(���H�����Hǅ���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HP���H��@���HǅH����   H��Hp���H��`���Hǅh����   H��HP���H��@���HǅH����   H�C0D�0H�C(� ��d���L�k@��oD ��   �   L����&��D����*���XSs L�s8H����  ���������L����pD ��   �   L���&��H����  L������M��M���`pD ��   �   L���p&��H��D��d�����  ��pD ��   �   L���I&��H���  �qD ��   �   L���)&��H��M���XJ  �gqD ��   �   L���&��H����O  ��qD ��   �   L����%��H����T  �rD ��   �   L����%��H����Y  L�������^ ������f�O f�������> ������H�) H������ǅ����+   �����H��H����H��0���Hǅ8���   H��0���H��H�������   0��Ɉ  �^  �����������������@�   ��L���)���Y�q ��y
���,�Dщ�l�������������	��  ��(��# ��)�������(��# ��)�������(v�# ��)�������(V�# ��)�������(6�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�p�����(��# ��)�`�����(��# ��)�P�����(v�# ��)�@�����(V�# ��)�0�����(6�# ��)� �����(�# ��)������(��# ��)� �����(��# ��)�����H��H`���H��0���Hǅ8����   H��0���H������Hǅ ����   H���# H�������(��# ��)� �����(t�# ��)�����H��H����H�����H��H������H��H������H��H��0�����"��L������L������������
��  ��(F�# ��)�������(&�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(f�# ��)�p�����(F�# ��)�`�����(&�# ��)�P�����(�# ��)�@�����(��# ��)�0�����(��# ��)� �����(��# ��)������(��# ��)� �����(f�# ��)�����H��H`���H�� ���Hǅ(����   H�� ���H������Hǅ ����   H�C�# H��P�����($�# ��)�@�����(�# ��)�0���H��H����H��H���H��H������H��H��0���H��H��0����I!����   L���< ��1Ʌ�H�Lc�L��A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L���� ��N��%O���Hǅ����oD Hǅ���   H�����L��L��� ��M�    L��L���d���x ��N���f��w f��L�����w ��H���H��w H��@���I��I��`�����   L���Y�����    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��N���L��L������N��-N���Hǅ ����oD Hǅ���   H�� ���L��L�����M�    L��L���|��f�9w f��L����(w ��H���H�w H��@���D��d���D��l���D��L�����  �gfff��d��������������   ����*�)���Xl ��Yl ��y
���,�Dщ�h�������������	��  ��(��# ��)�������(��# ��)�������(��# ��)�������(f�# ��)�������(F�# ��)�������(&�# ��)�������(�# ��)�������(��# ��)�p�����(��# ��)�`�����(��# ��)�P�����(��# ��)�@�����(f�# ��)�0�����(F�# ��)� �����(&�# ��)������(�# ��)� �����(��# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H���# H��������(��# ��)�������(��# ��)�p���H��H����H������H��H������H��H��p���H��H��0������L������L������������
��  ��(V�# ��)�������(6�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(v�# ��)�p�����(V�# ��)�`�����(6�# ��)�P�����(�# ��)�@�����(��# ��)�0�����(��# ��)� �����(��# ��)������(��# ��)� �����(v�# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H�S�# H��������(4�# ��)�������(�# ��)�����H��H����H������H��H������H��H�ư���H��H��0���������   L������1Ʌ�H�Lc�L��A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L���g��N��%O���Hǅ����pD Hǅ����   H������L��L���7��M�    L��L��������r ��N���f��r f��L�����r ��H���H��r H��@���I��I��`�����   L���������    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��P���L��L�����N��-P���Hǅ����?pD Hǅ����   H������L��L���O��M�    L��L�������(r �
  �gfffD�����������   D����z*�)���X�f ��Y�f ��y
���,�Dщ�h�������������	��  ��(�# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(n�# ��)�������(N�# ��)�������(.�# ��)�p�����(�# ��)�`�����(��# ��)�P�����(��# ��)�@�����(��# ��)�0�����(��# ��)� �����(n�# ��)������(N�# ��)� �����(.�# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H��# H�������(��# ��)� �����(��# ��)�����H��H����H�����H��H������H��H������H��H��0����Q��L������������
M����  ��(��# ��)�������(��# ��)�������(b�# ��)�������(B�# ��)�������("�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�p�����(��# ��)�`�����(��# ��)�P�����(b�# ��)�@�����(B�# ��)�0�����("�# ��)� �����(�# ��)������(��# ��)� �����(��# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H���# H��P�����(��# ��)�@�����(`�# ��)�0���H��H����H��H���H��H������H��H��0���H��H��0��������   L�����1Ʌ�H�Lc�A��   I���   MO�L��A��   M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L���-��N��%O���Hǅ����ppD Hǅ����   H������L��L������M�    L��L������n ��N���f�n f��L�����m ��H���H��m H��@���I��I��`�����   L��������    I�Lc�I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��O���L��L���K��N��%O���Hǅ�����pD Hǅ����   H������L��L�����M��    L��L�������Hm ��N���f�9m f��L����(m ��H���H�m H��@����7  ����*D����������@�   D����z*�)���X�a ��Y�a ��y
���,�Dщ�h�������������	��  ��(.�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(��# ��)�������(n�# ��)�������(N�# ��)�p�����(.�# ��)�`�����(�# ��)�P�����(��# ��)�@�����(��# ��)�0�����(��# ��)� �����(��# ��)������(n�# ��)� �����(N�# ��)�����H��H`���H��p���Hǅx����   H��p���H������Hǅ ����   H�+�# H��������(�# ��)�������(��# ��)�p���H��H����H������H��H������H��H��p���H��H��0�������L������������
M����  ��(��# ��)�������(��# ��)�������(��# ��)�������(b�# ��)�������(B�# ��)�������("�# ��)�������(�# ��)�������(��# ��)�p�����(��# ��)�`�����(��# ��)�P�����(��# ��)�@�����(b�# ��)�0�����(B�# ��)� �����("�# ��)������(�# ��)� �����(��# ��)�����H��H`���H��`���Hǅh����   H��`���H������Hǅ ����   H���# H��������(��# ��)�������(��# ��)�����H��H����H������H��H������H��H�ư���H��H��0����E����   L���8��1Ʌ�H�Lc�L��A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L������N��%O���HǅP����pD HǅX���   H��P���L��L�����M�    L��L���`���i ��N���f��h f��L�����h ��H���H��h H��@���I��I��`�����   L���U�����    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��P���L��L������N��-P���Hǅ@����pD HǅH���   H��@���L��L�����M�    L��L���x����(0h ��)�@���D��L���E��D��h���D��l���D��h���H���# H�������(��# ��)� �����(��# ��)�����L�{ L�����H��H������H��H��0������d �2��D������D��������d���������D������H��H����H��0���Hǅ8���    H�`     H��@���H�      H��H���H��P���HǅX���    Hǅ`���   Hǅh���   Hǅp���   ��(y�# ��)�������(Y�# ��)�������(9�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�������(��# ��)�����H��H0���H������H��H`���H�� ���Hǅ(����   H�� ���H������Hǅ�����   L������H��H�Ɛ���� �d I��I��0���L������H���# H��P�����(��# ��)�@�����(��# ��)�0�����(��# ��)� �����(�# ��)����H��H@���H�����I��Hǅ����   H�����H��(���Hǅ0����   H��H@���H�� ���I��Hǅ����   H�� ���H��H���HǅP����   H��H��������d L��������   L��M������A�ƾ�   L������A�DH��0���Lc�E1�M��L��IH�H��H���I��I)�I���L���   L�������AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L���:��K�<,H������Hǅ�����rD Hǅ����   H������L���
����   H��@����	�����    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L�����M�    L��L���z
��ǅ����    L������H��0���H������H������H�� H�D$H�D$    H�D$    H�$    H��H������H��H��p���1�E1�E1���D  H�� ��l�����tPH�nf_open1H�������l���H��H����H������Hǅ����   H������H��H��l����   0���n  H��H��@�����   I����	��A��H��H��@�����   ��	��A�DH��(���Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���	����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L���,
��K�|% H������Hǅ�����rD Hǅ����   H������L����	����   H��@����������    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L���	��M�    L��L���k��ǅ����    L������H��(���H������H������H�� H�D$H�D$    H�D$    H�$    H��H������H��H��t���1�E1�E1��C  H�� ��l�����tPH�nf_open2H�������l���H��H����H������Hǅ����   H������H��H��l����   0���l  H��H��@�����   I������A��H��H��@�����   ���A�DH�� ���Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���w����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L�����K�|% H������Hǅ�����rD Hǅ����   H������L��������   H��@����������    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L�����M�    L��L���\��ǅ����    L������H�� ���H������H������H�� H�D$H�D$    H�D$    H�$    H��H������H��H��x���1�E1�E1��A  H�� ��l�����tPH�nf_open3H�������l���H��H����H��p���Hǅx���   H��p���H��H��l����   0��j  H��H��@�����   I�����A��H��H��@�����   ���A�DH�����Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���h����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L�����K�|% H������Hǅ`����rD Hǅh���   H��`���L��������   H��@����������    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L�����M�    L��L���M��ǅ����    L��P���H�����H��X���H��P���H�� H�D$H�D$    H�D$    H�$    H��H������H��H��|���1�E1�E1��#?  H�� ��l�����tPH�nf_open4H������l���H��H���H��@���HǅH���   H��@���H��H��l����   0��h  H��H��@�����   I�����A��H��H��@�����   ���A�DH�����Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���Y����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L������K�|% H������Hǅ0����rD Hǅ8���   H��0���L��������   H��@����������    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L�����M�    L��L���>��ǅ����    L�� ���H�����H��(���H�� ���H�� H�D$H�D$    H�D$    H�$    H��H������H��H����1�E1�E1��4=  H�� ��l�����tPH�nf_open5H��0����l���H��H0���H�����Hǅ���   H�����H��H��l����   0��f  H��H��@�����   I�����A��H��H��@�����   ���A�DH�����Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���J����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L������K�|% H������Hǅ ����rD Hǅ���   H�� ���L�������   H��@���� �����    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L���r��M�    L��L���/ ��ǅ����    L������H�����H������H������H�� H�D$H�D$    H�D$    H�$    H��H������H��H����1�E1�E1��E;  H�� ��l�����tPH�nf_open6H��P����l���H��HP���H������Hǅ����   H������H��H��l����   0��d  H��H��@�����   I������A��H��H��@�����   �y���A�DH�� ���Lc�E1�M��L��IH�H��H���I��I)�I���L���   L���;�����AH�Lc�E9�MO�M)�    LH�A�   I��MN�M��LH�M)�LH�L��H��@���L�������K�|% H������Hǅ�����rD Hǅ����   H������L��������   H��@����������    I�Lc�M9�MN�M)�    LH�L�����L��H��H��@���L���c���M�    L��L��� ���ǅ����    L������H�� ���H������H������H�� H�D$H�D$    H�D$    H�$    H��H������H��H����1�E1�E1��V9  H�� ��l�����tPH�nf_open7H��p����l���H��Hp���H������Hǅ����   H������H��H��l����   0��xb  H�snowmeltH������H��H����H������Hǅ����   H������H��H��p���H��H�����   ������l�����L������L������L��@���L��8���te�WX �������GX ������H�2X H�������l���H��H����H������Hǅ����   H������H��H��l����   0��a  f��W f��������W ������H��H����H������Hǅ����   H������H��H��t���H��H�����   �2�����l�����te��W ��������W ������H��W H�������l���H��H����H��p���Hǅx���   H��p���H��H��l����   0���`  H�snowfallH�����H��H���H��`���Hǅh���   H��`���H��H��x���H��H�����   �r�����l�����te��V ��<�����V ��8���H��V H��0����l���H��H0���H��P���HǅX���   H��P���H��H��l����   0��`  ǅP���evapH��HP���H��@���HǅH���   H��@���H��H��|���H��H�����   ������l�����te�MV ��|����=V ��x���H�(V H��p����l���H��Hp���H��0���Hǅ8���   H��0���H��H��l����   0��Z_  ��U ��������U ������H��H����H�� ���Hǅ(���   H�� ���H��H�ǀ���H��H� ����   �������l�����te��U ��������U ������H�yU H�������l���H��H����H�����Hǅ���   H�����H��H��l����   0��^  �>U �������.U ������H��H����H�� ���Hǅ���   H�� ���H��H�Ǆ���H��H�¤����   �+�����l�����te��T ��������T ������H��T H�������l���H��H����H������Hǅ����   H������H��H��l����   0���]  ��T ������T �����H��H���H������Hǅ����   H������H��H�ǈ���H��H�¨����   �d�����l�����te�@T ��<����0T ��8���H�T H��0����l���H��H0���H������Hǅ����   H������H��H��l����   0��]  ��h�����P�����l�����T���ǅX���   ǅ\���   ǅ`���   ǅd���   ǅh���   H�C H�����l���L��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Hc1�H��HH�H������Hǅ����   H��P���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��`���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��p���H��H�Ɛ���H��H��p���H��H������I��I��0���E1���1  H����l�����tg��Q ������f��Q f������H��Q H�������l���H��H����H������Hǅ����   H������H��H��l����   0��Z  H��H��0���� �d � �d ������h�����������l���������ǅ����   ǅ����   ǅ����   ǅ����   ǅ����   H�C H���������L������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Hc1�H��HH�H�����Hǅ���   H������H��0���Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��t���H��H�Ɣ���H��H������H��H��0���I��I������E1��S0  H����l�����tg��O ������f��O f������H�tO H�������l���H��H����H������Hǅ����   H������H��H��l����   0��>X  H��H��0������d ���d �E�����h����������l��������ǅ���   ǅ���   ǅ ���   ǅ$���   ǅ(���   H�C H�����,���L��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Hc1�H��HH�H��h���Hǅp���   H�����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�� ���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��x���H��H�Ƙ���H��H��0���H��H������I��I������E1���.  H����l�����tg�9M ��Z���f�*M f��X���H�M H��P����l���H��HP���H������Hǅ����   H������H��H��l����   0���U  H��H��0���� �d � �d �������h�����p�����l�����t���ǅx���   ǅ|���   ǅ����   ǅ����   ǅ����   H�C H���������H��P���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Hc1�H��HH�H������Hǅ����   H��p���H������Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H������H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��|���H��H�Ɯ���H��H����H��H������I��I��P���E1��l-  H����l�����tg��J ������f��J f������H��J H�������l���H��H����H������Hǅ����   H������H��H��l����   0��WS  H��H��0������d ���d �^�����h�����������l���������ǅ����   ǅ����   ǅ����   ǅ����   ǅ����   H�C H���������H�CH������Hǅ����@   H�q     H�� ���H�    H�����Hǅ���    Hǅ���    Hǅ ���   Hc1�H��HH�H��(���Hǅ0���   H������H��P���HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H�ǀ���H��H�Ơ���H��H������H��H��P���I��I������E1���+  H����l�����tg�nH �����f�_H f�����H�IH H������l���H��H���H������Hǅ����   H������H��H��l����   0���P  H��H��0���� �d � �d �������h�����0�����l�����4���ǅ8���   ǅ<���   ǅ@���   ǅD���   ǅH���   H�C H�����L���H��X���H��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   Hc1�H��HH�H������Hǅ����   H��0���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��@���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H�Ǆ���H��H�Ƥ���H��H��P���H��H������I��I�����E1��*  H����l�����tg�F ��z���f��E f��x���H��E H��p����l���H��Hp���H��p���Hǅx���   H��p���H��H��l����   0��lN  H��H��0������d �� e �s�����h�����������l���������ǅ����   ǅ����   ǅ����   ǅ����   ǅ����   H�C H���������H�CH������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Hc1�H��HH�H������Hǅ����   H������H�����Hǅ���    H�q     H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H������H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H�ǈ���H��H�ƨ���H��H�°���H��H�����I��I��p���E1��)  H����l�����tg��C ������f��C f������H�~C H�������l���H��H����H��`���Hǅh���   H��`���H��H��l����   0���K  H��H��0���� e � e �����H��H��p����`�����l�����tY�C ������H�C H�������l���H��H����H��P���HǅX���	   H��P���H��H��l����	   0��mK  H��H��t����������l�����tY��B �����H��B H�� ����l���H��H ���H��@���HǅH���	   H��@���H��H��l����	   0���J  H��H��x����|�����l�����tY�TB ��(���H�?B H�� ����l���H��H ���H��0���Hǅ8���	   H��0���H��H��l����	   0��J  H��H��|����
�����l�����tY��A ��H���H��A H��@����l���H��H@���H�� ���Hǅ(���	   H�� ���H��H��l����	   0��J  H��H�ǀ���������l�����tY��A ��h���H�{A H��`����l���H��H`���H�����Hǅ���	   H�����H��H��l����	   0��I  H��H�Ǆ����&�����l�����tY�.A ������H�A H�������l���H��H����H�� ���Hǅ���	   H�� ���H��H��l����	   0��3I  H��H�ǈ���������l�����tY��@ ������H��@ H�������l���H��H����H������Hǅ����	   H������H��H��l����	   0���H  Hǅ����    Hǅx���    Hǅ����    Hǅp���    H�C LcM���L  M��I���1�A��H��X���H��P���L�SH�M  L��H����y����1�I����   1� ��}YL A��   ��|L ��}YA��   ��|��}YA��   ��|��Y��   ����Y��   ����}YA��   ��|��}YL ��|L ��}YL ��|L ��}YL ��|L ��YL ��L ��YL ��L ��}YL ��|L H��@H��H����$���L��H��yI��}YL ��|L ��}Y��|��}Y��|��Y����Y����}Y��|L��L��H)�H��|W��W���{*�����yYL� ��xL� ��yY���yY���x���x���Y������Y������yY���x�H��I9�~O��W���{*��{YL� ��{L� ��{Y���{���{Y���Y���{������Y������{Y���{�1�A���4  L��H��1�����I��}�A����������H�{�{  1������b})/ �b}%(/ �A���A������H�{f.�     f.�     f.�     �|)�p�����|| ALx��E����)����������cEK� �A|D �A|ALx��=����)������c=K� �A|,�A|ALx��-����)������A|(��c-K� �A|4��<Lx��-,��-��\���\���|��uK����}����U����K�P�Av��A����mK����X���eK�PALx��uK�P��|��|��|��|D ��X���(�������]K� ��}����}K� ��|l ��|l ��X���uK� ��U����UK� ��|L ��5X���(�������5K�p��|| �AE����X��cmK� ��EKƀ��|D ŝ�T �|(�p�����X��|(�������K��Bm-L �Bm-T �A5\��AeX��A5���5����eKڐ��X��cUK����\���|| ��EK� ��UK����|l ��|t ��uKȐ��m����uK� ��|L ��|D H��@H��H����������(��|)�H����   ��|l ��}�, ������)������A-����UK���|L ��|������uK� ��|<�A|�=���C=KҐ�A|�b}V, �%��B%-$�B%-,�A\�ŝ����}���-\��A|4�CK���Av��A%����EK� ��EK����X��|(�������MK����X���]K� ��X��cK����eX���eK� ��K%�eD  �A|,��|��|<��}���X���y���X����p�����X���}���X���y���X���X��c}���YX���y���X���X���}���X���y���X���X����p��������������������x���L��I9���   ��#+ H�K�    ��{L� ��.+ s��X�p������p���I�D�     ��{���.�* s��X������������I��    ��{���.�* s��X������������I��    ��{���.�s.��s\���{���.u* s��X�x������x���I��    ��{�H��L9��9���H���# H��������(s�# ��)�������(S�# ��)�������(3�# ��)�������(�# ��)�p���H��H����H������H��H��p���I��I��0�����e L����w�����H���# H�������(��# ��)� �����(��# ��)�������(q�# ��)�������(Q�# ��)�����H��Hx���H�����H��H�������@e L���b���H��# H��p�����(��# ��)�`�����(��# ��)�P�����(��# ��)�@�����(��# ��)�0���H��H����H��h���H��H��0���� e L�������H�f�# H�E���(J�# ��)E���(-�# ��)E���(�# ��)E���(��# ��)E�H��Hp���H�E�H��H�Ɛ��e L������H��H���A\A]A^A_H��]H��[Ë�L���������Ѓ���   )����������Y,( ��y
���,�Dщ�l�������������	��  ��(P�# ��)�������(0�# ��)�������(�# ��)�������(�# ��)�������(Ь# ��)�������(��# ��)�������(��# ��)�������(p�# ��)�p�����(P�# ��)�`�����(0�# ��)�P�����(�# ��)�@�����(�# ��)�0�����(Ы# ��)� �����(��# ��)������(��# ��)� �����(p�# ��)�����H��H`���H��0���Hǅ8����   H��0���H������Hǅ ����   H�M�# H�������(.�# ��)� �����(�# ��)�����H��H����H�����H��H������H��H������H��H��0�������L������������
��  ��(�# ��)�������(Ǭ# ��)�������(��# ��)�������(��# ��)�������(g�# ��)�������(G�# ��)�������('�# ��)�������(�# ��)�p�����(�# ��)�`�����(ǫ# ��)�P�����(��# ��)�@�����(��# ��)�0�����(g�# ��)� �����(G�# ��)������('�# ��)� �����(�# ��)�����H��H`���H�� ���Hǅ(����   H�� ���H������Hǅ ����   H��# H��P�����(ū# ��)�@�����(��# ��)�0���H��H����H��H���H��H������H��H��0���H��H��0����������   L�������1Ʌ�H�Lc�A��   I���   MO�L��A��   M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L���r���N��%O���Hǅ��� qD Hǅ���   H�����L��L���B���M�    L��L��������0 ��N���f��/ f��L�����/ ��H���H��/ H��@���I��I��`�����   L����������    I�Lc�I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��O���L��L������N��%O���Hǅ ���OqD Hǅ���   H�� ���L��L���`���M��    L��L�������=/ ��N���f�./ f��L����/ ��H���H�/ 鐶���gfff��L��������������   ��)����������Y�" ��y
���,�Dщ�l�������������	��  ��({�# ��)�������([�# ��)�������(;�# ��)�������(�# ��)�������(��# ��)�������(۩# ��)�������(��# ��)�������(��# ��)�p�����({�# ��)�`�����([�# ��)�P�����(;�# ��)�@�����(�# ��)�0�����(��# ��)� �����(ۨ# ��)������(��# ��)� �����(��# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H�x�# H��������(Y�# ��)�������(9�# ��)�p���H��H����H������H��H������H��H��p���H��H��0����>���L������������
��  ��(�# ��)�������(�# ��)�������(ҩ# ��)�������(��# ��)�������(��# ��)�������(r�# ��)�������(R�# ��)�������(2�# ��)�p�����(�# ��)�`�����(�# ��)�P�����(Ҩ# ��)�@�����(��# ��)�0�����(��# ��)� �����(r�# ��)������(R�# ��)� �����(2�# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H��# H��������(�# ��)�������(Ш# ��)�����H��H����H������H��H������H��H�ư���H��H��0���������   L������1Ʌ�H�Lc�L��A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L���#���N��%O���Hǅ����pqD Hǅ����   H������L��L�������M�    L��L�������+ ��N���f��* f��L�����* ��H���H��* H��@���I��I��`�����   L���������    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��P���L��L���;���N��-P���Hǅ�����qD Hǅ����   H������L��L������M�    L��L���������(0* ��)�@����_����gfff��L��������������   ��)����������Y� ��y
���,�Dщ�l�������������	��  ��(ç# ��)�������(��# ��)�������(��# ��)�������(c�# ��)�������(C�# ��)�������(#�# ��)�������(�# ��)�������(�# ��)�p�����(æ# ��)�`�����(��# ��)�P�����(��# ��)�@�����(c�# ��)�0�����(C�# ��)� �����(#�# ��)������(�# ��)� �����(�# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H���# H�������(��# ��)� �����(��# ��)�����H��H����H�����H��H������H��H������H��H��0�������L������������
��  ��(Z�# ��)�������(:�# ��)�������(�# ��)�������(��# ��)�������(ڦ# ��)�������(��# ��)�������(��# ��)�������(z�# ��)�p�����(Z�# ��)�`�����(:�# ��)�P�����(�# ��)�@�����(��# ��)�0�����(ڥ# ��)� �����(��# ��)������(��# ��)� �����(z�# ��)�����H��H`���H������Hǅ�����   H������H������Hǅ ����   H�W�# H��P�����(8�# ��)�@�����(�# ��)�0���H��H����H��H���H��H������H��H��0���H��H��0����]�����   L���P���1Ʌ�H�Lc�A��   I���   MO�L��A��   M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L�������N��%O���Hǅ�����qD Hǅ����   H������L��L������M�    L��L���r����*& ��N���f�& f��L����
& ��H���H��% H��@���I��I��`�����   L���g������    I�Lc�I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��O���L��L������N��%O���Hǅ�����qD Hǅ����   H������L��L�������M��    L��L�������`% ��N���f�Q% f��L����@% ��H���H�+% ��������*��d�������������@�   ����*�)���XB ��Y> ��y
���,�Dщ�h�������������	��  ��(�# ��)�������(ʤ# ��)�������(��# ��)�������(��# ��)�������(j�# ��)�������(J�# ��)�������(*�# ��)�������(
�# ��)�p�����(�# ��)�`�����(ʣ# ��)�P�����(��# ��)�@�����(��# ��)�0�����(j�# ��)� �����(J�# ��)������(*�# ��)� �����(
�# ��)�����H��H`���H��p���Hǅx����   H��p���H������Hǅ ����   H��# H��������(ȣ# ��)�������(��# ��)�p���H��H����H������H��H������H��H��p���H��H��0�������L������������
��  ��(��# ��)�������(a�# ��)�������(A�# ��)�������(!�# ��)�������(�# ��)�������(�# ��)�������(��# ��)�������(��# ��)�p�����(��# ��)�`�����(a�# ��)�P�����(A�# ��)�@�����(!�# ��)�0�����(�# ��)� �����(�# ��)������(��# ��)� �����(��# ��)�����H��H`���H��`���Hǅh����   H��`���H������Hǅ ����   H�~�# H��������(_�# ��)�������(?�# ��)�����H��H����H������H��H������H��H�ư���H��H��0���������   L�������1Ʌ�H�Lc�L��A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H��O���L������N��%O���HǅP��� rD HǅX���   H��P���L��L���b���M�    L��L�������'! ��N���f�! f��L����! ��H���H��  H��@���I��I��`�����   L���������    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H��P���L��L������N��-P���Hǅ@���OrD HǅH���   H��@���L��L���z���M�    L��L���7�����(O  麶��f.�     �L$ �    H�|$�H�L$�H�|$�0������L$ �    H�|$�H�L$�H�|$�0��p����L$ �    H�|$�H�L$�H�|$�0��P����L$ �    H�|$�H�L$�H�|$�0��0����L$ �    H�|$�H�L$�H�|$�0������L$ �    H�|$�H�L$�H�|$�0�������L$ �    H�|$�H�L$�H�|$�0������UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��	e �   0��;���A�E8��$�  H��H�Ǆ  ��	e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��	e �   0��;���A�E8��$�  H��H�Ǆ  ��	e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� 
e �   0��;���A�E8��$�  H��H�Ǆ  � 
e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@
e �   0��;���A�E8��$�  H��H�Ǆ  �@
e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L��蛿��H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��
e �   0��;���A�E8��$�  H��H�Ǆ  ��
e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L��蛾��H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��
e �   0��;���A�E8��$�  H��H�Ǆ  ��
e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L��蛽��H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� e �   0��;���A�E8��$�  H��H�Ǆ  � e ��  0������CtH��� ��)$A�D$tI�$�� ��)�$�  M�E $�$�  AH�$H��H���  0�L��L��蛼��H��H���[A\A]A^A_]�D  f.�     f.�     UH��AWAVSH���H��   H��H�    Lc2M��I���L�]L�UM����  1���W�A����  L��H��1�����I���)  1�L��H���������}�
 �f.�     ��$�x  ������v�������eK�P��eK�P��$��l ������eK�`��eK�`��X���X���l ��d@��  ������eK�P��eK�P��X���d@��d`������eK�P��eK�P��X���d`����   ��  ������eK�P��eK�P��X�����   ����   ������eK�P��eK� ��X�����   H���   H��H���H�������H��yD��}�	 f.�     ��������v�������UK�0��UK�0����X�H�� H��x���}���X���y���X�L��I9�~2��{	 �f.�     ���������T������X�H��L9�|�����W���*����^���.E	 �K  A�   A�A� H�K�# H�D$`��(.�# ��)D$P��(�# ��)D$@��(�# ��)D$0��(ԯ# ��)D$ ��(��# ��)D$��(��# ��)$H�t$8H�4$H��H�   �@e ��w賻��E���x  A��  �V  1�A��|kL��H��1�I��|01�D  ��������D ��D@��D`H��H��H���|�H��y%�     f.�     ������H�� H��x�D��Hc�L��H)�H��|Hc���W���˃�A9���   H�H��    ��   A�    A�A� H� �# H��$�  ��(�# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  ��(}�# ��)�$�  ��(\�# ��)�$�  ��(;�# ��)�$�  H��$�  H��H�ƀ  H��H�   �@e ��w�I�������e 0�H��L����w袹��H��H���[A^A_]��wÐSH��H���UUH�kH�l$H��AWAVAUATH��p  H������H������L������L������H������H��H����H��p���Hǅx���	   H��H`���H��P���HǅX���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H ���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��Hk���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H,���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��HP���H��@���HǅH����   H��H����H������Hǅ�����   L�s ��   L���u���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L������J��=O���H������Hǅp��� tD Hǅx���   H��p���L������H�{(��   �������    H�Lc�M9�MN�M)�LH�A�   I��ML�M��LH�M)�LH�L�����L��H�s(L��膶��M�Hǅ`���9tD Hǅh���   H��`���L��L���[���M��    L��L�������  ��N���f� f��L����  ��H���H�� H��@���L�s0��   L���������    H�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�M)�LH�H������L��L��諵��J��%����H������HǅP���OtD HǅX���   H��P���L���w���H�{8��   �y������    H�Lc�M9�MN�M)�LH�A�   I��MN�M��LH�1�M)�LH�L�����L��H�s8L������M�Hǅ@���TtD HǅH���   H��@���L��L������M��    L��L��謳����   H��H�ǀ����س��A��H��H��@�����   ��������    I�D�H������Hc�Lc�A9�LO�I��M)�    LH�E��DH�Mc�M9�MN�M)�LH�H��HH�H��H���I��I)�I���L��L��H��H��@���L���8���O�$<L��H��H�ƀ���L������M�    L��L���ܲ��H�u�# H��������(V�# ��)�������(6�# ��)�����L��0���H������H��8���H��0���H������H������H��H�Ɛ���H��H��������e 観����   I��I��@���L��菲��A�ƾ�   I��I�ǀ���L���u���D�H������Lc�M��L��    HH�H��H���I��I)�I���L���   L���7������    H�Lc�E9�MO�L��H��@���L������M)�    LH��   L����������    I�Lc�M9�MN�M)��    LH�O�l% L��H������L��豲��M�    L��L���n���ǅ����    L�� ���H������H��(���H�� ���H�� H�D$H�D$    H�D$    H�$    H��H�ƌ���H��H��p���1�E1�E1��  H�� ��l���L�sL�{��tPH�nf_open1H��p����l���H��Hp���H�����Hǅ���   H�����H��H��l����   0��  ǅ����densH��H����H�� ���Hǅ���   H�� ���H��H��p���H��H�����   �d�����l�����te�h �������X ������H�C H�������l���H��H����H������Hǅ����   H������H��H��l����   0��  ǅ����tempH��H����H������Hǅ����   H������H��H��p���H��H�����   諯����l�����L������te�� �������� ������H�� H�������l���H��H����H������Hǅ����   H������H��H��l����   0��E  ǅ���massH��H���H������Hǅ����   H������H��H��p���H��H�����   ������l�����te� ��<����� ��8���H�� H��0����l���H��H0���H������Hǅ����   H������H��H��l����   0��  �� ��T����� ��P���H��HP���H������Hǅ����   H������H��H��p���H��H�����   �$�����l�����te�` ��|����P ��x���H�; H��p����l���H��Hp���H������Hǅ����   H������H��H��l����   0���  ��
 ������f��
 f������H��H����H������Hǅ����   H������H��H��p���H��H� ����   �[�����l�����te��
 ��������
 ������H��
 H�������l���H��H����H��p���Hǅx���   H��p���H��H��l����   0���  �G
 �������7
 ������H��H����H��`���Hǅh���   H��`���H��H��p���H��H��l����   �ԭ����l�����te��	 ��������	 ������H��	 H�������l���H��H����H��P���HǅX���   H��P���H��H��l����   0��5  Hǅ@���    HǅH���    H��@���H��H��p���H��H��l���E1�L���Y�����l�����tg�[	 �����f�L	 f�����H�6	 H������l���H��H���H��0���Hǅ8���   H��0���H��H��l����   0��  A�E ��,���H������H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   IcE 1�H��HH�H��h���Hǅp���   Hǅ�����e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��p���H��H�Ɛ���H��H��0���H��H������I��I������E1��
  H����l�����tg�^ ��Z���f�O f��X���H�9 H��P����l���H��HP���H�� ���Hǅ(���   H�� ���H��H��l����   0��{  A�E ��l���H������H��p���Hǅx���@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   IcE 1�H��HH�H������Hǅ����   Hǅ�����e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��l���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��p���H��H�Ɣ���H��H��p���H��H������I��I��0���E1��c	  H����l�����tg�a ������f�R f������H�< H�������l���H��H����H�����Hǅ���   H�����H��H��l����   0��n  A�E ������H������H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   IcE 1�H��HH�H������Hǅ����   Hǅ��� e Hǅ���    H�q     H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H������H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ƙ���H��H�°���H��H�����I��I��p���E1��6  H����l�����tg�d ������f�U f������H�? H�������l���H��H����H�� ���Hǅ���   H�� ���H��H��l����   0��a
  A�E ������H������H������I��Hǅ����@   H�q     H�� ���H�    H�����Hǅ���    Hǅ���    Hǅ ���   IcE 1�H��HH�H��(���Hǅ0���   HǅP���@e HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��p���H��H�Ɯ���H��H������H��H��P���I��I������E1��  H����l�����tg�d �����f�U f�����H�? H������l���H��H���H������Hǅ����   H������H��H��l����   0��Q  A�E ��,���L��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   IcE 1�H��HH�H��h���Hǅp���   Hǅ�����e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��p���H��H�Ơ���H��H��0���H��H������I��I������E1���  H����l�����tg�n� ��Z���f�_� f��X���H�I� H��P����l���H��HP���H������Hǅ����   H������H��H��l����   0��K  H��H��p����̠����l�����tY��� ��x���H��� H��p����l���H��Hp���H������Hǅ����	   H������H��H��l����	   0���  McE ā{D����X�ā{D���   D9��  O��K��D��I��1�����   M��I��I��>M�I���1����� �H��H��H��H)���(�����R���B��@�����X���X�L��H)���G��X�����(�����R���XB���X���G���(�����R���XB���X���G���(�����R���XB���X���G�H��L9��n���L9�}<��	� �H��H��H��H)���(�����R���XB���X�L��H)���B�H��L9�|�H��H���A\A]A^A_H��]H��[Ðf.�     f.�     �L$ �    H�|$�H�L$�H�|$�0��П��UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� e �   0��;���A�E8��$�  H��H�Ǆ  � e ��  0������CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L��裝��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@e �   0��[���A�E8��$�  H��H�Ǆ  �@e ��  0��5����CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���Ü��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��e �   0��{���A�E8��$�  H��H�Ǆ  ��e ��  0��U����CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$��e �   0�蛜��A�E8��$�  H��H�Ǆ  ��e ��  0��u����CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� e �   0�軛��A�E8��$�  H��H�Ǆ  � e ��  0�蕛���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���#���H��H���[A\A]A^A_]ÐUH��AWAVSH���H��  H��I��H��H��H���   H��$�   HǄ$�      H��H��`H�L$PH�D$XP   �8 �\  H��H��PH�L$@I��H�D$HP   H�|$@�P   H��轙��H��# H��$`  ����(��# ��)�$P  ��(ӗ# ��)�$@  ��(��# ��)�$0  ��(��# ��)�$   ��(p�# ��)�$  ��(O�# ��)�$   ��(.�# ��)�$�   ��(�# ��)�$�   L�t$0H�D$8H�L$0H��$  H��$   L�|$ H�D$(P   H�D$ H��$X  HǄ$`  P   H��H���   H��H�  ��e 蠚��Ǆ$�       H��H�   H�D$H�D$   H�|$1��   0�����H��H���[A^A_]�f.�     f.�     f�     SH��H���UUH�kH�l$H��AWAVAUATH��	  H��p���H��`���L��x���L��h���H��X���H������H��H����H������Hǅ�����   L���   ��   L�����1Ʌ�H�Lc�A��   I���   MO�M)�LH��   I��M��LO�M��LH�M)�LH�H������L��L���`���N��-����Hǅ�����wD Hǅ����   H������L��L���0���M�    L��L��������� ������f�v� f�������e� ������H�P� H������Hǅ���   H��H��������   �ڗ��A��E��A�    EI�H���   ��   蹗��A��E�H���   ��   袗��B�D H��P���H�H��H���Mc�A9�LO�I��M)��    LH�I��A�   MN�M��LH�M)�LH�H��HH�H��H���I��I)�I���L��@���L��L��H��H������L������K�<.H��8���Hǅ�����wD Hǅ����   H������L���������   H���   �ߖ�����    H�Lc�M9�MN�M)�LH�A�   I��ML�M��LH�L��0���M)�LH�L�8���L��H���   L���x���M�Hǅ�����wD Hǅ����   H������L��L���M�����   H���   �L������    I�Lc�M9�MN�M)�    LH�A�   I��MN�M��LH�1�M)�LH�L�0���L��H���   L������M�Hǅ�����wD Hǅ����   H������L��L��躖��M��    L��L���w���H�       H�����H��H���H�����H������Hǅ ���    Hǅ(���    H��@���H������H��P���H������H������H��0���H��H���H��8���Hǅ@���    HǅP���    Hǅ`���    Hǅp���    Hǅx���    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ ���    Hǅ���    Hǅ ���    Hǅ0���    Hǅ8���    Hǅ@���    ����H��H�����0��(���I��I�ư�����e � e L���̔���@e ��e L��躔��H�Ö# H��������(��# ��)�������(��# ��)�����H��X���H������H��H�ư�����e L���g���H��# H�������(ї# ��)� �����(��# ��)�����H��`���H�����H��H�������@e L������H��# H��P�����(��# ��)�@�����(ޘ# ��)�0���H��h���H��H���H��H��0�����e L�������H�J�# H��������(+�# ��)�������(�# ��)�p���H��p���H������H��H��p����@e L���n���H�w�# H��������(X�# ��)�������(8�# ��)�����H��x���H������H��H�ư�����e L������H���# H�������(��# ��)� �����(e�# ��)�����H�CH�����H��H�������@ e L���˒��H�ԝ# H��P�����(��# ��)�@�����(��# ��)�0���H������H��H���H��H��0�����!e L���x����@#e ��$e L���f���H�o�# H��������(P�# ��)�������(0�# ��)�p���H�ChH������I��H��H��p�����$e L������H���# H��������(}�# ��)�������(]�# ��)�����H�CpH������I��H��H�ư����@&e L�������H�ɢ# H�������(��# ��)� �����(��# ��)�����H�CxH�����H��H��������'e L���p���H���# H��P�����(ڣ# ��)�@�����(��# ��)�0���H�C(H��H���H��H��0����@)e L��� ���H�)�# H��������(
�# ��)�������(�# ��)�p���H�C0H������H��H��p�����*e L���А��H�Y�# H��������(:�# ��)�������(�# ��)�����H�C H������H��H�ư����@,e L��耐����-e � /e L���n���H���# H�������(ب# ��)� �����(��# ��)�����H�C@H�����H��H�������@/e L������H�'�# H��P�����(�# ��)�@�����(�# ��)�0���H�C8H��H���H��H��0�����0e L���Ώ��H�W�# H��������(8�# ��)�������(�# ��)�p���H�CPH������H��H��p����@2e L���~���H���# H��������(h�# ��)�������(H�# ��)�����H�CHH������H��H�ư�����3e L���.���H���# H�������(��# ��)� �����(x�# ��)�����H�CXH�����H��H�������@5e L���ގ��H��# H��P�����(Ȯ# ��)�@�����(��# ��)�0���H�C`H��H���H��H��0�����6e L��莎���@8e ��9e L���|���L���   I�H������H�t�# H�E���(X�# ��)E���(;�# ��)�p���H��H����H�E�H��H��p�����9e L������H������I�I�GH������H���# H�E���(w�# ��)E���(Z�# ��)E�H��H����H�E�H��H�ư�@;e L���č��H������I�G��{E ��{^$��,�H�K�Hǅ���   H�       H��P���H������H��X���Hǅ`���    Hǅh���    Hǅp���    Hǅ����    Hǅ����    �P���H��H��P���0�������<e ��=e L��舍��H��H���A\A]A^A_H��]H��[�fD  f.�     SH��H���UUH�kH�l$H��AWAVAUATH��  M��L������I��L������H������H������H��������   L������A�ƾ�   L�������A�L*H��x���1҅ɸ    I�Lc�   I��M��LO�M��LH�M)�LH�Hc�H��HH�H��H���H��H)�H���H������H��Hǅ0����wD Hǅ8���   H��0���L���a�����   L���d������    H�Lc�M9�MN�M)�LH�   I��M��LO�M��LH�L��X���M)�LH�H������M�l L��H������L�������M�Hǅ ���xD Hǅ(���   H�� ���L��L���ˋ����   H�������ʊ�����    H�Lc�M9�MN�M)�LH�1ɸ   I��M��LO�M��LH�M)�LH�L�X���L��H������L���f���M�Hǅ���*xD Hǅ���   H�����L��L���;���M�    L��L�������H��# H�������(�# ��)� �����(ұ# ��)�����H������H�� ���H��x���H�����H�� ���H�����H�����H��H������H��H��0����@>e 車��HǅH���   ��   L������L��蜉��A�ƾ�   H������舉��A�L*H��h����ɸ    I�H��Lc�I��A�   MN�M���    LH�L��H���M)�L��LH�Hc�H��p���H��HH�H��H���I��I)�I���L��`���L��Hǅ����0xD Hǅ����   H������L���߉����   L���������    H�Lc�M9�MN�M)�LH�I��A�   MN�M��LH�L��P���M)�LH�H��H���N�,(L��H������L���w���M�Hǅ����@xD Hǅ����   H������L��L���L�����   H�������K������    I�Lc�M9�MN�M)�    LH�I��A�   MN�M��LH�1�M)�LH�L�P���L��H������L������M�Hǅ����ZxD Hǅ����   H������L��L��蹈��M�    L��L���v���H�       H��P���H��HH���H��X���I��Hǅ`���    Hǅh���    H��`���H������H��h���H������H������H��p���H��p���H��x���Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ���    Hǅ���    Hǅ ���    Hǅ0���    Hǅ@���    HǅP���    Hǅ`���    Hǅp���    Hǅx���    Hǅ����    �P���H��H��P���0��+�����?e ��@e I��I��0���L���φ��� Ae �@Be L��轆��H���# H��P�����(g�# ��)�@�����(G�# ��)�0���H������H��H���H��H��0�����Be L���j���H���# H�E���(��# ��)E���(z�# ��)�p���H������H�E�H��H��p���� De L��� ���H��# H�E���(ʹ# ��)E���(��# ��)E�H������H�E�H��H�ư��Ee L���܅��� Ge � He L���:���HǅH���   H�       H������L������Hǅ����    Hǅ����    Hǅ����    Hǅ����    Hǅ����    �����H��H�ǐ���0�����H��H���A\A]A^A_H��]H��[�fD  f.�     H�   @�z�?H�H�     `}@H�H�     ��@H�H�    � @H�H�   `�!	@I� �f.�     f.�     f.�     UH��AWAVAUATSH���H��   H��$�  I��L��$h  L��M��I��L��$�   ����)�$�   ���$  ��F��)�$p  ���$  H��H  H��$�  HǄ$�  @   H�`     H��$�  H�     H��$�  H��$�  HǄ$�  @   HǄ$�     HǄ$�     HǄ$�     ��(v�# ��)�$   ��(U�# ��)�$�  ��(4�# ��)�$�  H��H�  H��$�  H��H���  H��H��   ��He �����HcH��$�  H�EH����	  HcH��$�  H����	  H��H�   H��1�H��HH�H��H��$`  H��H��H��$X  ��(�$p  �����)�$@  ��(�$�   �����)�$�   H��H��H��$�  ��}���)�$�  ��}���)�$�  H�	H��$�   H�IH��$�   ���� L��H��$�  f.�     f.�     f.�     H�t$8�   1��^��$�  E1��    H��$�  H��$�  H9�H��HN�L)�E1�H���   H���H��H��E1�H���~C���$�  H��$�  L��$�  L��$h  �  f�     f.�     f.�     H��$�   ���$�  L��$�   H��$�  H��$�   H��$�  H�H��$�  H��$�   H�H��$�  H��$�  L�,L��$�  1�H��$�  L��$�  L��$h  �H��$�  ��{�H��āy�L��L��$�  ā{�H��$�  ��q���}���)�$@  ��X���}N� ��Y�ā{���q���)�$0  ��{���q���)�$   苁����(�$   ��u�$0  ��(�$�  ��\���Y���Y���(�$�  ��\�$@  �������Q�L��$�  ġ|)�4@  B�4�  L�$�  H��$�  ��{�H��$�  ��y�ā{���q���}���)�$`  ��X�$�  ��}k� ��Y�L�$�  ā{���q���)�$   ��{�H����q���)�$�  I��蚀��H��$�  ��(�$�  H��$�  ��(�$�  ��u�$   ��(�$�  ��\���Y���Y���\�$`  �������Q�I�ġ|)�4`  I�I��@L��$�  H��H�������L��$�   H��$�   L��$�  H��xI����   �    I��H��$�  H��$�  H���{���y�H��$�   H�H��$�   H���{���q���}���)�$�   ��X���}� ��Y���{���q���)L$p��{���q���)L$@L���]��I����(L$@��uL$p��(�$�  ��\���Y���Y���(�$�  ��\�$�   �������Q�ġ|�4@  M��L��$�   L��$h  H��$�  H��$�  ���$�  H9�I��LN�L��L)�L)�H��}H��$�  ��   fD  f.�     H��$�  L��H��$�  H��H�$�  H���{D� ��yL� ��)�$   H��$�  ��������)�$   ��X�$@  ��Y�� L�����$�  ��w��{�����$�  I����(�$�   ��\�$   ��Y���Y���(�$@  ��\�$   ������Q�ġx)��@  I��H��$�  M)�M9���(�$�  D��$�  ��   L��H��$�  H�$�  ��(�$�   H��$�  ��\����$  ��{D� ���$8  ��X�$p  ��Y�� L�����$�  ��w�Z|�����$�  I��D��$�  ��(�$�  H��$�  ��Y�$  ��Y���(�$p  ��\�$8  ������W���Q�ġ{��@  1�I��H��$�  }H��$�  ��  �     f.�     L��H��H��>L�H���M��I��M��I��1�L��I��I��H��$�   f.�     ����@  �ĸ  ��.�sB��A� ��.�s3E�B� ��A�_H��$�  H��H��$  I�\� H��$   ��(�����H  H���.�s@���� ��.�s1E�B�\ A�_H��$�  H��H��$  I�\� H��$   ��(�����P  H���.�sB���� ��.�s3E�B���A�_H��$�  H��H��$  I�\� H��$   ��(�����X  H���.�s@��G� ��.�s1E�B�\A�_H��$�  H��H��$  I�\� H��$   ��(�H�H��H9������M��� fD  f.�     f.�     H�H��L9�}U����@  ��.�s����� ��.�s�E�B�TA�WH��$�  H��H��$  I�T� H��$   ��(��@ H��$�  H��$�  H��$�  H��$�  H�$X  I��   H��   H��H;�$`  �Y���H�t$8H��H��$�  H9�� ���H��H��   ��Ie ��Je H����w�z����(�$�   ���$(  ��(�$p  ���$0  H��# H��$�  ��(��# ��)�$�  ��(խ# ��)�$�  ��(��# ��)�$p  ��(��# ��)�$`  ��(r�# ��)�$P  ��(Q�# ��)�$@  ��(0�# ��)�$0  ��(�# ��)�$   H��H(  H��$X  H��H0  H��$�  H��H��   �@Ke H���y����Me ��Ne H����x��A���$8  A�G��$<  H�߰# H��$@  ��(��# ��)�$0  ��(��# ��)�$   ��(}�# ��)�$  ��(\�# ��)�$   ��(;�# ��)�$�  ��(�# ��)�$�  ��(��# ��)�$�  ��(د# ��)�$�  H��H8  H��$�  H��H<  H��$8  H��H���  � Oe H���x��H�O�# H��$�  ��(/�# ��)�$�  ��(�# ��)�$�  ��(��# ��)�$�  ��(̱# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  ��(i�# ��)�$p  ��(H�# ��)�$`  H��H  H��$�  H��H   H��$�  H��H��`  �@Qe H���Fw����Se ��Te H���4w���@Ue �@Ve H���"w��H��H���[A\A]A^A_]�fD  f.�     SH��H���UUH�kH�l$H��AWAVAUATH��  H��h���H��`���L��p���L��x���H��8���H�C`HcH��P���1�H��HI�H��   H���I��I)�I���L��H�CXHc��*���=� ��^�H�CPLcM��I���A�B�H�KpL�ChH�{0H��L�[�a  ���Y  ��$���H��@���L��L��@���H�������  H��0���M��I��1�A��L��0���L��p���L��x���L��h���L��`����y  1�H��H������H��L��M��L������H�{ L�[@L�S8L�C�     f.�     H��8�����L�   ��������H�S����   H��@�����L�   ��������   ��|LA�   ����H�S(����   ��|LA�   ������|L A��   ��|LA�   ������|A��   ��L�   ������|A��   ��|L A�    ���� L��������|)A��   H�� H�������L������H������M��I��L������H��@���H9�L��L��@���L��8���H��L�[ H�S��	  �     f.�     ��{D���W���W������H��H�{���H����{D���W�����T� ��{���{D���W������H�K(�����{D���W�����������{D���W������H�K8�����{D���W������H�K@�����D������L��������{�H��H9��;����	  ��*���Z� ��^�E1�A��H�C@I����  �x(�L��0���H��H������H��H������H��H��H��������}���)�����H�H������H�RH������H��0����z�� �{�� L��L������M��L�K8I��L�{(L�c M��H�CH��x���H��p���H��x���L��h���L��`���H��@���H��8��� f.�     f.�     H��H���H��P���L��X���L��`���H��h���H��p�������)�0�����B��)������ ��)�������@��)�������|��)������A|S��|��)������A|b�|/�|w�|>��FH��x�����1���`H��0���H������f�f.�     �G��W���*���rX�ŲY���Z�ţ\���}���Y�0�����}���(�����������yt� L�āyt� ��}�H������L�āyt� H�7��yt� ��Y�������(������������Y���y4�āy4���}�āy4���y4���Y������|)��������Y���y4�āy4���}�āy4���y4���Y������|)��������Y���y4�āy4���}�āy4���y4�ŕY��|)��������Y���y4�āy4���}�āy4���y4�ŅY���(��������Y���y4�āy4���}�āy4���y4���Y���ݨ����ġy���}�ġy����H��H9��d���H������I�I�I�I�I�I�H�H��H���H�� H��P���H�� L��X���I�� L��`���I�� H��h���H�� H��p���H�� H��x��� H�������0���L������L������L�[L��0���L�{@H�{0��(�����H��0����x)�L��L)�H��H��H���I���d  L��0���L��H��H������M�H�C L�M��H�C(L�$H�L�C8I�4J�<9L��H������M��H�KI����x���)� ���H��H����x���)����L��8�����x���)� �����x���)�����L��@�����x���)������Ax4�L��`����Ax,��Ax�L��h����Ax<��Ax�L��x�����x,���x$�L��p�����x<���x4�O�<�����)�����1��     f.�     �Q��W���*���X!� ��Y���Z��x(����� ��\��{�ŹY� �������x(������b���Ay$�J�	�Ay$�ŹY������x)�������Y���y���y���Y��x)�������Y���y���y���Y��x)�������Y���������QY���(�������Y���������AY���(�������Y�������ŹY� �����(���������x)�J�1��yL� I���yL� H��L9������L��0���L������I����(������H������L��0���L��0���M��L�C8L�{ M��L������H��@���L9�L��p���L��L��x���L��h���L��`���H�C0L��@���L�[(�8  M��H������H��H�KH�49M��M�<M��L��M��M�48L�<8M�9H�C@L�I��H��0���L��H�C������������L��������H��8���������������L����������{����������{L����������{L� ��������A{T��{��{l�H��x����{<��{d�H��p�����$���l�I��1�H������H���     f.�     f.�     �z��W���*���X5a� ��Y��JZ���� ��s\��KY������Ax(���������B���{4���Y�������(��x)��������������Y���{���Y������x)��⩩���Y���{�ţY��x)��⑩���Y���{�ŃY��x)��♩���Y���{���Y��x)���ѩ���Y���{���Y�������������bɩ���(�J�<�{�H��H9�����M��L��x���L��h���L��`���L��@���H��@���L��`���L��h���L��x���L��p���L��0���E��H��0���D��A)�L�ChL������H��I����$���H�KpH�{I��L����r  L��H��@���A��D���  �B  L��@���E1Ƀ�L�[(I��L��8�����  H��@���ġcYD��H��`���ġcYT��H��h���ġcYt��H��x���ġcYd��H��p���ġcYl����}���}���}���}���}�H��H��@���H�KĢ}l��Ă}t��H��I��H��0���I���L��H��I����   H�K@I��H�S8I��H�C0H�K H�S�    ��42��1��|3��0��|<7��|$2��|,4��t2 ��L1 ��|D3 ��T0 ��||7 ��|d2 ��|l4 ��t2@��L1@��|D3@��T0@��||7@��|d2@��|l4@��t2`��L1`��|D3`��T0`��||7`��|d2`��|l4`H��H��H����=���L��8���H��yiH�K@I��H�S8H�C0H�K L�{fD  f.�     f.�     ��|47��1��|3��0��<2��|$2��|,4H�� H��x�L��8���H��0���Ic�H��H)�H����  L��H�{0L�k �s  A��Ic�M��H��H��8���H�t��0�L��I����w��f��H��@���ġ{D��L�{XIc��*ʉ�A����^����X���H�H�K H�<�H��H��X���0��f��H��`���ġ{D��Ic��W���*ʉ�A����^����`���H�H�K(H�<�H��H��`���0��9f��H��h���ġ{D��Ic��W���*ʉ�A����^����h���H�H�K0H�<�H��H��h���0���e��H��x���ġ{D��Ic��W���*ʉ�A����^����p���H�H�K8H�<�H��H��p���0��e��H��p���ġ{D��Ic��W���*ʉ�A����^����x���H�H�K@H�<�H��H��x���0��ge��IcD��Ic�I�<�H�CJ�t��0��Ie��L�ChL���S  H��@���H��Ic�H�ā{D��H�C���H��@���ġcYD�����M��L�k ��xD� H��`���ġcYD�������x�H��h���ġcYD�����H�{0���H��x���ġcYD�����H�s8���H��p���ġcYD�����H�C@���H�sġ{D��L����x�A��D9�H����   E�H��@�����YD��I�D��Ic�H�sH����{D� H��`�����YD����{�H��h�����YD�����H��x�����YD��H��p�����YL��H�C8���H�C@���H�CH�D��I��H�CXD�0D��H�KpA� ��9���  �kxD ��   �   L�kxL��A����w�d��E��Ic�H������I��I��L��H���H��H���M��I���I)�H�KHH��H�����H��P���)�H��P�����*����� ��^����(���H��H���I��I���L)�H�� ���H��H�����H��?H�������������H��Hc�H�������  �qxD ��   �   L��E��L��������w�c��L������E��H����  1�A���   ��  L��H������������������H�sH�� ����     ��XT ��XL@��XD`��$��  ��  ��X�H��H��H9�|�H���{�(���y*�    f.�     f.�     ��XH�� H��x���X���X���X���}���X���y���X�Ic�L9���  Ic�L���f.�     f.�     f.�     ��X�H��H9�|���  �^xD ��   �   L�{xL����w�b��H���z  �dxD ��   �   L����w�]b��H���W  Ic�H�K��D����{L��L��������Y�� ���� �������X�� L�kH��{D����  1�A���   �=  L��H������������������L�SH�� ���f�f.�     ��mXT ��uXL@��}XD`��|$A��  A��  ��X�H��H��H9�|�H��L�kHL�����(���y@ f.�     ��eXH�� H��x���X���X���X���}���X���y���X�Ic�L9���  Ic�L��f.�     f.�     f.�     ��{X�H��H9�|��[  L������Ic�H�K��D����Xy� ��y� ����x� L�kH��{L���  ��W�A��L�S��   L��H��1�����I��|w1�L��H��������f.�     f.�     f.�     ��}X
��|L
 ��uXL
@��X���|�
�   ��uX�
�   ��uXL
`��X�H���   H��H���H���H��y fD  f.�     ��}X
H�� H��x���}���X���y���X���W���X�L��Ic�H9�L�kHL�����(���~(Ic� f.�     f.�     ��{X�H��H9�|���Y���X.� ��.� ���-� A��  ��   1�A��L������������   L��H����}�1�I��|A1�fD  f.�     f.�     ����L ��L@��L`H��H��H���|�H��y%f�f.�     f.�     ��H�� H��x�D��L������Ic�Hc�H)�H��|Hc������σ�H��P���L�����A9�~bH�H�SH����UL�������������Ic�H��H�ƈ���0�M��E����w�]��E�����(���M��H��P���L�����L�������������  ���  ��W���ѩ�� ���� ����� 1�@ ��{�H��L9�|���  ��W�A��H�s��   L��H��1�����I��|n1�L��H��������     f.�     f.�     ��X��L ��XL@��X�����   ��X��   ��XL`��X�H���   H��H���H���H��y%f�f.�     f.�     ��XH�� H��x���}���X���y���X���W���X�L��Ic�H9��{�(���~/Icɐf.�     f.�     f.�     ��X�H��H9�|�H�� ���ŻY�1�A���   ��   L��H������������������ ��eX\ ��mXT@��uXL`��|(,A��  A��  ��X�H��H��H9�|�H��y-f�     f.�     f.�     ��]X$H�� H��x���X���X���X���}���X���y���X�Ic�L9��  Ic�L��f.�     f.�     f.�     ��sX�H��H9�|���   ��W�A����   L��H��1�����I��|X1�L��H���������uX��|(T ��|(��   ��mXT@��X���eX��   ��mXT`��X�H���   H��H���H���H��y fD  f.�     ��uXH�� H��x���}���X���y���X���W���X�L��Ic�H9�~Ic�fD  ��sX�H��H9�|�I��ŻY���Ys� ��s� �������Xn� A��  }FL��1�A��M��L�����D�������  H��H��H����}�1�H��H���"  H�SH�D  L�������������Ic�H��H�ƀ���0�H�{HE����w�Y��E���{�(���M��H��P���L�����L�����D������d  H��������tnA����  1����� @ ��WɅ�u%I��1���W�@ f.�     ��X�H��L9�|���ѩr� ��(����u� ��{�H��L9�|��  1�L��H���I����F� f�f.�     f.�     M��1�����L���������������     I�<��X$��ǀ  ��X\� ��XT�@���  ��XL�`H��H��L9�|�H��H��H���y/�     f.�     f.�     I�<��X$�H�� H��x�M����X���X���X���}���X���y���X�H9�����|.I��H��fD  f.�     f.�     ��X�H��L9�|���ѩ2� ��(����5� ��{�H��L9�������?  1�H�SH��
��L
 ��L
@��L
`H��H��H���|�H��y%f�f.�     f.�     ��
H�� H��x��Ic�Hc�H)�H��}	L�������L������Hc����H�SH��ʃ�H��P���L�����A9�~H�H�KH�������  D��<���L��H���L��H���H������1�L��H��H��X������ ��� ��� E1�H������f.�     ��W�E����  H������1҅��  ��Wۃ�<�����   ����L��H��X���I��|h����H������L��H��X����     ��eX;��|d; ��]Xd;@��X���|�;�   ��]X�;�   ��]Xd;`��X�H���   H��H���H���H��y fD  f.�     ��eX;H�� H��x���}���X���y���X���W���X�L��H9������  I����X�H��L9�|���   @ f.�     ����H��X���������������f�     I���X4��Ȁ  ��Xl� ��Xd�@���  ��X\�`H��H��L9�|�H��y I���X4�H�� H��x���X���X���X���}���X���y���X�L9�����|DI��L��f�     f.�     f.�     ��X�H��L9�|�f�     f.�     ŻY���W�E����  H������1҅���   ��W䃽<�����   ����L��H��X���I��|^����H������L��H��X���@ ��X$>��l> ��Xl>@��X����>�   ��X�>�   ��Xl>`��X�H���   H��H���H���H��y%f�f.�     f.�     ��X$>H�� H��x���}���X���y���X���W���X�L��H9�������   D  H�
��X$�H��L9�|���   f.�     ����H��X���������������f�     H���X<��Ȁ  ��Xt� ��Xl�@���  ��Xd�`H��H��L9�|�H��y H���X<�H�� H��x���X���X���X���}���X���y���X�L9�����|$L�ȐH���X$�H��L9�|�D  f.�     ŻY���Y�������X���{�I��H��L9�������-  H��H���I��I��1�I��I���I����� 1��    f.�     M������H��L��I��|X����L��H��L����uX:��|T: ��mXT:@��X���|�:�   ��mX�:�   ��mXT:`��X�H���   H��H���H���H��y fD  f.�     ��uX:H�� H��x�M����}���X���y���X���W���X�L9�����|I��L����X�H��L9�|���ѩ� ��(������ ��{�H��H��L9�����H�@�# H��������(!�# ��)�������(�# ��)�����H�CpH������H��H�Ɛ���I��I��������Ve L����w�Q��H�#�# H��������(�# ��)�������(�# ��)�����H�CXH������H��H������� Xe L���:Q��H�CH�����Hǅ���@   H�`     H�� ���I��H�     H��(���I��H��0���Hǅ8���@   Hǅ@���   HǅH���
   HǅP���   ��(��# ��)�������(��# ��)�������(a�# ��)�p���H��H���H������H��H��p����@Ye L���rP��H������H������Hǅ����@   L������L������H������Hǅ����@   Hǅ����   Hǅ����
   Hǅ����   ��(0�# ��)�0�����(�# ��)� �����(�# ��)����H��H����H�� ���H��H�������Ze L����O��H�CHH��P���HǅX���@   L��`���L��h���H��p���Hǅx���@   H�E�   H�E�
   H�E�   ��(˓# ��)E���(��# ��)E���(��# ��)E�H��HP���H�E�H��H�ư��[e L���+O��H��H���A\A]A^A_H��]H��[�f�     f.�     UH��AWAVAUATSH���H��   I��I��L�L$L��H�t$H�EHc H�M(�	)�H�U0HcH��Lc�I�J�D��H�D$ H�D�# H�D$`��('�# ��)D$P��(	�# ��)D$@H��H�� H�D$XH��H��@I��I�ǀ   � ee L���YN��K�D��H�D$(H�(�# H��$�  ��(�# ��)�$�  ��(�# ��)�$�  H��H��(H��$�  H��H�ƀ  �@fe L����M��K�D��H�D$0H��# H��$�  ��(�# ��)�$�  ��(Ɲ# ��)�$�  H��H��0H��$�  H��H���  ��ge L���M��H�D$J�D��L�d$I�$L�uH�D$0I���p� ��]D$ ����he ��ie L���TM���@je �@ke L���BM��H�ۡ# H��$`  ��(��# ��)�$P  ��(��# ��)�$@  ��(y�# ��)�$0  ��(X�# ��)�$   ��(7�# ��)�$  ��(�# ��)�$   H��$8  H��H��   ��ke L���L��H�G�# H��$�  ��('�# ��)�$�  ��(�# ��)�$�  ��(�# ��)�$�  ��(Ģ# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  L��$�  H��H�ƀ  ��me L���L��H���# H��$`  ��(��# ��)�$P  ��(r�# ��)�$@  ��(Q�# ��)�$0  ��(0�# ��)�$   ��(�# ��)�$  ��(�# ��)�$   L��$8  H��H��   ��oe L���K��H�D$(H�D$8H��# H��$�  ��(��# ��)�$�  ��(ԥ# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  ��(q�# ��)�$�  ��(P�# ��)�$�  H��H��8H��$�  H��H�ƀ  ��qe L����J����se ��te L����J���@ue �@ve L���J��H��H���[A\A]A^A_]Ðf.�     UAWAVAUATSH��(����Y�� Lc&H��$�   ġ{L��L��$�   ā{D����(���{ H��$�   ġ{D��L��$�   K�D��    L��$�   K�D��    H��$�   H�t$xL�T$pH�l$`�   D9��n  ��D$��L$��E ��{Y���� ��^���Y� ��� ������\$L�$H���I����L$N�L� O�l� O��O��N�4�H��$�   N��D��I���sY�� H�$����Y�� 1Ƀ�}��d$��|$H�t$xL���  L�l$ M��I��?M�I���1��{b� ��d$��|$@ f.�     H��H��L��H)�H�|$ H)�L��H)�L��H)�L��H)�M��I)�Ż����aK�P��\���Y���Y���Y���X���]���\���Y���Y�Ż����aK�p��Y���X���]���n���~���O���O�H�@�    H�@�    H�B�    H�B�    H�C�    H�C�    ��X���{o���X���{g�H��L9��(���H�t$xL��L�l$ L9���   I���{O� �    Ż����aK�P��\���Y���Y���Y���X���]�H��H��H��H)���z�L��H)���J�L��H)�L��H)�L��H)���X�L��H)�H�B�    H�F�    H�G�    ��c�H��L9��{���L��H�    H��([A\A]A^A_]��wÐf.�     f.�     f.�     UAWAVAUATSH��  ����s� ���r� ��{ ��Ym� Lc.H��$�  H��$�  L��$�  H��$�  L��$�  H��$�  M���_  �� ��Ť ��W����M��I�����YΓ ����b}��b}������}����$�  ����b}�N�4�����L��H���H��H�������k  ��)�$0  ���$H  ��)�$P  A����  ��)�$�  ��)�$�  ��)�$   ��)�$  L��$(  E1�A����  L��$�  H��$h  L��$�  ��)�$p  I��H��$�  I��H��$(  L��$�  �|�$�  �|�$`  �|�$@  f�f.�     ġ|L6����$@  BL6�ġuYD3���(�BL3�ġ|D7�BL7�ā|D7�CL7���}΢ ��Y���}ɢ ��X���^��������$`  ��}�� ��Y���}�� ��X���}�� ��Y����$   ����I��H���E������$   ��}}� ��Y����$   ��(��WE����}�� ��(���X���^���}Q� ��(����$   ��Y���\���Y���Y�$`  ��}.� ��Y�����$@  ��(���}-� ��}5� ��ͨ���}-� �������Y���Y���Y���\�$`  ��}%�� ��}-�� ��ը���Y���X���Y�$@  ��^����$�  ��^���Q�ă�L5���Y�CL5����� ���$�  �D�����$�  ���$�  �ED�����$�  H��L��H��$�  L��$�  ��Y�$�  ��}l� �������}ޠ ��]�����ā|D7���}̠ ġ|D1�BL1�I���I�������L��$�  ��(�$p  L��$�  H��$h  L��$�  H��$h  I��H��M9���  H���  ���L  ��)�$0  ���$H  ��)�$P  A���:  ��)�$�  ��)�$�  ��)�$   ��)�$  H��$h  1�A����  L��$�  L��$�  ��)�$p  H��H��$�  H��L��I��M���|�$�  �|�$`  �|�$@  �     ġ|L6����$   BL6�āuYD5���(�CL5�ā|D7�CL7�ġ|D3�BL3���}n� ��Y���}i� ��X���^��������$   ��}Q� ��Y���}L� ��X���}G� ��Y����$�  �����'B������$�  ��}#� ��Y����$�  ��(���A����}L� ��(���X���^���}�� ��(����$�  ��Y���\���Y���Y�$   ��}Ԟ ��Y�����$   ��(���}-�� ��}5�� ��ͨ���}-�� �������Y���Y���Y���\�$`  ��}%�� ��}-�� ��ը���Y���X���Y�$@  ��^����$�  ��^���Q�ă�L4���Y�CL4������ ���$�  ��@�����$�  ���$�  ��@�����$�  H��$�  ��Y�$�  ��}� �������}�� ��]�����ġ|D3���}�� H��$�  ġ|D1�BL1�I���H�������H��$�  ��(�$p  L��I��M��L��$�  L��$�  I9��D  L��I��H��$�  H���'  A���		  1���]�� A����   H��H����}�H����   ��(���}Ԝ @ ā|\1�C�1`���ġeY\6�B�6`���ġ|\7�B�7`���ā|L7�C�7`���ġ|T1�B�1`���ġ|\6�āeY\1�ġ|\7�ā|L7�ġ|T1�I���H��H����{�����(���(�H��y,ġ|T6�āmYT1�ġ|T7�ā|L7���}� ġ|L1�H����(�I9�|6L��H)���L����qYL����L�������xL����(z� ��L��H��H9���  I)�ġ{L��āsYL��ġ{L���  L��L)���\����)�$�  H��$(  ��YD����xD����xD����Y!� ��)�$p  ��X � ��^���y���)�$�  ��Y� ��y���X� ��Y� ��)�$p  ��w��;����)�$�  ��y�$p  ��;����(�� ��X���^���\���(�$�  ��Y� ��Y���(�� ��Y���Y�$�  ��Y� ��y�$�  ��(-� ��٨-� ��٨-�� ��Y���Y���Y���X�$�  ��(� ����� ��Y���X���Y�$�  ��^���(�$   ��^���Q�H��$�  ��yL����Y���W ��)�$@  �:����)�$`  ��(�$@  ��:����(�$p  ��Y�$`  ��(�$  ����� ��]�� ��y���xD����(l� ��D��H��I��H��$h  L9�L��H��$�  L��$(  ��  M)�ġ{L�����$�  āsYD��ġ{D��ā{D�����$�  ��Yz� ��Xz� ��Yz� �����)�$p  I��I��H��L����w�:����y����� ��X����$�  ��X� �����Y
� �����^���y���Y� ��\���Y���%� ��Y���Y���Y� ��� ���$�  ��ɩ�� ��(����-�� ��Y���Y���Y���\�$0  ��� ���ߘ ��Y���X���Y�$P  ��^����$H  ��^���Y� ��Q�H��$�  ġ{YD����ʘ ��W����$8  �:�����$X  ���$8  �:��I��H��L��L����(�$p  ��Y�$X  ��c� ��  L��H)���\����)�$   ��aYD��L����D��I����xD����Y_� ��)�$p  ��X^� ��^���y���)�$  ��YS� ��y���XU� ��Y]� ��)�$�   H��$0  ��w�48����)�$�   ��y�$�   �8����(3� ��X���^���\���(�$�   ��Y&� ��Y���(*� ��Y���Y�$  ��Y%� ��y�$   ��(-"� ��٨-)� ��٨-0� ��Y���Y���Y���X�$�  ��(#� ���*� ��Y���X���Y�$�  ��^���(�$   ��^���Q���yL����Y���W � ��)�$�   ��6����)�$�   ��(�$�   �+7��H��$0  ��(�$p  ��Y�$�   ��(�$  ����� ��]�� ��y���xD����(�� H��$�  ��D��H��H��H��$h  H9�H��$�  L��I����  I)�ġ{L�����$(  āsYD��ġ{D��ā{D�����$�   ��Y�� ��X�� ��Y�� �����)�$p  H��I��L����w�G6����y���є ��X����$�   ��XT� �����Y@� �����^���y���YR� ��\���Y���%J� ��Y���Y���YB� ��J� ���$(  ��ɩ0� ��(����-3� ��Y���Y���Y���\�$0  ��&� ���� ��Y���X���Y�$P  ��^����$H  ��^���Y;� ��Q�ā{YD����� ��W����$�   �F6�����$�   ���$�   ��6��I��L��H��H��$�  ��(�$p  ��Y�$�   ���� �������]� ā{D��H�     8��J�D������{Y����{H�     8��H�H�Ĉ  [A\A]A^A_]��w�f�     f.�     f.�     UH��AWAVAUATSH���H���  H�T$@L��H�t$HHcL���   H���   L���   L���   L���   L���   H��L���   H���   H����  L�L$ =  �-  H�T$E1҃�H���   I��H���   I����   H���H��H��1�H��|M1��f.�     ������1��3��D1 ��D3 ��D1@��D3@��D1`��D3`H��H��H���|�H��y3�    f.�     f.�     ������1��3H�� H��x�A��Ic�H�t$H)�H��H���   I��L�L$ L���   H���   |Ic���W����������A��H���   H�D$D9�~qIc�H��    H��    �\���e 0�H��I���4�����e 0�H��L���4��H���   H���   I��H���   I��H���   I��L�L$ L���   H���   I�    H�EH�     H�E H�     H�E(H�     H�E0H�     H�E8H�     H�EH�     H�E@H�     H�EHH�     H�EPH�     H�EXH�     H�E`H�     H�EhH�     H�EpH�     H�ExH�     H���   H�     H���   H�     H�    I�$    E� H�D$H�D�����A��D$0H���   D�8D�����A��D$,H���   ��ؙ��A�E A��H�D$@i�3��șA��A��D$<�șA��H���   ��D$8�ș����D$4I��I�ƀ   ��ve ��we L����w��2���@xe �@ye L����2����ye ��ze L����2��D�d$PD�|$T�\$XH��H��PH��$�  HǄ$�      H�`     H��$�  I��H�      H��$�  H��H��$�  HǄ$�      HǄ$�     HǄ$�     HǄ$�     ��(@�# ��)�$   ��(�# ��)�$  ��(��# ��)�$   ��(ݕ# ��)�$�  ��(��# ��)�$�  H��H�  H��$  H��H���  �@{e L����1���D$0�D$`�D$,�D$dD�l$hH��H��`H��$@  HǄ$H      L��$P  H��$X  H��$`  HǄ$h      HǄ$p     HǄ$x     HǄ$�     ��(�# ��)�$�  ��(�# ��)�$�  ��(ɖ# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  H��H@  H��$�  H��H�Ơ  � }e L����0���D$<�D$p�D$8�D$t�D$4�D$xH��H��pH��$   HǄ$      L��$  H��$  H��$   HǄ$(      HǄ$0     HǄ$8     HǄ$@     ��(ӗ# ��)�$�  ��(��# ��)�$�  ��(��# ��)�$�  ��(p�# ��)�$p  ��(O�# ��)�$`  H��H   H��$�  H��H��`  ��~e L����/�����e ���e L����/��H��H���[A\A]A^A_]�fD  f.�     f.�     H�    I�     I�    H�D$H�     H�D$H�     H�    ����    SH��H���UUH�kH�l$H��AWAVAUATH��0'  L������Hc7H��H���H��  L��  L��  H��L���   H���   ���  H�������?��*����������� ��W������)� �����n���y�������)�P���Ic H��H���H������Hc	H��(���H����v  I���<  �����L�d�L��H��H������L�t�L��H��H�H��H��H�� ���H��H��=L�D�L��H��H�H��H��H��H��=L�T�L��H��H�H��H��H�|�H��H��H��=H�H��H��H��H��=H�H��H��H�H�|�H��H��H��=H�H��H��H��H��=H�H��H��H�H�L�H��H��H��=H�H��H��H��H��h���H��H��=H�H�C(���������H���H������H���H������H���H������I���L������I���L������H�Cx�� ���H���A���l���A��������(���w�X+��H�C� H��������*	H�K �{9�x)�����Ł���(;� �i^��x)�`�����?� ����>� �����(= ��_���)�������%&� ����%%� ��%� ����$� ��$� ����#� H�K8�{1�x)� ���H�K@�{!�x)� ��������_���)�0��������_���)� ��������_���)�������ӌ �{_��{_��b}���_��{_���}��A}_��|)�0�����}���_���)������}���_���)������}���_���)��������������cY���Y� ��*���)������[Y��{�����I���L�����I���L�� �����f� ŀW���Y�� ��sY����������Y}� ��{Y����������Y����p�����W���*���Yd� ��sY����h����bZ��x)�0��������)�@�����ׇ ��k\���)������-�w ��Y���)�������X������)�@�����}���)�0�����}���)�������}���)�����ŋY���)�������Y���}���)�������}���)�p������H�����}���)����������)�����H�K����)�������}���)�0���1������)�������}���)�p�����}���)������}���)�������}���)�P�����n���y�������)�0�����Y����������y���)�������{���)� ��������)�0��������)� ��������)���������)� ���śYG� �������ŃX?� ������������)�������{���)������{���)�������}���)�������)�������)�������)�P�����)� �����)�������)�`�����)�@�����)�������)�p�����)�0�����)�������)�������)�������)�P�����)�������)�p�����)�0�����)�������)�p�����)�0�����)�������)�������)�0�����)������)�������)�p�����)�������)�������)�������)�p�����)������)�������)�������)�P�����)�������)�0���L������L���   f.�     f.�     f.�     H��x����A1�H������������A8W�H��   M��L���   L���   M���    Mc1M����   ���H����{�h����`zD ��   �   I����w��'��L��H���H������H����   �fzD ��   �   L����w�'��H��L���   �)  �zzD ��   �   L����w�'��H����  ��zD ��   �   L����w�\'��H����  A��L���   L������H��������  E1���(������\  f.�     f.�     ���H����{�h���L���   �+  @ ��zD ��   �   L����w��&��H��L���   �  ��zD ��   �   L����w�&��H���s  A���I  E1�L���   ��(�����L�������  �f.�     �mzD ��   �   L����w�V&��H���]  �szD ��   �   L����w�3&��H����  A���0  E1�L���   ��(�����L��������  �     f.�     A����  E1�L���   ��(�����L��������	  @ f.�     f.�     A���  E1�L���   ��(�p���L��������  @ f.�     f.�     A��L���   L������H��������  E1���(������  �    f.�     A����  E1�L���   ��(�0���L��������  @ f.�     f.�     A���  E1�L���   ��(����L��������  @ f.�     f.�     E1�I��I����(�����D  f.�     ā|L5 ��)�P���CL5x��}J� ����O�$3��P���tp��u-$��Y�p�����)�������}�� ��^���(�������sY��� ��^���}���\�M���4$����(�����M����)�������Y�������Y�0�����)�������v�������P�����   ��)�������(�p�����}Y$��}n� ��^�AL$x��(�������sYM����T� ��^���}���\��#����(�������(�������uK� ��)�������Y������Y������(�������eK� �f.�     f.�     M����(�����L������L���   ��(�������(�P�����\���Y���(�0����������]�ā|D5 I�� I��M���F���L���   L������Mc1L��L)�H����  ��)�����H���   ġx���)��������O� ��P���tyĂq-���Y�������)�`�����(9� ��^���(�������sY�� ��^������\�M����w� ����(�`���M����)�0�����Y�������Y�������)�p�����v���W���P�����  ��)�P�����(�����āyY���(�� ��^���(�������sYM�����~ ��^������\���w�����(�0�����(�P�����qK� ��)�0�����Y�������Y�������(�p�����qK� �  �    f.�     f.�     E1�L������I��L���   ��(�����@ ā|L5 ��)�����CL5x��}�} ����O�$3��P���tp��u-$��Y�p�����)�������}�} ��^���(�������sY���} ��^���}���\�M���� ����(�����M����)�������Y�0�����Y�0�����)�������v�������P�����   ��)�P�����(�p�����}Y$��}} ��^�AL$x��(�������sYM�����| ��^���}���\��B ����(�������(�P�����uK� ��)�������Y������Y������(�������eK� ��(������f�f.�     M����(�������(�����L������L���   ��(�������\���Y���(�0����������]�ā|D5 I�� I��M���F���L���   L������Mc1L��L)�H���d  ��)�����H���   ġx���)�@�������~ ��P���tyĂq-���Y�������)�0�����(�~ ��^���(�������sY���{ ��^������\�M����w�����(�0���M����)� �����Y�0�����Y�������)�P�����v���W���P����J  ��)������(�����āyY���(L~ ��^���(�������sYM���� { ��^������\���w�,����(� �����(������qK� ��)� �����Y�������Y�������(�P�����qK� ��  �    f.�     f.�     E1�L������I��L���   ��(�����@ ā|L5 ��)����CL5x��}�z ����O�$3��P���tp��u-$��Y�p�����)�������}:z ��^���(�������sY��)z ��^���}���\�M���t����(�����M����)�p�����Y�������Y�0�����)�������v�������P�����   ��)�p�����(�p�����}Y$��}�y ��^�AL$x��(�������sYM�����y ��^���}���\�������(�p�����(�p�����uK� ��)�p�����Y������Y������(�������eK� ��(������f�f.�     M����(�������(�����L������L���   ��(������\���Y���(�0����������]�ā|D5 I�� I��M���F���L���   L������Mc1L��L)�H����  ��)�����H���   ġx���)�P�������{ ��P���tyĂq-���Y�������)�@�����(y{ ��^���(�������sY��Px ��^������\�M����w�Y����(�@���M����)�������Y�������Y�������)�@�����v���W���P�����  ��)� �����(�����āyY���(�z ��^���(�������sYM�����w ��^������\���w������(�������(� �����qK� ��)�������Y� �����Y�������(�@�����qK� �F  �    f.�     f.�     E1�L������I����(�p����f.�     ā|L5 ��)�����CL5x��}*w ����H���   N�$1��P���to��u-$��Y�p�����)�������}�v ��^���(�������Y	���v ��^���}���\�I�������(�����L����)�0�����Y�0�����Y�0�����)�p�����v�������P�����   ��)�������(�p�����}Y$��}Hv ��^�AL$x��(�������Y	��2v ��^���}���\������(�0�����(�������uK� ��)�0�����Y������Y������(�p�����eK� ��(������f.�     ��(�������(�p���L���   ��(�������\���Y���(�0����������]�ā|D5 I�� I���S���L���   L���   L������L������Mc1L��L)�H����  ��)�p���H���   ġx���)�������.x ��P���tyĂq-���Y�������)�������(x ��^���(�������sY���t ��^������\�M����w������(�����M����)�������Y�0�����Y�������)�������v���W���P�����  ��)�������(�����āyY���(�w ��^���(�������sYM����_t ��^������\���w�k����(�������(�������qK� ��)�������Y� �����Y�������(�������qK� �E  fD  f.�     f.�     A����  1�L���   ��(�����L�������  D  f.�     f.�     E1�L������I��L���   ��(�����@ ā|L5 ��)�0���CL5x��}�s ����O�$3��P���tp��u-$��Y�p�����)�p�����}:s ��^���(�������sY��)s ��^���}���\�M���t����(�p���M����)�P�����Y�0�����Y�0�����)�������v�������P�����   ��)�0�����(�p�����}Y$��}�r ��^�AL$x��(�������sYM�����r ��^���}���\�������(�P�����(�0�����uK� ��)�P�����Y������Y������(�������eK� ��(������f�f.�     M����(�������(�����L������L���   ��(�0�����\���Y���(�0����������]�ā|D5 I�� I��M���F���L���   L������Mc1L��L)�H���D  ��)�����H���   ġx���)���������t ��P���tyĂq-���Y�������)� �����(yt ��^���(�������sY��Pq ��^������\�M����w�Y����(� ���M����)�p�����Y�0�����Y�������)�������v���W���P����*  ��)� �����(�����āyY���(�s ��^���(�������sYM�����p ��^������\���w������(�p�����(� �����qK� ��)�p�����Y�������Y�������(�������qK� �  �    f.�     f.�     E1�I��I����(�����D  f.�     ā|L5 ��)����CL5x��}*p ����O�$3��P���tp��u-$��Y�p�����)�P�����}�o ��^���(�������sY���o ��^���}���\�M�������(�P���M����)�p�����Y�������Y�0�����)�������v�������P�����   ��)������(�p�����}Y$��}No ��^�AL$x��(�������sYM����4o ��^���}���\������(�p�����(������uK� ��)�p�����Y������Y������(�������eK� ��(������f�f.�     M����(�������(�����L������L���   ��(������\���Y���(�0����������]�ā|D5 I�� I��M���F���L���   L������Mc1L��L)�H����  ��)�����H���   ġx���)��������/q ��P���tyĂq-���Y�������)������(q ��^���(�������sY���m ��^������\�M����w������(����M����)�������Y�������Y�������)�������v���W���P�����  ��)�������(�����āyY���(�p ��^���(�������sYM����`m ��^������\���w�l����(�������(�������qK� ��)�������Y� �����Y�������(�������qK� �F  �    f.�     f.�     E1�L������I����(�0����f.�     ā|L5 ��)�����CL5x��}�l ����H���   N�$1��P���to��u-$��Y�p�����)�������}sl ��^���(�������Y	��cl ��^���}���\�I�������(�����L����)������Y�0�����Y�0�����)�0�����v�������P�����   ��)�������(�p�����}Y$��}�k ��^�AL$x��(�������Y	���k ��^���}���\�� ����(������(�������uK� ��)������Y������Y������(�0�����eK� ��(������f.�     ��(�������(�0���L���   ��(�������\���Y���(�0����������]�ā|D5 I�� I���S���L���   L���   L������L������Mc1L��L)�H����  ��)�0���H���   ġx���)���������m ��P���tyĂq-���Y�������)�������(�m ��^���(�������sY���j ��^������\�M����w�����(�����M����)�`�����Y�0�����Y�������)�������v���W���P�����  ��)�������(�����āyY���(+m ��^���(�������sYM�����i ��^������\���w�����(�`�����(�������qK� ��)�`�����Y� �����Y�������(�������qK� �%  fD  f.�     f.�     E1�L������I����(�����f.�     ā|L5 ��)�����CL5x��}ji ����H���   N�$1��P���to��u-$��Y�p�����)�������}i ��^���(�������Y	��i ��^���}���\�I���N����(�����L����)�������Y�0�����Y�0�����)������v�������P�����   ��)�������(�p�����}Y$��}�h ��^�AL$x��(�������Y	��rh ��^���}���\�������(�������(�������uK� ��)�������Y������Y������(������eK� ��(������f.�     ��(�������(����L���   ��(�������\���Y���(�0����������]�ā|D5 I�� I���S���L���   L���   L������L������Mc1L��L)�H���c  ��)����H���   ġx���)��������nj ��P���tyĂq-���Y�������)�������(Xj ��^���(�������sY��/g ��^������\�M����w�8����(�����M����)�P�����Y�0�����Y�������)�������v���W���P����I  ��)�������(�����āyY���(�i ��^���(�������sYM�����f ��^������\���w�����(�P�����(�������qK� ��)�P�����Y� �����Y�������(�������qK� ��
  fD  f.�     f.�     M����(�p�����)�p�����(������(�������\���Y���(� ���������]�H���   ġx�I��M����(�����L������M9���   ��)�������(�����ākY���Zi ��^���kYM�����e ��^���X�L���   ā{����x�����w�A�����x�����.ye ���h�����
  ���p����
  f�f.�     f.�     ��)������s  f�     f.�     M����(�P�����(�����L������H���   ��(�@�����)�P�����(������\���Y���(� ���������]�ġx�I��M����)�����M9���  ��(�����ākY���:h ��^���kYM����vd ��^���X�L���   ā{L� ���h�����w� �����h����;  f�     f.�     M����(�@�����(�����L������H���   ��(�P�����)�@�����(������\���Y���(� ���������]�ġx�I��M����)�����M9��  ��(�����ākY���Zg ��^���kYM�����c ��^���X�L���   ā{L� ���p�����w�@�����p�����.xc ���h���L���   �{  ��������n  @ f.�     E1�L������I����(������f.�     ā|L5 ��)�p���CL5x��}
c ����H���   N�$1��P���to��u-$��Y�p�����)�P�����}�b ��^���(�������Y	���b ��^���}���\�I��������(�P���L����)�0�����Y�0�����Y�0�����)�������v�������P�����   ��)�0�����(�p�����}Y$��}(b ��^�AL$x��(�������Y	��b ��^���}���\��`����(�0�����(�0�����uK� ��)�0�����Y������Y������(�������eK� ��(������f.�     ��(�������(�����L���   ��(�p�����\���Y���(�0����������]�ā|D5 I�� I���S���L���   L���   L������H������Mc1L��H)�H���i  ��)�����H���   �����)��������d ��P���t��q-���Y�������)�@�����(�c ��^���(�������sY���` ��^������\�I��M����w������(�@���M��L����)�������Y�0�����Y�������)�0�����v���W���P����D  ��)�0�����(�������yY���(fc ��^���(�������sYI��M����7` ��^������\���w�C����(�������(�0�����qK� ��)�������Y� �����Y�������(�0�����qK� ��  �     f.�     M����(�������(�p���L������H���   ��(������)�������(������\���Y���(� ���������]�ġx�I��M��M9�~|��)�p�����(�����ākY����b ��^���kYM����:_ ��^���X�L���   ā{����(�����w�������(����?  @ f.�     f.�     ��)�p����3  f�     f.�     M����(�������(�����L������H���   ��(�������)�������(������\���Y���(� ���������]�ġx�I��M����)�����M9���  ��(�����ākY����a ��^���kYM����6^ ��^���X�L���   ā{L� ��������w�� ����������.^ �������L���   �  ���p����  @ f.�     M����(�������(�����L������H���   ��(�������)�������(������\���Y���(� ���������]�ġx�I��M����)�����M9���  ��(�����ākY����` ��^���kYM����6] ��^���X�L���   ā{L� ��� �����w�������� �����.] ���h���L���   v������� f.�     f.�     ��Y��x(�����Ż\���Y�0�������Ż]�ā{D� L���   �  f�f.�     M����(�������(�0���L������H���   ��(�������)�������(������\���Y���(� ���������]�ġx�I��M��M9�~|��)�0�����(�����ākY����_ ��^���kYM�����[ ��^���X�L���   ā{����������w�������������   @ f.�     f.�     ��)�0�����  f�     f.�     M����(�������(����L������H���   ��(�������)�������(������\���Y���(� ���������]�ġx�I��M��M9���   ��)������(�����ākY����^ ��^���kYM�����Z ��^���X�L���   ā{����������w�������������.�Z �������v���������Y��x(�����Ż\���Y�0�������Ż]�ā{�L���   M��L���   L���   �x(�`�����(�0���L��������  @ ��)�����  f�     f.�     I��M����(�0�����)�0�����(������(�������\���Y���(� ���������]�H���   ġx�I��M����(�����L������L��I9���   M����)�������(�������kY���Q] ��^���kYM�����Y ��^���X�H���   ����������I����w�6������������.nY �������v���������Y��x(�����Ż\���Y�0�������Ż]�H���   ġ{�L���   M��L���   L���   �x(�`�����(�0���M���fD  f.�     f.�     ��)�����@ f.�     f.�     L���   L���   L���   �x(������x(�`�����(�0�����l����0  H��H������H��H���H��H��H)�H���H������H��H�BH�� ���L�DH�����H�DH��`���H������L�|H������H�|H������H�tH������L�dIc�H���L�r�m  H��h���H��p���L��x���E1����  H��H���H�����E1�H��H�ڐf.�     f.�     f.�     H��(���ā|+��)�����CL+x��}�W ��Y���}�W ��X���^���)�������}�W ��Y�������}�W ��Y���)�����H���   ġ|(��)�����BL(x��}dW ��X���}_W ��Y��N���H��(���L���   ��}�V ��X���^���\��|(�����ŽY���}.W ��(���Y���Y�������}W ��Y���(�������(���}-W ��}5W ��ͨ���}-W ��Ũ���Y�ŽY���Y���\�������}%�V ��}-�V ��ը���Y���X���}TV ��}CV ��(������������Y���^�H���   ġ|(BL(x��Y���^���Y�������(���(�P��������ā|.CL.x��Y���������Y ā|/CL/xā|,CL,x��Y�����H��p���ġ|(BL(x��(�p��������H��x���ġ|(BL(xH��h���ġ|(BL(xI�� H�������L������L�����IcH��L)�H�����  āx���)�������YV ��XV ��^���)�p�����YV H��������w������)�P���H���   ġx���)�`�����X�U ��Y�U �����H������L��������(�U ��X���^���(�P�����Y�U ��\���Y���(%�U ��Y���Y�p�����Y�U ��(-�U ��(�`�����ɨ-�U ��ɨ-�U ��Y���Y���Y���X�@�����(%�U ���%�U ��Y���(1W ��(�������٨0W ��X���Y���^�H���   ġx���Y���^���Y�@�����(���(���������āx)���Y� �����W�T āx�āx���Y�0���H��p���ġx���(��������H��x���ġx�H��h���ġx�I��H��L9�L���   L��x���H��p���H��h�����  ā{���� �����YnS ��XnS ��^����������YbS ��w������������H���   ġ{����������X=S ��Y=S ����H��h���H��p���L��x���L������L���   ��]R ��X���^���\����������Y�R ��Y���%�R ��Y���Y�������Y�R ��-�R �����������-�R ��(���ѩ5�R ��Y���Y���Y���\�������%�R ���%�R ��Y���X���.R ��� ������%R ��Y���^���Y�����H���   ġ{���Y���^���X���(� ��������ā{���Y��������R ��W�ā{�ā{���Y�����ġ{���(��������ā{�ġ{�H��p���L��x���H��h���Mc)ā{T�����������Y�Q ��X�Q ��^����������Y�Q ��w������������H���   ġ{D�����������XaQ ��YaQ �����L������H���   ���P ��X���^���\����������Y2Q ��Y���%.Q ��Y���Y�������Y"Q ��-*Q �����������-Q ��(���ѩ5Q ��Y���Y���Y��x(�������C\���%Q ���%�P ��Y���X���bP ����������YP ��Y���^���Y�����H���   ġ{T����Y���^���X���Y� �����X�ā{\�����? ��Y���%�P ��W�ā{T��ā{T����Y������Y�L��p���ā{T����\�H��x���ġ{D��H��h���ġ{T��D����I��I����  1�A���  I���L��H���H��H��1�H����   1�H��`�����|D
A�
@  ��Y��|LA�@  ����D��|LA�@  ��YL��  ��X���D�@  ��|D
(��YD ��|L(����D(��|L(��YL0��X���D(H��@H��H����`���H��y;��|D
��Y��|L����D��|L��YL��X�H��`�����D
 Ic	H��H)�H���x=H��H����x���Y���x�������x�H����Y���X�H��`������H���H9�~<��{D���Y���{L���T�������{D���YD���X�H��`�����D�I��Ic	��{D����YD����{L����T��������{D��H���   H�������Y���X�H��`�����D����{��{��YO��{X ���������I�ø<  �����H�T�H��H��H�H�4�   H���H��I��I)�I���L�����L ��{^��X���sYM�jH�����{LM�D��Y���{BI���ʃ���x(�`�����(�0�����  H��H��H�   H��1�H��HH�H��H��h���H��A�   H��H������E1�L��f�f.�     f.�     L9�M��LN�M)�1�I���4  L��L��H��H��>L�H���1��    f.�     N�/ā{D�C�̀   ā{�ā{T�C�΀   ā{\�C�π   �������^�ā{D�C�Ȑ   ā{L�Ă��L�ā{D���^�ā{L�Ă��L�ā{T���^�ā{D�ā{L�ā{D� ā{T� Ă�D� ��^�H��H9�ā{D� �@���I��H�������Af.�     f.�     J�4/��{D���{����D���{L���^���{D�H��L9�|�L9�M��LN�M)�1�I���F  I���L��H��1�H�����   L������H��H���1���}�J fD  I�<M�>I�8ġ|�B��  L��H����|�A�ɸ  ������^��ȸ  ��)�=0���L��H�����L��H�����������=������^���)�=P���H��@H��H����o���H��`���H������H��H���L������H��y\I�4<H�����L������I��L��H��H��������I�8ġ|�I�>�����L��L��������}�I ��^����=0���L��L9�M��LN�L��L)�H)�H��|3J�/��x�H����x�������(�J ��^���)��0���H��M)�I9�~-J�/��{D���{����D���#I ��^�����0���1�I���  L��H��H��>L�H���1�D  f.�     J�/��{D�A�Ā   ��{���T��   ���������0�����������Y���{D�A�   ��{L�����L���Y��8�����{D���{L�����L���Y��@�����{L� ����L� ��Y��H�����{D���{L� H��H9��G����0D  J�/��{D���{����D���Y��0�����{D�H��L9�|�I��   I��   H��H;�h���H�����������L������H���   H���   H��������Ic��{L�����D����D�����K ��]���D���   9�L���   H���J  L�4�I��I��I���1�����   M��M��I��I��>M�I���1���;K f�     f.�     H��H)���{T��A��X�����{\��A��X����������]�H��H��L��H)���F��X�����{T�����D����]���F���{T�����D����]���F���{T�����D����]���F�H��L9��f���M��L9�L��}[H����xJ fD  f.�     H��H)���{T�����D����]�H��H��L��H)���A�H��L9�|�H��f.�     H�AH�H������H�`L���   I��L���   L���   ��l����Y  H��H�� ���H��H���H��H��H)�H���H������H��H�BH�� ���H�|H�����H�tI��I��IcH���  H������E1����  H���H��H��E1�H���W  H��`���H������E1�L�� ���L���   �    H��@���H���   ġ|0��}�E ��Y���}�E ��X���}�E ��Y���}%tE ��(���}%_E ��ݨ�L�� ���ā|7BL0xCL7x��^���}�E ��Y���)������c�����}jE ��Y���)�����ā|4��)�����CL4x��}1E ��X���},E ��Y�������}jD ��X���^���\���(�������Y���}	E ��(���Y���Y�������(�������(���}�D ��(���}�D ������}�D ��ͨ���Y���Y���Y���\�������}�D ��(���}�D ������Y���X�L������ā|D5 CL5xH���   ġ|L0 ��}2D ��Y���}-D ��X���^���}$D ��Y���}�C ��}%�C ��ݨ�ā|L7 ��}!D ��Y���)�����������}�C ��Y���)�p���ā|L4 ��)�������}�C ��X���}�C ��Y�����H��@�����}�B ��X���^���\���(�p�����Y���}�C ��Y���(�������(���}%�C ��}-�C ��ը���Y������}%|C ��ͨ���Y���Y���Y���\�������}bC ��}%aC ��ݨ���Y���X�ā|D5 I��@H��H�������L���   �x(������x(�`�����(�0���H������L�� ���H��`���H��x
I���  f�H��`���H���   ġ|0��}{B ��Y���}vB ��X���^���)�������}5B ��}$B �����ā|4��}HB ��Y�I���D�����}KB ��Y���)�����H���   ġ|0��)�������}B ��X���}B ��Y������L����(�0����x(�`����x(�������}/A ��X���^���}�A ��(�������Y���\���Y���Y�������}�A ��Y���}�A ��}%�A ��(�������ը���}�A ��ը���Y���Y���Y���\�������}�A ��}%�A �������Y���X�H������ġ|0L��`���L������M��Mcm L��L)�H��}L�� ����^  f�     f.�     H���   ġx���YkA ��XsA ��(�C �����C ��^���)�����āx)�L�� �����YQA I����w�V�����)�����H���   ġx���)�������X1A ��Y9A �$���L����(�0����x(�`����x(�������(!A ��X���^���(�������YA ��\���Y���(A ��Y���Y�������YA ��(%A ��(�������Ѩ%A ��Ѩ%#A ��Y���Y���Y���X�@�����(A ���A ��Y���X�H������ġx�I��L���   M9�L������L���   L�� ���H�������W  ā{���YS? ��XS? ��^����������? ����? ā{���Y0? M��M��I����w�����������H���   ġ{����������X? ��Y? �e���L��M��L�� �����(�0����x(�`����x(�����M��H������L���   L���   ��> ��X���^���\����������Y�> ��Y����> ��Y���Y�������Y�> ��%�> �����������%�> ��C\���٩=�> ��Y���Y���Y����> ��ѩp> ��Y���X�ġ{�H�    H���   H��������Ic	��{\D����YD��L���   ��{^D����R> ��W���D���   )ȅ���  E1����9  I��I���I���M��I��1�I����   L�� ���1���(_A I��I��D  f.�     f.�     I�<��L����T����P���M�T ��\���|T��A�����I�4��Y���T����^����������I���L���������L����\Lϸ��uYLʸ��^Lθ����H���I��I�����Lȸ�j���L��L���   L��L�� ���M��y6I���D����\D��H���YD��I���^D�����e@ H�:��D��L������Ic	H��L)�H���x3H��L)���xD����y\D����YD����y^D����W�< ��D��I��H��H���L9�~FL)���{D����{\D����YD����{^D�����< ��W���D��f�f.�     f.�     I�CI�Ic	�ȃ��H���   ��  1�����  H��H��H���H��H��1�H���  1���}2+ �     f.�     ��|��Y��|\A��  ��d������X���^���  ��|\A�   ��Y���^���\��  ��\��Y���Y���^���|TA�   ��|\ ��Yd ��|l(��X���ݘl(��^���]Yd(��^���l(��\l ��Y���Y���^���\���|L��|L(��\���|L(H��@H��H�������H��yi��|��Y��|T��X�����T��^���uYL��^���T��\��Y���}�) ��Y���^���|L��\���|D�f.�     Ic	H��H)�H���xf��x���Y�H��H����x���X�������^���qY���(� �����^������\���Y���Y�< ��^���x���\���x�H��H��H9�~a��{���Y���{T���\�������X���^���sYL���(�������^���T���\���Y ) ��Y���^���{XD���{D�H������H�`L���   H��   H������{�H���   ���H���   �{,�H���   �{<�H���   �{4�H���   ���H���   ��ȃ�l���H���   H���   u>ŃX���%"< ��.�r>��%D8 ��.�r0f.�     f.�     f.�     �AX���W��A W���W���)�0���Ic�A(W�H����  1��A(W҃��  H���H��H��1�����H����   1�����D  ����^�x  �x  ��|,A�x  ��\���|��|
A�
x  ��\���t ��^t ��|
��\���|d ��\���\���|t ��|D
 ��\���|D
 H��@H��H����j���H��y.����^��|,��\���|��|
��\���|
��\���}���X���y��[X�H��IcH��H)�H��|<�����^���x$���\���y���X��+\���x���x���\���x�H��H9�~.�����^���{$���\���{���{���\���{��+\���X���\�Mcġ{XD��ġ{D����%6 ��(�ġS^d��ŋY����������S^��Y�āXl����X�ā{l����W���.�s$��Y��"f.�     f.�     f.�     ��Y���\�ā{T��ā{T����^�ġ{d����^䃽l���tR���x����������L���   L���   ��(�������(������{�h����   f.�     f.�     I����\�ġ{D����\�ā{D��E���~  ���x�����������{������{�����L��N�4�O�<�K�L� H���   N��N��H���   N��1��{�4 �{%�8 f�H��H��H��H)���F���.�8 ��  Ż\�ųY���x4 ��\���^����8 ����8 H��H��L��H)���~���Y�L��H)���F�ŻY���\���]�L��H)���f���5W4 ��\���Y���w4 ��٩v4 ��Y���Y5"8 ��.�w
��.���  ��.�s|��X�H��H��L��H)���~�L��H)���XN���Y5�7 ��N������Y5�7 �����^�H��H)���F���y���X�L��H)���F��e@ f.�     f.�     ��X�H��H��L��H)���~�L��H)���XN���N���^�H��H)���F�L��H)�H�   `fq@H�F���(�H��H��H��H)�Ż\F�ųY����2 ��\���^���7 ����7 ��Y���\���]���\�L��H)���V���.�s;��\���\�H��H��L��H)���V���X���   @ f.�     f.�     ��\���.���   ��X�H��H��L��H)���N���\��   f�     f.�     ��X���Yl6 �����Y-h6 �����^���X�H��H��L��H)���V�H��H)���F�L��H)���N�L��H)���XF���    ��X�H��H��L��H)���F���W�H��H��L��H)�L��H)���N���R��⹫�H��H)�Ż\B�ųY�ţ\���^��x)���٩�5 ��Y���]�L��H)���b���.�s��\���\�H��H��L��H)���R���X�H��L9������H�CP�� L������Mc1�I��L���   ��(�������(������{�h����{�������  I���L��H��1�I����   �{�h����{�p����{�x����{�����H���   I��H���   H���   H���   I����������  f�     f.�     ���x����������L���   L���   ��(�������(������{�h���L������L���  �    f.�     f.�     ��������{�h����{�p����{�x����{�����1���}-0 �b}-\0 �b}5K0 H���   I��H���   H���   H���   I��H���   �     f.�     ��Lx��|$ALx������������������M-,Lx��(�����e\�ŭY���Y��b}=3 ��UY����������|�A]X��C-K�p�|�bE-Lx�A]X�����mK�p���BE-��}=3 �=Y��b}%(3 �A]Y��A^��UX��A=^�������mK�`�A|d ALx�CK�p�A|L �-X��C-K�`�AuX���uK�p���A%����]K�p��|$�|����=K�`��m^ALx��5K�`��|T ��}K�`����\���]K�`��|��|D ��l ť������������M-T ��(������\���Y���Y���mY�����������| �}X���EK� ��| �bm-L �A}X��|T �C-K� �|L �Bm-T �A=^��A|\ �C%K� ��}�1 �=Y���}�1 �}Y��A^��AUX���UK� ��}K%�yD  ������]X���5K� �EX���EK� ��| ��t ��=K� �A|T ��M^t ��-K� ��|t ��|D ��\���}K� ��l ��UK� ��T ��|D H��@H��H����a����{�h�����������������H��x$��������|   f.�     f.�     ��|��$���������b}=�, ��]������H���   ��M-,
��}=�, �x)��b}�, �bݨ�Ņ\�ŽY���Y���}=�0 ��Y����������|
�AuX��C=K�p�|
��������bE-�AuX��|�C%K�p�AUX��|�BE-$�����A5^��C-K�`�A|\ �C%K�p�A|\ �=X��b}-0 �A5Y��b}50 �AuY��A^��C=K�`�A]X���]K�p��$��uK%�yD p��|�|
�x(���5K�`��m^�|��%K�`��|T ��]K�`����\���uK�`��|L����(�������(������{������{�x����{�p����{�h����{�����IcL���   L���   L���   L���   L��H���   H9��.  �{�����fD  ��{���{���.+ u��   ��W���.���   �����%�* ��\���Y���-	+ ���-+ ��Y���Y%�. ��.�v>��X������[X���{���s^����H�   `fq@I����\���{��X@ ��X������{X$���{$���{$���Y-]. ��YM. ��������^������y���X���{�I��    H��H9�� ����&f.�     f.�     f.�     �{�������W���.�u{���H���@ ��W������u%��X�ŃX����H�����X�������\���Y�- ��\�x������������\�������X���+X�Ż\�������X��sX�Ic��{D����(� �����.���   I�L��I��H�L��H����G ��Y���{���YT����sY\�������{���sYd����{$���sYL����{�H������I��*���Y������Y�������{l����{D����T����{\����{d����{L���@A�Ic��{D����.�������   ��L����sYT����\����{d��������X���^���{T����L����{XD����^���{D����L����{D����{XD����{D����{D����{XD����{D����{D����{XD����{D��I�D��    H�D��    I�D��    H�D��    I�D��    I�D��    I�D��    I�D��    �@�A�Ic	H���	  ��(�������  1���}�+ H������ ����  ��]L ��]�H��@H��x���}���]���y���]����������.���  I��1�H��H8������   �(  H��H��X���H���H��H��1�H���}  1�f�     f.�     f.�     ��|�3@  A�3x  ��|3��|�5@  A�5x  ��|D5 ��|�6@  A�6x  ��|6���7@  �7x  ��7��|�2@  A�2x  ��|2��|�7@  A�7x  ��|7��|�0@  A�0x  ��|0��|�4@  A�4x  ��|4��|�3`  ��|D3 ��|�5`  ��|D5 ��|�6`  ��|D6 ���7`  ��D7 ��|�2`  ��|D2 ��|�7`  ��|D7 ��|�0`  ��|D0 ��|�4`  ��|D4 H��@H��H��������H��y��|�3@  ��|3��|�5@  ��|D5 ��|�6@  ��|6���7@  ��7��|�2@  ��|2��|�7@  ��|7��|�0@  ��|0��|�4@  ��|4H��X���H)�H��6�����   ��x��@  ��x�L������@  �����x��@  ��x�����@  �����x��@  ��x���x��@  ��x���x��@  ��x���x��@  ��x�H���f�L��H9�~`I���@  I��H���@  H��I���@  I��H���@  H��I���@  I��I���@  I��I���@  I��I���@  I�Ը8���AA�H�H��H��H������f�     f.�     f.�     ������|A�x  ���x  ��|A�x  ���x  ��|
A�
x  ��|A�x  ��|A�x  ��|A�x  ��|D ��D ��|D ��D ��|D
 ��|D ��|D ��|D H��@H���D���IǄ�@      HǄ�@      IǄ�@      HǄ�@      IǄ�@      IǄ�@      IǄ�@      IǄ�@      IcH=�   �  ����  L��E1����  I��I���M��I��E1�I���  L������E1�fD  f.�     f.�     H���   N�9��|D��A��X�����|��   M��O�\= ��|D��A��X�����|��   M��O�48��|D��A��X�����|��   H���   N�$9��|D��A��X�����|��   H���   J�9��D����X�������   H���   J�4:��D����X�������   H���   J�:��D����X�������   H���   J�<?��D����X�����|L����|���  ��|L����|���  M����|L����|���  M����|L����|���  ��L�������  ��L�������  ��L�������  ����   ��D�������  I���I��I����S���H���   L���   H���   H���   L���   L���   L������M����   J�>��D������   K�;��D������   K�>��D������   J�9��D������   J�?��D������   K�:��D������   K�8��D������   K�'��D������   L������H��L)�H��}I���   fD  ��D������  ��xD����x��  ��xD����x��  I����xD����x��  ��D������  ��xD����x��  ��xD����x��  ��xD����x��  I��L9�~nH��L)�H�T��H���  I�T��I���  I�T��I���  I�T��I���  H�T��H���  I�T��I���  I�T��I���  I�T��I���  �@dA���}�   ��|��?������|
��6��|��|��|��}� ��|$�    H������f�     f.�     ��|A�x  ��<�x  ������|A�x  ��4�x  ��|A�x  ��|A�x  ��| A� x  ��|A�x  ��|\ ��| ��|T ��t ��|T ��|D ��|T  ��|L H��@H���D���M��Mc��� ā{YD��H���   ġ{D���   D9��3  O��H���   J�<�D��I��1�����   M��M��I��I��>M�I���1���5 �f.�     f.�     H��H��L��H)���YQ���Y��@���������X�H��H)���B��X�����YQ���Y�������X���B���YQ���Y�������X���B���YQ���Y�������X���B�H��L9��r���M��L9�}J��� f�H��H��L��H)���YQ���Y�������X�H��H)���A�H��L9�|� f.�     L��h���E��H��   M����  M��A��  �6  1�A��H���   M����   H������H��H��1�H��|O1�f�     ������|
����|D
 ��D ��|D
@��D@��|D
`��D`H��H��H���|�H������H��yPfD  f.�     ������|
��H�� H��x�H��������f.�     f.�     f.�     Hc�L��H)�H��|Hc���W���x���ƃ�A9���   Hc�I��    H���   H��    �}   f.�     f.�     � �e 0�L��L��M���{�h����{�������w�(���� �e 0�H���   L�������{������{�h���L���   M��L���   H��   M��H�����H��H�������W�H;�(��������H��h�������   1҃���(�P���|wH������H��H��1�H��|;1���|A��  ��\���|��|D ��\���|D H��@H��H���|�H������H��y��|��\���|H������@ H��H)�H��|��(� �����yX���x�H��H9�~��{���\������{��{�����Mc��W���W�M��L��H����  1���W���W�A����  L��H���H��H��1�����H��}������   f.�     �Ax(�1�������}$ ��}� ��}% ����H���   I��H���   H���   f.�     f.�     f.�     ��,��  ������M-<��\���Y���| �E���B=-L ��ݨ���}K�`��\���UY���ݨ���}Kŀ��X��  ��XL H��@H��H���|��Ax(�H��yXH���   ����}G ����H���   ��e-$��}-� ��\���Y���}%# �������}K�0H���   ��X��}���X���y���X���}���X���y���X�I9�~kH���   @ f.�     f.�     H���   �����.� w)��Z ��\�H���   ��Y���� ������(���X�H��L9�|��{�h���H���   H��8  H������I��8  M��H������I��$8  H�� ����{������������� ���H�(%# H�� �����(	%# ��)������(�$# ��)� �����(�$# ��)�������(�$# ��)�������(�$# ��)�������(i$# ��)�������(I$# ��)�������()$# ��)�������(	$# ��)�������(�## ��)�������(�## ��)�p���H��H����H������H��H����H��0���Hǅ8���@   H�`     H��@���H�     H��H���H��P���HǅX���@   Hǅ`���   Hǅh���   Hǅp���   H��H0���H������L������M��H��H���H������Hǅ����@   H������H������H������Hǅ����@   Hǅ����   Hǅ����   Hǅ����   H��H����H�����H��H��p���H��H������� �e ��w����H��x���H����W�L9�����H�CxH�     H�CH��������� H�CX���h����� H�C`H�     H�ChH�     H�CpH�     H���   �  L������A�? �S  ���&  I��I��8  H�����E1�H��8  H�����fD  f.�     1�A�F������A����   H������H��H��1�H��|I1�H���   f�f.�     ����  ��\�����D ��\���D H��@H��H���|�H������H��y H���   ����\���H������@ L��H)�H��|H���   ��(� �����X����H��I9�~H���   �����\�������Ic1���W���W�H���[  H���H��H��1�����H��}�����    1�������}I ��}  ��}%? ����H���   H���    f.�     ��|l A��  ������M-<��\���Y���|| �E���b=-L ��ݨ���}K�`��\���UY���ݨ���}Kŀ��X��  ��XL H��@H��H���|�H��yS��|T ��}� ����H���   ��e-$
��}-2 ��\���Y���}%i �������}K�0H���   ��X
��}���X���y���X���}���X���y���X�H��IcH9�~m�f.�     f.�     f.�     ��{T� ��.� w)��� ��\�H���   ��Y���� ������(�H���   ��X�H��H9�|�H�����H������H�����H������H���   H��8  H�� ���Hǅ���    ��������� ���H�o!# H��������(P!# ��)�������(0!# ��)�������(!# ��)�p�����(� # ��)�`�����(� # ��)�P�����(� # ��)�@�����(� # ��)�0�����(p # ��)� �����(P # ��)������(0 # ��)� �����( # ��)�����H��H����H��(���H��H����H������Hǅ����@   H�`     H������H�     H������H������Hǅ����@   Hǅ����   Hǅ����   Hǅ����   H��H����H��@���L��p���H��H���H�����Hǅ���@   H�� ���H��(���H��0���Hǅ8���@   Hǅ@���   HǅH���   HǅP���   H��H���H������H��H������H��H������� �e ��w�d���I��L;�H�����(�P���������5  ���;  L������	  I��I��8  H������E1�H��8  H������I��L��H����    L��1�A�O������A��|H���   ��(� �����X �� �   I9�I��~H���   �����\�������H������H������H������H������H���   H��8  H�� ���Hǅ���    Hǅ���    Hǅ ���    H��"# H��������(f"# ��)�������(F"# ��)�������(&"# ��)�p�����("# ��)�`�����(�!# ��)�P�����(�!# ��)�@�����(�!# ��)�0�����(�!# ��)� �����(f!# ��)������(F!# ��)� �����(&!# ��)�����H��H����H��(���H��H����H������Hǅ����@   H�`     H������H�     H������H������Hǅ����@   Hǅ����   Hǅ����   Hǅ����   H��H����H��@���L��p���H��H���H�����Hǅ���@   H�� ���H��(���H��0���Hǅ8���@   Hǅ@���   HǅH���   HǅP���   H��H���H������H��H������H��H������� �e �}���I��M9��q����Z
  I��8  H��8  H�����E1�M��I��H���   H��8  H�� ���D  f.�     f.�     A�D$������Mc1���W���W�I���m  L��H���H��H��1�����H��}%����H���   ��   �f.�     f.�     1�������}� ��}� ��}%� ����H���   H���    f.�     ��|l A��  ������M-<��\���Y���|| �E���b=-L ��ݨ���}K�`��\���UY���ݨ���}Kŀ��X��  ��XL H��@H��H���|�H��yL��|T ��} ������e-$��}-� ��\���Y���}%� �������}K�0H���   ��X��}���X���y���X���}���X���y���X�I9�~SH���   D  ��{T� ��.� w"��? ��\���Y���v ������(�H���   ��X�H��L9�|�L������H�����H������H�� ���H�� ���Hǅ���    ��������� ���H�$# H�� �����(# ��)������(�# ��)� �����(�# ��)�������(�# ��)�������(�# ��)�������(e# ��)�������(E# ��)�������(%# ��)�������(# ��)�������(�# ��)�������(�# ��)�p���H��H����H������H��H����H��0���Hǅ8���@   H�`     H��@���H�     H��H���H��P���HǅX���@   Hǅ`���   Hǅh���   Hǅp���   H��H0���H������L������H��H���H������Hǅ����@   H������H������H������Hǅ����@   Hǅ����   Hǅ����   Hǅ����   H��H����H�����H��H��p���H��H������� �e ��w����I��L;�H����������  I��8  H������H��8  H������E1�M��H���   L��8  L��H��� f.�     A�D$������H������H������H������H������L�� ���Hǅ���    Hǅ���    Hǅ ���    H��# H�� �����({# ��)������([# ��)� �����(;# ��)�������(# ��)�������(�# ��)�������(�# ��)�������(�# ��)�������(�# ��)�������({# ��)�������([# ��)�������(;# ��)�p���H��H����H������H��H����H��0���Hǅ8���@   H�`     H��@���H�     H��H���H��P���HǅX���@   Hǅ`���   Hǅh���   Hǅp���   H��H0���H������L������H��H���H�E�H�E�@   H�M�H�U�H�E�H�E�@   H�E�   H�E�   H�E�   H��H���H�����H��H��p���H��H������� �e 诨��I��M9�������  I��H��h���L������M)�I��8  H������H��8  H������L��H��H������E1�M��I�� f.�     f.�     A�M1�������L��I��|_1�L��H���   �f.�     f.�     f.�     ����  ��\�����D ��\���D H��@H��H���|�H��yH���   ����\���L��I��|"H���   ��(� ���ġyX�ġx�H������H9�h���~H���   �����\�������H������H������H������H������H���   H��8  H�� ���Hǅ���    Hǅ���    Hǅ ���    H�<# H��������(# ��)�������(�# ��)�������(�# ��)�p�����(�# ��)�`�����(�# ��)�P�����(}# ��)�@�����(]# ��)�0�����(=# ��)� �����(# ��)������(�# ��)� �����(�# ��)�����H��H����H��(���H��H����H������Hǅ����@   H�`     H������H�     H������H������Hǅ����@   Hǅ����   Hǅ����   Hǅ����   H��H����H��@���H������H��p���H��H���H�����Hǅ���@   H�� ���H��(���H��0���Hǅ8���@   Hǅ@���   HǅH���   HǅP���   H��H���H������H��H������H��H������� �e ��w�*���I��L;�H�����(�P��������H�CXH�     H��H���A\A]A^A_H��]H��[��w�f�f.�     SH��H���UUH�kH�l$H��AWAVAUATH��P  H������L������H������Hc1�H��HI�H�IH��   H���I��I)�I���M�t� H��N�)L��Hc:H�C(L�[ L�cL�{H����  L������1�������  L������H��H���H��H����}�E1�H��}��)E���)�����H�������  ��)�����H������H������E1�M����)E�f�H�M�H�Cġ|8��}! ��Y���} ��X���} ��Y���}%�  ��(���}%�  ��ݨ�ā|)>BL8xCL>x��^���} ��Y���)�P����У����}�  ��Y���)�����ā|<��)�0���CL<x��}�  ��X���}�  ��Y�舣����}�� ��X���^���\���(�������Y���}v  ��(���Y���Y�P�����(�0�����(���}a  ��(���}\  ������}V  ��ͨ���Y���Y���Y���\U���}?  ��(���}:  ������Y���X�L������ā|D= CL=xH�Cġ|L8 ��}�� ��Y���}�� ��X���^���}�� ��Y���}Z� ��}%I� ��ݨ�ā|)L> ��}�� ��Y���)�p����c�����}j� ��Y���)�����ā|L< ��)������}6� ��X���}1� ��Y�� ���H�M���}k� ��X���^���\���(�������Y���}
� ��Y���(������(���}%� ��}- � ��ը���Y�p�����}%�� ��ͨ���Y���Y���Y���\U���}�� ��}%�� ��ݨ���Y���X�ā|D= I��@H��H�������M��H������H���U  H�Cġ|8��}.� ��Y���})� ��X���^���)�������}�� ��}�� �����ā|T= ��}�� ��Y�I���������}�� ��Y���)�����ā|<��)�������}�� ��X���}�� ��Y�贠��L����} � ��X���^���}�� ��(�������Y���\���Y���Y�������}�� ��Y���}�� ��}%�� ��(�������ը���}~� ��ը���Y���Y���Y���\U���}g� ��}%f� �������Y���X�H������ġ|8L������H��������(�������V� ��W�H��H)�H��}��)�����L�{L�[ L�������X  ��)�����L�{��)�p�����x���Y(� ��X0� ��(H� ����O� ��^���)�P�����x)\� ��Y� H������H��������w������)� ���H��������x���)�@�����X�� ��Y�� �ܜ��H������H��������(�� ��X���^���(� �����Y�� ��\���Y���(�� ��Y���Y�� ��Y�P�����(%�� ��(�@�����Ѩ%�� ��Ѩ%�� ��Y���Y���Y����p�����X���(�� ����� ��Y���X�L��������x)�H��L�[ H9��?  ��{����8������ ���� ��{D� ��{����h�����Y� ��X� ��Y� ���H������L������H��������w講��L������H������L�[ ��y���$� ��X����8�����X�� �����Y�� �����^���y���\���Y�� ��Y���%�� ��Y���Y���Y�� ���� ���h�����ɩ�� ��\��������5�� ��Y���Y���Y����� ��ѩo� ��Y���X�H��������{�I�    H������Hc H��������D����{\D����{YD����{^D����O� ��W���{D���   )�����  1�I��I������+  L��H���I��I��1�H����   L������H������1���(�� I��D  f.�     I�<ġ|L��ġu\L��I�B��P���ġ|T��B�����I���Y�ġ|T����^�B���������I�4ġ|L��B�����ġ|L��ġu\LǸġuYL��ġu^L������H���I��I���ġ|LƸ�d���L��H������L������M��y6I���D����\D��I���YD��I���^D������� J�2��D��H��H)�H���x5H��H)���xD����y\D����yYD����y^D����W�� ��xD��H��I9�~5H��H)���{D����{\D����{YD����{^D������ ��W���{D��I�GI�������x  H��������W���* H��H���1�����  H��H���H��H����}�1�H���  1���}A� ���|��eY$��|lA��  ��X���ݘl��^�A��  ��|lA�   ��Y���^���|lA��  ��U\,��Y���Y���^���|dA�   ��|l ��UYt ��||(��X���͘|(��^���MYt(��^���||(��E\| ��Y���Y���^���\���|\��|\(��\���|\(H��@H��H����
���H��ya��|��mY��|d��X����d��^���eY\��^���|\��e\��Y���}�� ��Y���^���|T��\���|LH)�H���xf��x���qY�H��H����x���X�������^���iYT� �����^���x���a\���Y���Y�� ��^���x���\���x�H��H9�~]��{���sY���{\���{d�������X���^���kYT���^���{T���k\���Y1� ��Y���^���{XD���{D�H��H���A\A]A^A_H��]H��[��w�@ f.�     f.�     SH��H���UUH�kH�l$H��AWAVAUATH��P  H������H������Hc1�H��HI�Hk�XH��H���I��I)�I���M��H��H��J�H������H�I�4�H�� ���H��H��J�4H��N�H��N�,L����*1��{(���� ��\�H�C@��8��� ��W�Lc:A�G�H�S(L�s ����  ��)�`�����)����L������L��H������� ��Y���)�������)�0�����Y���)�������)� ���L��H���H��p�����X���)�����E1�A����  �������L������L��p���I���L�����H������I��E1�L�����L������L�s0L�k(D  f.�     H�C ġ|8��)U�BL8x��}z� ��Y���}u� ��X���^���)E���}g� ��Y��f�����}m� ��Y���)�P���ā|L= ��)�p���CL=x��}3� ��X���}.� ��Y�������}l� ��X���^���\��|(�P���ŽY���}� ��(���Y���}� ��Y���YM���(�p�����(���}-�� ��}5�� ��ͨ���}-�� ��Ũ���Y�ŽY���Y���}������\���}%�� ��}-�� ��ը���Y���X���}/� ��}� ��(e��������Y���^�ā|>��Y���^���}�������Y���(���}� ��������H��H���ġ|)8BL8x��}�������Y����'� H�� ���ġ|8BL8xH�����ġ|)8BL8x��}�������Y�H������ġ|)8BL8x��}�0��������H������ġ|8BL8xH������ġ|8CL>xBL8xI�� I�������L�����L�s H������L������L�������������L�����H������L������L��L)�H�����  āx���)�P�����Y�� ��X�� ��^���)�@�����Y�� ���������w蚑����)� ���H�C(ġx���)�0�����Xx� ��Y�� �k������������({� ��X���^���\���(� �����Yo� ��Y���(%s� ��Y���Yw� ��(-� ��(�0�������-~� ��Y�@�������-}� ��Y���Y���Y����`�����X���(%m� ���%t� ��Y���X���(�� ��(�P�������� ��Y���^�H�K0ġx���Y���^����������Y���(� �����������H��H���ġx)����������Y���W� H�� ���ġx�ġx�āx)T� ���������Y�H������ġx)���(�0�����������L������āx�L������āx)�I���%L��������(�0�����(� ���L������H������L9�p���L��H�����(����H�S(��  ā{�������ġ{���������Y�� ��X�� ��Y�� ���L��������)�0�����)� ����������I��L�����M��H��������w�?�����(����L��H���H������M��L�����L��H�K0L������L�s ��y��{�� ��kX���������X� �����Y� �x(���a���^���y�Ż\���Y� ��Y���-� ��Y���Y���Y� ��%� ���������%�� ��(���٩5�� ��Y����������Y���Y���(� �����\���%�� ���%�� ��Y���X���<� �Ⱪ;� ��Y���^���Y�ġ{���Y���^���X���Y���kX�ā{����� ��Y���%�� ��W���(�0���H�� ���ġ{�ā{T� ��Y���Y�ġ{�Ż\�ā{�ā{���)�������������)� ���L��H�����)�0���L�����L������ā{D�����x���ġ{L�����������Y�� ��X�� ��Y�� ���M��I��M��M����w�1���M��M��M��L��������y��{�� ��kX����x�����X.� �����Y� ��(���a���^���y�Ż\���Y#� ��Y���-� ��Y���Y���Y� ��%� �����������%� ��(���٩5	� ��Y���Y���Y���(���\������%�� ���%�� ��Y���X���X� ��ѩW� ��Y���^���Y�����H�C0ġ{T����Y���^���X���Y� �����kX�H��H���ġ{\������ ��Y���%�� ��W�H�� ���ġ{T��H�����ġ{T����Y�0�����Y�ā{T��Ż\�ā{D��ā{T����(�D�����H��H��I��H��������  E1�M��I���A���  M��I���L��H��1�I����   1Ґf.�     f.�     ��|TA�@  ��mY��|\A�@  ���T��|\A�@  ��eY\A��  ��X���T�@  ��|T(��mYT ��|\(���T(��|\(��eY\0��X���T(H��@H��H����[���H��y3��|T��mY��|\���T��|\��eY\��X���TL��L)�H���x9L��H����xT� āiY���x�������x�I��āaY���X����M9�~9ā{T�ākY�ā{\�ā{d�����ā{T�ākYT���X�ġ{T�L������Ic H�S��YT��āsYL����ѩ���X�ġ{D����{E ��{X��{��sYN��{�������A$H�� L�D$H�CH�D$H�$0�L��M��L��������w�9   H�� I�FI�H��H���A\A]A^A_H��]H��[�f�f.�     f.�     SH��H���UUH�kH�l$H��AWAVAUATH��  H������H��x���McE1�M��MI�L��H��H��H���I��I)�I���L����d� ��^����Y������ā{�O����Y��{M��I���L��h���H�s H�KD�؃����  L��X���L��H�   L��`���H��1�H��HH�H��H��p���A�   E1��f.�     f.�     f.�     H������M9�M��MN�M)�H��x���N��N�$�J��1�I���9  L��H��H��>L�H���1�fD  f.�     f.�     N�,>��D���   ā{���{T�A��   �������{D�A���   ��^�ā{D�C��   ��L�����L���{D���^�ā{D���L�����L���{D���^�ā{D���L� ����L� ��{D� ��^�ā{D� H��H9��E���L��h����H�     f.�     f.�     J�>��D���{����D���{L���^���{D�H��L9�|�M9�M��MN�M)�1�I���  I���M��I��1�I�����   1���}� �f.�     f.�     M�,	H�
L�$ā|L� C���  L��H����|�A���  ������^���  ��)�����L��H�����L��H����|T� ������H�����^���)�����H��@I��I����l���L��h���M��y:H�
L��H�����I�	ġ|�H��������};� ��^��������L��M9�M��MN�L��L)�H)�H��|2J�8��x�H�����������(y� ��^���)�Ő���H��M)�I9�~,J�8��D���{����D����� ��^����Ő���H������N�$�J�4�1�I���  M��I��I��>M�I���1�f�f.�     f.�     J�8��D��ƀ   ��{���{T�A�Ā   ��������Ő����������Y���{D�A�ʐ   ��L�����L���Y�Ř�����{D���L�����L���Y�Š�����L� ����L� ��Y�Ũ�����{D���{L� H��L9��J����2�     J�8��D���{����D���Y�Ő�����{D�H��L9�|�I��   I��   H������H��H;�p��������L��`���L��X���H�KH�s Hc��D��ā{L��Ă�D��ā{D����� ��]�ā{D���   D9��)  K��1�A����   M��I��I��>M�I���1����� f�f.�     L��H)���{T��A��X�����{\��A��X����������]�H��H��H��H)���A��X�����{T�����D����]���A���{T�����D����]���A���{T�����D����]���A�H��L9��f���L9�}S��!� D  f.�     f.�     L��H)���{T�����D����]�H��H��H��H)���F�H��L9�|�H��H���A\A]A^A_H��]H��[��w�f�     f.�     UAWAVAUATSH���  I��H��M��Hc.L��$  L��$  H���v)  ��{���$�  ��w�ԁ�����$�  ��W���*���$�  ��{��)�$�  ��{$��)�$�  ��zD ��   �   L��薃������ ���$�  ����� ��)�$�  ���� ����� ��(���)�$P  ��u� ���t� ��)�$@  ��}����$�  �����)�$0  �����)�$p  ��}����$�  �����)�$`  ��}����$�  ���$�  ��Y�� ��Z���)�$P  ��}����$`  �����)�$p  ���$�  ��Yh� �����)�$�  ��}����$   ��YM� �����)�$�  ��}����$   I��I���H����(�$�  �����)�$P  ��}����$�  ��(�$�  ��}����$@  �����)�$`  �K  ��zD ��   �   L����w�����H��L��$  �   ��zD ��   �   L����w�Ł������ ���$�  ����� ��)�$0  �����)�$p  ��}�H����  ���$@  ��zD ��   �   L����w�d���1�H����  �����$@  �  H��$�  ��}�� ��_����$   ���$�  ��_����$   L��H��1����$   �     f.�     ��|L- ���$�  AL-x��}�� ����M�4,��P���tq��u-��Y�$�  ���$`  ��}X� ��^���(�$�  ��sY$��E� ��^���}���\�蓀�����$`  ���$   ��Y�$   ��Y�$   ���$�  ��v�������P�����   ���$@  ���$�  ��}Y��}�� ��^�ANx��(�$�  ��sY$���� ��^���}���\�� ������$   ���$@  ��uK� ���$   ��Y�$   ��Y�$   ���$�  ��eK� �f�f.�     ���$�  ���$@  ���$�  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H���  ��xD� ��)�$   ����� ��P���tz��q-���Y�$P  ��)�$0  ��(�� ��^���(�$�  ��sY$���� ��^������\���w�|����(�$0  ��(���(�$p  ��_M� ��Y���Y�$�  ��)�$p  ��)�$  ��v���W���P����  ��)�$   ��(�$P  ��yY���(� ��^���(�$�  ��sY$���� ��^������\���w��{����(�$p  ��(�$   ��qK� ��(�$0  ��_�� ��Y���Y�$�  ��(�$  ��YK� ��(�$p  ��(�$`  �  ��zD ��   �   L����w�}��H��L��$  �5  ��zD ��   �   L����w�z}��1�H���P  ���  H��$�  ��}�� ���$�  ��_����$�  ���$�  ��_����$�  L��H��1����$�  �     ��|L- ���$@  AL-x��}�� ����M�4,��P���tq��u-��Y�$�  ���$�
  ��}x� ��^���(�$�  ��sY$��e� ��^���}���\��|�����$�
  ���$�  ��Y�$�  ��Y�$   ���$�  ��v�������P�����   ���$@
  ���$�  ��}Y��}�� ��^�ANx��(�$�  ��sY$���� ��^���}���\�� |�����$�  ���$@
  ��uK� ���$�  ��Y�$�  ��Y�$   ���$�  ��eK� �f�f.�     ���$�  ���$@  ���$@  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H����  ��xD� ��)�$�  ����� ��P���tz��q-���Y�$P  ��)�$�  ��(�� ��^���(�$�  ��sY$���� ��^������\���w�x����(�$�  ��(���(�$`  ��_m� ��Y���Y�$�  ��)�$P  ��)�$�  ��v���W���P�����  ��)�$�  ��(�$P  ��yY���(:� ��^���(�$�  ��sY$��� ��^������\���w�x����(�$P  ��(�$�  ��qK� ��(�$0  ��_�� ��Y���Y�$�  ��(�$�  ��YK� ��(�$p  ��(�$`  �U  ��zD ��   �   L����w��y��H���:  ��zD ��   �   L����w�y��H����  1ۃ��  H��$�  ��}�� ���$�  ��_����$�  ���$�  ��_����$`  L��H��1����$`	  fD  f.�     ��|L- ���$�
  AL-x��}�� ����M�4,��P���tq��u-��Y�$�  ���$�  ��}�� ��^���(�$�  ��sY$���� ��^���}���\���x�����$�  ���$`	  ��Y�$�  ��Y�$   ���$   ��v�������P�����   ���$`  ���$�  ��}Y��}
� ��^�ANx��(�$�  ��sY$���� ��^���}���\��@x�����$`	  ���$`  ��uK� ���$`	  ��Y�$`  ��Y�$   ���$   ��eK� �f�f.�     ���$   ���$@  ���$�
  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H����  ��xD� ��)�$�  ���� ��P���tz��q-���Y�$P  ��)�$�   ��(�� ��^���(�$�  ��sY$���� ��^������\���w��t����(�$�   ��(���(�$`  ��_�� ��Y���Y�$�  ��)�$  ��)�$  ��v���W���P�����  ��)�$�   ��(�$P  ��yY���(Z� ��^���(�$�  ��sY$��/� ��^������\���w�;t����(�$  ��(�$�   ��qK� ��(�$0  ��_�� ��Y���Y�$�  ��(�$  ��YK� ��(�$p  ��(�$`  �U  1ۃ��  H��$�  ��};� ��_����$�  ���$�  ��_����$�  L��H��1����$�  �f.�     f.�     f.�     ��|L- ���$`  AL-x��}I� ����M�4,��P���tq��u-��Y�$�  ���$�
  ��}�� ��^���(�$�  ��sY$���� ��^���}���\��3u�����$�
  ���$�  ��Y�$�  ��Y�$   ���$�  ��v�������P�����   ���$`
  ���$�  ��}Y��}j� ��^�ANx��(�$�  ��sY$��R� ��^���}���\��t�����$�  ���$`
  ��uK� ���$�  ��Y�$�  ��Y�$   ���$�  ��eK� �f�f.�     ���$�  ���$`  ���$@  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H���  ��xD� ��)�$�  ���r� ��P���tz��q-���Y�$P  ��)�$�  ��(Z� ��^���(�$�  ��sY$��/� ��^������\���w�;q����(�$�  ��(���(�$p  ��_�� ��Y���Y�$�  ��)�$`  ��)�$�  ��v���W���P����  ��)�$�  ��(�$P  ��yY���(�� ��^���(�$�  ��sY$���� ��^������\���w�p����(�$`  ��(�$�  ��qK� ��(�$p  ��_B� ��Y���Y�$�  ��(�$�  ��YK� ��(�$p  ��(�$`  �u  1ۃ��  H��$�  ��}�� ���$�  ��_����$   ���$�  ��_����$   L��H��1����$�  f�f.�     f.�     ��|L- ���$�  AL-x��}�� ����M�4,��P���tq��u-��Y�$�  ���$ 
  ��}X� ��^���(�$�  ��sY$��E� ��^���}���\��q�����$ 
  ���$�  ��Y�$   ��Y�$   ���$   ��v�������P�����   ���$ 
  ���$�  ��}Y��}�� ��^�ANx��(�$�  ��sY$���� ��^���}���\�� q�����$�  ���$ 
  ��uK� ���$�  ��Y�$   ��Y�$   ���$   ��eK� �f�f.�     ���$   ���$@  ���$�  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H���  ��xD� ��)�$�  ����� ��P���tz��q-���Y�$P  ��)�$�  ��(�� ��^���(�$�  ��sY$���� ��^������\���w�m����(�$�  ��(���(�$`  ��_M� ��Y���Y�$�  ��)�$@  ��)�$�  ��v���W���P����  ��)�$p  ��(�$P  ��yY���(� ��^���(�$�  ��sY$���� ��^������\���w��l����(�$@  ��(�$p  ��qK� ��(�$p  ��_�� ��Y���Y�$�  ��(�$�  ��YK� ��(�$p  ��(�$`  �u  �����$@  �  H��$�  ��}�� ��_����$�  ���$�  ��_����$@  L��H��1����$@	  @ f.�     f.�     ��|L- ���$�
  AL-x��}	� ����M�4,��P���tq��u-��Y�$�  ���$�  ��}�� ��^���(�$�  ��sY$���� ��^���}���\���m�����$�  ���$@	  ��Y�$�  ��Y�$   ���$   ��v�������P�����   ���$@  ���$�  ��}Y��}*� ��^�ANx��(�$�  ��sY$��� ��^���}���\��`m�����$@	  ���$@  ��uK� ���$@	  ��Y�$@  ��Y�$   ���$   ��eK� �f�f.�     ���$   ���$@  ���$�
  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H����  ��xD� ��)�$�  ���2� ��P���tz��q-���Y�$P  ��)�$�   ��(� ��^���(�$�  ��sY$���� ��^������\���w��i����(�$�   ��(���(�$p  ��_�� ��Y���Y�$�  ��)�$   ��)�$   ��v���W���P�����  ��)�$�   ��(�$P  ��yY���(z� ��^���(�$�  ��sY$��O� ��^������\���w�[i����(�$   ��(�$�   ��qK� ��(�$p  ��_� ��Y���Y�$�  ��(�$   ��YK� ��(�$p  ��(�$`  �=  ���  H��$�  ��}]� ���$�  ��_����$   ���$�  ��_����$�  L��H��1����$ 	  @ f.�     f.�     ��|L- ���$�	  AL-x��}i� ����M�4,��P���tq��u-��Y�$�  ���$   ��}� ��^���(�$�  ��sY$��� ��^���}���\��Sj�����$   ���$ 	  ��Y�$   ��Y�$   ���$�	  ��v�������P�����   ���$�  ���$�  ��}Y��}�� ��^�ANx��(�$�  ��sY$��r� ��^���}���\���i�����$ 	  ���$�  ��uK� ���$ 	  ��Y�$�  ��Y�$   ���$�	  ��eK� �f�f.�     ���$�	  ���$�	  ���$@  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H����  ��xD� ��)�$@  ����� ��P���tt��q-���Y�$P  ��)L$`��(}� ��^���(�$�  ��sY$��R� ��^������\���w�^f����(L$`��(���(�$`  ��_� ��Y���Y�$�  ��)�$�   ��)�$   ��v���W���P�����  ��)D$@��(�$P  ��yY���(�� ��^���(�$�  ��sY$���� ��^������\���w��e����(�$�   ��(T$@��qK� ��(�$p  ��_n� ��Y���Y�$�  ��(�$   ��YK� ��(�$p  ��(�$`  ��(�$@  �q  1ۃ��  H��$�  ��}�� ���$�  ��_����$   ���$�  ��_����$�  L��H��1����$ 	  D  f.�     f.�     ��|L- ���$�	  AL-x��}�� ����M�4,��P���tq��u-��Y�$�  ���$   ��}x� ��^���(�$�  ��sY$��e� ��^���}���\��f�����$   ���$ 	  ��Y�$   ��Y�$   ���$�	  ��v�������P�����   ���$�  ���$�  ��}Y��}�� ��^�ANx��(�$�  ��sY$���� ��^���}���\�� f�����$ 	  ���$�  ��uK� ���$ 	  ��Y�$�  ��Y�$   ���$�	  ��eK� �f�f.�     ���$�	  ���$�	  ���$@  ��\���Y����$`  �������]���|D- H�� H���V���L��H��$�  H��H)�H����	  ��xD� ��)�$P  ����� ��P���tt��q-���Y�$P  ��)L$p��(�� ��^���(�$�  ��sY$���� ��^������\���w�b����(L$p��(���(�$`  ��_s� ��Y���Y�$�  ��)�$�   ��)�$0  ��v���W���P�����  ��)D$P��(�$P  ��yY���(C� ��^���(�$�  ��sY$��� ��^������\���w�$b����(�$�   ��(T$P��qK� ��(�$p  ��_�� ��Y���Y�$�  ��(�$0  ��YK� ��(�$p  ��(�$`  ��(�$P  �q  E1����  H��$�  ��}� ���$�  ��_����$�  ���$�  ��_����$�  L��H��1����$`  @ f.�     f.�     ��|L- ���$�  AL-x��})� ����M�4,��P���tq��u-��Y�$�  ���$   ��}ؿ ��^���(�$�  ��sY$��ſ ��^���}���\��c�����$   ���$`  ��Y�$�  ��Y�$   ���$�  ��v�������P�����   ���$�  ���$�  ��}Y��}J� ��^�ANx��(�$�  ��sY$��2� ��^���}���\��b�����$`  ���$�  ��uK� ���$`  ��Y�$�  ��Y�$   ���$�  ��eK� �f�f.�     ���$�  ���$�  ���$@  ��\���Y����$`  �������]���|D- H�� H���V���M��H��$�  H��L)�H���|  āxD� ��)�$�   ���R� ��P���ttĂq-���Y�$P  ��)L$��(=� ��^���(�$�  ��sY$��� ��^������\���w�_����(L$��(���(�$`  ��_�� ��Y���Y�$�  ��)T$ ��)�$�   ��v���W���P�����  ��)$��(�$P  āyY���(�� ��^���(�$�  ��sY$��|� ��^������\���w�^����(L$ ��($��qK� ��(�$p  ��_6� ��Y���Y�$�  ��(�$�   ��YK� ��(�$p  ��(�$`  ��  ��(�$p  ��(�$`  ��(�$  ��(�$   ��\���Y��������]���xD� H��H9���  ��(�$�  ��kY���c� ��^���kY$���� ��^���X���{L� ���$  ��w�Q_�����$  ��.%�� �V  �  ��(�$p  ��(�$`  ��(�$�  ��(�$�  ��\���Y��������]���xD� H��H9���  ��(�$�  ��kY����� ��^���kY$��� ��^���X���{L� ���$�  ��w�^�����$�  ��  ��(�$p  ��(�$`  ��(�$�  ��(�$�  ��\���Y��������]���xD� H��H9��8  ��(�$�  ��kY���� ��^���kY$��S� ��^���X���{L� ���$�  ��w�^�����$�  �  ��(�$p  ��(�$`  ��(�$�  ��(�$�  ��\���Y��������]���xD� H��H9���  ��(�$�  ��kY���u� ��^���kY$���� ��^���X���{L� ���$�  ��w�c]�����$�  ��  ��(�$p  ��(�$`  ��(�$  ��(�$�  ��\���Y��������]���xD� H��H9���  ��(�$�  ��kY���ս ��^���kY$��� ��^���X���{L� ���$�  ��w��\�����$�  ��.%�� �H  ��Y�$�  ��Yӹ ��(�$�  �C  ��(�$p  ��(�$`  ��(�$   ��(�$�  ��\���Y��������]���xD� H��H9��0  ��(�$�  ��kY���� ��^���kY$��K� ��^���X���{L� ���$�  ��w��[�����$�  ��.%2� �d  ��Y�$�  ��Y� ��(�$0  �{  ��(�$p  ��(�$`  ��(�$@  ��(�$   ��\���Y��������]���xD� H��H9��h  ��(�$�  ��kY���E� ��^���kY$���� ��^���X���{L� ���$`  ��w�3[�����$`  �   ��(�$p  ��(�$`  ��(�$P  ��(�$0  ��\���Y��������]���xD� H��H9���   ��(�$�  ��kY����� ��^���kY$��� ��^���X���{L� ���$h  ��w�Z�����$h  ��.%ʷ v��Y�$�  ��Y�� ��(�$P  ���Y�$�  ��Y�� ��(�$@  ��_
� ��Y���(�$P  ��(�$�  ��\���Y�������]���{D� H���  [A\A]A^A_]��w���(�$p  ��(�$`  ��(�$�   ��(�$�   ��\���Y��������]�āxD� I��L9�~���(�$�  ākY����� ��^���kY$���� ��^���X�ā{L� ��L$8��w�tY����d$8��.%�� v��Y�$�  ��Y�� ��(�$P  ���Y�$�  ��Yw� ��(�$@  ��_� ��Y���(�$P  ��(�$�  ��\���Y�������]�ā{D� �����f�     f.�     UAWAVAUATSH��  H��$�   I��L��$�   L��$�   H��$   H��$�   HcH��$�  ��\��H��$�  ��T��H��$�  �{L��H��$�  �{D��H��$�  ��D��H��$�  ��d�����$  H��$0  �L��$P  L��$H  H��$@  L��$8  ��H��$  H��$   L��$�  L��$�  L��$�  L��$�  H��$�  H��$�  H��$X  uLųX���{)��:���� ��.�r��ʴ ��.�s#��$�   ��X����$  ��
ųX���{	�)��$�   ��kX�HǄ$      HǄ$      �A0W�L��I�     Mc$M����  L��M��1���W�A���  L��H���I��I��E1�����H����   E1�����fD  ġ|ġu^B�x  B�x  ā|$C�x  ��\�ā|ġ|B�x  ��\�ġ|l ġU^l ġ|��\�ā|\ ��\���\�ā|l ġ|L ��\�ġ|L I��@I��I����d���M��y0ġ|ġu^ā|$��\�ā|ġ|��\�ġ|��\�L��$8  ��}���X���y���X�L��H)�H��|:�����^���x$���\���y���X���\���x������\����H��I9�~,�����^���{$���\���{������\������\�M����{I��H��$�  H��$�  H��$X  ŻX���\�ġsX\��ġ{\����sXM ��{M ��M ġqL����(%�� ��^���Y�����(�āCXL����y�ŻY���X�ā{L����/�x(���W���.�s��Y����Y���\�ā{L��L��$`  H��$�   ��ز ��W���{H��$h  ��ā{T����iU ��^���W���*��y���Y� ġyD����Y�H��$�  ����(����$  ��^���W�H��$�  ����$�   H��$X  �  H��$�   L��$�   H��$   H��H��$�   	L��$p  AH��$`  �$  �$  L��$�  A
A
L��$  A$L��$   AM M H��$`  H�D$HH�D$HH��$�   H�l$@L�l$8L�d$0L�T$(H�\$ H�T$L�|$H��H  H�D$H��H  H�$0�M��H��H��$�   L��$`  �{�$�   ���$�   ���$�   �{�$�   ���$�   ���$�   ��w�  H��$   H��$�   H��$`  H��$�   	AM AA$AM L�l$L�d$L�t$H�$0�L��$�   L��$`  L��L��M��I���
  ���$�   ���$�   �{�$�   ���$�   ���$�   L��$�  L��$�  L��$8  H��$X  �{�$�   M��H��$  ��W�A�;uųX���{M ������\E ��Y<� H��$�  ����\���;X��X���CXH��$x  ����\���X���X��� Ic$L��$P  ��{��y.D����   H��$�   H�D$8H��$�   H�D$0H��$�   H�D$(H��$(  H�D$ H��$   H�D$H�t$L�T$H��$  H�$0�M��H��$   H��L��L��L��$�  L��$�  H����w��  M��H��Ic$��{D����{��Y� ��.�w_H��$(  H�D$H��$   H�D$H�t$L�$0�L��H��$   I��L��H��L��$�  L��$  H��$�  ��w�7  I��L��A�<$�	  ��   1���}ȱ H�������    f.�     ����  ��]L ��]�H��@H��x���}���]���y���]���{M ��X� ��.�w_H��$(  H�D$H��$   H�D$H�t$L�$0�M��H��$   H��L��H��L��$�  L��$  H��$�  ��w�^  M��H��A�<$�   iH��$(  H�D$(H��$   H�D$ H�t$L�T$H��$  H�D$H��$�  H�$0�L��H��$   L��L��I��L��$�  L����w�  I��Mc$ā{D����Yr� H��$  ġ{D���   D9��  K��J��D��I��1�����   L��H��H��>L�H���1���!� �H��H��H��H)���YQ���Y��@���������X�H��H)���B��X�����YQ���Y�������X���B���YQ���Y�������X���B���YQ���Y�������X���B�H��H9��r���L9�}@���� D  H��H��H��H)���YR���Z�������X�H��H)���B�H��L9�|�H��  [A\A]A^A_]��w�@ f.�     f.�     UAWAVAUATSLc*H�L$Pġ{D��H�T$8��
��\�ġ{D��H�T$Xġs^D��H�\$`ġ{L����\�ġ{D��L��$�   H�D$xL�|$pH�t$hL�d$HL�t$@M����  J�<�J��J�,�O�$�N��O�<���{�{I� �3^�1��{S� f�     f.�     H��H��H��H)���N���.� ��  ��\�ūY����� ��\���^��� � ����� H��H��H��H)���~���Y�L��H)���^��{Y�ţ\���]�L��H)���N���%ש ��\���Y���5�� ���5�� ��Y���Y5�� ��{&��.�w
��.���  ��.�sw��X�H��H��H��H)���~���kX3��{3L��H)���Xv���Y-\� ��v������Y5S� �����^�H��H)���^���y���X�L��H)���N��`f�     ��X�H��H��H��H)���~���KX��{L��H)���XN���N���^�H��H)���N�L��H)�H�   `fq@H�F���(�H��H��H��H)���\N�ūY���s� ��\���^���{� ���z� ��Y�ţ\���]���\�H��H)���V���.�s1��\���\�H��H��H��H)���V���X��9�    f.�     ��\���.���   ��X�H��H��H��H)���V���\���{�   @ f.�     ��Y� ��X������Y5� �����^���X�H��H��L��H)�I�    ��N�H��H)���V�H��H)���n���[X��{L��H)���XN���N��1�    f.�     ��X�H��H��H��H)���N�I�     H��H��H��H)�L��H)���V���Z������H��H)���\J�ūY�ų\���^��x)���٩� ��Y���]�H��H)���Z���.�s$��\���\�H��H��H��H)���Z���kX��{H��L9��q�����{H�D$x��{��W���.�u{�� I�    [A\A]A^A_]��w�f�     f.�     f.�     Hc��H�T$��\L����^���*� ��\���^���2� ���1� ����(���{T����Y�H�L$��\���������]�H�L$��L����{ ��.�s��\���\���L����X���{ ���\���.�v��X���L����\���{ ��w���X���D��I�     ��w�f.�     f.�     f.�     AVSPLcH�T$8H�t$0L�t$(H�D$ M��M���   1�A����  I���L��H��1�����I���!  1�������}-D� ��}5k� �    f.�     f.�     �|Lx�A|ALx�����Am���-���A%���B-ALx�A|(���}� �bͨ��AU\��A%Y��A%Y���}=�� �5Y��A=���A%���A|,�A=X��CK��A|,�b%-<Lx�A=X�����mKװ���B%-<ŵX��A^��A����mK����Lx��eK߰����}%A� �Y���}%,� ��(�ŽY���]^��AX��CK��ŭX���-K���|$�c=K%@zD ��|�A|,�CK���A^,��ALx��eK������uX���uKҰŵX���]K����|��-\���-K������D �A|L ���������5���A]����--d ��mK����X��A|(���}�� �bͨ��AU\��A%Y�ŽY��]Y���}����]���A|\ �A}X��C%K�@�A|\ �b]-l ��mK�@��}X��|l ��K�@��\ �B]-l �A^���}�� �Y���}� �}Y��A^��|| �CK�@�A5X��C5K�@��}K<%@zD @ŭ��ŽX���mK�@ŽX���eK�@��%X��c%K�@�A|T ��T ��K�@�|l ��m^T ��K�@��T ��| ��E\���EK�@�A|L ��5K�@��|\ ��T H��@H��H�������H���Z  ��|������������}� ����������E-,��}5,� �b}� �b����\���MY���Y���}5�� ��Y����������A|�AmX��C=K�`�A|�bM-�AmX��|�C%K�`�AUX��|�BM-$�����A5^��C-K�p�|�C%K�`�|�=X��b}-Z� �A5Y��b}5D� �AmY��A^��C=K�p�AeX���eK�`��|��mK�`��$�A|�C5K�p�A=^�|�C%K�p�|��X���}K�`��X���}K�p��eK�p��|��\���]K�p����}���X���y���X���X��L��I9��.  �    �����{���W���.��  ��.� u��   ��{���i� ��\���Y���%�� ���%�� ��Y���Y4� ��.�v^��X���{���X����X������s^����H�   `fq@I����\�����yf.�     f.�     f.�     ��X���{���X����X������{���Y%�� ��Y�� ��������^������y���X���{�H��    H��L9������H��[A^��w�@ f.�     UAWAVAUATSH��   L��H��$�   ��*�� � ��^���{H��$  ����
��
H��$�   ��
L��$  ����	��{	H��$�   �{H��$�   ��H��$�   ��L��$(  ������{��H��$   �� L��$H  ��{.L��$@  �ASXH��$X  ��2H��$P  �KXL��$x  ��{} L��$�  ��{,$H��$�  ��X+H��$p  ��X;I�    L��$�  ��{7H�     H��$�  ��X0H�    H��$�  ��H�     H��$�  ��XH�     H��$0  ����] ��] L��$  �B��L��$   ����H��$8  ���� �� �A{�{
��{} I��L����{7�A{�x)���{��{M ��+��(ŋD�6�A��H�΅�H��$h  H��$`  ��  ��r*���%�� ��^���Z���W���[*�L��$  ��sY.��5� ��W���{6��Y���{��sY��{��Y���{��YU ��U ��\���sY(��{(H��$8  I����sY$��X���X���X���{$���$�   ��^���H��$�   H��$�   AAA	A
AM A$�$�   H��$H  
	H��$�   H��H��$X  
H��$�  H�t$xH��$�  H�t$pH�\$hH��$�  H�t$`L�l$XL�|$PH��$x  H�\$HH�T$@H�L$8H�|$0H�D$(H��H�   H�D$ L�d$H�l$L�\$L�$0�H��$�   H��$�   H��$�   L���Z  H�Ĩ   [A\A]A^A_]��w�fD  f.�     UAWAVAUATSH��xM��M��H�|$HHcH��$   H�T$pH��$  H�T$hH��$  H�T$`H��$  H�T$XH��$   H�T$PH��$�   L��$�   H��L��$�   L��$�   L��$�   H��$�   H�T$@H��$�   L��$�   L��$�   H��$�   H���q  H�t$0=  ��  H�T$81҃��  H���H��H��1�H����   1�f.�     ������1��|D5 ��|7��3��|4��|6��D1 ��|D5 ��|D7 ��D3 ��|D4 ��|D6 ��D1@��|D5@��|D7@��D3@��|D4@��|D6@��D1`��|D5`��|D7`��D3`��|D4`��|D6`H��H��H����O���H��y;f.�     ������1��|D5 ��|7��3��|4��|6H�� H��xω�Hc�H�t$8H)�H��|-Hc���W������xD� ��x������x���xƃ�H�t$0H�D$89���   Hc�H��    I�D�     I��    H��    I��    I��    �   �@�e 0�H�|$(H��H�T$8L�D$ L�L$L�T$L�\$�g;���@�e 0�L��I��H�l$8H���M;���@�e 0�L��H���;;���@�e 0�H��H���);���@�e 0�L��H��H��L���;���@�e 0�L��H����:��L�\$L�T$L�L$L�D$ H�|$(H�t$0Lc6M���J  M��M��M��L��A��  ��  1�A��H�t$@�  L��H���H��H��1�H����   1�������D ����|��|D ��|����D ��D ��|D ��|D ��|D ��D ��D@��D@��|D@��|D@��|D@��D@��D`��D`��|D`��|D`��|D`��D`H��H��H����R���H��y> f.�     ������D ����|��|D ��|��H�� H��x�Hc�L��H)�H��|KHc���W���D� ��D� ��������x���x���xD� ��xD� ��x���x���˃�A9���   H�H�D�     H��    I��    I�D�     I��    H��    �{�@�e 0�H�|$(H��L����w��8���@�e 0�H�|$@L����8���@�e 0�L��L����8���@�e 0�L��L���8���@�e 0�L��L���8���@�e 0�H��L���8��H�|$(H�D$HLc0M��L�|$pL�d$hL�l$`H�\$XH�l$P��  A��  ��  1�A���)  L��H���H��H��1�H����   1�f�     f.�     f.�     ��������D ����|D ��|��|��D ��D ��D ��|D ��|D ��|D ��D@��D@��D@��|D@��|D@��|D@��D`��D`��D`��|D`��|D`��|D`H��H��H����R���H��y> f.�     ��������D ����|D ��|��|H�� H��x�Hc�L��H)�H��|JHc���W���������D� ��D� ��������xD� ��xD� ��x���x���xσ�A9�~4H�H��    H�D�     H��    I�D�     I��    I��    H��x[A\A]A^A_]��wþ@�e 0�L����w�i6���@�e 0�H��L���W6���@�e 0�H��L���E6���@�e 0�L��L���36���@�e 0�L��L���!6���@�e 0�L��L��H��x[A\A]A^A_]�6���    f.�     HcI�T��I��H�T��H��L�L$��^� ��{YL����{���{YT����{�L�T$��{Y\����{�H�|$��Yd����$�H�L$ ��YD�����H�T$0��**H�T$@��*2��Y�H�T$8��W���*2��^�H�T$(��l����{L����{T����{\����d����D������fD  f.�     f.�     Hc��D����{YL����T����{\��������X����L�T$��{T����kXT�������^���yL����L����D����{T��L�\$��{D����{XD����{D����{D����{XD����{D��H�|$��D����XD����D��I�D��    H�D��    I�D��    H�D��    I�D��    I�D��    H�D��    H�L$ H�D��    �ȉ�f�     UAWAVAUATSLc&E��$8���L�T$PL�\$HH�D$@H�|$8E���  E1�M��I��8���A���   ��  M��I���L��H��1�I���i  1�D  f.�     f.�     ��|�(@  A�(x  ��|(���)@  �)x  ��)���/@  �/x  ��/���*@  �*x  ��*���(@  �(x  ��(��|�)@  A�)x  ��|)��|�+@  A�+x  ��|+��|�*@  A�*x  ��|*��|�(`  ��|D( ���)`  ��D) ���/`  ��D/ ���*`  ��D* ���(`  ��D( ��|�)`  ��|D) ��|�+`  ��|D+ ��|�*`  ��|D* H��@H��H��������H��yx��|�(@  ��|(���)@  ��)���/@  ��/���*@  ��*���(@  ��(��|�)@  ��|)��|�+@  ��|+��|�*@  ��|*M)�I��6�����   āx��@  āx�ġx��@  ġx�ġx��@  ġx�ġx��@  ġx�ġx��@  ġx�āx��@  āx�āx��@  āx�āx��@  āx�I��M9�~`K���@  K��J���@  J��J���@  J��J���@  J��J���@  J��K���@  K��K���@  K��K���@  K��D�>D��Hc�H��H��H������D  ������|A�x  ���x  ���x  ���x  ���x  ��|A�x  ��|A�x  ��|A�x  ��|D ��D ��D ��D ��D ��|D ��|D ��|D H��@H���J���IǄ�@      HǄ�@      HǄ�@      HǄ�@      HǄ�@      IǄ�@      IǄ�@      IǄ�@      [A\A]A^A_]��w��     f.�     f.�     UAWAVAUATSI��L�L$�HcH�|$`H�l$XL�d$PL�T$HL�\$@L�t$8H����  E1����  I��I���L��H��1�I����  L�l$�H�t$�L�|$�H�L$�1�f�     M�$��|D��A��X�����|��   M�,��|D��A��X�����|��   M���|D��A��X�����|��   I�4��D����X�������   M���|D��A��X�����|��   H�L$PL���|D��A��X�����|��   L��L�t ��|D��A��X�����|��   L�<��|D��A��X�����|L����|���  ��|L����|���  ��|L����|���  L�L$���L�������  ��|L����|���  L�T$H��|L����|���  L�\$@H�t$PI����|L����|���  I����|��   ��|D����|���  H���H��H����X���H�L$�L�|$�H�t$�L�l$�H����   I���D������   I���D������   I���D������   I���D������   I���D������   I���D������   H�T ��D������   H�;��D������   H��L)�H����   ��xD����x��  ��xD����x��  ��xD����x��  ��xD����x��  ��xD����x��  ��xD����x��  ��D������  ��D������  I��L9�~nH��L)�I�T��I���  I�T��I���  I�T��I���  I�T��I���  I�T��I���  I�T��I���  H�T��H���  H�T��H���  �@d���{��	��Y���}���}�   ��|��}���|������|$$��}���|��|#��|��e ��}%� ��'�    H������f�f.�     ��|A�x  ��| A� x  ������|,A�x  ��|A�x  ��|,A�x  ��|A�x  ��l �x  ��$�x  ��|D ��|L  ��|l ��|T ��|l ��|\ ��l ��d H��@H���C���[A\A]A^A_]��w�f.�     f.�     �SH��H���UUH�kH�l$H��AWAVAUATH��  H��P���H��@���L��X���L��H���H��`���H��H���H�E�H�E�   H��H���H�E�H�E�   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H ���H�����Hǅ���   H��H����H������Hǅ����	   H��H����H��p���Hǅx���   H��Hl���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H-���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��Hl���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H,���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����	   H��HP���H��@���HǅH����   L�{0��   L����&��1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��O���L��L���'��N��-O���Hǅ ����zD Hǅ(���   H�� ���L��L���_'��M��    L��L���&���̉ ��N���f��� f��L������ ��H���H��� H��@���H��H��@�����   I���&��A��H�{ ��   � &��F�$8L�s(��   L����%��B�D H��8���Lc�M��L��    HH�E��DH�Mc�E9�MO�H��H���H��H)�H���H��h���H��L��L���&��M)��    LH��   L���|%�����    H�Lc�M9�MN�M)�LH�A�   I��MN�M��LH�L��0���M)�LH�H��h���M�<L��H�s(L���&��M�Hǅ����zD Hǅ���   H�����L��L����%����   H�{ ��$�����    I�Lc�M9�MN�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L�0���L��H�s L���%��M�Hǅ ����zD Hǅ���   H�� ���L��L���\%��M�    L��L���$���E�    H��h���H������H��8���H������H������H��0H�D$ H�D$    H�D$    H�D$    H�$    H��H���H��H��<���1�E1�E1���  H��0��l���L�sL�{��tR�Q� ������H�<� H������H��H����H������Hǅ����	   H������H��H��l����	   0��]����� �������� ������H��H����H������Hǅ����   H������H��H��<���H��H��8���A�   L��`���L���:#����l�����t`��� ������f��� f������H�� H������H��H����H������Hǅ����   H������H��H��l����   0�萇��ǅ���dens�E�   H��H���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I��p���E1��\#��H��P��l�����t`��� ��:���f��� f��8���H�m� H��0���H��H0���H������Hǅ����   H������H��H��l����   0��n���ǅP���temp�E�   H��HP���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I��t���E1��:"��H��P��l�����t`��� ��z���f�q� f��x���H�[� H��p���H��Hp���H������Hǅ����   H������H��H��l����   0��L���ǅ����mass�E�   H��H����H��p���Hǅx���   H��p���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I��x���E1��!��H��P��l�����t`�n� ������f�_� f������H�I� H������H��H����H��`���Hǅh���   H��`���H��H��l����   0��*����� �������� �������E�   H��H����H��P���HǅX���   H��P���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I��|���E1�����H��P��l�����t`�N� ������f�?� f������H�)� H������H��H����H��@���HǅH���   H��@���H��H��l����   0�������� �����f�� f������E�   H��H���H��0���Hǅ8���   H��0���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I������E1����H��P��l�����t`�,� ��:���f�� f��8���H�� H��0���H��H0���H�� ���Hǅ(���   H�� ���H��H��l����   0��ȁ��ǅP���year�E�   H��HP���H�����Hǅ���   H�����H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��<���H��H���H��H��8���I��I������E1����H��P��l�����t`�� ��z���f�� f��x���H�� H��p���H��Hp���H�� ���Hǅ���   H�� ���H��H��l����   0�覀��H��H��<���1�1�1�E1��  ��l�����tR�� ������H�� H������H��H����H������Hǅ����	   H������H��H��l����	   0��2���H��@���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic$1�H��HH�H������Hǅ����   H��H�$    H��H��<���H��H��p���H��H�°���1�E1�E1��@  H����l�����t`��~ �����f��~ f�����H��~ H�����H��H���H������Hǅ����   H������H��H��l����   0����H��H���H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Ic$1�H��HH�H��h���Hǅp���   H��H�$    H��H��<���H��H��t���H��H��0���1�E1�E1���  H����l�����t`��} ������f��} f������H�q} H������H��H����H������Hǅ����   H������H��H��l����   0���}��H��P���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic$1�H��HH�H������Hǅ����   H��H�$    H��H��<���H��H��x���H��H�°���1�E1�E1��  H����l�����t`��| �����f�w| f�����H�a| H�����H��H���H������Hǅ����   H������H��H��l����   0���|��H��X���H��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Ic$1�H��HH�H��h���Hǅp���   H��H�$    H��H��<���H��H��|���H��H��0���1�E1�E1��   H����l�����t`�v{ ������f�g{ f������H�Q{ H������H��H����H������Hǅ����   H������H��H��l����   0��{��L������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Ic$1�H��HH�H������Hǅ����   H��H�$    H��H��<���H��H�ƀ���H��H�°���1�E1�E1���  H����l�����t`�mz �����f�^z f�����H�Hz H�����H��H���H������Hǅ����   H������H��H��l����   0��z��L��0���Hǅ8���@   H�q     H��@���H�    H��H���HǅP���    HǅX���    Hǅ`���   Ic$1�H��HH�H��h���Hǅp���   H��H�$    H��H��<���H��H�Ƅ���H��H��0���1�E1�E1��n  H����l�����tU�dy �E�f�Xy f�E�H�Ey H�E�H��H���H������Hǅ����   H������H��H��l����   0��y��H��H��<��������l�����tDH�nf_closeH�E�H��H���H������Hǅ����   H������H��H��l����   0��.y��H��H���A\A]A^A_H��]H��[�f�f.�     f.�     �L$(�    H�|$�H�L$�H�|$�0����H��(1���W���)$�    f.�     0Ƀ<� uH���H��|��t0�H��(�n���D$    �D$    �D$    �D$$   L$L$ L$L$$H��H��H��H�� H��H��I��I��$0��'��H��(�f�UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$�@�e �   0�L���1���C8��$�  H��H�Ǆ  �@�e ��  0����L�AA$�$�  AH��H���  0�L��L��L�����H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���q���C8��$�  H��H�Ǆ  ���e ��  0��L��L�AA$�$�  AH��H���  0�L��L��L������H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L������C8��$�  H��H�Ǆ  ���e ��  0����L�AA$�$�  AH��H���  0�L��L��L���7��H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$� �e �   0�L�������C8��$�  H��H�Ǆ  � �e ��  0�����L�AA$�$�  AH��H���  0�L��L��L���w��H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$�@�e �   0�L���1���C8��$�  H��H�Ǆ  �@�e ��  0����L�AA$�$�  AH��H���  0�L��L��L�����H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���q���C8��$�  H��H�Ǆ  ���e ��  0��L��L�AA$�$�  AH��H���  0�L��L��L������H��H���[A\A^A_]� f.�     f.�     SH��H���UUH�kH�l$H��AWAVAUATH��0  H������H������L������L������H�����H������H��H���H�E�H�E�   H��H���H��p���Hǅx���   H��H`���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����	   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��HM���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��HL���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��HL���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���	   H��H����H������Hǅ�����   L�{P��   L���7
��1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H������L��L����
��N��-����Hǅ�����{D Hǅ����   H������L��L���
��M��    L��L���d	���4n ������f�%n f�������n ������H��m H������I��I��������   L���Y	��A��H�{@��   �H	��F�$8H�{H��   �6	��B�D H������Lc�M��L��    HH�E��DH�E1�Mc�E9�MO�H��H���H��H)�H���H�� ���H��L��L����	��M)�MH��   H�{H�������    H�Lc�M9�MN�M)�LH�A�	   I��	MN�M��LH�L������M)�LH�H�� ���M�<L��H�sHL���`	��M�Hǅ����|D Hǅ����	   H������L��L���5	����   H�{@�7�����    I�Lc�M9�MN�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L�����L��H�s@L������M�Hǅ����!|D Hǅ����   H������L��L�����M��    L��L���e���E�    H�� ���H������H������H������H������H��0H�D$ H�D$    H�D$    H�D$    H�$    H��H�ƬH��H������1�E1�E1��Q  H��0�����L�c(L�k L�sL�{��tY��k ��X���H��k H��P�������H��HP���H������Hǅ����	   H������H��H������	   0��k���ik ��t����Yk ��p���H��Hp���H��p���Hǅx���   H��p���H��H������H��H������A�   H������z���������tg�k ������f��j f������H��j H����������H��H����H��`���Hǅh���   H��`���H��H������   0���j��H�constantH�������E�   H��H����H��P���HǅX���   H��P���H��H������H��H�°H��H������A�   ��  �������tg�Kj ������f�<j f������H�&j H����������H��H����H��@���HǅH���   H��@���H��H������   0���i��ǅ����dens�E�   H��H����H��0���Hǅ8���   H��0���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�´H��H������I��I�����E1�����H��P�������tg�2i �����f�#i f�����H�i H���������H��H���H�� ���Hǅ(���   H�� ���H��H������   0���h��ǅ0���temp�E�   H��H0���H�����Hǅ���   H�����H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�¸H��H������I��I�����E1����H��P�������tg�h ��Z���f�
h f��X���H��g H��P�������H��HP���H�� ���Hǅ���   H�� ���H��H������   0��g��ǅp���mass�E�   H��Hp���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�¼H��H������I��I�����E1��r��H��P�������tg� g ������f��f f������H��f H����������H��H����H������Hǅ����   H������H��H������   0��}f����f ��������f �������E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I�����E1��;��H��P�������tg��e ������f��e f������H��e H����������H��H����H������Hǅ����   H������H��H������   0��Fe���ue ������f�fe f�������E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I�� ���E1����H��P�������tg��d �����f��d f�����H��d H���������H��H���H������Hǅ����   H������H��H������   0��d��ǅ0���year�E�   H��H0���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��$���E1������H��P�������tg��c ��Z���f��c f��X���H�rc H��P�������H��HP���H������Hǅ����   H������H��H������   0���b��H�refreezeH��p����E�   H��Hp���H��p���Hǅx���   H��p���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��(���E1�����H��P�������tg�wb ������f�hb f������H�Rb H����������H��H����H��`���Hǅh���   H��`���H��H������   0��a��fǅ����dz�E�   H��H����H��P���HǅX���   H��P���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��,���E1�����H��P�������tg�_a ������f�Pa f������H�:a H����������H��H����H��@���HǅH���   H��@���H��H������   0��`��ǅ����drho�E�   H��H����H��0���Hǅ8���   H��0���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��0���E1��X���H��P�������tg�F` �����f�7` f�����H�!` H���������H��H���H�� ���Hǅ(���   H�� ���H��H������   0��c_����_ ��2���f��_ f��0����E�   H��H0���H�����Hǅ���   H�����H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��4���E1�����H��P�������tY�_ ��X���H�_ H��P�������H��HP���H�� ���Hǅ���   H�� ���H��H������   0��8^��ǅp���time�E�   H��Hp���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I��8���E1�����H��P�������tY�^ ������H��] H����������H��H����H������Hǅ����   H������H��H������   0��]��H��H������1�1�1�E1��  �������tY��] ������H��] H����������H��H����H������Hǅ����	   H������H��H������	   0��\��H������H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   H�����Hc 1�H��HH�H�����Hǅ���   H��H�$    H��H������H��H�����H��H������1�E1�E1��*  H���������tg��\ ��:���f��\ f��8���H�s\ H��0�������H��H0���H������Hǅ����   H������H��H������   0��u[��L��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   L�����Ic1�H��HH�H������Hǅ����   H��H�$    H��H������H��H�����H��H��P���1�E1�E1���  H���������tg��[ ������f�s[ f������H�][ H����������H��H����H������Hǅ����   H������H��H������   0��OZ��H������H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   H��H�$    H��H������H��H�����H��H������1�E1�E1��^  H���������tg�lZ ��:���f�]Z f��8���H�GZ H��0�������H��H0���H������Hǅ����   H������H��H������   0��)Y��L��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   Ic1�H��HH�H������Hǅ����   H��H�$    H��H������H��H�����H��H��P���1�E1�E1���
  H���������tg�]Y ������f�NY f������H�8Y H����������H��H����H������Hǅ����   H������H��H������   0��
X��L������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   H��H�$    H��H������H��H�� ���H��H������1�E1�E1��
  H���������tg�NX ��:���f�?X f��8���H�)X H��0�������H��H0���H������Hǅ����   H������H��H������   0���V��L��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   Ic1�H��HH�H������Hǅ����   H��H�$    H��H������H��H��$���H��H��P���1�E1�E1��A
  H���������tg�?W ������f�0W f������H�W H����������H��H����H��p���Hǅx���   H��p���H��H������   0���U��H�C0H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   H��H�$    H��H������H��H��(���H��H������1�E1�E1���	  H���������tg�,V ��:���f�V f��8���H�V H��0�������H��H0���H��`���Hǅh���   H��`���H��H������   0��T��H�C8H��P���HǅX���@   H�q     H��`���H�    H��h���Hǅp���    Hǅx���    Hǅ����   Ic1�H��HH�H������Hǅ����   H��H�$    H��H������H��H��,���H��H��P���1�E1�E1��{	  H���������tg�U ������f�
U f������H��T H����������H��H����H��P���HǅX���   H��P���H��H������   0��S��H������H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ ���   Ic1�H��HH�H�����Hǅ���   H��H�$    H��H������H��H��0���H��H������1�E1�E1��	  H���������tg�T ��:���f��S f��8���H��S H��0�������H��H0���H��@���HǅH���   H��@���H��H������   0��`R��H��H������H��H��4���1�H�������>	  �������tY�~S ��X���H�iS H��P�������H��HP���H��0���Hǅ8���   H��0���H��H������   0���Q��H��H������H��H��8���1�H�������9	  �������tY�	S ��x���H��R H��p�������H��Hp���H�� ���Hǅ(���   H�� ���H��H������   0��VQ��H��H������������������tKH�nf_closeH�E�����H��H���H�����Hǅ���   H�����H��H������   0���P��H��H���A\A]A^A_H��]H��[�fD  f.�     f.�     �L$(�    H�|$�H�L$�H�|$�0��`���H�constantH��   H�t$�H�D$�   H�t$�A�   0��=��� f.�     H��(1���W���)$�    f.�     0Ƀ<� uH���H��|��t0�H��(������D$    �D$    �D$    �D$$   L$L$ L$L$$H��H��H��H�� H��H��I��I��$0�����H��(�f�UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$� �e �   0�L�������C8��$�  H��H�Ǆ  � �e ��  0�����L�AA$�$�  AH��H���  0�L��L��L���7���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$�@�e �   0�L��������C8��$�  H��H�Ǆ  �@�e ��  0������L�AA$�$�  AH��H���  0�L��L��L���w���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���1����C8��$�  H��H�Ǆ  ���e ��  0�����L�AA$�$�  AH��H���  0�L��L��L������H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���q����C8��$�  H��H�Ǆ  ���e ��  0��L���L�AA$�$�  AH��H���  0�L��L��L�������H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$� �e �   0�L�������C8��$�  H��H�Ǆ  � �e ��  0�����L�AA$�$�  AH��H���  0�L��L��L���7���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$�@�e �   0�L��������C8��$�  H��H�Ǆ  �@�e ��  0������L�AA$�$�  AH��H���  0�L��L��L���w���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���1����C8��$�  H��H�Ǆ  ���e ��  0�����L�AA$�$�  AH��H���  0�L��L��L������H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L���q����C8��$�  H��H�Ǆ  ���e ��  0��L���L�AA$�$�  AH��H���  0�L��L��L�������H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$� �e �   0�L�������C8��$�  H��H�Ǆ  � �e ��  0�����L�AA$�$�  AH��H���  0�L��L��L���7���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H��  I��I��H��L�$$�@�e �   0�L�������A���$�  A$�$�  H��H���  0�H��L��L�������H��H���[A\A^A_]�f�UH��AWAVATSH���H��  I��I��H��L�$$���e �   0�L���y���A���$�  A$�$�  H��H���  0�H��L��L���R���H��H���[A\A^A_]�f�UAWAVAUATSI��M��L�D$�I�ҋ��>H���{��Z�H��$�   ��D��H��$�   Hc1H)�H�p��{��Z������   H�,vH�(��{$��Z������  H����{ ��Z�����T  H�<�H�8H�\$H����Z�����  H�hH�\$P����Z������  Hk�H�L�l$8��{E ��Z������  H��L�\$@��{��Z�����t  L�L$X��{��Z�H��H�����<  L�D$`��{ ��Z�H�x����  H�L$h��Hk�H���Z������  H�L$p����Z�H�������  L�t$x��{��Z�Hk�H�����\	  H��$�   ����Z�Lk�I�,����$
  H�<H�<L��$�   ��{��Z������
  Hk�H��H��$�   ��E ��Z�H�4�����  H��$�   ��H�4��Z�����|  H��$�   ����Z�H�X����D  I�    I�$    I�E     I�    H�D$HH�     H�D$�H�     H�D$PH�     I�    I�     I�    H�    I�    H�E     H�    [A\A]A^A_]�D  f.�     SH��H���UUH�kH�l$H��AWAVAUATH��   E���A��L��p����H�U�H�s�6L�shL�C`L�{XL�KPL�[H9�L�k8�E�   �U�~H�U�)��U��u��E�Hc?I��I���H����   ��*���A ��^���W���r*���^���Z�1҃�|yL��H����}�1�I��|M1�@ f.�     f.�     ��uYT A��  ��|T ��uYT ��|T H��@H��H���|�L��H��y��uYL ��|L L��H��H)�H��|�����qYL� ��xL� H��H9�~��{YD� ��{D� L��h���H��x���D�U�H�M�D)�H�M��AH�1�H��HI�H��   H���H��H)�H��H�u�H���Hc�H�E�H��H���H�M�I��I�������M�L��H���H�E�Ic�H�E�H�ԅ��O  1�H�M�����   H�E�H��H��H��H�{H�� ��   @ f.�     Ic���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ǘ   ��}���)F@��Z�Ǹ   ��Z���   ��}���)F`H��H��   H�� H���|�H�E�H��yF�    f.�     f.�     Ic���ZD����ZL���}���H�� H��@H��x�H�E�H9E�|Ic�H�H�s��ZD����)�H��H9E�|'H�M�H�sH��@ ��D����Z����H��L9�|�L�e�D��H��p���Lc I��2�}� ��   LcU��E�Lc�1�����   L��H��H��=L�H���1�������   C�D��M�D�C�D��M�D�C�D��M�D�C�D��M�D�C�D��M�D�C�D��M�D�C�D��M�D�C�D��M�H��H9�|��@ f.�     ��C�D��M�H��L9�|�H�E�H��H��H���H��H)�H���H�ă}� �>  1�H�u�����   H�U�I��I��H��H�{ H�� |}�Hc���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ט   ��}���)F@��Z�׸   ��Z���   ��}���)F`H��H��   I�� I���|�H�U�M��yF�    f.�     f.�     Hc���ZD����ZL���}���H�� H��@I��x�H�U�H9U�|Hc�H�H�{ ��ZD����)�H��H9U�|'H�u�H�{ H�4�H�}���D����Z����H��H9�|�D�]�E����   HcU�Mc�1�A����   E��M��I��I��=M�I���1��     f.�     �4����   A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�H��L9�|�E����f.�     �4�A�t��L�H��L9�|�H�E�H��H��H���H��H)�H���H�ă}� H�{8A��L�m��S  1�H�M�����   H�M�I��I��H��H�� }	H�K(�   H�K( f.�     Ic���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ј   ��}���)F@��Z�Ѹ   ��Z���   ��}���)F`H��H��   I�� I���|�H�U�M��yF�    f.�     f.�     Ic���ZD����ZL���}���H�� H��@I��x�H�U�H9U�|Ic�H�H�K(��ZD����)�H��H9U�|'H�M�H�s(H�4�@ ��D����Z����H��L9�|�E����   HcU�Mc�1�A����   E��M��I��I��=M�I���1�f�f.�     f.�     �4����   A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�H��L9�|�E��L�m���    �4�A�t��L�H��L9�|�H�E�H��H��H���H��H)�H���H�ă}� �^  1�H�M����  H�M�I��I��H��H�� }	H�K0�   H�K0@ f.�     f.�     Ic���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ј   ��}���)F@��Z�Ѹ   ��Z���   ��}���)F`H��H��   I�� I���|�H�U�M��yF�    f.�     f.�     Ic���ZD����ZL���}���H�� H��@I��x�H�U�H9U�|Ic�H�H�K0��ZD����)�H��H9U�|'H�M�H�s0H�4�@ ��D����Z����H��L9�|�E����   HcU�Mc�1�A����   E��M��I��I��=M�I���1�f�f.�     f.�     �4����   A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�H��L9�|�E����f.�     �4�A�t��L�H��L9�|�H�E�H��H��H���H��H)�H���H�ă}� �^  1�H�M����  H�M�I��I��H��H�� }H���   H��fD  f.�     f.�     Ic���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ј   ��}���)F@��Z�Ѹ   ��Z���   ��}���)F`H��H��   I�� I���|�H�U�M��yF�    f.�     f.�     Ic���ZD����ZL���}���H�� H��@I��x�H�U�H9U�|Ic�H���ZD����)�H��H9U�|+H�M�H�4�f�f.�     ��D����Z����H��L9�|�E��E����   HcU�Mc�1�A����   M��I��I��=M�I���1�f�f.�     f.�     �4����   A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�H��L9�|��@ f.�     �4�A�t��L�H��L9�|�H�M�H��H���H��H)�H���H�ă}� L�Sp�]  1�H�M����  H�M�I��I��H�� I��H��H�K@��   D  f.�     f.�     Ic���ZD����ZL���}���)��ZD�8��ZL�X��}���)F ��ZD�x��Z�ј   ��}���)F@��Z�Ѹ   ��Z���   ��}���)F`H��H��   I�� I���|�M��yF�f.�     f.�     f.�     Ic���ZD����ZL���}���H�� H��@I��x�L��H9U�|Ic�H�H�K@��ZD����)�H��H9U�|(H�K@H�u�H�4�D  ��D����Z����H��L9�|�E����   HcU�Mc�1�A����   M��I��I��=M�I���1�D  f.�     f.�     �4����   A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�t�A�t��L�H��L9�|��@ f.�     �4�A�t��L�H��L9�|�H��x�������   ��  ��   1����|   H��h���H��H��1�H��|;1�fD  f.�     ��������D ��D@��D`H��H��H���|�H��y%�     f.�     ������H�� H��x��Hc�H��H)�H��|Hc���W���σ�9�~H�H��    ����e 0���w����H��H���A\A]A^A_H��]H��[��w�f�     f.�     SH��H���UUH�kH�l$H��AWAVAUATH��   H�U�L�M�L�E�H�u�Lc1�M��L��HH�LcM��II�H��H��H�T�H���H��H)�H���H�u�L�4�H�H��H��N�$2H��L�,�H�IM�<�H�{0L�K(M����  A��  �\  1�A��L�K(��   L��H���H��H��1�H����   1�f.�     ������|��|��|��|D ��|D ��|D ��|D ��|D ��|D@��|D@��|D@��|D@��|D`��|D`��|D`��|D`H��H��H���|�H��yC�     f.�     f.�     ������|��|��|��|D H�� H��x�Hc�L��H)�H��|#Hc���W���x���x���x���xD� ��A9���   H�I��    I��    I��    I�D�     �c���e 0�L��L��L�E�L�U���������e 0�L��H�U��������e 0�L��H�U��������e 0�L��H�U�����L�U�H�{0L�E�L�K(L�E�L��H���H�E�E���<  H�E�=  }31�H�Mȃ�H�{0L�K(��   H�U�H��H��1�H��}5H�U�H�s8�|0�H�}�H�s8H�U�L�U���w�6���L�U�H�{0L�K(��   1�H�U�H�s8�f.�     ����)
��D ��)D
 ��D@��)D
@��D`��)D
`H��H��H���|�H��y3�f.�     f.�     f.�     ����
H�� H��x�H�U�H�E�Hc�H�U�H)�H��|Hc�H�S8���H�U���ʃ�H�M�9�~H�H�K8H��H�U�H��H�E�D� H�E���0E��L�]���  ��f) ��^Ƹ   ��W�D�E� f.�     ��W���*���Y�Ic���{T����X���.��|   ��Y�Ic�H�K��d��Hc����d����{d��H�K ��d�����d����{d��H�K8��^T����{\�����\����{\��A���   �    f.�     f.�     ��.�uz�v��\���Y�Ic�H�K��l��Hc���٩l����{l��H�K ��l����٩l����{l��H�K8��^d����{l����٩l����{l����\���{T����X��k�     ��Y�Ic�H�K��d��Hc����d����{d��H�K ��d�����d����{d��L��L�K8��k^T��I����{\�����\����{\��A������(�D9��n���D�E���W�H�Eȅ�H�s8�  H�E�=�   ��   H�E�H��H���H���H)�H�E�H��1������������������     f.�     f.�     ��XT ��XL@��XD`��$��  ��  ��X�H��H��H9�|�H��y ��XH�� H��x���X���X���X���}���X���y���X�H�E�H;E��1  H�E��f.�     f.�     f.�     ��X�H��H;E�|���   1���W�H�Mȃ���   H�U�H��H��1�����H��|i1�H�U�H�������f�f.�     f.�     ��X��L ��XL@��X�����   ��X��   ��XL`��X�H���   H��H���H���H��y%f�f.�     f.�     ��XH�� H��x���}���X���y���X���W���X�H�E�H9E�~�f.�     ��X�H��H;E�|���+ ��.��  H�Eȅ��  H�E�=  ��   1�H�Mȃ���   H�U�H��H��1�H��|D1�f�����|)��D ��|)D ��D@��|)D@��D`��|)D`H��H��H���|�H��y0�    f.�     f.�     ����|H�� H��x�H�E�Hc�H�U�H)�H��|Hc������x˃�H�M�9�~CH�H��I���70�L��H�U�L�U�D�E���u���w�D�����u�D�E�H�s8L�]�L�U�H�{0E����   �   ��W�f�     Ic���{T����X���*���Y���.�sDIc���^L����T��Hc����T����{T��A���   �     f.�     f.�     ��.�uz�6��\�Ic���^\����d��Hc����d����{d����\���{T����X��$�Ic���^L����T��Hc����T����{T��A������(�D9��.���H�E�� �H�M��9E���  L�[@M��I��I��=M�I���H�KHc	H��21�Hc�A����   1�fD  ��{���Z�A���   ��zD��H���{D���Z���zD��H���{D���Z���zD��H���{D���Z���zD����{D� ��Z�H���zD��H���{D�(��Z���zD��H���{D�0��Z���zD��H���{D�8��Z���zD��H�H��L9��C�������{���Z���zD��H�H��L9�|�1�Hc�A��}	L�cH��   1�L�cHf.�     ��{���Z�A���   ��zD��H���{D���Z���zD��H���{D���Z���zD��H���{D���Z���zD����{D� ��Z�H���zD��H���{D�(��Z���zD��H���{D�0��Z���zD��H���{D�8��Z���zD��H�H��L9��C�������{���Z���zD��H�H��L9�|�1�Hc�A��}	L�{P��   1�L�{Pf.�     ��{���Z�A���   ��zD��H���{D���Z���zD��H���{D���Z���zD��H���{D���Z���zD����{D� ��Z�H���zD��H���{D�(��Z���zD��H���{D�0��Z���zD��H���{D�8��Z���zD��H�H��L9��C�������{���Z���zD��H�H��L9�|�1�H�A����   1�H�SXfD  f.�     ��{D� ��Z�A���   ��D��H���{D���Z���D��H���{D���Z���D��H���{D���Z���D����{D� ��Z�H���D��H���{D�(��Z���D��H���{D�0��Z���D��H���{D�8��Z���D��H�H��L9��J���L9�}!H�SX���{D� ��Z���D��H�H��L9�|�H�Eȅ���   H�E�=  ��   1�H�Mȃ�|sH�U�H��H��1�H��|51�f.�     ��������D ��D@��D`H��H��H���|�H��y%�     f.�     ������H�� H��x��Hc�H�U�H)�H��|Hc���W���σ�H�M�9�~H�H��    ����e 0�H�U���w����H��H���A\A]A^A_H��]H��[��w�f.�     SH��H���UUH�kH�l$H��AWAVAUATH��D  I��H������L��p���H������H��h���H��H@���H��0���Hǅ8���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����	   H��Hp���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H ���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��H0���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hp���H��`���Hǅh���   H��H@���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��HP���H��@���HǅH���   H��H0���H�� ���Hǅ(���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hp���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H ���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��Hl���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H+���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HH���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H	���H������Hǅ����   H��H����H������Hǅ����   H��Hj���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H,���H�����Hǅ���   H��H ���H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HK���H��0���Hǅ8���   H��H ���H�����Hǅ���   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��Hl���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H����H������Hǅ�����   ��   L����w�'���1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H������L��L���Ƕ��N��-����Hǅ�����}D Hǅ����   H������L��L��藶��M��    L��L���T���� ������f�� f�������� ������H�� H������������)�������)�������)�������)������)�0�����)�P����    ǅp���    ǅt���    ��)�����H������@ ������)�������)�����H��@H��x�ǅ����    ǅ����    H��H�ǰ���� �e ��  0���w辴��H��H�ǀ�����   I��藴��A�Ǿ�   H������胴��F�48��   L������L���k���B�D0H��`���Lc�1�M��L��HH�E��DH�Mc�E9�MO�H��H���H��H)�H���H��x���H��L��L������M)��    LH��   L����������    H�Lc�M9�MN�M)�LH�   I��M��LO�M��LH�L��X���M)�LH�H��x���M�<L��H������L��蒴��M�Hǅ�����}D Hǅ����   H������L��L���g�����   H�������f������    I�Lc�M9�MN�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L�X���L��H������L�������M�Hǅ�����}D Hǅ����   H������L��L���Գ��M�    L��L��葲��ǅL���    H��x���H��p���H��`���H��x���H��p���H��0H�D$ H�D$    H�D$    H�D$    H�$    H��H��L���I��I��,���1�E1�E1�L���V  H��0��l���ǅp���timeL��h���A���2��l���H��Hp���H��`���Hǅh���   H��`���H��H��l���H��H��,���A�   L���	�����������tj�c ������f�T f������H�> H�����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w�U����,���������fǅ����zsH��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   ǅP���   H��H����H��@���HǅH���   H��@���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��P���H��H������I��I��,���E1���w蹰��H��P��������L��p���tj�� ��:���f�� f��8���H�� H��0��������H��H0���H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�����,�����L���ǅP���viceH��HL���H��p���Hǅx���    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   ǅT���   H��HP���H�� ���Hǅ(���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��T���H��H��p���I��I������E1���w�����H��P��������tj� ������f� f������H�� H�����������H��H����H�����Hǅ���   H�����H��H�Ǭ����   0���w������,���������ǅ����vaccH��H����H�����Hǅ���    H�q     H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   ǅX���   H��H����H�� ���Hǅ���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��X���H��H�����I��I������E1���w�H���H��P��������tj�n ��z���f�_ f��x���H�I H��p��������H��Hp���H������Hǅ����   H������H��H�Ǭ����   0���w�0����,����������� ������f�� f������H��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   ǅ\���   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��\���H��H������I��I������E1���w胫��H��P��������tj�� �����f�� f�����H�� H����������H��H���H������Hǅ����   H������H��H�Ǭ����   0���w�k����,�����,����H ��4����8 ��0���H��H,���H��P���HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   ǅ`���   H��H0���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��`���H��H��P���I��I��L���E1���w�����H��P��������tj� ������f�� f������H�� H�����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�����,����������� �������� ������H��H����H������Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   ǅd���   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��d���H��H������I��I�����E1���w�����H��P��������tj�S ��Z���f�D f��X���H�. H��P��������H��HP���H������Hǅ����   H������H��H�Ǭ����   0���w������,�����l���ǅp���vsubH��Hl���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   ǅh���   H��Hp���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��h���H��H������I��I������E1���w�H���H��P��������tj�� ������f�� f������H�� H�����������H��H����H��p���Hǅx���   H��p���H��H�Ǭ����   0���w�0
����,��������ǅ���vsndH��H���H��0���Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   ǅl���   H��H���H��`���Hǅh���   H��`���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��l���H��H��0���I��I������E1���w蓤��H��P��������tj�	 ������f��
 f������H��
 H�����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w�{����,���������f��
 f��������
 ������H��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   ǅp���   H��H����H��@���HǅH���   H��@���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��p���H��H������I��I��l���E1���w�΢��H��P��������tj�\	 ��:���f�M	 f��8���H�7	 H��0��������H��H0���H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�����,�����L���f�� f��T����� ��P���H��HL���H��p���Hǅx���    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   ǅt���   H��HP���H�� ���Hǅ(���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��t���H��H��p���I��I��4���E1���w�	���H��P��������t\�� ������H�� H�����������H��H����H�����Hǅ���   H�����H��H�Ǭ����   0���w������,����������O ������f�@ f�������/ ������H��H����H�����Hǅ���    H�q     H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   ǅx���   H��H����H�� ���Hǅ���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��x���H��H�����I��I������E1���w�F���H��P��������t\� ��x���H�� H��p��������H��Hp���H�����Hǅ����   H�����H��H�Ǭ����   0���w�<����,���������f�� f�������� ������H��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   ǅ|���   H��H����H�����Hǅ���   H�����H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H��|���H��H������I��I������E1���w菝��H��P��������t\�c �����H�N H����������H��H���H��п��Hǅؿ��   H��п��H��H�Ǭ����   0���w�����,�����,���H�refreezeH��0���H��H,���H��P���HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   �E�   H��H0���H������Hǅȿ��   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H��P���I��I������E1���w����H��P��������t\�� ������H�� H�����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�������,���������ǅ����rainH��H����H������Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   �E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H������I��I��T���E1���w�F���H��P��������t\�: ��X���H�% H��P��������H��HP���H������Hǅ����   H������H��H�Ǭ����   0���w�<�����,�����l���H�surfmeltH��p���H��Hl���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��Hp���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H������I��I�����E1���w螘��H��P��������t\���  ������H���  H�����������H��H����H��p���Hǅx���   H��p���H��H�Ǭ����   0���w������,���������B�  ������2�  �����H��H���H��0���Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   �E�   H��H���H��`���Hǅh���   H��`���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H��0���I��I������E1���w����H��P��������t\��  ������H���  H�����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w�������,������������  ������f���  f���������  ������H��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   �E�   H��H����H��@���HǅH���   H��@���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H������I��I������E1���w�2���H��P��������t\�f�  ��8���H�Q�  H��0��������H��H0���H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�(�����,�����L���ǅP���Rho0H��HL���H��p���Hǅx���    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��HP���H�� ���Hǅ(���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��,���H��H�H��H��p���I��I��t���E1���w葓��H��P��������t\���  ������H���  H�����������H��H����H�����Hǅ���   H�����H��H�Ǭ����   0���w�������  �������}�  ������H�h�  H�������E�  �|H��H����H�� ���Hǅ���   H�� ���H��H��,���H��H��,���H��H���A�   ��w�?�����������tj�#�  ��$�����  �� �����(��  ��)���������H��H���H�����Hǅ����   H�����H��H�Ǭ����   0���w�������  ��L������  ��H���H���  H��@����E�  �|H��H@���H�����Hǅ���   H�����H��H��,���H��H������H��H���A�   ��w�S�����������tj�g�  ��t����W�  ��p�����(9�  ��)�`��������H��H`���H��о��Hǅؾ��   H��о��H��H�Ǭ����   0���w������  ���������  ������H���  H�������E�  �|H��H����H������HǅȾ��   H������H��H��,���H��H�Ƽ���H��H���A�   ��w�g�����������tj���  ���������  ��������(}�  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w������Q�  �������A�  ������H�,�  H�������E�  �|H��H����H������Hǅ����   H������H��H��,���H��H�Ƅ���H��H���A�   ��w�{�����������tj���  ��������  �������(��  ��)� ��������H��H ���H������Hǅ����   H������H��H�Ǭ����   0���w��������  ��<������  ��8���H�p�  H��0����E�  �|H��H0���H������Hǅ����   H������H��H��,���H��H��L���H��H���A�   ��w菍����������tj�3�  ��d����#�  ��`�����(�  ��)�P��������H��HP���H��p���Hǅx���   H��p���H��H�Ǭ����   0���w��������  ���������  ������H���  H�������E�  �|H��H����H��`���Hǅh���   H��`���H��H��,���H��H�����H��H���A�   ��w裌����������tj�w�  �������g�  ��������(I�  ��)����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w�������  ��������  ������H���  H�������E�  �|H��H����H��@���HǅH���   H��@���H��H��,���H��H������H��H���A�   ��w跋����������tj���  ��������  �� �����(��  ��)����������H��H����H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�����a�  ��,����Q�  ��(���H�<�  H�� ����E�  �|H��H ���H�� ���Hǅ(���   H�� ���H��H��,���H��H�Ƥ���H��H���A�   ��w�ˊ����������tj���  ��T������  ��P�����(��  ��)�@��������H��H@���H�����Hǅ���   H�����H��H�Ǭ����   0���w�'������  ��|������  ��x���H���  H��p����E�  �|H��Hp���H�� ���Hǅ���   H�� ���H��H��,���H��H��l���H��H���A�   ��w�߉����������tj�C�  �������3�  ��������(�  ��)����������H��H����H�����Hǅ����   H�����H��H�Ǭ����   0���w�;������  ���������  ������H���  H�������E�  �|H��H����H�����Hǅ���   H�����H��H��,���H��H��4���H��H���A�   ��w������������tlf���  f�������u�  ��������(W�  ��)����������H��H����H��н��Hǅؽ��   H��н��H��H�Ǭ����   0���w�M����+�  �������  �����H��  H������E�  �|H��H���H������HǅȽ��   H������H��H��,���H��H������H��H���A�   ��w������������tlf���  f��D������  ��@�����(��  ��)�0��������H��H0���H������Hǅ����   H������H��H�Ǭ����   0���w�_����m�  ��l����]�  ��h���H�H�  H��`����E�  �|H��H`���H������Hǅ����   H������H��H��,���H��H������H��H���A�   ��w������������tlf�
�  f���������  ��������(��  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�q������  ���������  ������H���  H�������E�  �|H��H����H������Hǅ����   H������H��H��,���H��H�ƌ���H��H���A�   ��w�)�����������tlf�L�  f�������;�  ��������(�  ��)����������H��H����H��p���Hǅx���   H��p���H��H�Ǭ����   0���w�������  ��������  �����H���  H�� ����E�  �|H��H ���H��`���Hǅh���   H��`���H��H��,���H��H��T���H��H���A�   ��w�;�����������tlf���  f��4����}�  ��0�����(_�  ��)� ��������H��H ���H��P���HǅX���   H��P���H��H�Ǭ����   0���w�����3�  ��\����#�  ��X���H��  H��P����E�  �|H��HP���H��@���HǅH���   H��@���H��H��,���H��H�����H��H���A�   ��w�M�����������tlf���  f���������  ��������(��  ��)�p��������H��Hp���H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�����u�  �������e�  ������H�P�  H�������E�  �|H��H����H�� ���Hǅ(���   H�� ���H��H��,���H��H������H��H���A�   ��w�_�����������tlf��  f��������  ��������(��  ��)����������H��H����H�����Hǅ���   H�����H��H�Ǭ����   0���w�������  ���������  ������H���  H�������E�  �|H��H����H�� ���Hǅ���   H�� ���H��H��,���H��H�Ƭ���H��H���A�   ��w�q�����������tlf�T�  f��$����C�  �� �����(%�  ��)���������H��H���H�����Hǅ����   H�����H��H�Ǭ����   0���w��������  ��L������  ��H���H���  H��@����E�  �|H��H@���H�����Hǅ���   H�����H��H��,���H��H��t���H��H���A�   ��w胁����������tlf���  f��t������  ��p�����(g�  ��)�`��������H��H`���H��м��Hǅؼ��   H��м��H��H�Ǭ����   0���w�����H��H��,���1�1�1�E1���w�b&  ��������t\��  ������H���  H�����������H��H����H������Hǅȼ��	   H������H��H�Ǭ����	   0���w�\���A���2������L������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   A���2Hc�1�H��HH�H������Hǅ����   Hǅ��� �e Hǅ���    H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H������H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��,���H��H�°���H��H�����I��I��p���E1���w�%  H����������t^�7�  ��������(�  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�_���A���2������Ic�H2Hc�1�H��HH�H� H)�I����   H������Hǅ����    H�q     H�� ���H�     H�����Hǅ���    Hǅ���    Hǅ ���   H��(���Hǅ0���   HǅP���@�e HǅX���    H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H������H��H������H��H��P���I��I������E1���w��#  H����������tj�:�  �����f�+�  f�����H��  H����������H��H���H������Hǅ����   H������H��H�Ǭ����   0���w�D���A���2��,���Ic�H2Hc�1�H��HH�H�@H)�I����  H��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   H��h���Hǅp���   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��,���H��H�Ƽ���H��H��0���H��H������I��I������E1���w�"  H����������t^�=�  ��P�����(�  ��)�@��������H��H@���H������Hǅ����   H������H��H�Ǭ����   0���w�5���A���2��l���Ic�H2Hc�1�H��HH�H��    H)�I���X  H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H��l���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��,���H��H�Ƅ���H��H��p���H��H������I��I��0���E1���w�n!  H����������t^�J�  ��������(,�  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�"���A���2������Ic�H2Hc�1�H��HH�H��H)�I���   H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ��� �e Hǅ���    H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�¬���H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��L���H��H�°���H��H�����I��I��p���E1���w�?   H����������t^�[�  ��������(=�  ��)����������H��H����H��p���Hǅx���   H��p���H��H�Ǭ����   0���w����A���2������Ic�H2Hc�1�H��HH�H� H�RH)�I����  H������Hǅ����    H�q     H�� ���H�     H�����Hǅ���    Hǅ���    Hǅ ���   H��(���Hǅ0���   HǅP���@�e HǅX���    H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H�����H��H������H��H��P���I��I������E1���w�  H����������t^�h�  �������(J�  ��)� ��������H��H ���H��`���Hǅh���   H��`���H��H�Ǭ����   0���w� ���A���2��,���Ic�H2Hc�1�H��HH�Hk�H)�I����  H��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   H��h���Hǅp���   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��,���H��H������H��H��0���H��H������I��I������E1���w��  H����������t^�y�  ��P�����([�  ��)�@��������H��H@���H��P���HǅX���   H��P���H��H�Ǭ����   0���w�����A���2��l���Ic�H2Hc�1�H��HH�H��    H)�I���x  H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H��l���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��,���H��H�Ƥ���H��H��p���H��H������I��I��0���E1���w�  H����������t^���  ��������(h�  ��)����������H��H����H��@���HǅH���   H��@���H��H�Ǭ����   0���w�����A���2������Ic�H2Hc�1�H��HH�H��H)�I���@  H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ��� �e Hǅ���    H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�¬���H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��l���H��H�°���H��H�����I��I��p���E1���w�{  H����������t^���  ��������(y�  ��)����������H��H����H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�����A���2������Ic�H2Hc�1�H��HH�H� H��H)�I���  H������Hǅ����    H�q     H�� ���H�     H�����Hǅ���    Hǅ���    Hǅ ���   H��(���Hǅ0���   HǅP���@�e HǅX���    H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��4���H��H������H��H��P���I��I������E1���w�H  H����������t`f���  f�������(��  ��)� ��������H��H ���H�� ���Hǅ(���   H�� ���H��H�Ǭ����   0���w����A���2��,���Ic�H2Hc�1�H��HH�Hk�H)�I����  H��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   H��h���Hǅp���   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��,���H��H������H��H��0���H��H������I��I������E1���w�  H����������t`f���  f��P�����(��  ��)�@��������H��H@���H�����Hǅ���   H�����H��H�Ǭ����   0���w����A���2��l���Ic�H2Hc�1�H��HH�H��    H�RH)�I����  H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H��l���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��,���H��H������H��H��p���H��H������I��I��0���E1���w��  H����������t`f���  f��������(��  ��)����������H��H����H�� ���Hǅ���   H�� ���H��H�Ǭ����   0���w����A���2������Ic�H2Hc�1�H��HH�Hk�H)�I���`	  H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ��� �e Hǅ���    H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�¬���H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H�ƌ���H��H�°���H��H�����I��I��p���E1���w�  H����������t`f���  f��������(��  ��)����������H��H����H�����Hǅ����   H�����H��H�Ǭ����   0���w����A���2������Ic�H2Hc�1�H��HH�Hk�H)�I���(
  H������Hǅ����    H�q     H�� ���H�     H�����Hǅ���    Hǅ���    Hǅ ���   H��(���Hǅ0���   HǅP���@�e HǅX���    H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��T���H��H������H��H��P���I��I������E1���w�|  H����������t`f���  f�������(��  ��)� ��������H��H ���H�����Hǅ���   H�����H��H�Ǭ����   0���w�n���A���2��,���Ic�H2Hc�1�H��HH�H��H�RH)�I����
  H��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   H��h���Hǅp���   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H��,���H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��,���H��H�����H��H��0���H��H������I��I������E1���w�G  H����������t`f���  f��P�����(��  ��)�@��������H��H@���H��л��Hǅػ��   H��л��H��H�Ǭ����   0���w�Y���A���2��l���Ic�H2Hc�1�H��HH�H��H��H)�I����  H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ������e Hǅ����    H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H��l���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H��,���H��H������H��H��p���H��H������I��I��0���E1���w�  H����������t`f���  f��������(��  ��)����������H��H����H������HǅȻ��   H������H��H�Ǭ����   0���w�E���A���2������Ic�H2Hc�1�H��HH�Hk�H)�I����  H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   H������Hǅ����   Hǅ��� �e Hǅ���    H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�¬���H��p���Hǅx���    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H�Ƭ���H��H�°���H��H�����I��I��p���E1���w��  H����������t`f���  f��������(��  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�4���A���2������Ic�H2Hc�1�H��HH�H� H��H)�I���H  H������Hǅ����    H�q     H�� ���H�     H�����Hǅ���    Hǅ���    Hǅ ���   H��(���Hǅ0���   HǅP���@�e HǅX���    H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��,���H��H��t���H��H������H��H��P���I��I������E1���w�  H����������t`f��  f�������(��  ��)� ��������H��H ���H������Hǅ����   H������H��H�Ǭ����   0���w����H��H��,�����w�[����������tSH�nf_closeH��0��������H��H0���H������Hǅ����   H������H��H�Ǭ����   0���w����H��H���A\A]A^A_H��]H��[��wÐf.�     f.�     �L$(�    H�|$�H�L$�H�|$�0�� [��H��(1���W���)$�    f.�     0Ƀ<� uH���H��|��t0�H��(��Z���D$    �D$    �D$    �D$$   L$L$ L$L$$H��H��H��H�� H��H��I��I��$0��[��H��(�f�UH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��Z��A�E8��$�  H��H�Ǆ  ���e ��  0��uZ���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L����X��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��Y��A�E8��$�  H��H�Ǆ  ���e ��  0��Y���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���X��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� �e �   0���X��A�E8��$�  H��H�Ǆ  � �e ��  0��X���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���#W��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@�e �   0���W��A�E8��$�  H��H�Ǆ  �@�e ��  0���W���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���CV��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��W��A�E8��$�  H��H�Ǆ  ���e ��  0���V���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���cU��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��;V��A�E8��$�  H��H�Ǆ  ���e ��  0��V���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���T��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� �e �   0��[U��A�E8��$�  H��H�Ǆ  � �e ��  0��5U���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���S��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@�e �   0��{T��A�E8��$�  H��H�Ǆ  �@�e ��  0��UT���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L����R��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��S��A�E8��$�  H��H�Ǆ  ���e ��  0��uS���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L����Q��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��R��A�E8��$�  H��H�Ǆ  ���e ��  0��R���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���Q��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� �e �   0���Q��A�E8��$�  H��H�Ǆ  � �e ��  0��Q���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���#P��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@�e �   0���P��A�E8��$�  H��H�Ǆ  �@�e ��  0���P���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���CO��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��P��A�E8��$�  H��H�Ǆ  ���e ��  0���O���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���cN��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��;O��A�E8��$�  H��H�Ǆ  ���e ��  0��O���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���M��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$� �e �   0��[N��A�E8��$�  H��H�Ǆ  � �e ��  0��5N���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���L��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$�@�e �   0��{M��A�E8��$�  H��H�Ǆ  �@�e ��  0��UM���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L����K��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��L��A�E8��$�  H��H�Ǆ  ���e ��  0��uL���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L����J��H��H���[A\A]A^A_]ÐUH��AWAVAUATSH���H�� !  H��I��M��I��I��I�E0   H�C0   I�D$0   H�<$���e �   0��K��A�E8��$�  H��H�Ǆ  ���e ��  0��K���CtH�� �$A�D$tI�$� ��$�  M�E $�$�  AH�$H��H���  0�L��L���J��H��H���[A\A]A^A_]ÐSH��H���UUH�kH�l$H��AWAVAUATH��p'  H������H������L������L������H������H������H��H���H�E�H�E�   H��H���H��p���Hǅx���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H8���H�� ���Hǅ(���	   H��H ���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H`���H��P���HǅX���   H��H@���H��0���Hǅ8���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��Hp���H��`���Hǅh���   H��HP���H��@���HǅH���   H��H0���H�� ���Hǅ(���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H+���H�����Hǅ���   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��HX���H��@���HǅH���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HL���H��0���Hǅ8���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��Hx���H��`���Hǅh���   H��H[���H��@���HǅH���   H��H0���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H0���H�� ���Hǅ(����   L�{0��   L����w��F��1Ʌ�H�Lc�A��   I���   MO�M)�LH�A�   I��MN�M��LH�M)�LH�H��/���L��L���G��N��-/���Hǅ@��� �D HǅH���   H��@���L��L���PG��M��    L��L���F�����  ��.���f�޳  f��,����ͳ  ��(���H���  H�� ���������)�P�����)�p�����)�������)�������)�������)������    ǅ���    ǅ���    ��)�P���L�{(L�c H�������f.�     f.�     ������)�P�����)�p���H��@H��x�ǅ0���    ǅ4���    H��H��P����@�e ��  0���w�^E��H��H�� �����   �:E��A�ƾ�   L���*E��F�$0��   L���E��B�D H������M��Lc�E1�M��L��IH�E��EH�Mc�E9�MO�H��H���H��H)�H���H������H��H��H�� ���L���E��M)�MH��   L���D�����    H�Lc�M9�MN�M)�LH�   I��M��LO�M��LH�L������M)�LH�H������M�4L��H�s(L���=E��M�Hǅ0���'�D Hǅ8���   H��0���L��L���E����   H�{ �D�����    I�Lc�M9�MN�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L�����L��H�s L���D��M�Hǅ ���+�D Hǅ(���   H�� ���L��L���D��M��    L��L���BC���E�    H������H�����H������H�����H�����H��0H�D$ H�D$    H�D$    H�D$    H�$    H��H�ƬI��I������1�E1�E1�L����  H��0������ǅ ���timeL������A���2������H��H ���H�� ���Hǅ���   H�� ���H��H������H��H������A�   L����B����L�����tj�J�  ��*���f�;�  f��(���H�%�  H�� ����L���H��H ���H������Hǅ����   H������H��H��L����   0���w������  ��D����կ  ��@���H��H@���H������Hǅ����   H������H��H������H��H������A�   L������L����w��A����L�����tj���  ��j���f�q�  f��h���H�[�  H��`����L���H��H`���H������Hǅ����   H������H��H��L����   0���w�2���������������������������ǅ����densH��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�°H��H������I��I������E1���w�@��H��P��L�����tj�խ  ��
���f�ƭ  f�����H���  H�� ����L���H��H ���H������Hǅ����   H������H��H��L����   0���w�w�����������(�����������,���ǅ0���tempH��H(���H��P���HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   �E�   H��H0���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�´H��H��P���I��I������E1���w��>��H��P��L�����tj�*�  ������f��  f������H��  H�������L���H��H����H������Hǅ����   H������H��H��L����   0���w輢��������������������������ǅ����yearH��H����H������Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   �E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�¸H��H������I��I��X���E1���w�=��H��P��L�����tj��  ��J���f�p�  f��H���H�Z�  H��@����L���H��H@���H��p���Hǅx���   H��p���H��H��L����   0���w������������h�����������l���� �  ��r���f��  f��p���H��Hh���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��Hp���H��`���Hǅh���   H��`���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H�¼H��H������I��I�� ���E1���w�N;��H��P��L�����tj�Ĩ  ������f���  f������H���  H�������L���H��H����H��P���HǅX���   H��P���H��H��L����   0���w�6��������������������������G�  ������7�  �����H��H���H��0���Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   �E�   H��H���H��@���HǅH���   H��@���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H��0���I��I������E1���w�9��H��P��L�����tj��  ������f���  f������H��  H�������L���H��H����H��0���Hǅ8���   H��0���H��H��L����   0���w�m���������������������������ǅ����dRhoH��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   �E�   H��H����H�� ���Hǅ(���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H������H��H���H��H������I��I������E1���w��7��H��P��L�����tj�`�  ��*���f�Q�  f��(���H�;�  H�� ����L���H��H ���H�����Hǅ���   H�����H��H��L����   0���w貛����  ��L������  ��H���H��  H��@����E�  �|H��H@���H�� ���Hǅ���   H�� ���H��H������H��H������H��H���A�   ��w�j5����L�����tj���  ��t������  ��p�����(p�  ��)�`����L���H��H`���H������Hǅ����   H������H��H��L����   0���w�ƚ���D�  �������4�  ������H��  H�������E�  �|H��H����H������Hǅ����   H������H��H������H��H�Ɛ���H��H���A�   ��w�~4����L�����tj��  �������ң  ��������(��  ��)������L���H��H����H������Hǅ����   H������H��H��L����   0���w�ڙ�����  �������x�  ������H�c�  H�������E�  �|H��H����H������Hǅ����   H������H��H������H��H��X���H��H���A�   ��w�3����L�����tj�&�  �������  �������(��  ��)� ����L���H��H ���H������Hǅ����   H������H��H��L����   0���w�����̢  ��<������  ��8���H���  H��0����E�  �|H��H0���H������Hǅ����   H������H��H������H��H�� ���H��H���A�   ��w�2����L�����tj�j�  ��d����Z�  ��`�����(<�  ��)�P����L���H��HP���H������Hǅ����   H������H��H��L����   0���w������  ������� �  ������H��  H�������E�  �|H��H����H������Hǅ����   H������H��H������H��H������H��H���A�   ��w�1����L�����tj���  ���������  ��������(��  ��)������L���H��H����H��p���Hǅx���   H��p���H��H��L����   0���w�����T�  �������D�  ������H�/�  H�������E�  �|H��H����H��`���Hǅh���   H��`���H��H������H��H�ư���H��H���A�   ��w��0����L�����tj��  �������  �� �����(Ġ  ��)������L���H��H����H��P���HǅX���   H��P���H��H��L����   0���w�*���H��H������1�1�1�E1���w�/  ��L�����t\�o�  ��(���H�Z�  H�� ����L���H��H ���H��@���HǅH���	   H��@���H��H��L����	   0���w評��A���2��H���A���L���H������H��P���HǅX���    H�q     H��`���H�     H��h���Hǅp���    Hǅx���    Hǅ����   A���2H�1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������Hǅ����@�e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ ���   Hǅ���   Hǅ���   H��H���H��0���Hǅ8���    H��@���H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H��H�$    H��H������H��H������H��H��P���H��H������I��I��0���E1���w�  H����L�����tR��(F�  ��)������L���H��H����H��0���Hǅ8���   H��0���H��H��L����   0���w�|���A���2������A�������H������H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   A���2H�1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H�� ���H�����Hǅ0�����e Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H������H��H�Ɛ���H��H�°���H��H��0���I��I������E1���w�  H����L�����tR��()�  ��)������L���H��H����H�� ���Hǅ(���   H�� ���H��H��L����   0���w�O���A���2�����A������H�CH�����Hǅ���    H�q     H�� ���H�     H��(���Hǅ0���    Hǅ8���    Hǅ@���   A���2H�1�H��HH�H��H���HǅP���   HǅX���   IcH��HH�H��`���H��h���Hǅ������e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�����H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H������H��H��X���H��H�����H��H������I��I������E1���w�
  H����L�����tR��(�  ��)�@����L���H��H@���H�����Hǅ���   H�����H��H��L����   0���w�%���A���2��h���A���l���H������H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   A���2H�1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������Hǅ���� �e Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��h���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H������H��H�� ���H��H��p���H��H������I��I��P���E1���w�x	  H����L�����tR��(�  ��)������L���H��H����H�� ���Hǅ���   H�� ���H��H��L����   0���w�����A���2������A�������H������H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ ���   A���2H�1�H��HH�H�����Hǅ���   Hǅ���   IcH��HH�H�� ���H��(���HǅP���@�e HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H������H��H������H��H������H��H��P���I��I������E1���w�k  H����L�����tR��(Օ  ��)� ����L���H��H ���H������Hǅ����   H������H��H��L����   0���w�ˊ��A���2��(���A���,���H�CH��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   A���2H�1�H��HH�H��h���Hǅp���   Hǅx���   IcH��HH�H������H������Hǅ������e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��(���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H������H��H�ư���H��H��0���H��H������I��I�����E1���w�a  H����L�����tR��(��  ��)�p����L���H��Hp���H������Hǅ����   H������H��H��L����   0���w衈��H��H��������w�#����L�����tNH�nf_closeH�E��L���H��H���H������Hǅ����   H������H��H��L����   0���w�7���H��H���A\A]A^A_H��]H��[��w��     f.�     f.�     �L$(�    H�|$�H�L$�H�|$�0��"��H��(1���W���)$�    f.�     0Ƀ<� uH���H��|��t0�H��(�n"���D$    �D$    �D$    �D$$   L$L$ L$L$$H��H��H��H�� H��H��I��I��$0��'#��H��(�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���e �   0��"���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���e ��  0���!��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���$ ��H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �e �   0��� ���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �e ��  0�� ��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L�����H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$�@�e �   0������C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�e ��  0����A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���e �   0�����C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���e ��  0��h��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���e �   0�����C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���e ��  0��H��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L�����H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �e �   0��s���C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �e ��  0��(��A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L�����H��H���[A\A]A^A_]�f�SH��H���UUH�kH�l$H��AWAVAUATH���$  H��8���H�� ���L��H���L��@���H��X���H��P���H��H���H�E�H�E�   H��Hp���H��`���Hǅh���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��HP���H��@���HǅH���   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H ���H������Hǅ����	   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��H0���H�� ���Hǅ(���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H(���H�����Hǅ���   H��H����H������Hǅ����   H��H����H��p���Hǅx���   H��H`���H��P���HǅX���   H��H����H������Hǅ����   H��H����H������Hǅ����   H��HM���H��0���Hǅ8���   H��H���H�� ���Hǅ���   H��H����H������Hǅ����   H��Hx���H��`���Hǅh���   H��H���H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H����H������Hǅ����   H��H|���H��`���Hǅh���   H��H����H������Hǅ�����   HcH��(���1�H��HH�H�����H��H��H���H��H)�H���H��h���H��L�{(��   L����w������    I�Lc�A��   I���   MO�M)�    LH�A�   I��MN�M��LH�1�M)�LH�H������L��L���&��N��-����Hǅ����@�D Hǅ����   H������L��L������M��    L��L��L��h������H�����M�<Ɗ��  ������f���  f���������  ������H���  H������������)�������)�������)�������)������)�0�����)�P����    ǅp���    ǅt���    ��)�����L�k L�cH������f�     ������)�������)�����H��@H��x�ǅ����    ǅ����    H��H�ǰ������e ��  0���w����H�� ����� ��{��Y�c  ��{H��(���������x  ��X�H��H���1���|{H��H��H��>H�H���1�fD  f.�     f.�     ��{L�A�ƀ  ��X���{L���X���{L���X���{L� ��X�H��H9�|����{L���X�H��H9�|�H��H����~"��E�L��H��H��H���0���w�u���   1�����   H��H���H��H����}Ⱦ   H��|A�   f�     f.�     ��|7��|L7 ��|L7@��|L7`H��H��H���|�H��y"�     f.�     ��|7H�� H��x�)ǃ�|Hc������xL���Hc�H9�~	H���{D�L��0���H��H�ǀ�����   I����w�&��A�Ǿ�   L�����F�$8��   L�����B�D H�����M��Lc�1�M��L��HH�E��DH�Mc�E9�MO�H��H���H��H)�H���H��`���H��L��L�����M)��    LH��   L��������    H�Lc�M9�MN�M)�LH�A�
   I��
MN�M��LH�L�����M)�LH�H��`���M�<L��H�s L���.��M�Hǅ����h�D Hǅ����
   H������L��L�������   H�{������    I�Lc�M9�MN�M)��    LH�A�   I��MN�M��LH�1�M)�LH�L����L��H�sL�����M�Hǅ����r�D Hǅ����   H������L��L���v��M��    L��L���3���E�    H��`���H��p���H�����H��x���H��p���H��0H�D$ H�D$    H�D$    H�D$    H�$    H��H�ƬI��I��0���1�E1�E1�L���<  H��0��\���ǅ`���timeL��P���A�$��2��\���H��H`���H��`���Hǅh���   H��`���H��H��\���H��H��0���A�   L�������������tj���  ������f�s�  f������H�]�  H�����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w��u����  ��������  ������H��H����H��@���HǅH���   H��@���H��H��0���H��H������A�   L��X���L����w������������L��h���tj���  ������f���  f������H���  H�����������H��H����H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w�u����0���������������������ǅ����densH��H����H�����Hǅ���    H�q     H�� ���H�
     H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   �E�   H��H����H�� ���Hǅ(���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H�°H��H�����I��I��0���E1���w�x��H��P��������tj�  ��j���f��~  f��h���H��~  H��`��������H��H`���H�����Hǅ���   H�����H��H�Ǭ����   0���w�`s����0���������������������ǅ����tempH��H����H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��H����H�� ���Hǅ���   H�� ���H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H�´H��H������I��I������E1���w���H��P��������tj�[}  ��
���f�L}  f�����H�6}  H�� ��������H��H ���H������Hǅ����   H������H��H�Ǭ����   0���w�q����0�����(�����������,�����|  ��2���f��|  f��0���H��H(���H��P���HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   �E�   H��H0���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H�¸H��H��P���I��I������E1���w����H��P��������tj��{  ������f��{  f������H�{{  H�����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w��o���������������/{  �������{  ������H��H����H������Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   �E�   H��H����H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H�¼H��H������I��I������E1���w�5
��H��P��������tj��y  ��Z���f��y  f��X���H��y  H��P��������H��HP���H������Hǅ����   H������H��H�Ǭ����   0���w�n����������l���fǅp���dzH��Hl���H������Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   �E�   H��Hp���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H���H��H������I��I��P���E1���w���H��P��������tj�Ux  ������f�Fx  f������H�0x  H�����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�ol����0�������������������H�refreezeH�����H��H���H��0���Hǅ8���    H�q     H��@���H�
     H��H���HǅP���    HǅX���    Hǅ`���   Hǅh���   Hǅp���   �E�   H��H���H������Hǅ����   H������H��PH�D$@   H�D$8    H�D$0    H�D$(    H�D$     H�D$    H�D$    H�D$    H�$    H��H��0���H��H���H��H��0���I��I�����E1���w����H��P��������tj��v  ������f��v  f������H�~v  H�����������H��H����H��p���Hǅx���   H��p���H��H�Ǭ����   0���w�j���Kv  �������;v  ������H�&v  H�������E�  �|H��H����H��`���Hǅh���   H��`���H��H��0���H��H��0���H��H���A�   ��w�e����������tj��u  ��������u  ��������(�u  ��)����������H��H����H��P���HǅX���   H��P���H��H�Ǭ����   0���w��i����u  ������u  �����H�ju  H�� ����E�  �|H��H ���H��@���HǅH���   H��@���H��H��0���H��H������H��H���A�   ��w�y����������tj�-u  ��4����u  ��0�����(�t  ��)� ��������H��H ���H��0���Hǅ8���   H��0���H��H�Ǭ����   0���w��h����t  ��\�����t  ��X���H��t  H��P����E�  �|H��HP���H�� ���Hǅ(���   H�� ���H��H��0���H��H������H��H���A�   ��w�����������tj�qt  �������at  ��������(Ct  ��)�p��������H��Hp���H�����Hǅ���   H�����H��H�Ǭ����   0���w��g���t  �������t  ������H��s  H�������E�  �|H��H����H�� ���Hǅ���   H�� ���H��H��0���H��H�ƈ���H��H���A�   ��w�����������tj��s  ��������s  ��������(�s  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w��f��H��H��0���1�1�1�E1���w�b  ��������t\�2s  ������H�s  H�����������H��H����H������Hǅ����	   H������H��H�Ǭ����	   0���w�|f��L�����Hǅ���@   H�q     H�� ���H�    H��(���Hǅ0���    Hǅ8���    Hǅ@���   Ic1�H��HH�H��H���HǅP���   H��H�$    H��H��0���H��H�ƈ���H��H�����1�E1�E1���w��  H����������th�/r  ��|����r  ��x���H�
r  H��p��������H��Hp���H������Hǅ����   H������H��H�Ǭ����   0���w�Ye��H��0���H������Hǅ����@   H�q     H������H�    H������Hǅ����    Hǅ����    Hǅ����   Mc71�M��L��HH�H������Hǅ����   H��H�$    H��H��0���H��H��P���H��H����1�E1�E1���w�  H����������th�q  �������q  ������H��p  H�����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�,d��A�$��2�����D�����H��8���H�����Hǅ���    H�q     H�� ���H�     H��(���Hǅ0���    Hǅ8���    Hǅ@���   A�$��2H�1�H��HH�H��H���HǅP���   HǅX���   IcH��HH�H��`���H��h���Hǅ������e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H�����H������Hǅ����    H�� ���H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��H�$    H��H��0���H��H��0���H��H�����H��H������I��I������E1���w��	  H����������t`f��n  f��P�����(�n  ��)�@��������H��H@���H������Hǅ����   H������H��H�Ǭ����   0���w��a��A�$��2��h���A���l���H��@���H��p���Hǅx���    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ����   A�$��2H�1�H��HH�H������Hǅ����   Hǅ����   IcH��HH�H������H������Hǅ������e Hǅ����    H�q     H�� ���H�
     H�����Hǅ���    Hǅ���    Hǅ ���   Hǅ(���   Hǅ0���   H��h���H��P���HǅX���    H��`���H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��0���H��H������H��H��p���H��H������I��I��P���E1���w��  H����������t`f��l  f��������(�l  ��)����������H��H����H������Hǅ����   H������H��H�Ǭ����   0���w�_��A�$��2������A�������H��H���H������Hǅ����    H�q     H������H�     H������Hǅ����    Hǅ����    Hǅ ���   A�$��2H�1�H��HH�H�����Hǅ���   Hǅ���   IcH��HH�H�� ���H��(���HǅP��� �e HǅX���    H�q     H��`���H�
     H��h���Hǅp���    Hǅx���    Hǅ����   Hǅ����   Hǅ����   H������H������Hǅ����    H������H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��H�$    H��H��0���H��H������H��H������H��H��P���I��I������E1���w��  H����������t`f��j  f�������(�j  ��)� ��������H��H ���H������Hǅ����   H������H��H�Ǭ����   0���w�w]��A�$��2��(���A���,���H�CH��0���Hǅ8���    H�q     H��@���H�     H��H���HǅP���    HǅX���    Hǅ`���   A�$��2H�1�H��HH�H��h���Hǅp���   Hǅx���   IcH��HH�H������H������Hǅ����@�e Hǅ����    H�q     H������H�
     H������Hǅ����    Hǅ����    Hǅ����   Hǅ����   Hǅ����   H��(���H�����Hǅ���    H�� ���H��(���Hǅ0���    Hǅ8���    Hǅ@���   HǅH���   HǅP���   H��H�$    H��H��0���H��H�����H��H��0���H��H������I��I�����E1���w�  H����������t`f��h  f��p�����(wh  ��)�`��������H��H`���H������Hǅ����   H������H��H�Ǭ����   0���w�=[��H��H��0�����w������������tNH�nf_closeH�E������H��H���H��p���Hǅx���   H��p���H��H�Ǭ����   0���w��Z��H��H���A\A]A^A_H��]H��[��w�@ f.�     f.�     �L$(�    H�|$�H�L$�H�|$�0��@���H��(1���W���)$�    f.�     0Ƀ<� uH���H��|��t0�H��(�����D$    �D$    �D$    �D$$   L$L$ L$L$$H��H��H��H�� H��H��I��I��$0������H��(�f�UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L��������C8��$�  H��H�Ǆ  ���e ��  0�����L�AA$�$�  AH��H���  0�L��L��L���W���H��H���[A\A^A_]� f.�     f.�     UH��AWAVATSH���H�� !  H��I��I��H�C0   L�$$���e �   0�L�������C8��$�  H��H�Ǆ  ���e ��  0������L�AA$�$�  AH��H���  0�L��L��L������H��H���[A\A^A_]� f.�     f.�     UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$� �e �   0��3����C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  � �e ��  0������A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���D���H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$�@�e �   0������C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  �@�e ��  0������A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L���$���H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���e �   0�������C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���e ��  0�����A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L������H��H���[A\A]A^A_]�f�UH��AWAVAUATSH���H��@!  I��H��M��I��I��H�C0   H�CH   I�E0   I�D$0   H�<$���e �   0�������C8��$�  �CP��$�  ��$�  ��$�  ��$�  ��$�  H��H�Ǩ  ���e ��  0�����A�EtI�E ��$�@�D$A�D$tI�$���$�  �@��$�  L�$�$�  AH�$H��H���  0�L��L�������H��H���[A\A]A^A_]�f�SH��H���UUH�kH�l$H��AWAVAUATH��0  L��h���L������H�����H������H��x���H�C0HcH������Hǅ����    H�x     H�� ���H�     H�����Hǅx���    H������H������Hǅ����    H�� ���H�����Hǅx���    H������H������Hǅ����    H�� ���H�����Hǅx���    H������H������Hǅ����    H�� ���H�����Hǅx���    H������H������Hǅ����    H�� ���H�����Hǅx���    H������H������Hǅ����    H�� ���H�����H�CLc01�M��LH�H��I��LH�Ik�8H��H���H���   Lc9M��LH�L������Hǅ���    Hǅ���    Hǅ����    Hǅ����    Hǅ���    Hǅ���    Hǅ����    Hǅ����    Hǅ���    Hǅ���    Hǅ����    Hǅ����    Hǅ���    Hǅ���    Hǅ����    Hǅ����    Hǅ���    Hǅ���    Hǅ����    Hǅ����    Hǅ���    Hǅ���    K�D� H��J�D�H���I��I)�I���L������L��H�s8������ ����������������������������H��H��H����H�D$H��H����H�$H��H�Ǭ���H��H�� ���H��H������I��I������I��I������0������H��K�7M��M��I��K�H��0���L��H��I�8H������K��L��H��L��H���I�L��H��L�<L��H��H��J�?H��(���I�L��H��H�<O�Lm K��H������I��H�
H�9H�� ���I�9H��H���Ik�8L�K�<vL��H��I�M��H�����M�L������H�H�
H�<H�H��K��H��K��H�����H������M��L��  L���   L�chH�{ �#  L������H��8���H������L������H�� ���L������H��H���L������L������L������A��  }T1�A��H������L������L�������F  L��H���H��H��1�H���%  L��(���L�� ���H��8����S  ���e 0�L��M���>������e 0�L��L���,������e 0�L��8���L��L���������e 0�H�����L����������e 0�H��(���L����������e 0�H��0���L����������e 0�L��L���������e 0�H�� ���L���������e 0�H�����L������M��H������L������L������H�� ���L������H������L������L��  L��L�chL���   �q  H�� ���1�H��0���L��(���L�� ���H��8���L�����L�����f���������|����|��|D ����|)��|��|��D ��|D ��D ��|D ��|D ��D ��|)D ��|D ��|D ��D@��|D@��D@��|D@��|D@��D@��|)D@��|D@��|D@��D`��|D`��D`��|D`��|D`��D`��|)D`��|D`��|D`H��H��H�������L������L������H�� ���H����   H�� ���L������L������L��0���L�����H�����fD  f.�     f.�     ��������|������|D ��|��|��|��|H�� H��x�L������L������H�� ���H��8���L�� ���L��(���H������Hc�L��H)�H��L������M��I��|cHc���W������x�H��8������H��������H��(������H��0��������x�H�� ������H�������ȃ�H�� ���H������H�� ���A9�L������L������H��H���H��8���L��  H��L�chL���   H������H�{ ��   H�� ���H�H��H���H��    H������H��    H��    H�����H��    H��(���H��    H��0���H��    H������H��    H�� ���H��    H�����H��    H������H�{ H��8���L������H������L������H������L������H�� ���L������L������	A	A
A
AL������Hǅ�����   H��H������H��   H�D$hH������H�D$`H�|$XH�CH�D$PH��h���H�D$HH������H�D$@L�\$8L�T$0H�t$(H�T$ L�|$L�L$H�L$L�4$H�D$p�   0�H���   H���   H���   L���   L���   H���   M����w�Gn��H�Ā   I�H��p���H�On! H��������(0n! ��)�������(n! ��)�p���H��Hp���H������H��H��p���H��H�°������e �����L�s0AM�������H�S@
A$�p���H�CpH��H�CXI��AL��I��L���   M��AL���   A� ���L��0���A$L�����AL��H���AM H������H��(���H��  H��p���Hǅx����   H��p���H��`H�D$PH�� ���H�D$HH�t$@H�|$8L�l$0L�t$(L�d$ H��H ���H�D$L�|$L�\$L�$H�D$X�   H��H�Ƭ���0�M��L��I��I��p�������H��`L������Hǅ����@   H�`     H������H�     H������L������Hǅ����@   Hǅ����   Hǅ����
   Hǅ����   ��(�m! ��)�0�����(|m! ��)� �����(\m! ��)����H��H����H�� ���H��H��������e H��H�°��������AL�������H��x���
L���   A
L�C`AA$L������AH��8���	AM AH�����H�{pH��0H�|$(H�D$ L�t$L�l$H�L$L�<$0�I��L��H��H�Ƭ���M��M��L������H������H��0H�{H��D �H��H�������  H��(���H������L�CAL�SA
A�����L������A
L�cXA$L�[pAL���   AL��H�SP
L�KhA	�p���H�Kx	� ���H��  H��`���Hǅh����   H��`���H��   H��$  H�C(H��$  H�����H��$   H���   H��$�   H������H��$�   H������H��$�   H������H��$�   H������H��$�   H������H��$�   H������H��$�   H�� ���H��$�   H�����H��$�   H�� ���H��$�   H��(���H��$�   L�����L��$�   H������H��$�   L��H���L��$�   H��0���H��$�   L��8���L��$�   H������H�D$xH��H����H�D$pH��H����H�D$hH��H����H�D$`H��H����H�D$XH��H����H�D$PH��H ���H�D$HH�L$@H��Hp���H�D$8L�L$0H�T$(H�t$ L�D$L�\$L�d$L�$HǄ$  �   I��I������0�H������H�SL�c0M��H�K��w�>���H��   H������L��H�� ���M��H��(���L�����M��M����  H��(���H��H�°���� �e M��� �e ��w�D���AL�������L��H���AA$L��8���AM L�����AH��(���I��L��0���AH��  H��P���HǅX����   H��P���H��  H��@���HǅH����   H��@���H��   H��0���Hǅ8����   H��0���H��   H�� ���Hǅ(����   H�� ���H��PH�t$(H�T$ H�L$H�D$L�D$L�$H�D$H�   H�D$@�   H�D$8�   H�D$0�   H��H�Ƭ���0�L��M��M��L����#��H��PH��������H������L��H�� ���M��H��(���L�����M��L�s0�E  I��I��H��������  ��   E1ۃ�L����   H��(���H��H��1�H��|E1���}�:  �    f.�     ��|7��|D7 ��|D7@��|D7`H��H��H���|�H��y"��}�:  f�     ��|7H�� H��x�H��(���A��Ic�H��H)�H��L�����|Ic���(�:  ��x�A��D9�H��(���H�� ���L��8���L��~JMc�H�     8��K���7� �e 0�L����w�)���H��(���H�� ���L��8���L��L��L�����H��H�� ���L�����H������H��H���L�����L��8���H��(���AH������H�S(
	I��A� ������������������������������������������������H���   H��Hx���H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H���   H��$�   H���   H��$�   H���   H��$�   H��H����H��$�   H��H����H�D$xH��H����H�D$pH��H����H�D$hH��H����H�D$`H��H����H�D$XH��H����H�D$PH��H����H�D$HH��H����H�D$@H��H����H�D$8H��H����H�D$0H��H����H�D$(H��H����H�D$ H��H����H�D$H��H����H�D$H��H����H�D$H��H����H�$M��I��I�� ���0�L��M����w����H���   A�����L��H���AL������A
AM A$H��(���AH��   H�����Hǅ����   H�����H��  H�� ���Hǅ����   H�� ���H��  H������Hǅ�����   H������H��@H�T$ H�L$H�D$L�|$H�4$H�D$8�   H�D$0�   H�D$(�   H��H�Ƭ���0�L��L��M��M��L��菮��H��@H��������H�� ���H!�I��H��H�� ���Hǅ ���   ��������2��<���H�H��`���1�H��HH�H��0���H��(���Hǅ0���   Hǅ8���   Hǅ@���   H��H���H�      H��0���L������L��8���Hǅ����    Hǅ����    H������H��0���1�1�E1�0��&���0�L������H������L!�H��H������Hǅ����   ��������2Lc�M��    LH�L��p���L������Hǅ����   Hǅ����   H���   Lc8M��LH�L������L������L������I�      L��p���L��p���L��x���Hǅ����    Hǅ����    H������H��p���1�1�E1�0��I���0�L���/���H������I��������L!�H��H������Hǅ����   L������Hǅ����   Hǅ����   L������L������L������H��Hp���H������Hǅ����    Hǅ����    H������H��H�ǰ���1�1�E1�0�����0�H��H��p����}���H�� ���L!�H��H�� ���Hǅ ���   L��(���Hǅ0���   Hǅ8���   L��@���L��H���L������I��I������L������Hǅ����    Hǅ����    H������H��H������1�1�E1�0������0�L�������H�� ���I��������L!�H��H�� ���Hǅ ���   L��(���Hǅ0���   Hǅ8���   L��@���L��H���L��0���H��H����H��8���Hǅ����    Hǅ����    H������H��H��0���1�1�E1�0��J���0�H��H�������)���H������L!�H��H������Hǅ����   L������Hǅ����   Hǅ����   L������L������L��p���I��I��p���L��x���Hǅ����    Hǅ����    H������H��H��p���1�1�E1�0�����0�L������H�� ���L!�H��H�� ���Hǅ ���   L��(���Hǅ0���   Hǅ8���   L��@���L��H���I�      L������I��I������L������Hǅ����    Hǅ����    H������H��H�ǰ���1�1�E1�0������0�L�������H������L!�M��H��H������Hǅ����   ��������2Lc�M��    LH�L������Hǅ����   Hǅ����   H������H������I��L������L������L��p���L������Hǅp���    Hǅx���    H��p���H������1�1�E1�0��4���0�L������H�� ���L!�H��H�� ���Hǅ ���   L��(���Hǅ0���   Hǅ8���   L��@���L��H���I�      L��0���I��I������L��8���Hǅ`���    Hǅh���    H��`���H��H��0���1�1�E1�0�����0�L���n���H�� ���I��������L!�H��H�� ���Hǅ ���   L��(���Hǅ0���   Hǅ8���   L��@���L��H���L��p���I��I������L��x���HǅP���    HǅX���    H��P���H��H��p���1�1�E1�0������0�L�������L��H#�����H��H������Hǅ����   L������Hǅ����   Hǅ����   H������H������L������L������I��I��p���L������Hǅ@���    HǅH���    H��@���H��H�ǰ���1�1�E1�0��3���0�L������H��`������  H��0���H�	H����<�����   H����   H��H���H��H��1�H������H���N1���}B  f.�     f.�     f.�     ��)9��)D9 ��)D9@��)D9`H��H�� H���|�H��y%��}�A   f.�     ��)9H�� H��x�H��H)�H��|��y�A  ���H��Hc�H9�~6�f.�     H���  �|��Hc�H9�|��H�������@�e 0���w�e���L������L��H��H��p���H��H��H����  L��I��  �_  M����  1�I���X  L��H���H��H��1�L��p���L��p���L������L������H��p���H������H�����   1���}�@  f.�     f.�     f.�     ��|)��|)��|)��|)��)��)��|)D ��|)D ��|)D ��|)D ��)D ��)D ��|)D@��|)D@��|)D@��|)D@��)D@��)D@��|)D`��|)D`��|)D`��|)D`��)D`��)D`H��H�� H����T���H��y;��}�?  fD  ��|)��|)��|)��|)��)��)H�� H��x�Hc�L��H)�H��|WHc�H��p�����y�?  ���H��p������H���������H���������H��p������H�����������Hc�I9��  L��p���L��p���L������H������H��p���H������@ f.�     f.�     H�A��  �|A��  �|A��  �|��  �|��  �|��  �|��Hc�L9�|��   H��p����@�e 0�L����w����H��p����@�e 0�L������H�������@�e 0�L������H�������@�e 0�L���w���H��p����@�e 0�L���a���H�������@�e 0�L���K���M��L��H��L������L��H��H���B  M��I��  ��  M���(  1�I����   L��H���H��H��1�L��p���L������H������H��p���H�����   1���}�=  f�f.�     f.�     ��|)��|)��)��)��|)D ��|)D ��)D ��)D ��|)D@��|)D@��)D@��)D@��|)D`��|)D`��)D`��)D`H��H�� H���|�H��y)��} =  ��|)��|)��)��)H�� H��x�Hc�L��H)�H��|NHc�H��p�����y�<  ������H������������H������������H��p��������Hc�I9���   L��p���H������H������H��p����    f.�     f.�     H�A��  �|��  �|��  �|��  �|��Hc�L9�|��[H��p����@�e 0�L����w�,���H�������@�e 0�L������H�������@�e 0�L��� ���H��p����@�e 0�L�������H��H�°������e ���e ��w����H�CHc H������ǅ����   H����  1�L�{0L������L���   L��H���f�H��@���L��H���D�nD������A�����AH�KX	L�CpAH���   I��A	A$A
H��8���
H��  H��0���Hǅ8����   H��0���H��0H�D$H�T$L�T$L�$$H�D$ �   H��H�Ƭ���0�L��L����w��E��H��0H�����D� A����   D��,���A�����AL�CPAL���   A	H�� ���	H�CI��AL��8���A
L��H���AH��0���H�����H�Sp
H��@H�T$0H�D$(H�|$ L�t$L�T$L�\$H�$H��H��,���H��H�¬���0�L��L������L����w�1��H��@A��L�cM��M����   D��,���A�����A$L�� ���A	H��8���
H��H���H��0���H�spH�� H�t$H�|$H�D$H�$H��H��,���H��H�¬���0�L��L��M����w�y#��M��H�� ����D����������i�@ D��)�D)�Hc�H��@���H�H������H��H��p���H���L���   ��   H��V! H�E���(�V! ��)E���(�V! ��)�p�����(`V! ��)�`�����(@V! ��)�P���H��H����H��h���H��H ���H�E�H��H��P���H��H�°���� �e ��w�z���L�CD��,���A�����A$L�K(A	H�sh�p���H�Sp
AH�Kx	� �������������������������H��   H��H����H��$  H��H����H��$  H��H����H��$   H��H����H��$�   H��H����H��$�   H���   H��$�   H�����H��$�   H�����H��$�   H�� ���H��$�   H��(���H��$�   H�����H��$�   H������H��$�   H��H���H��$�   H��0���H��$�   H��8���H��$�   L������L��$�   H������H��$�   H������H��$�   H������H��$�   H������H�D$xH������H�D$pH�� ���H�D$hH��H����H�D$`H��H ���H�D$XH��H���H�D$PH��H���H�D$HH��H����H�D$@H��H���H�D$8H��H����H�D$0H��H ���H�D$(H�L$ L�t$H�T$H��Hp���H�D$H�4$H��H��,���H��H�¬���0�L��L��M����w�l��H��   Hǅ����    Hǅ����    Hǅ����    Lc�����M���0  E1���W���W�A��H��H���L��0���H��(����  H�Sp��
��Y!  ��^���}�M��I���L��H��1���}�������}-%  �������Ґf.�     f.�     f.�     ��<1�18  ������P���t�BM-2��\���EY��������eK�`��X0�08  ��uXL5 A�58  H�� H��x���}���X���y���X���}���X���y���X���}���X���y���X���X������������L��L)�H����   ġx����j4  ��P҅�tRĂi-$�H�Sp��*�����\���Y���%   ��^������Y���W���YK� ��y���X���X������������ġx���y���X���X�āxT� ��y���X���X�I��M9�~Cġ{���.�#  w&H�Kp����\�ākY���^���X������������ġ{X�āsXL� ��������������L��@���E�l$D��,��������A� �������������������������� ���������������������������������H������H���   H��H������H��$�   H��$�   H��Hx���H��$�   H��Hp���H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H��$�   H��H����H�D$xH��H����H�D$pH��H����H�D$hH��H����H�D$`H��H����H�D$XH��H����H�D$PH��H����H�D$HH��H����H�D$@H��H����H�D$8H��H����H�D$0H��H����H�D$(H��H����H�D$ H��H ���H�D$H��H���H�D$H��H���H�D$H��H����H�$H��H��,���H��H������H��H�� ���I��I������I��I�����0�L����w��}��L��H���   D��������D��)�H�H��H)�H���L��H���L���   �Z  D��,���A�����A�����L���   AA
H��8���H��(���
L�����AH������	L�����A$H��p���H��p���H������H������L��p���L��p���L������L������H��p���H������H��pH�D$`H�|$XL�\$PL�T$HL�L$@L�D$8L�d$0L���   H�L$(L�|$ H������H�T$H�t$H��H���H�D$L�4$H��H�Ƭ���H��H��,���I��I������I��I������0�H�{0I��I����w�@���H��@���L��H���H��pD��������D��)�H�H��H)�H����,  D��,���L�������H���   	L���   A�����A
L��8���AL��(���AM H�� ���
H��@���H��0���H��p���H������H������H��p���L��p���L������L������H��p���H��PH�D$HL�T$@L�L$8L�|$0H�t$(H�T$ L�l$L�\$H��H���H�D$H��H����H�$H��H�Ƭ���H��H��,���I��I������0�I����w�����H��@���L��H���H��PH��H;����������H��������������L��H���H����������   1���H�����||1�1�H�{(+��n���y�����H��(���H��H��H��|@1��     f.�     ��X��  ����XL ��L H��@H��H���|�H��y
��X��H��H��H)�H��|"1�H�{(+��W���*������X����H��H9�~!1�H�S(+
��W���*�H�������X����I��I�ư������e ���e L����w����� �e � �e L�������H�O! H�E���(�N! ��)E���(�N! ��)E�H��  H�� ���I��Hǅ(����   H�� ���H�E�H�E��   H��H�ư���e L��蕻�������H������H��   H�����I��Hǅ����   H�����L�� ���M��Hǅ����   H�� ���H��  H������I��Hǅ�����   H������L������H��H�D$�   H�$�   H��H������A��   0��o���H�������H���   H��p���H��p���H������H������H��p���H������H��p���H��p���L������L������L��p���L������L������Hǅ�����   L������L������Hǅ�����   H������L������Hǅ�����   H������H��@H�|$ H�D$L�|$L�\$L�$H�D$8�   H�D$0�   H�D$(�   H��H������0���m��H��@�����H���   H���   
H��p���H������H������H��p���H��p���L������L������L��p���L������L��H��@���Hǅ�����   L������L������Hǅ�����   H������L������Hǅ�����   H������H��@H�|$H�D$L�\$L�$H�D$0�   H�D$(�   H�D$ �   H��H������0��F���H��@�����L�{0AL�������H��H���L������AH������L��8���A	L�����AL��(���AH�����	H�� ���H��0���H��@���H������Hǅ�����   L������L��p���Hǅx����   L��p���L��`���Hǅh����   H��`���H��`H�T$@L�T$8L�D$0H�|$(H�T$(H������H�D$ H�L$H�D$H��X���L�|$H�D$H��P���L�\$L�$H�$H������H�D$X�   H�D$P�   H�D$H�   H��H�Ǆ���H��H�¬���0�I��M��L������M��L��H���L��芧��H��`���e ���e I��I�ư���L���J����@�e �@�e L���8���H�{H������A$AM H������H������L��P���A$AL��X���AM H���   I��AH���   I��A	H���   I��A
H���   I��AH���   H��
H���   H��	H��   H������H�D$pH������H�D$hH������H�D$`H������H�D$XH������H�D$PH�� ���H�D$HH�L$@H�T$8L�\$0L�T$(L�L$ L�D$L�l$L�|$L�$$0�L������L������H�S0H�������`w��H�Ā   ���e ���e L�������@�e �@�e L���յ�����e ���e L���õ��H�      H������H��H����H������H��Hp���H�� ���H��H����H�����H��Hp���H�����H��H����H�����H��Hp���H�� ���H��H����H��(���H��Hp���H��0���H��H����H��8���H��Hp���H��@���H��H����H��H���HǅP���    HǅX���    H��P���H��H������1�1�E1�0��Ͳ��H��H���A\A]A^A_H��]H��[�f.�     f.�     �                               ��fffff.�     H�l$�L�|$�H�-�C  L�=�C  L�d$�L�l$�L�t$�H�\$�H��8L)�A��I��H��I��賱��H��t1�@ L��L��D��A��H��H9�r�H�\$H�l$L�d$L�l$ L�t$(L�|$0H��8��    UH��S� �d H��H�KC  H���tD  H����H�H���u�H��[]�f�UH���   H��]������                    MALLOC_MMAP_MAX_= MALLOC_TRIM_THRESHOLD_= CRAY_MALLOPT_OFF= ���K       ����?    o����*>��L>  �>���=    8~A   ��6?     �p@      �?      �    8~A    8~A                                ANT27           /FM_Data/INPUT/ANT27_averages/  /scratch/ms/nl/ _ANT27_79-16_ave.ncXPEN055      /FM_Data/INPUT/XPEN055_averages//scratch/ms/nl/ _XPEN055_79-16_ave.ncFGRN11     /data/input/era_files/averages/ /scratch/ms/nl/ _FGRN11_60-79_ave.ncFGRN055     /data/input/era055_files/averages/      /scratch/ms/nl/         _FGRN055_60-80_ave.ncPAT055     /FM_Data/INPUT/PAT055_averages/ /scratch/ms/nl/ _PAT055_79-12_ave.ncXDML055                     /FM_Data/INPUT/XDML055_averages//scratch/ms/nl/ _XDML055_79-15_ave.ncASE055     /FM_Data/INPUT/ASE055_averages/ /scratch/ms/nl/ _ASE055_79-15_ave.ncDMIS055     /FM_Data/INPUT/DMIS055_averages//scratch/ms/nl/ _DMIS055_79-17_ave.nc   no valid domainANT27    ../lsm_ANT27.ncLSM      nf_inq_varidLat nf_inq_varid11Lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lon  ../ism_ANT27.ncISM      nf_inq_varid    nf_get_var_ism  nf_close1       nf_close6XPEN055        ../Height_latlon_XPEN055.ncmask2d       nf_inq_varidlat nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lon  iceshelves      nf_inq_varid    nf_get_var_ism  nf_close1FGRN11                 ../mask/FGRN11_Masks_wholedomain.ncicemask      nf_inq_varid_lsmlat     nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lonFGRN055           ../mask/FGRN055_Masks.ncIcemask_GR              nf_inq_varid_lsmlat     nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lonPAT055            ../lsm_PAT055.ncnf_inq_varidlat nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lonXDML055   ../lsm_XDML055.ncism    nf_inq_varidlat nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lonASE055            ../Masks_ASE055.ncLSM   nf_inq_varidISM nf_inq_varidlat nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lonDMIS055           ../Masks_DMIS055.ncLSM  nf_inq_varidISM nf_inq_varidlat nf_inq_varid11lon       nf_inq_varid12  nf_get_var_lsm  nf_get_var_lat  nf_get_var_lon  no valid domain snowmeltsnowmeltprecipff10mtskinevapsndiv       nf_inq_varidprecip      nf_inq_varidff10m       nf_inq_varidtskin       nf_inq_varid    nf_inq_varidsndiv       nf_inq_varid    nf_get_var_avemelt              nf_get_var_aveacc               nf_get_varavewind               nf_get_varavetsurf              nf_get_varavesubl               nf_get_varavesndiv      nf_close1       nf_close2       nf_close3       nf_close4       nf_close5       nf_close7ANT27          /FM_Data/INPUT/ANT27_files/     /scratch/ms/nl/.nc      _ANT27_79-16_pXPEN055   /FM_Data/INPUT/XPEN055_files/   /scratch/ms/nl/.nc              _XPEN055_79-16_pFGRN11          /data/input/era_files/files/    /scratch/ms/nl/.nc      _FGRN11_60-16_pFGRN055  /data/input/era055_files/files/ /scratch/ms/nl/.nc              _FGRN055_57-20_pPAT055          /FM_Data/INPUT/PAT055_files/    /scratch/ms/nl/.nc      _PAT055_79-12_pXDML055  /FM_Data/INPUT/XDML055_files/   /scratch/ms/nl/.nc              _XDML055_79-15_pASE055          /FM_Data/INPUT/ASE055_files/    /scratch/ms/nl/.nc      _ASE055_79-15_pDMIS055  /FM_Data/INPUT/DMIS055_files/   /scratch/ms/nl/.nc              _DMIS055_79-17_pno valid domain snowmeltprecip  snowfallevaptskinsndivff10m     nf_inq_varid1precip     nf_inq_varid2   nf_inq_varid3   nf_inq_varid4tskin      nf_inq_varid5sndiv      nf_inq_varid6ff10m      nf_inq_varid7   nf_get_var1     nf_get_var2     nf_get_var3     nf_get_var4     nf_get_var5     nf_get_var6     nf_get_var7     nf_close1       nf_close2       nf_close3       nf_close4       nf_close5       nf_close6       nf_close7       /FM_Data/INPUT/ini_files//      /scratch/ms/nl/_ini_.nc nf_inq_varid1   nf_inq_varid2   nf_inq_varid3depth      nf_inq_varid4lwcnf_inq_varid5layer      nf_inq_dimid1   nf_inq_dim1     nf_get_var1     nf_get_var2     nf_get_var3     nf_get_var4     nf_get_var5     nf_close1            8�@����F߁?      �?   ���@   ज़�?    \_X@   `fq�   �p=@   ���v@   `fq@     L�@     ��@��A����?��Xr��?     0�@   ��|@     c@Nb�)q>   `fq@     8��   ��c?   �G�d@    �Xw�      |�   @�z��   ��#@���_�?tw_f;��@    2 �   ����>   �t��?   �"� @    �m?      $@      $�����F߁?����F߁?       �       �     8��     8��   ��c?   ��c?   �G�d@   �G�d@    �Xw�    �Xw�      |�      |�   @�z��   @�z��      �?      �?   ��#@   ��#@���_�?���_�?tw_f;��@tw_f;��@   ����>   ����>    2 �    2 �   �t��?   �t��?    �m?    �m?   �"� @   �"� @      $@      $@   `fq@   `fq@      $�      $�  �?  ��               �       �       �       �/data/ms_files/ /scratch/ms/nl/ model_settings__.txt    /perm/ms/nl/            /code/DATA/input_settings_.txt  /perm/ms/nl/    /code/DATA/input_settings_.txtFGRN11FGRN055FGRN11FGRN055              �?      �?     0�@     0�@     L�@     L�@   ��|@   ��|@     c@     c@      �?      �?Nb�)a>   ��Sÿ   ����?   �uq{?    S�?    ��ɿ   ����?    ��ҿ   `��@      �?UUUUUU�?      �     L��   `f
q@   ��ư>     ��@   �/�?   �rh�?	uU8;'�>    �ZA   ��c�@n���O�      �     p�@A�����Q?I13���K                       �       �       �       �                                       �       �       �       �       �       �       �       �                                FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055FGRN11FGRN055   ����>/data/output/era055/    /scratch/ms/nl/_ini_.nc nf_createlayer  nf_def_dim1     nf_def_var1     nf_def_var2     nf_def_var3depthnf_def_var4lwc  nf_def_var5     nf_def_var6     nf_enddef       nf_put_var1     nf_put_var2     nf_put_var3     nf_put_var4     nf_put_var5     nf_put_var6     /data/output/restart/   /scratch/ms/nl/ _restart_.nc    nf_createlayer  nf_def_dim1     nf_def_dim2     nf_def_var1     nf_def_var2     nf_def_var3depthnf_def_var4lwc  nf_def_var5     nf_def_var6     nf_def_var7dz   nf_def_var8     nf_def_var9kUL  nf_def_var10    nf_def_var11    nf_enddef       nf_put_var1     nf_put_var2     nf_put_var3     nf_put_var4     nf_put_var5     nf_put_var6     nf_put_var7     nf_put_var8     nf_put_var9     nf_put_var10    nf_put_var11    constant/data/output/era055/    /scratch/ms/nl/_1D_.nc  nf_def_dim1zs   nf_def_var1     nf_def_var2     nf_def_var3vfc  nf_def_var4vmeltnf_def_var5vbouynf_def_var6     nf_def_var7     nf_def_var8vtotal       nf_def_var9Runoff       nf_def_var10FirnAir     nf_def_var11TotLwc      nf_def_var12    nf_def_var13    nf_def_var14    nf_def_var15solin       nf_def_var16icemass     nf_def_var17    nf_def_var18    missing_value   nf_def_att_miss_val_1   missing_value           nf_def_att_miss_val_2   missing_value           nf_def_att_miss_val_3   missing_value           nf_def_att_miss_val_4   missing_value           nf_def_att_miss_val_5   missing_value           nf_def_att_miss_val_6   missing_value           nf_def_att_miss_val_7   missing_value           nf_def_att_miss_val_8   missing_value           nf_def_att_miss_val_9   missing_value           nf_def_att_miss_val_10  missing_value           nf_def_att_miss_val_11  missing_value           nf_def_att_miss_val_12  missing_value           nf_def_att_miss_val_13  missing_value           nf_def_att_miss_val_14  missing_value           nf_def_att_miss_val_15  missing_value           nf_def_att_miss_val_16  missing_value           nf_def_att_miss_val_17  missing_value           nf_def_att_miss_val_18  nf_enddef               nf_put_var_speed1       nf_put_var2             nf_put_var_speed3               nf_put_var_speed4               nf_put_var_speed5               nf_put_var_speed6               nf_put_var_speed7               nf_put_var_speed8               nf_put_var_speed9               nf_put_var_speed10              nf_put_var_speed11              nf_put_var_speed12              nf_put_var_speed13              nf_put_var_speed14              nf_put_var_speed15              nf_put_var_speed16              nf_put_var_speed17              nf_put_var_speed18              /data/output/era055/    /scratch/ms/nl/_2D_.nc  nf_def_dim2layernf_def_dim3     nf_def_var1     nf_def_var2     nf_def_var3lwc  nf_def_var4depthnf_def_var5     nf_def_var6     missing_value   nf_def_att_miss_val_1   missing_value           nf_def_att_miss_val_2   missing_value           nf_def_att_miss_val_3   missing_value           nf_def_att_miss_val_4   missing_value           nf_def_att_miss_val_5   missing_value           nf_def_att_miss_val_6   nf_enddef               nf_put_var_grid1nf_put_var_grid2nf_put_var_grid3nf_put_var_grid4nf_put_var_grid5nf_put_var_grid6/data/output/era055/    /scratch/ms/nl/ _2Ddetail_.nc   nf_def_dim4layernf_def_dim5     nf_def_var1     nf_def_var2lwc  nf_def_var3depthnf_def_var4dz   nf_def_var5     nf_def_var1     missing_value           nf_def_att_miss_val_1   missing_value           nf_def_att_miss_val_2   missing_value           nf_def_att_miss_val_3   missing_value           nf_def_att_miss_val_4   nf_enddef       nf_put_var1_1   nf_put_var2_1           nf_put_var_detail1              nf_put_var_detail2              nf_put_var_detail3              nf_put_var_detail6                   p�@     p�@  �|;p  �   �����  |����  ����  <���  <���4  �?��\  �?��t  �@���  �A���  �B���  �B��  �C��,  D��D  <E��l  F���  �F���  H���  <H���  \I��$  <J��L  K��t  <K���  \L���  <M���  N��	  <N��	  \O��D	  <P��l	  Q���	  <Q���	  \R���	  <S���	  T��$
  <T��<
  \U��d
  <V���
  W���
  <X���
  \X���
  |Y��  \Z��D  <[��l  \\���  |\���  �\���  �\���  �\���  �\��  ]��$  <^��L  \_��t  |`���  �a���  �b���  �c��  ����<  ����T  ����l  �����  ����  <����  \����  |����  |���  |���4  |���\  |����  |����  |����  |����  ����  \���D  |���\  \����  <����  ����  �����  ����$  ����D  ���l  \	���  �	���  ����  �?���  \D��$  �G��\  �X���  \_���  �_���  ����  ���$  ���L  ��t  �5���  �@���  |E��  \F��,  �L��L  \P���  |X���  \Y���  \Z���  �^��  <d��L  �y��t  �y���  �z���  \{���  |���  �|��  �}��D  \~��l  ���  ����  <����  |����  ���  ܣ��,  ����T  \���|  ����  ܦ���  �����  \���  ���D  ܩ��l  \����  ܪ���  �����  |���  |���<  \2��d  |2��|  3���  �3���  �4���  �5��  �6��4  |7��\  \8���  <9���  :���  �:���  �;��$  �<��L  �=��t  |>���  \?���  <@���  A��  �A��<  �B��d  �j���  �j���  �k���  �l���  �m��  �n��4  p��\  <q���  \r���  <����  \����  ����  ����,  |���T  ����|  �����  ܝ���  �����  |���  ����4  T����         zR x�  $      ����   FJw� ?;*3$"       D   ����    A�DI       d   ����              |   ����<   A�C7 $   �   (����   A�CT�����    $   �    ���e�   S�CO����         �   X9��            $     `9��   A�CT�����    $   ,  X:���    A�CR����      $   T  ;���    A�CR����         |  �;��            $   �  �;��   A�CT�����       �  �<��            $   �  �<��   A�CT�����    $   �  �=���    A�CR����      $   $  �>���    A�CR����      $   L  8?��   A�CT�����       t  0@��            $   �  8@��   A�CT�����    $   �  0A���    A�CR����      $   �  �A���    A�CR����           �B��            $     �B��   A�CT�����    $   D  �C���    A�CR����      $   l  XD���    A�CR����         �  E��            $   �  E��   A�CT�����    $   �  F���    A�CR����      $   �  �F���    A�CR����         $  �G��            $   <  �G��   A�CT�����    $   d  �H���    A�CR����      $   �  8I���    A�CR����         �  �I��            $   �  �I��   A�CT�����    $   �  �J���    A�CR����      $     �K���    A�CR����      $   D  `L��   A�CT�����       l  XM��            $   �  `M��   A�CT�����    $   �  XN���    A�CR����      $   �  O���    A�CR����      $   �  �O��   A�CT�����       $  �P��               <  �P��               T  �P��               l  �P��               �  �P��               �  �P��            $   �  �P��   A�CT�����    $   �  �Q��   A�CT�����    $     �R��   A�CT�����    $   ,  �S��   A�CT�����    $   T  �T��   A�CT�����    $   |  �U��   A�CT�����    $   �  �V���d   S�CO����         �  X���               �  `���               �  h���                 p���               ,  x���               D  ����               \  ����            $   t  �����    A�CT�����    $   �  h����    A�CT�����    $   �  @����    A�CT�����    $   �  ����    A�CT�����    $   	  ����    A�CT�����    $   <	  ȿ���    A�CT�����    $   d	  �����    A�CT�����       �	  x���_   A�CP���$   �	  ����k   S�CO����         �	  ���            $   �	  ����    A�CT�����    $   
  �����    A�CT�����    $   <
  �����    A�CT�����    $   d
  @����    A�CT�����    $   �
  �����    A�CT�����       �
  �����   A�CP���$   �
  p����   S�CO����      $   �
  (����   S�CO����         $  ����B           $   <  ���0   A�CT�����    $   d   ���'   S�CO����      $   �  �.��U   A�CT�����    4   �  03��!   ABB B(B0A8D`������      4   �  86���   ABB B(B0A8G�	������     $   $   G���   A�CT�����       L  �M��9           $   d  �M��t�   S�CO����      $   �  ����   S�CO����      $   �  �����   S�CO����      $   �  ����   S�CO����      4     �����*   ABB B(B0A8G�������     4   <  H#���
   ABB B(B0A8G�������     ,   t  �-���   ABB B(B0A8������    �  `2���              �  (3��2   BAA ��  4   �  H9���   ABB B(B0A8G�������     4     �<��   ABB B(B0A8D�������        L  �D���              d  �E���           ,   |  hF���   ABB B(B0A8������ ,   �  �J��+   ABB B(B0A8������ $   �  �O���   S�CO����           `e��                 he���    D0    $   4  �e���    A�CR����      $   \  �f���    A�CR����      $   �   g���    A�CR����      $   �  �g���    A�CR����      $   �  Ph���    A�CR����      $   �  �h���    A�CR����      $   $  �i���"   S�CO����         L  X���               d  `���3              |  �����    D0    $   �  ����    A�CR����      $   �  �����    A�CR����      $   �  @����    A�CR����      $     ؎���    A�CR����      $   4  p����    A�CR����      $   \  ����    A�CR����      $   �  �����    A�CR����      $   �  8����    A�CR����      $   �  Б���    A�CR����      $   �  h���~    A�CR����      $   $  ����~    A�CR����      ,   L  ����   ABB B(B0A8������ $   |  �����   S�CO����      $   �  `����   S�CO����      $   �  8����d   S�CO����         �  ���                 ����    D0    $   $  ����    A�CT�����    $   L  8���    A�CT�����    $   t  ����    A�CT�����    $   �  ����    A�CT�����    $   �  `���    A�CT�����    $   �  ���    A�CT�����    $     ����    A�CT�����    $   <  ����    A�CT�����    $   d  @ ���    A�CT�����    $   �  � ���    A�CT�����    $   �  �!���    A�CT�����    $   �  h"���    A�CT�����    $      #���    A�CT�����    $   ,  �#���    A�CT�����    $   T  �$���    A�CT�����    $   |  H%���    A�CT�����    $   �   &���    A�CT�����    $   �  �&���    A�CT�����    $   �  p'���'   S�CO����           HO��               4  PO���    D0    $   L  �O��   A�CT�����    $   t  �P��   A�CT�����    $   �  �Q��   A�CT�����    $   �  �R��   A�CT�����    $   �  �S��   A�CT�����    $     �T��   A�CT�����    $   <  �U���%   S�CO����         d  `{��               |  h{���    D0    $   �  �{���    A�CR����      $   �  �|���    A�CR����      $   �   }��   A�CT�����    $     ~��   A�CT�����    $   4  ��   A�CT�����    $   \  ���   A�CT�����    $   �   ���K;   S�CO����         �  X���           $   �  P����    J��f@����               p@      @             ��������        ��������                                     �              �              �                           �             �             �             �                          )             H             �             �             �             �                                       &             J                          @            �dD             �d     !                     �d                          X@            �	@            �@     
       �                                          @�d            8                           �@            �@                  	              ���o    �@     ���o           ���o    *@                                                                                                             (�d                     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@     @     @     &@     6@     F@     V@     f@     v@     �@     �@     �@     �@     �@     �@     �@     �@                                                                             @       h                                                                                 @       h                                                                                 @       h                                                                                 @       h                                                                                 @       h                                                         @       h                                                                                 @       h                                                                                 @       h                                                                                 @       h                                                                                 @       h                                                         @       h                                                         @       h                                                         @       h                                                         @       h                                                         @       h                                                         @       h                                                                                                                                                                                                                                                                                                                 �            �  ��d     $                               ------------------------------------                                                                                                                                                                                                                                                                                            �            �   �d     $                               ----- FIRN DENSIFICATION MODEL -----                                                                                                                                                                                                                                                                                            �            �  ��d     $                               ------------------------------------                                                                                                                                                                                                                                                                                            �            �   �d                                    TEST                                                                                                                                                                                                                                                                                                                            �            �  ��d                                    Read all averaged values                                                                                                                                                                                                                                                                                                        �            �   �d                                                                                                                                                                                                                                                                                                                                                                    � 	           �  ��d                        ��e                                                                     ------ Point number:                                                                                                                                                                                                                                                                                                            �            �  ��d                                                 �   �d     
                                                                                           Run for Lon:                                                    and Lat:                                                                                                                                                                                                                                                                                                                       �            �   �d                                                  �  @�d                                                                                                 Gridpoint:                                                     ,                                                                                                                                                                                                                                                                                                                               � 	           �   �d                        `�e                                                                      Number of spin up times:                                                                                                                                                                                                                                                                                                       � 	           �  ��d     #                   �e                                                                      Grounded (0) or Floating (1) ice:                                                                                                                                                                                                                                                                                              � 	           �  ��d     &                    �e                                                                      Implicit (1) or Explicit (2) scheme:                                                                                                                                                                                                                                                                                           �            �   �d     $                               ------------------------------------                                                                                                                                                                                                                                                                                            �            �  ��d                                                                                                                                                                                                                                                                                                                                                                    �            �   �d     '                               Got all variables from the NetCDF files                                                                                                                                                                                                                                                                                         �            �  ��d     '                               ---------------------------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       �            �                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        (I1)                                                            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����(I1)            `    TMF      �           4         ����(I2)            `    TMF      �           4         ����                 @                                   ��d     ��d                                                                                                                                                                                                              �                                                       @                                    �d     ��d                                                                                                                                                                                                              �                                                       @                                   @�d     0�d                                                                                                                                                                                                              �                                                       @                                   ��d     p�d                                                                                                                                                                                                              �                                                       @                                   ��d     ��d                                                                                                                                                                                                              �                                                       @                                    �d     ��d                                                                                                                                                                                                              �                                                       @                                   @�d     0�d                                                                                                                                                                                                              �                                                       @                                   ��d     p�d                                                                                                                                                                                                              �                                                       @                                   ��d     ��d                                                                                                                                                                                                              �                                                       @                                    �d     ��d                                                                                                                                                                                                              �                                                       @                                   @�d     0�d                                                                                                                                                                                                              �                                                       @                                   ��d     p�d                                                                                                                                                                                                              �                                                       @                                   ��d     ��d                                                                                                                                                                                                              �                                                       @                                    �d     ��d                                                                                                                                                                                                              �                                                       @                                   @�d     0�d                                                                                                                                                                                                              �                                                       @                                   ��d     p�d                                                                                                                                                                                                              �                                                                                                                                                                                                                                                                                                                      �                                                                                                                                                                                                                                                                                                                      �                                                  �                                                                                                                                                                                                                                                                                                                              � 	           �                            �                                                                                                                                                                                                                                                                                                                                          �            �  @�d                                    Read snowmelt...                                                                                                                                                                                                                                                                                                                �            �  ��d                                    Read precipitation...                                                                                                                                                                                                                                                                                                           �            �  @�d                                    Read snowfall...                                                                                                                                                                                                                                                                                                                �            �  ��d                                    Read sublimation...                                                                                                                                                                                                                                                                                                             �            �  @�d                                    Read skin temperature...                                                                                                                                                                                                                                                                                                        �            �  � e                                    Read snow drift...                                                                                                                                                                                                                                                                                                              �            �  @e                                    Read wind speed...                                                                                                                                                                                                                                                                                                              � 	           �   e     	                                                                                          wegPsol:                                                                                                                                                                                                                                                                                                                        � 	           �  �e     	                                                                                          wegPliq:                                                                                                                                                                                                                                                                                                                        � 	           �  �e     	                                                                                          wegPtot:                                                                                                                                                                                                                                                                                                                        � 	           �  @	e     	                                                                                          wegMelt:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 �            �  �e     (                                            �   e                                    No melt, implicit scheme!    Melt sum =                          mm yr-1                                                                                                                                                                                                                                                                                                                        �            �  �e     '                                            �   e                                    Melt, explicit scheme!      Melt sum =                           mm yr-1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   �            �                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     �            �  @e                      �                            �  �e                      �                                                                          netCDF error (                                                  ):                                                                     �e                                                                                                                                                                                                                                                                                                                      �                                                                  @e                                                                                                                                                                                                                                                                                                                      �                                                                  �e                                                                                                                                                                                                                                                                                                                      �                                                             @e                                                                                                                                                                                                                                                                                                                      �                                                             �e                                                                                                                                                                                                                                                                                                                      �                                                             @e                                                                                                                                                                                                                                                                                                                      �                                                             �e                                                                                                                                                                                                                                                                                                                      �                                                             @!e                                                                                                                                                                                                                                                                                                                      �                                                             �"e                                                                                                                                                                                                                                                                                                                      �                                                             @$e                                                                                                                                                                                                                                                                                                                      �                                                                  �%e                                                                                                                                                                                                                                                                                                                      �                                                            @'e                                                                                                                                                                                                                                                                                                                      �                                                            �(e                                                                                                                                                                                                                                                                                                                      �                                                            @*e                                                                                                                                                                                                                                                                                                                      �                                                             �+e                                                                                                                                                                                                                                                                                                                      �                                                             @-e                                                                                                                                                                                                                                                                                                                      �                                                             �.e                                                                                                                                                                                                                                                                                                                      �                                                                  @0e                                                                                                                                                                                                                                                                                                                      �                                                             �1e                                                                                                                                                                                                                                                                                                                      �                                                             @3e                                                                                                                                                                                                                                                                                                                      �                                                             �4e                                                                                                                                                                                                                                                                                                                      �                                                             @6e                                                                                                                                                                                                                                                                                                                      �                                                             �7e                                                                                                                                                                                                                                                                                                                      �                                                            @9e                                                                                                                                                                                                                                                                                                                      �                                                                  �:e                                                                                                                                                                                                                                                                                                                      �                                                            @<e                                                                                                                                                                                                                                                                                                                      �                                                                                                                                                                                                                                                                                                                     �            �   >e                                    testest                                                                                                                                                                                                                                                                                                                         �            �                                                 �@e                                                                                                                                                                                                                                                                                                                      �                                                                   Be                                                                                                                                                                                                                                                                                                                      �                                                                  �Ce                                                                                                                                                                                                                                                                                                                      �                                                              Ee                                                                                                                                                                                                                                                                                                                      �                                                             �Fe                                                                                                                                                                                                                                                                                                                      �                                                                                                                                                                                                                                                                                                                      �            �  @He                                    Read input settings                                                                                                                                                                                                                                                                                                             �                                                                                                                                                                                                                                                                                                                        �            �   Ke                                                                                                                                                                                                                                                                                                                                                                    �            �   Me                                                 �  @Me     
                                                                                           Lon:                                                            and Lat:                                                                                                                                                                                                                                                                                                                       �            �  �Ne     $                               ------------------------------------                                                                                                                                                                                                                                                                                            �            �  �Pe                                                  �   Qe                                                                                                 Closest gridpoint:                                             ,                                                                                                                                                                                                                                                                                                                               �            �   Se                                                 �  @Se     
                                                                                                    with Lon:                                              and Lat:                                                                                                                                                                                                                                                                                                                       �            �  �Te     A                               -----------------------------------------------------------------                                                                                                                                                                                                                                                                                                                               �            �  �Ve                                                                                                                                                                                                                                                                                                                                                                    �                                                                                                                                                                                                                                                                                                                      �                                                                                                                                                                                                                                                                                                                      �                                                                                                                                                                                                                                                                                                                        �                                                                                                                                                                                                                                                                                                                        �                                                        (A19,1X,F8.3,1X,A1)                                             `    TMF                                     	                            4         ����                                                                                                                                                                                                                                                                                                                                (A19,1X,F8.3,1X,A13)                                                            `    TMF                                     	                            4         ����                                                                                                                                                                                                                                                                                                                                (A19,1X,F8.3,1X,A5)                                                             `    TMF                                     	                            4         ����                                                                                                                                                                                                                                                                                                                                (A19,1X,F8.3,1X,A13)                                                            `    TMF                                     	                            4         ����                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                �                                                                                                                                                                                                                                                                                                                     �                                                                                                                                                                                                                                                                                                                     �                                                                                                                                                                                                                                                                                                                     �            �   je                                                                                                                                                                                                                                                                                                                                                                    �            �  �ke     $                               ------------------------------------                                                                 @]e      ]e                                                                                                                                                                                                              �            �  @me                                                 �  �me                                     Average Tsurf:                                                 K                                                                                                    @_e     �^e                                                                                                                                                                                                              �            �  @oe                                                 �  �oe                                     Average Acc:                                                   mm w.e. yr-1                                                                                         @ae     �`e                                                                                                                                                                                                              �            �  @qe                                                 �  �qe                                     Average Wind:                                                  m s-1                                                                                                @ce     �be                                                                                                                                                                                                              �            �  @se                                                 �  �se                                     Total Melt:                                                    mm w.e. yr-1                                                                                                                                                                                                                                                                                                                    �            �   ue     $                               ------------------------------------                                                                                                                                                                                                                                                                                            �            �  �ve                                                                                                                                                                                                                                                                                                                                                                    �            �   xe                                                                                                                                                                                                                                                                                                                                                                    �            �  �ye                                    Output variables                                                                                                                                                                                                                                                                                                                �            �   {e                                    Prof, Speed, Detail                                                                                                                                                                                                                                                                                                             � 
           �  �|e     
                                                                                             writein...                                                                                                                                                                                                                                                                                                                      � 
           �  �~e                                                                                                  numOutput...                                                                                                                                                                                                                                                                                                                    � 
           �  @�e     	                                                                                             output...                                                                                                                                                                                                                                                                                                                       �            �  ��e                                                                                                                                                                                                                                                                                                                                                                    �            �  ��e                                                                                                                                                        After spin-up #                                                                                                                                                                                                                                                                                                                 �            �  ��e                                                                                                                                                        After spin-up #                                                                                                                                                                                                                                                                                                                 �            �  ��e                                                                                                                                                        After spin-up #                                                                                                                                                                                                                                                                                                                 �            �  ��e                                                                                                                                                        After spin-up #                                                                                                                                                                                                                                                                                                                 �            �  ��e                                                                                                                                                        After spin-up #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 8��                                                          �|                                                                                                                                                                                                                                                                                                                            �                                                                                                                                                                                                                                                                                                                     �                                                                                                                                                                                                                                                                                                                        �            �  @�e                                    Load from ini-file                                                                                                                                                                                                                                                                                                              �            �  ��e                                    Start of time loop                                                                                                                                                                                                                                                                                                              � 	                                                                                                                                                                                                                                                                                                                                                                                �            �  ��e                                    End of Time Loop                                                                                                                                                                                                                                                                                                                �            �  @�e                                                                                                                                                                                                                                                                                                                                                                    �            �                                                                                                                                                                                                                                                                                                          �            �   �e                                    Written output data to files                                                                                                                                                                                                                                                                                                    �            �  ��e                                                                                                                                                                                                                                                                                                                                                                    �            �   �e                                    Cleared all variables                                                                                                                                                                                                                                                                                                           �            �  ��e                                                                                                                                                                                                                                                                                                                                                                    �            �   �e     #                               ___________________________________     GCC: (SUSE Linux) 4.3.4 [gcc-4_3-branch revision 152973] GCC: (GNU) 4.8.1 20130531 (Cray Inc.) Cray Fortran : Version 8.5.8 GCC: (GNU) 4.9.0 20140422 (Cray Inc.)       ,             @     *                       ,            @                            ,    `       �@     �                      ,    �       �3@     cJ                     ,    %       @~A     �j                      ,    �        �A     �                     ,    �       ��B     �:                     ,    6       �(D     K;                      ,    �        dD     �                       !    w   �   o   _IO_stdin_used     6    �  0  M   __libc_csu_fini i   __libc_csu_init     s            @     >@     ../sysdeps/x86_64/elf/start.S /usr/src/packages/BUILD/glibc-2.11.3/csu GNU AS 2.23.1 ��       ~   
      @@     @@     a   g   ;   .   )   i   N   int u   p   X   �   	 eD     W    Z    _   �   `  �              �   �  < @            �V   �   �V        .    d    �   �   �@     �3@     mainprogram_firnmodel.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �Y    �   �  �3@     #~A     openNetCDF.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �X    �   t  @~A     �A     ini_model.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �X    �   �    �A     ��B     time_loop.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �U    �   �,  ��B     ~(D     output.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �Y    	  �A  �(D     �cD     subprogram.f90 /perm/ms/nl/rumb/code/Program GNU AS 2.26.0 �,     ~          dD     �dD     �D  u   J  �?   )   int :  � dD     dD     w�  jdD     �dD     ?   �   �  iF   x   Q  i�   �   �  i�   �   $  z�   	    
i {4       �   �   p   4     �      F   �   �    �     ;�   )  =�     %   %  $ >  $ >  $ >  4 :;I?
  & I   %U  .:;'@�B  4 :;I  $ >    %    %    %    %    %    %   %  $ >   :;I  $ >  . ?:;'@
  .?:;'@   :;I  4 :;I  	U  
4 :;I   I  & I  I  !   '   I  4 :;I?<   ]    4   �      ../sysdeps/x86_64/elf  start.S     	@     � .>!>L$ uvx[ #       �       init.c     m    O   �      ../../../cray-gcc-4.8.1/libgcc/config/i386  crtfastmath.c     	 @     <� Z �   0   �       mainprogram_firnmodel.f90      	�@     "!!(]�$6\�%6[�&6Z�'6Y�)6#Ti,�T.0XPJ0HP�1�DLq4=�K?5?�J?68/I?78/H?88/G?98/F?:?�EB;80Cz=6/B[>8/A[?8/@[� 8/�[� ?��[48t��[� 8��_� ?[��� �!$ba�]qwrg�$!��!!!!&��� �!&�h��2�-��� �XvXe���t �   %   �       openNetCDF.f90      	�3@     	�w��qX�#������6XG(#XX�pe(G�me(G:je(G.ge(G:de(G.ae(G:^e"(Gt<;[e%(Gn��#���������~G�#�~<S�-(S�.jR�//Q#/-QX0�P.0#PX1tO#1-OX2�N52#NX3tM#3-MX4�L54#LX7t.�WPH58#HX9��F5:#FX;��D5<#DX� @�� (@�� j��� /�#� ���� ��X� ��.� #�X� t��WP�5� #�X� ��'� #�X� tK�+� (��� (��� j��� /�#� -�X� ��2� #�X� t�#� -�X� ��5� #�X� t�#� -�X� ��5� #�X� t-�WP�5� #�X� ���5� #�X� ���5� #�X� t�%� -�X� ��'� #�X� t�WP�5� #�X� ��+� #�X� ���� (��� j��� /�/� -�X� ��$� #�X� t�#� -�X� ��5� #�X� t�#� -�X� ��5� #�X� ���WP�5� #�X� ���5� #�X� ���5� #�X�t��~��(�~��j�~��/�~%�-�~X���~$�#�~X�t�~#�-�~X���~5�#�~X�t�~#�-�~X���~5�#�~X����WP�~5�#�~X����~5�#�~X����~5�#�~X�t��~��(�~��j�~��/�~.�=�~.�#�~X�t�~#�-�~X���~5�#�~X�t�~#�-�~X���~5�#�~X����WP�~5�#�~X����~5�#�~X����~5�#�~X�t��~��(�~��j�~��/�~#�-�~X���~.�#�~X�t�~#�-�~X���~5�#�~X�t�~#�-�~X���~5�#�~X����WP�~5�#�~X����~5�#�~X����~5�#�~X�t��� [�~��(�~��j�~��/�~#�-�~X���~2�#�~X�t�~#�-�~X���~'�#�~X�t�~#�-�~X���~5�#�~X�t�~#�-�~X���~5�#�~X�t-�WP�~5�#�~X����~5�#�~X����~5�#�~X�t�WY�~��(�~��j�~��/�~#�-�~X���~2�#�~X�t�~#�-�~X���~'�#�~X�t�~#�-�~X���~5�#�~X�t�~#�-�~X���~5�#�~X�t-�WP�~5�#�~X����~5�#�~X����~5�#�~X�t�W	T�~��,�~��E���@Z�~V�(��~Z�q�~��<�~��,�~j�j�~��<�~��,�~j�j�}��<�}��,�}j�j�}��<�}��,�}j�j�}��<�}��,�}j�j�}��2�}��@�}@�#�}X���}#���}����}X���}'�#�}X���}����}����}X���}'�#�}X���}����}����}X���}'�#�}X���}.���}����}'�#�}X���}����}����}X���}'�#�}X��'�WS�}+�#�}X��'�WS�})�#�}X��'�WS�})�#�}X��'�WS�}+�#�}X��'�WS�})�#�}X��'�WS�}+�#�}X�<��}'�#�}X��u�}'�#�}X��u�}'�#�}X��u�}'�#�}X��u�}'�#�}X��u�}'�#�}X���.)���0$��#��#��}6	������	����	�	����	����	����	�������������������������(��qt�2-'�#��'�N� #�X�*��m�+ol�+ahj(Vgj(U"H0���a�+o`� +a\j$(V[j%(��V0���U�++hT�,+aPp0(VOd1(L�V0���I�7+hH�8+aDj<(VCj=(h�� �e��� +(�W� o���� ,�j� j�%� 9��� ,�j� j�%� 9��� ,�j� j�%� 9��� ,�j� j�%�9���,�j�j�~%�9�~��,�~j�j�~%�9�~��,�~j�j�~%�/�~��=�~V�#�~X�t�~#���~����~X���~:�#�~X�t�~����~����~:�#�~X�t�~.���~����~:�#�~X�t�~����~����~X���~:�#�~X�t�~����~����~X���~:�#�~X�t�~����~����~X���~:�#�~X�t,-�WP�~<�#�~X��,-�WP�~<�#�~X��,-�WP�~<�#�~X��,-�WP�~<�#�~X��,1�WP�~<�#�~X��,-�WP�~<�#�~X��,1�WP�~<�#�~X����~.�#�~X�tK�~.�#�~X�tK�~.�#�~X�tK�~.�#�~X�tK�~.�#�~X�tK�~.�#�~X�tK�~.�#�~X�t����7��gKK!!E��׻�z��׻���z��gegg��z��g�gWgY�w.�q9��"?ؽ��lfX��cijifuh�p�h�2Ĥ��30W�yf�cwui�wht#L�������p�XUijf%����ijjf�y�ot���m��q.������xfu���h��uu`�~�/������ +h��� +a�p� (V�d� (E$������ +h��� +a�j� (V�j� (($������ +h��� +a�p� (V�d� (E"H0������ +h��� +a�j� (V�j� (*�������������� �
fuNV�%I�I(I#I��6H�dY&K�O�ug��ug�Q$�qe0qZ(PpJfpffp00Yp�fp1(-�[n�>��z��z9^z:^k\rj%/g.=f:#fXte.=dA#dXtc.=b:#bXta�-aX �`: #`X!t_#!-_X"�^:"#^X%t[�%-[X&�Z:&#ZX(IX<(#XX+t�T<,#TX-t�R<.#RX/t�P<0#PX1t�N<2#NX3t�L<4#LX7�H.8#HX<tH!׭�/�H�4�V-������Pu�#t�Z�s<� h   $   �       ini_model.f90      	@~A     qsg(GIY�YYr�,�r_(�r_(rQ*�!"SSSSSPS!SSPPPP!PPPPPP�`\Zy"(Ny|%yn(yh(yAE����Y,gxN%Y'gx9(Y0gx0(xM*��!SJD"r(������4�t�5�$��=�>�)�>/&=�=��=�=��=�=��=�=�9,9��=�=y7�$�K�!��!""|fk<�0�'&%#u,g$sK�rL�qM�pN�oO��pC�����/��z%p733//%v.======z������Ku�cY�K�Kuz<#�///�E�X�����/q����KK�Kz%"��!!6G'�ggg"sUH4�����YYYYyX�ztrCYYgYggzfgguguuztgguguuztgguguuxtOgYgYYgxf>&FFFFFz1H>/�u���x��=��u/!����y��h	��m�JWx�H�Es"e)W&�Yj�*���Os�xUg�)��>��le)W.�Y�%(��dMrt!;��L�]T���+	/xaag�"!!����!"�%�[�����C\YKK�yt;f���d#�"xtj�gggxff��f��%����y<dJ"f���xt��ggg����yX
�v1&�o3�o�E4�$^z��z��'y84z��O��cXVl:4�$^z�Bz��'y84z��Rx�c�{fX�d�0�d�Jwct�JwcfY��c��Gx#lz���wJ	y�%z��OvMc(�I7zf	�y�'z��Oo�(x#lz�"�wJ	y�%z��RvMc�#�F7zf	�y�'z��R��d��X�2%�Wgegege7WLW��C����������������uu��g�=/�#!!���#,uu��vuM �   $   �       time_loop.f90      	 �A     ,��Lc<4�M
CxJ�s�Xs-h
Js��v-*W
�v�)�
&v�
 sJ^�Rs��$�s��
�n�@�X�$��v�6
�
�#Wueueue4WMW��6:@�f�+�LV>9+8Yh�d<������z��S1!��?�+/NV>c'AY��s�JS5!�C��SL'�[�+/GV>9:Y��ISD-��s<tJ_5!���S�/t'`x�<p.3z�Jl��-My�"@01-��w�z��	�y��j4%�}۶�v�)lz�&"�;2-��v+Hlz�zJ�&�z;�$"��vJ&dS6gyJ������Mu?J��s�6�1�5u9r��4J$Bg"?�';�'63UI�x1<,)T<��$�0y;)z�$t����%/��(j [m�,j5wJ	�yt�y��z��t���"Y'u��jfKk��<iwJ	�y�&z�Pt��Z�-uןi<�H<jvJ
�x�&y��s����-u/��`�g3�9Y@$$<v.3mf#�O��I7�2)9��=�.�=3�h'�s4�/�g�0�G�K�.�)����aX,$q�V$�aS^�)�q�\0�^tDaXaX"$q�V$�aS^�)�q�\0�^tDYXaX*$q�V$�aS^�)�q�\0�^rW�aX9$q�V$�aS^�)�q�\0�^rWXaX0$q�V$�aS^�)�q�\0�^r^�aX2$q�V$�aS^�)�q�\0�^rUXaX2$q�V$�aS^�)�q�V0�[r`XaX3$q�V$�aS^�)�q�V0�[raXaX2$q�V$�aS^�)�q�V-�Wh^tXf(tf�7^tXf(tC^tXf(tC^tXf(tC^tXf(tf��C^tXf(tf��C^tXf(tC^tXf(tf��r���w9^t f(tf<�r���04�����ho�HI�K Y�CtsJJ��Ih,�eHv�vh�hIh�ἝhhҒ�hIh�'Y�$���i	�w�H�MUw>;��k(.�}�sH�׭�/�H�8�-y �?�9JDJLJJfJ� .����vt�N0sJ��K�JhNp����Y��Np����yCKb��Iwqeu\���#�xtjb\x���LKY�m��Kz+k</)���K�J_r Kk��tJ�t"S�$*�/��h�lzJ	XUzfK#����e�w��M�p��
tyfF/!��sgg�ջ��F��"dbm)mx�{�of�&eg��h_yXY	fwJ^zf�G�c��eYw����u�z��)W�n�����t)���'�
����.ɻ�Y���u��O�]�P�Pz�Qr�+�˒YK�wX�	Xg�$
JvZ����GKh�h���%���uhYgzfguvguztguvguztguvguyt�uhYgyf�uhYgyfK�����z�4-v�x�y�z)Q�YgugSgguuuagguuuagguuuzf$�YgugzX��׻zX���������v�w�x�y�{�gYugaggguuoggguuoggguuzt$�gYugzf�׻�׻zf�������x ��v�w�x�y*��K!!��-uuughM�/��dg�K�"������גM���YYYY��y�����xg�������x���������x򑻻������h)������y�ugggguuyt��������* y�k���#��y<Y�Y�=y�B/////=y.��y���������_�]�Y���YYYy%��/����y�uuuuuugxf� �   !   �       output.f90      	��B     �re(Y�\�\j\n#fn:(�n%fn1(n�l'#lXth�-h��g5#gXtc��b5#bXta� �`5 #`X!t_�"�^5"#^X#t]($�\5$#\X%t[*&�Z5&#ZX'tY�(�X5(#XX+=T',#TX/t�P50#PX1t�N52#NX3t�L54#LX5t�J56#JX7t�H58#HX<t�C*=#CX� ���� /O.�������
��	re(V�\�\e�n�fn:(�n%fn1(n��l.#lXth�-hf g<#gXtf�Ee<#eXta� �`< #`X!t_�"�^<"#^X#t]�$�\<$#\X%t[(&�Z<&#ZX'tY*(�X<(#XX)tW�*�V<*#VX+tU�,�T<,#TX-tS�.�R<.#RX/tQ�0�P<0#PX2tN*3�M.3#MX4tL�5�K.5#KX8=G.9#GX<t�C<=#CX>t�A<?#AX� t��<� #�X� t��<� #�X� t��<� #�X� t��<� #�X� ���<� #�X� ���<� #�X� t��<� #�X� ��.� #�X� ��.� #�X� ���� /�2�@������������tZ���Y��������"�"I*�Lu���uuuuuuu�w�X؟(��������.5o<�j<�uX<u<�Lggquuuquuuquuuq9ggq&ggqY����v�w��vL&$&u�#�K&u��u�&��	<h��@���$�K�yt�Ch�������%�sh(GMQ�l�(lk(l�}h&-lf�g<#gX�d�md��cC#cX�b�mb��a<#aX �`� m`�!�_<!#_X"�^�"m^�#�]<##]X$�\�$m\�%�[<%#[X&�Z�&mZ�'�Y<'#YX(�X�(mX�)�W<)#WX*�V�*mV�+�U<+#UX,�T�,mT�-�S<-#SX.�R�.mR�/�Q./#QX0�P�0mP�1�O.1#OX2�N�2mN�3�M.3#MX4�L�4mL�5�K.5#KX6�J�6mJ�7�I.7#IX8�H�8mH�9�G.9#GX:�F�:mF�;�E.;#EX<�D�<mD�=�C.=#CX>�B�>mB�?�A.?#AX� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��<� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ��6� ���� ��>� #�X� ���.� #�X� ��y�S�0� #�X� ����S�<� #�X� ����S�0� #�X� ����S�0� #�X� ����S�0� #�X� ����S�0� #�X� ����S�0�#�X�����S�~0�#�~X�����S�~0�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X�����S�~2�#�~X��u�~%�2�~0� ���������������������qh(GMj�j�(je(j�zf&-jf�e<#eX�d��d��d�c<#cX �`� m`�"�^<"#^X#�]�#m]�%�[<%#[X&�Z�&mZ�(�X<(#XX)�W�)mW�+�U<+#UX,�T�,mT�.�R<.#RX/�Q�/mQ�1�O<1#OX4�L64�L�5�K<5#KX6�J66�J�7�I<7#IX8�H68�H�9�G<9#GX:�F6:�F�;�E<;#EX<�D6<�D�=�C<=#CX>�B6>�B�?�A<?#AX� ���.� #�X� �I�WS�$� #�X� �I�WS�$� #�X� ���WS�$� #�X� �I�WS�$� #�X� �I�WS�$� #�X� ���WS�$� #�X� �u��� 2�7,���������d�2or(��4M^>�=�K�-5Gv��$Y���%e=���a�(ae(a�#z]'$A\<$#\X%�[�%�[�%�[�&ZC&#ZX)�W�)mW�+�U<+#UX,�T�,mT�.�R<.#RX/�Q�/mQ�1�O<1#OX2�N�2mN�3�M<3#MX4�L�4mL�5�K<5#KX6�J�6mJ�8�H<8#HX;�E6;�E�<�D<<#DX=�C6=�C�>�B<>#BX?�A6?�A�� �@<� #@X� ��6� ���� ��<� #�X� ���.� #�X� ���:� #�X� ���:� #�X� �;��S�2� #�X� �W��S�2� #�X� �W��S�2� #�X� ����S�2� #�X� �u��� 2�3&�������� i   %   �       subprogram.f90      	�(D     %Lx<42d�2a�2�,XW,aX"*w��rL+تL��H��LH�LHL�gJ%~k�تw�x�yߤz�{y�|x��gYguYggxfguguuguuxtguguuguuxtguguuguuxt}gYYugggxfcg����g�x�u�������OX5��`F�=����Y�}� �K�T� ����@� ���� 4�{��J� 4��[� A��\� 7��[� A��\� 7�zD��t� 7��~� 4�xD��t�7x�	��~\�7��~g�7��Xuuuuo9gggYSuuuugauuuugauuuuga,gggYS1�����Yuuuuџ��uuo(su�rv�qw�px�oy�Xuuq1gYUuugcuugcuugc"gYU6�Yuu��uuq%su�rv�qw#7XJ<P�~F�W�	�/$�����>I�m��==:��i=G�6g���/���C�K$�~,�/��~t���~w���~_���~w���[�!!�~ c��'�~��X �    Z   �      /usr/lib64/gcc/x86_64-suse-linux/4.3/include  elf-init.c    stddef.h     	 dD     �^��o�to<foJ<[�ǒ#          ���� x�              dD            $       dD     �       J��f@����/usr/src/packages/BUILD/glibc-2.11.3/csu long unsigned int short unsigned int short int _IO_stdin_used unsigned char long int GNU C 4.3.4 [gcc-4_3-branch revision 152973] mxcsr GNU C 4.8.1 20130531 (Cray Inc.) -mlong-double-80 -mlong-double-80 -mfxsr -msse -mtune=generic -march=x86-64 -g -g -g -O2 -O2 -O2 -fbuilding-libgcc -fno-stack-protector -fpic ../../../cray-gcc-4.8.1/libgcc/config/i386/crtfastmath.c /tmp/peint/cray-gcc/4.8.1/BUILD/snos_objdir/x86_64-suse-linux/libgcc set_fast_math envp argc __libc_csu_init elf-init.c __init_array_start size __init_array_end __libc_csu_fini size_t argv @     @     	 �t�
@�!�@     @      �t                       @        w@       �        w�                        U        UU       �        ^                       U        TU       �        ]                       U        QU       �        \                \       {        S                 @     @                            (       U       v       I       M       @       C                                                           @                   @                   <@                   X@                   �@                   �	@                   *@                   �@                  	 �@                  
 �@                   @                    @                    @                   �dD                    eD                   ��D                   ��D                   �d                   �d                    �d                   �d                    �d                   (�d                   8�d                   @�d                   ��d                   @�e                                                                                                                                 !                      "                      #                      $                      %                 �!   ��                     ��                    ��                     @@                   �dD             '    ��                2      �d             @     �d             N      �d             [     `@             ]     �@             p     �@             �     @�e            �     H�e            �     @@             '    ��                �     �d             �     �D             �      �d             �     �dD                 ��                �    ��                �      @               ��                    �d               ��                #    @�d            4    @�d     (       E    ��d            V    ��d     (       g    @�d            x    @�d     (       �    ��e     �      �    ��d            �    ��d     (       �    @�d            �    @�d     (       �    ��d            �    ��d     (           @�d                @�d     H       +     �d     �       =     �d            O    @�d     �       a    @�d            s    ��d            �    ��d     H       �    @�d            �    @�d     H       �     �d            �     �d     H       �    ��d            �    ��d     (           @�d                @�d     (       '    ��d            9    ��d     (       K    @�d            ]    @�d     (       o    ��d     (       �     �d     (       �    ��d     (       �     �d            �    ��d            �     �d            �    ��d            �    ��d            �     �d                 �d                 @�d            2     �d             D    ��d     (       V    ��d     (       h     �d     (       z    ��d            �     �d     (       �    ��d     (          ��                �    �eD            �    fD            �    efD            �    �fD            �    gD                dgD                �gD            +    hD            =    hhD            O    �eD            `    �eD            q     fD            �     fD             �    @fD            �    PfD            �    pfD            �    �fD            �    �fD            �    �fD     "       �    �fD                 gD                  gD            2    @gD            D    PgD            V    �gD             h    �gD            z    �gD            �    �gD            �    �gD            �     hD            �     hD             �    @hD            �    PhD            �    whD            
    QiD                )jD            .    �jD            @    �kD            R    lD            d    �lD            w    >mD            �    �mD            �    �hD            �    �hD            �    �hD            �    �hD            �    �hD            �    �hD            	    �hD                �hD            -    �hD            ?    �hD            Q     iD            c    iD            u    iD            �    (iD            �    8iD     	       �    HiD     	       �    `iD            �    {iD            �    �iD            �    �iD            	    �iD            	    �iD            )	    �iD            ;	    �iD            M	    �iD            _	    �iD            q	    �iD     
       �	     jD            �	    jD            �	     jD     	       �	    @jD     #       �	    cjD            �	    pjD            �	    �jD            
    �jD            
    �jD            %
    �jD            7
    �jD            I
    �jD            [
    �jD            m
    �jD            
    kD     
       �
     kD            �
    0kD            �
    8kD            �
    FkD            �
    PkD            �
    `kD            �
    pkD                �kD            !    �kD            3    �kD            E    �kD            W    �kD            i    �kD            {    �kD            �    �kD            �    �kD            �    lD            �     lD            �    1lD            �    8lD            �    DlD                HlD                 VlD            3    `lD            F    plD            Y    �lD            l    �lD            #    @�e                �lD            �    �lD            �    �lD            �    �lD            �    �lD            �    �lD            �    �lD                �lD                 mD            *    mD            =     mD            P    0mD            c    PmD            v    cmD            �    hmD            �    tmD            �    xmD            �    �mD            �    �mD            �    �mD            �    �mD                �mD            !    �mD            4    �mD            G    �mD            o    ��d     (       4    ��d            Z    �mD            m     nD            �    nD            �    nD            �    nD            �    nD            �     nD            �    ,nD            �    8nD                DnD                PnD            +    \nD            >    hnD            Q    xnD            d    �nD            w    �nD            �    �nD            �    �nD            �    �nD            �     oD            �     oD            �    @oD            �    XoD     	           hoD     	       "    xoD     	       5    �oD     	       H    �oD     	       [    �oD     	       E    ��d            V     �d            �    @�d            g    ��d            x    ��d            �     �d            �    @�d            �    ��d            �    ��d            �     �d            �    @�d            �    ��d            �    ��d            �     �d            �    @�d                ��d                ��d            �     �d            =    @�d            +    ��d            �    ��d            �     �d            a    @�d            O    ��d                ��d                  �d            s    @�d            �    ��d            2    ��d            �     �d            �    @�d            D    ��d            �    ��d            �     �d            n    �oD            �    pD            �    `pD            �    �pD            �    qD            �    gqD            �    �qD            �    rD                prD            V    ��d            �    ��d     (       �     �d            h     �d     (           �oD            ,    �oD            ?    �oD            R    �oD                @�d                @�d     (       z    ��d            '    ��d     (       e    pD            x    0pD            �    ?pD            �    PpD            9    ��d            �    ��d     (       K     �d            ]     �d     (       �    ppD            �    �pD            �    �pD            �    �pD            �    @�d            �    @�d     (           ��d            !    ��d     (       3    �pD            F    �pD            Y    �pD            l     qD                ��d     (       �    ��d            �     �d     p       �     �d            �    ��d     H       �    ��d            �    �rD            �    �rD                �rD            $    �rD            7    �rD            J    �rD            ]    �rD            p    �rD            �    �rD            �    �rD            �    �rD            �    �rD            �    �rD            �     sD            �    sD                sD                %sD            .    0sD            A    @sD            T     �d            f     �d     (       x    PsD            �    ��d            �    ��d     (       �    `sD            �     �d            �     �d     (       �    psD            �    ��d                ��d     (           �sD            0     �d            B     �d     (       T    �sD            g    ��d            y    � e     (       �    �sD            �     e            �     e     (       �    �sD     	       �    �sD     	       �    �sD     	       �    �sD     	           �sD     	       !     tD     	       4    tD     	       G    �e     H       Y    �e            k    @e     H       ~    @e            �     e     H       �     e            �    �e     H       �    �e            �    ��d            �    ��d     (             �d                 �d     (       $     qD            7    @qD            J    OqD            ]    XqD            p    @�d            �    @�d     (       �    ��d            �    ��d     (       �    pqD            �    �qD            �    �qD            �    �qD                ��d                ��d     (       (     �d            :     �d     (       L    �qD            _    �qD            r    �qD            �    rD            �    @�d            �    @�d     (       �    ��d            �    ��d     (       �     rD            �    @rD                OrD                `rD            ,    �	e            ?    �	e            R     
e            e    @
e            x    �
e            �    �
e            �     e            �    @e     h       �    @e            �    @e     h       �    @e            �    ��e                 tD            #    9tD            6    @tD            I    OtD            \    TtD            o    �e     (       �    �e            �    XtD            �    htD            �    xtD            �    �tD            �    �tD            �    �tD                �tD                �tD            -    �tD            @    �tD            S    �tD            f    �tD            y    �tD            �    uD            �    uD            �    (uD     	       �     e            �    @e            �    �e            �    �e                 e            $    �e     �       7    �e            J    @�d            \    ��d            n    @�d            �    ��d            �    @�d            �    � e            �    @e            �     e            �    �e            �    �e                 @	e                �e     (       &     e            9    �e     (       L     e            _    @e            r    �e               ��                �    �wD            O    �wD            `    �wD            q    �wD            �    �wD            #    �e            o     e            E    @e            �    �e            �     e     (       g    �e            �    �e     (       �    @e            �     e     (       �    �e            �    �e     (       �    @e            �      e     (           �e            �    �!e     (       =    @ e            O     #e     (       �    �!e                @#e            s    �$e            �     &e     (       �    �$e            �    �'e     (       �    @&e            �     )e     (       �    �'e                �*e     (       �    @)e            '     ,e     (           �*e            K    �-e     (       9    @,e            ]    �-e            �     /e            �    �0e     (           @/e                 2e     (       �    �0e            �    �3e     (       p    @2e                 5e     (       �    �3e            �    �6e     (       (    @5e            �     8e     (       �    �6e            �    @8e            �    �9e            �     ;e     (       �    �9e            J    �<e     (       T    @;e            �    �<e            �    �=e     (       �    �wD            �    xD            �    *xD            �    @?e     (       �    @>e            �    0xD            �    @xD            �    ZxD            n    �?e                �@e            �     Ae            B    @Be            y    �Ce     (       �    �Be            �    @Ee     (       �     De            G    �Fe     (       �    �Ee            �     Ge            ~     He     (       �    �Ie     0       �    �He            �    �Ie            �    �Je     (            @Le     �       �    @Ke            R    �Me            e    �Ne     (       �     Pe     �       �     Oe                @Re     �       �    @Qe            �    �Se            9    �Te     (       �    @Ue            o    @Ve     (       �    kxD            �    qxD            �    ^xD            �    dxD            �    �We     (       �    �Ve                 Ye     (       �     Xe            $    @Ze     0       7    @Ye            r    �[e     0       _    �Ze            �    �\e     0       �    �[e            �     fe     (       �     ee            �    @ge     (       �    @fe            �    �he     (       
    �ge                �he            0    �ie     (       C    @je            V    @ke     (       i    �le     h       |    �ke            �    �ne     h       �    �me            �    �pe     h       �    �oe            �    �re     h       �    �qe                �se                �te     (       '    @ue            :    @ve     (       M    ��e            `    �ve            s    �we     (       �    @xe            �    @ye     (       �    �ye            �    �ze     (       �    @|e     P       �    @{e            �     ~e     P             }e                 �e     P       1     �~e            D     ��e            W     ��e     (       4    �e            V    @e            x    �e            �    @e            �    �e            �    @e                �e            +    @!e            a    �"e                 @$e            2    �%e            D    @'e            V    �(e            h    @*e            z    �+e            �    @-e            �    �.e            !    @0e                 �1e            �    @3e                �4e            :    @6e            �    �7e                @9e            �    �:e            f    @<e            \     >e            �    �@e            0     Be            g    �Ce            �     Ee            Y    �Fe            k    @He            �     Ke            ,     Me            ?    @Me            x    �Ne     (       �    �Pe            �     Qe            &     Se            �    @Se            L    �Te     H       �    �Ve            j      je            }     �ke     (       �     @me            �     �me            �     @oe            �     �oe            �     @qe            �     �qe            !    @se            !    �se            (!     ue     (       ;!    �ve            N!     xe            a!    �ye            t!     {e            �!    �|e            �!    �~e            �!    @�e            �!    ��e               ��                �    `zD            O    fzD            �    zzD            �    �zD            �    �zD            �    �zD            `    mzD            q    szD            #     �e            o     �e     �       4     �e            �     �e     �       V     �e            �     �e     �       �     �e            �     �e     �       x     �e            �     �e     �       �     �e            �    �zD            �    �zD            �    �zD            �    �zD            �    �zD                �zD            �    �zD            �    �zD            �    @�e            E    ��e            g    ��e            �    ��e            �    ��e            �    ��e               ��                �    �zD            O    �zD            `    �zD            q    �zD            �     {D     	       �    	{D            �    {D            �     {D            �    0{D            �    @{D            �    K{D            �    P{D            �    [{D            �    `{D            �    p{D                �{D     	       �    �{D                 �{D            2    �{D            D    �{D                �{D            V    �{D            #    @�e            4    ��e            o    ��e            E     �e            V    @�e            �    ��e            h    �{D            z    |D                |D     	       �    !|D            �    (|D     	       �    1|D            +    8|D            �    H|D            �    X|D            �    h|D            =    x|D            �    �|D            �    �|D            �    �|D            �    �|D            �    �|D            �    �|D            	    �|D                �|D            -    �|D            ?    �|D            Q    �|D            c    }D     	       u    }D            �    (}D            �    8}D            �    H}D            
    X}D            �    h}D            �    x}D            �    �}D            �    �}D            	    �}D            	    �}D            g     �e            x    @�e            �    ��e            �    ��e            �     �e            �    @�e            �    ��e            �    ��e            �     �e            �    @�e            �    ��e            �    ��e                ��e            ;	    �}D            M	    �}D                 �e            _	    �}D            q	    �}D            �	     ~D            �	    ~D                 ~D            �	    0~D            �	    ;~D            �	    @~D            �	    K~D            
    P~D            
    [~D            %
    `~D            7
    p~D            I
    �~D            [
    �~D            .    �~D            m
    �~D            
    �~D            �
    �~D            �
    �~D            �
    �~D            �
    �~D            �
    �~D            �
     D            �
    D                D            @    (D            !    4D            3    @D            E    PD            W    `D            i    pD            {    �D            �    �D            �    �D            �    �D            R    �D            �     �D            �    �D            �    0�D            �    H�D                `�D                 x�D            3    ��D            F    ��D            Y    ��D            l    ؀D            d    ��D                �D            �     �D            �    8�D            �    P�D            �    h�D            �    ��D            �    ��D                ��D                ȁD            *    ��D            =    ��D            P    �D            w    (�D            c    @�D            v    X�D            �    p�D            �    ��D            �    ��D            �    ��D     	       �    ЂD            �    �D            �     �D                 �D            !    @�D            4    `�D            �    ��D            G    ��D            Z    ��D            m    ��D            �     �D            �     �D            �    @�D            �    `�D            �    ��D            �    ��D            �    ��D                ��D            �    ��e            =    ��e            +     �e            �    @�e            �    ��e            a    ��e            O     �e                @�e                 ��e            s    ��e            �     �e            2    @�e            �    ��e            �    ��e            D     �e            �    @�e            �    ��e            V    ��e                 �D            +    �D            �    @�e            >    '�D            Q    +�D            d    0�D            w    ;�D            �    @�D            �    P�D            �    `�D            �    p�D            �    {�D            �    ��D            �    ��D                ��D            "    ��D            5    ��D            H    ��D            [    ؅D            n    ��D                �D            ,     �D            ?    8�D            R    P�D            �    h�D            e    ��D            x    ��D            �    ��D            �    ȆD     	       �    ��D            �    ��D            �     �D            �    �D            �     �D            �    0�D            �    ��e            h     �e                @�e                ��e            z    ��e            '     �e            3    @�D            F    X�D            9    ��e            Y    h�D     
       l    r�D            �    x�D            $    ��D            7    ��D            J    ��D            ]    ��D            �    ��D            �    ��D            �    ÇD            �    ȇD            �    ؇D            L    �D            _    ��D            r    �D            �    (�D            �    @�D            �    X�D            �    p�D                ��D                ��D                ��D     	       �    ȈD            �    ؈D                ��D            $    �D            7    0�D            J    P�D            �    ��e            K    ��e            ]     �e            �    @�e            �    ��e                ��e            �    �|D            )	    �}D            �	    ~D            �    ӇD               ��                #    ��e            V    ��e     (       E    ��e            g    ��e     0       �    ��e            x     �e            �     �e     (       4     �e            o    @�e            �    ��e            �    ��e     (       �     �e     H       �     �e            �    ��e            �    ��e     (            �e                 �e     (       +    ��e     (       =    ��e            �    ��e            �    ��e     (       O    @�e                @�e     (       s    ��e            �    ��e     (       �    @�e            �    @�e     (       �    ��e            �    ��e     (       �    @�e            �    ��e            �    ��e            �    @�e            a     �e                  ��e            2     �e            D    ��e            V     �e     (       �!   ��                    ��d            �!   ��                     ��                �!      �d             "     �d             "    @�d             *"    (�d             3"     ��D             F"    @�C           {"    ��@             �"     7C     �      �"    �@     �      �"    �zA     �       �"    �RA             #    �"D     �       <#    �C     �       r#    ��B     �       #                     �#   (�e             &    @             �#                      �#                     �#     �@     �       �#                     $    ��@     �       9$    `'D           o$                     !     �dD            �$     dD            �$    �yA     �       �$    ��@     �       �$     �B     �       
%    ��@     �       @%  !  @�e            `%    @�@             �%     �C     �       �%                     �%    ��C     �       &    ��d             &                     "&  !   �e            >&    `�@           s&     �C     �       �&    (�e             �&    ��@             �&                     �&    ��A     �'      �&    `�C           '    �/C     �       T'  !   �e            n'    �C     �       �'    ��@           �'     XA     �       (    ��@           F(    ��C     �       M%    ��C     �%      {(    @&D           �(    ��@     �       �(    �!D     �       )     �@           A)  !  ��e            `)    �2C     �       �!    H�e             �)    ��@     �       �)     �d     P      ,<    ��A     U      �)                     *                      *                     4*    `�C     �       i*    �.C     �       �*  !  @�d            �*     $D           �*     ,C     �       +     RA             6+     �C     �       k+    ��@           �+                     �+    @�e            �+    �3C     ~       �+                     /&    ��B     �      ,    ��C     �       7,    �!D             [,     ZA     _      i,     �C     �       �,    ��C     �       �,                     �,    p@     <      �,    ��C     �       -    �e           ,-    �C             P-     TA     �       �-                     �-     VA     �       �-    ��e     �       �-    @RA             .    @�@           :.    �+C     3       c,    `�B     �      _.     1C     �       �.    ��@           �.    (�e             �.    ��B     2      �.    @~A     �      �.    @#D     �       %/    @�e     H      =/    ��B     �
      G/    ��A     9       S/    ��@     �       �/    ��C     �       �/                     �/                      �/     �C     �       0    �RA             <0     %D           r0                     y0                     �0                     �0    ��C     �       �0     �B     �      �0    ��@           1                     1    `�@             W4    `|A     �      81    ��@           n1    �,C     �       �1     YA     �       �1    ��@             �1    `C             2     C     �       U2     �C           �2     �B     �       �2                     �2    ��@     �       �2    �{A     �       3    dD     �       3    ��@             :3                     �*    `�@     �d      Q3    `RA             s3                     �@     WC     �d      {3    ��A     B       �3    ��@           �3     .C     �       �3    �wA             4    ��@           J4  !  @e            c4    ��@           �4    `C     �       �4    ��C     �       3    @             5                     5    �C     �       N5    ��B     �      \5    @�C     �       �5                     N)    �C     �"      �5    ��@     �       o@    `^A     k      �5                     �5                     �5    ��C           06    @�@             Q6     �A     �      a6    ��B     �*      k6    �xA     �       �6     �A     �      �6    �C     �       �6    ��@             �6                     �6     �@           +7    `-C     �       a7                     �7     WA     �       �7    ��@           �7     �A     t�      �7     UA     �       a'    `�C     �'      *8    �1C     �       `8    �RA             �8    ��@             �8    ��@             �    �@     �      �8                     �8    ��C     �       9   ��  @             �;    �3@     e�      9    `�@           L9    ��@           �9                     �9     xA     �       �9                     �9                     :                     
:    `�@     �       ?:   �d             &     ��d             L:                     _:    @�A     0      i:    �(D     K;      u:                      �:     C     �       %>                     �:   Ȭd             �:                     �:    ��C     �       ;                     ;    �RA             @;    ��B     +      K;    ��@             m;                     s;  !  ��d            �;                     �;    ��d     0      �;    ��@     �       �;    @�C     �       !<     ]e     �      ;<    �sB     �      I<     �e     H      a<    @�@             �<    ��@           �<                     �<                     �<     �d     @      �<    ��@     �       =    ��C           F=    ��A     !      P=    `�C             t=    ��@           �=    `0C     �       �=     SA     �       >                     ->    �+C             Q>                     b>    ��@             �>                     �>    `3C     ~       �>    `4C     �      �>     HC     �      �>                     �>                     �>    ��@     �       &?    �cD             8?    ��@     �       n?     �A     �      |?                     �?    ��C     �       �?    ��B           �?     eD            �?    ��@     �       @     �@             -@    ��C           b@  !  @e            }@    ��@     �       �@     �e     @       �@    ��B     �      �@  !  ��e            �@     �@             A    ��B           A                     &A    `�C     �        initfini.c call_gmon_start _real_fini crtstuff.c __CTOR_LIST__ __DTOR_LIST__ __JCR_LIST__ deregister_tm_clones __do_global_dtors_aux completed.6298 dtor_idx.6300 frame_dummy __CTOR_END__ __FRAME_END__ __JCR_END__ __do_global_ctors_aux crtfastmath.c set_fast_math no_mmap.c p The Cpu Module __STATIC_LOCAL_0 __STATIC_LOCAL_1 __STATIC_LOCAL_3 __STATIC_LOCAL_4 __STATIC_LOCAL_6 __STATIC_LOCAL_7 $data$mainprogram_ __STATIC_LOCAL_9 __STATIC_LOCAL_10 __STATIC_LOCAL_12 __STATIC_LOCAL_13 __STATIC_LOCAL_15 __STATIC_LOCAL_16 __STATIC_LOCAL_18 __STATIC_LOCAL_19 __STATIC_LOCAL_22 __STATIC_LOCAL_21 __STATIC_LOCAL_26 __STATIC_LOCAL_25 __STATIC_LOCAL_29 __STATIC_LOCAL_30 __STATIC_LOCAL_32 __STATIC_LOCAL_33 __STATIC_LOCAL_35 __STATIC_LOCAL_36 __STATIC_LOCAL_38 __STATIC_LOCAL_39 __STATIC_LOCAL_41 __STATIC_LOCAL_42 __STATIC_LOCAL_44 __STATIC_LOCAL_45 __STATIC_LOCAL_47 __STATIC_LOCAL_48 __STATIC_LOCAL_2 __STATIC_LOCAL_5 __STATIC_LOCAL_8 __STATIC_LOCAL_11 __STATIC_LOCAL_14 __STATIC_LOCAL_17 __STATIC_LOCAL_20 __STATIC_LOCAL_23 __STATIC_LOCAL_24 __STATIC_LOCAL_27 __STATIC_LOCAL_28 __STATIC_LOCAL_31 __STATIC_LOCAL_34 __STATIC_LOCAL_37 __STATIC_LOCAL_40 __STATIC_LOCAL_43 __STATIC_LOCAL_46 __STATIC_LOCAL_49 __CONSTADDR_GV_0 __CONSTADDR_GV_4 __CONSTADDR_GV_8 __CONSTADDR_GV_12 __CONSTADDR_GV_16 __CONSTADDR_GV_20 __CONSTADDR_GV_24 __CONSTADDR_GV_28 __CONSTADDR_GV_32 __CONSTADDR_GV_1 __CONSTADDR_GV_2 __CONSTADDR_GV_3 __CONSTADDR_GV_5 __CONSTADDR_GV_6 __CONSTADDR_GV_7 __CONSTADDR_GV_9 __CONSTADDR_GV_10 __CONSTADDR_GV_11 __CONSTADDR_GV_13 __CONSTADDR_GV_14 __CONSTADDR_GV_15 __CONSTADDR_GV_17 __CONSTADDR_GV_18 __CONSTADDR_GV_19 __CONSTADDR_GV_21 __CONSTADDR_GV_22 __CONSTADDR_GV_23 __CONSTADDR_GV_25 __CONSTADDR_GV_26 __CONSTADDR_GV_27 __CONSTADDR_GV_29 __CONSTADDR_GV_30 __CONSTADDR_GV_31 __CONSTADDR_GV_33 __CONSTADDR_GV_50 __CONSTADDR_GV_65 __CONSTADDR_GV_76 __CONSTADDR_GV_87 __CONSTADDR_GV_97 __CONSTADDR_GV_108 __CONSTADDR_GV_121 __CONSTADDR_GV_134 __CONSTADDR_GV_34 __CONSTADDR_GV_35 __CONSTADDR_GV_36 __CONSTADDR_GV_37 __CONSTADDR_GV_38 __CONSTADDR_GV_39 __CONSTADDR_GV_40 __CONSTADDR_GV_41 __CONSTADDR_GV_42 __CONSTADDR_GV_43 __CONSTADDR_GV_44 __CONSTADDR_GV_45 __CONSTADDR_GV_46 __CONSTADDR_GV_47 __CONSTADDR_GV_48 __CONSTADDR_GV_49 __CONSTADDR_GV_51 __CONSTADDR_GV_52 __CONSTADDR_GV_53 __CONSTADDR_GV_54 __CONSTADDR_GV_55 __CONSTADDR_GV_56 __CONSTADDR_GV_57 __CONSTADDR_GV_58 __CONSTADDR_GV_59 __CONSTADDR_GV_60 __CONSTADDR_GV_61 __CONSTADDR_GV_62 __CONSTADDR_GV_63 __CONSTADDR_GV_64 __CONSTADDR_GV_66 __CONSTADDR_GV_67 __CONSTADDR_GV_68 __CONSTADDR_GV_69 __CONSTADDR_GV_70 __CONSTADDR_GV_71 __CONSTADDR_GV_72 __CONSTADDR_GV_73 __CONSTADDR_GV_74 __CONSTADDR_GV_75 __CONSTADDR_GV_77 __CONSTADDR_GV_78 __CONSTADDR_GV_79 __CONSTADDR_GV_80 __CONSTADDR_GV_81 __CONSTADDR_GV_82 __CONSTADDR_GV_83 __CONSTADDR_GV_84 __CONSTADDR_GV_85 __CONSTADDR_GV_86 __CONSTADDR_GV_88 __CONSTADDR_GV_89 __CONSTADDR_GV_90 __CONSTADDR_GV_91 __CONSTADDR_GV_92 __CONSTADDR_GV_93 __CONSTADDR_GV_94 __CONSTADDR_GV_95 __CONSTADDR_GV_96 __CONSTADDR_GV_98 __CONSTADDR_GV_99 __CONSTADDR_GV_100 __CONSTADDR_GV_101 __CONSTADDR_GV_102 __CONSTADDR_GV_103 __CONSTADDR_GV_104 __CONSTADDR_GV_105 __CONSTADDR_GV_106 __CONSTADDR_GV_107 __CONSTADDR_GV_109 __CONSTADDR_GV_110 __CONSTADDR_GV_111 __CONSTADDR_GV_112 __CONSTADDR_GV_113 __CONSTADDR_GV_114 __CONSTADDR_GV_115 __CONSTADDR_GV_116 __CONSTADDR_GV_117 __CONSTADDR_GV_118 __CONSTADDR_GV_119 __CONSTADDR_GV_120 __CONSTADDR_GV_122 __CONSTADDR_GV_123 __CONSTADDR_GV_124 __CONSTADDR_GV_125 __CONSTADDR_GV_126 __CONSTADDR_GV_127 __CONSTADDR_GV_128 __CONSTADDR_GV_129 __CONSTADDR_GV_130 __CONSTADDR_GV_131 __CONSTADDR_GV_132 __CONSTADDR_GV_133 __CONSTADDR_GV_135 __CONSTADDR_GV_136 __CONSTADDR_GV_137 __CONSTADDR_GV_138 __CONSTADDR_GV_139 __CONSTADDR_GV_140 __CONSTADDR_GV_141 __CONSTADDR_GV_142 __CONSTADDR_GV_143 __CONSTADDR_GV_144 __CONSTADDR_GV_145 __CONSTADDR_GV_146 __CONSTADDR_GV_147 __CONSTADDR_GV_148 __CONSTADDR_GV_149 __CONSTADDR_GV_150 __CONSTADDR_GV_151 __CONSTADDR_GV_152 __CONSTADDR_GV_153 __CONSTADDR_GV_154 __CONSTADDR_GV_155 __CONSTADDR_GV_156 __CONSTADDR_GV_157 __CONSTADDR_GV_158 __CONSTADDR_GV_159 __CONSTADDR_GV_160 __CONSTADDR_GV_161 __CONSTADDR_GV_162 __CONSTADDR_GV_163 __CONSTADDR_GV_164 __CONSTADDR_GV_169 __CONSTADDR_GV_174 __CONSTADDR_GV_179 __CONSTADDR_GV_184 __CONSTADDR_GV_189 __CONSTADDR_GV_194 __CONSTADDR_GV_199 __CONSTADDR_GV_204 __CONSTADDR_GV_165 __CONSTADDR_GV_166 __CONSTADDR_GV_167 __CONSTADDR_GV_168 __CONSTADDR_GV_170 __CONSTADDR_GV_171 __CONSTADDR_GV_172 __CONSTADDR_GV_173 __CONSTADDR_GV_175 __CONSTADDR_GV_176 __CONSTADDR_GV_177 __CONSTADDR_GV_178 __STATIC_LOCAL_50 __STATIC_LOCAL_51 __STATIC_LOCAL_52 __CONSTADDR_GV_180 __CONSTADDR_GV_181 __CONSTADDR_GV_182 __CONSTADDR_GV_183 __STATIC_LOCAL_70 __STATIC_LOCAL_69 __STATIC_LOCAL_72 __STATIC_LOCAL_71 __STATIC_LOCAL_74 __STATIC_LOCAL_73 __CONSTADDR_GV_205 __CONSTADDR_GV_206 __CONSTADDR_GV_207 __CONSTADDR_GV_208 __CONSTADDR_GV_209 __CONSTADDR_GV_210 __CONSTADDR_GV_211 __CONSTADDR_GV_212 __CONSTADDR_GV_213 __CONSTADDR_GV_214 __CONSTADDR_GV_215 __CONSTADDR_GV_216 __CONSTADDR_GV_217 __CONSTADDR_GV_218 __CONSTADDR_GV_219 __CONSTADDR_GV_220 __CONSTADDR_GV_221 __CONSTADDR_GV_222 __CONSTADDR_GV_223 __STATIC_LOCAL_75 __STATIC_LOCAL_76 __CONSTADDR_GV_224 __STATIC_LOCAL_78 __STATIC_LOCAL_79 __CONSTADDR_GV_225 __STATIC_LOCAL_81 __STATIC_LOCAL_82 __CONSTADDR_GV_226 __STATIC_LOCAL_84 __STATIC_LOCAL_85 __CONSTADDR_GV_227 __STATIC_LOCAL_87 __STATIC_LOCAL_88 __CONSTADDR_GV_228 __STATIC_LOCAL_90 __STATIC_LOCAL_91 __CONSTADDR_GV_229 __STATIC_LOCAL_93 __STATIC_LOCAL_94 __CONSTADDR_GV_230 __CONSTADDR_GV_231 __CONSTADDR_GV_232 __CONSTADDR_GV_233 __CONSTADDR_GV_234 __CONSTADDR_GV_235 __CONSTADDR_GV_236 __STATIC_LOCAL_97 __STATIC_LOCAL_96 __STATIC_LOCAL_100 __STATIC_LOCAL_99 __STATIC_LOCAL_103 __STATIC_LOCAL_102 __STATIC_LOCAL_106 __STATIC_LOCAL_105 __STATIC_LOCAL_53 __STATIC_LOCAL_54 __STATIC_LOCAL_55 __STATIC_LOCAL_56 __CONSTADDR_GV_185 __CONSTADDR_GV_186 __CONSTADDR_GV_187 __CONSTADDR_GV_188 __STATIC_LOCAL_57 __STATIC_LOCAL_58 __STATIC_LOCAL_59 __STATIC_LOCAL_60 __CONSTADDR_GV_190 __CONSTADDR_GV_191 __CONSTADDR_GV_192 __CONSTADDR_GV_193 __STATIC_LOCAL_61 __STATIC_LOCAL_62 __STATIC_LOCAL_63 __STATIC_LOCAL_64 __CONSTADDR_GV_195 __CONSTADDR_GV_196 __CONSTADDR_GV_197 __CONSTADDR_GV_198 __STATIC_LOCAL_65 __STATIC_LOCAL_66 __STATIC_LOCAL_67 __STATIC_LOCAL_68 __CONSTADDR_GV_200 __CONSTADDR_GV_201 __CONSTADDR_GV_202 __CONSTADDR_GV_203 __STATIC_LOCAL_108 __STATIC_LOCAL_109 __STATIC_LOCAL_110 __STATIC_LOCAL_111 __STATIC_LOCAL_112 __STATIC_LOCAL_113 __STATIC_LOCAL_114 __STATIC_LOCAL_117 __STATIC_LOCAL_116 __STATIC_LOCAL_121 __STATIC_LOCAL_120 __STATIC_LOCAL_115 __CONSTADDR_GV_237 __CONSTADDR_GV_238 __CONSTADDR_GV_239 __CONSTADDR_GV_240 __CONSTADDR_GV_241 __STATIC_LOCAL_125 __STATIC_LOCAL_124 __CONSTADDR_GV_242 __CONSTADDR_GV_243 __CONSTADDR_GV_244 __CONSTADDR_GV_245 __CONSTADDR_GV_246 __CONSTADDR_GV_247 __CONSTADDR_GV_248 __CONSTADDR_GV_249 __CONSTADDR_GV_250 __CONSTADDR_GV_251 __CONSTADDR_GV_252 __CONSTADDR_GV_253 __CONSTADDR_GV_254 __CONSTADDR_GV_255 __CONSTADDR_GV_256 __CONSTADDR_GV_257 __STATIC_LOCAL_126 __STATIC_LOCAL_127 __STATIC_LOCAL_128 __STATIC_LOCAL_129 __STATIC_LOCAL_130 __STATIC_LOCAL_132 __STATIC_LOCAL_131 __STATIC_LOCAL_77 __STATIC_LOCAL_80 __STATIC_LOCAL_83 __STATIC_LOCAL_86 __STATIC_LOCAL_89 __STATIC_LOCAL_92 __STATIC_LOCAL_95 __STATIC_LOCAL_98 __STATIC_LOCAL_101 __STATIC_LOCAL_104 __STATIC_LOCAL_107 __STATIC_LOCAL_118 __STATIC_LOCAL_119 __STATIC_LOCAL_122 __STATIC_LOCAL_123 __STATIC_LOCAL_133 __STATIC_LOCAL_134 __STATIC_LOCAL_136 __STATIC_LOCAL_135 __STATIC_LOCAL_138 __STATIC_LOCAL_137 __STATIC_LOCAL_140 __STATIC_LOCAL_139 __STATIC_LOCAL_142 __STATIC_LOCAL_141 __STATIC_LOCAL_143 __STATIC_LOCAL_144 __STATIC_LOCAL_146 __STATIC_LOCAL_147 __STATIC_LOCAL_150 __STATIC_LOCAL_149 __STATIC_LOCAL_154 __STATIC_LOCAL_153 __STATIC_LOCAL_158 __STATIC_LOCAL_157 __STATIC_LOCAL_162 __STATIC_LOCAL_161 __STATIC_LOCAL_165 __STATIC_LOCAL_166 __STATIC_LOCAL_168 __STATIC_LOCAL_169 __STATIC_LOCAL_171 __STATIC_LOCAL_172 __STATIC_LOCAL_173 __STATIC_LOCAL_175 __STATIC_LOCAL_176 __STATIC_LOCAL_178 __STATIC_LOCAL_179 __STATIC_LOCAL_182 __STATIC_LOCAL_181 __STATIC_LOCAL_185 __STATIC_LOCAL_184 __STATIC_LOCAL_188 __STATIC_LOCAL_187 __STATIC_LOCAL_190 __STATIC_LOCAL_191 __STATIC_LOCAL_145 __STATIC_LOCAL_148 __STATIC_LOCAL_151 __STATIC_LOCAL_152 __STATIC_LOCAL_155 __STATIC_LOCAL_156 __STATIC_LOCAL_159 __STATIC_LOCAL_160 __STATIC_LOCAL_163 __STATIC_LOCAL_164 __STATIC_LOCAL_167 __STATIC_LOCAL_170 __STATIC_LOCAL_174 __STATIC_LOCAL_177 __STATIC_LOCAL_180 __STATIC_LOCAL_183 __STATIC_LOCAL_186 __STATIC_LOCAL_189 __STATIC_LOCAL_192 atpSigHCommData.c elf-init.c __init_array_end __init_array_start _GLOBAL_OFFSET_TABLE_ _DYNAMIC __GNU_EH_FRAME_HDR nf90_put_var_2d_fourbytereal$netcdf_clone_116140_87_ nf90_open$netcdf_clone_115290_24_ to_out_2d_ main nf90_get_var_1d_eightbytereal$netcdf_clone_115290_68_ nf90_open$netcdf_clone_115290_56_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_100_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_14_ mergelayers_ __ALLOCATE __TMC_END__ _ITM_deregisterTMCloneTable _COS_V nf90_get_var_2d_eightbytereal$netcdf_clone_115290_35_ __atpHandlerInstall nf90_get_var_2d_eightbytereal$netcdf_clone_115290_15_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_105_ nf_put_vara_real_ __libc_csu_fini nf90_get_var_1d_eightbytereal$netcdf_clone_115290_67_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_26_ splitlayers_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_14_ $data_init_5$save_out_2ddetail_ nf90_open$netcdf_clone_115290_41_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_72_ strncmp@@GLIBC_2.2.5 nf90_put_var_1d_fourbytereal$netcdf_clone_116140_77_ __data_start nf_get_vara_double_ $data_init_1$write_initial_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_2_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_64_ __bss_start nf90_open$netcdf_clone_115290_20_ _DEALLOC interpol_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_88_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_33_ $data_init_4$save_out_2d_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_10_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_37_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_62_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_21_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_68_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_104_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_31_ nf90_enddef$netcdf_clone_116140_99_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_47_ $data_init_2$save_out_restart_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_37_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_27_ $data_init$getnetcdfaverage_ nf90_put_att_one_fourbytereal$netcdf_ __gmon_start__ nf90_strerror$netcdf_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_62_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_32_ $data_init_2$getvarnetcdf_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_102_ nf90_enddef$netcdf_clone_116140_28_ nf90_open$netcdf_clone_115290_50_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_65_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_25_ __cray2_EXP_V _ATP_VERSION_2.1.0.1 nf90_put_var_fourbyteint$netcdf_clone_116140_39_ __cray2_ALOG nf90_put_var_1d_fourbytereal$netcdf_clone_116140_67_ nf90_create$netcdf_clone_116140_92_ checkformelt_ nf90_enddef$netcdf_clone_116140_59_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_66_ _CLOSE __no_mmap_for_malloc nf90_put_var_1d_fourbytereal$netcdf_clone_116140_61_ $data_init$ini_fromfile_ nf90_create$netcdf_clone_116140_40_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_58_ mallopt nf90_get_var_1d_eightbytereal$netcdf_clone_115290_60_ $data_init$save_out_2ddetail_ nf90_open$netcdf_clone_115290_51_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_49_ nf90_def_dim$netcdf_clone_116140_16_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_35_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_17_ _edata lwrefreeze_ model_settings_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_101_ $data_init$save_out_2d_ vertgrid_ spinup_var_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_23_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_76_ nf90_close$netcdf_ _Jv_RegisterClasses nf90_put_var_1d_fourbytereal$netcdf_clone_116140_73_ nf90_open$netcdf_clone_115290_54_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_103_ _STOP2 __cray_dcopy_HSW nf_create_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_75_ speedcomp_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_11_ nf_enddef_ nf90_open$netcdf_clone_115290_42_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_46_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_29_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_63_ nf90_open$netcdf_clone_115290_16_ nf90_create$netcdf_clone_116140_1_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_13_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_86_ lwcontent_ memset nf90_get_var_2d_eightbytereal$netcdf_clone_115290_18_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_69_ __libc_csu_init nf90_open$netcdf_clone_115290_28_ nf90_inq_varid$netcdf_ nf90_open$netcdf_clone_115290_52_ getarg_ constants_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_45_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_31_ nf90_open$netcdf_clone_115290_64_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_13_ $data_init_4$handle_err_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_29_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_12_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_60_ nf_put_vara_double_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_11_ deletelayers_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_63_ _F90_LEN_TRIM_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_19_ __cray_sset_HSW nf90_def_dim$netcdf_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_89_ nf90_open$netcdf_clone_115290_5_ input_settings_ densific_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_66_ ini_temp_ nf90_enddef$netcdf_clone_116140_8_ nf90_open$netcdf_clone_115290_38_ _COS nf90_get_var_2d_eightbytereal$netcdf_clone_115290_48_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_30_ nf90_inquire_dimension$netcdf_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_61_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_32_ spinup_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_59_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_36_ nf90_open$netcdf_clone_115290_55_ nf90_open$netcdf_clone_115290_12_ nf90_open$netcdf_clone_115290_43_ nf_open_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_69_ __executable_start nf90_get_var_2d_eightbytereal$netcdf_clone_115290_6_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_44_ __libc_start_main@@GLIBC_2.2.5 nf90_get_var_1d_eightbytereal$netcdf_clone_115290_65_ __cray_dset_HSW nf90_def_var_manydims$netcdf_ _FRF nf90_get_var_2d_eightbytereal$netcdf_clone_115290_4_ __DTOR_END__ __pgas_register_dv findgrid_ subprogram_ _ITM_registerTMCloneTable nf90_put_var_1d_eightbytereal$netcdf_clone_116140_9_ __dso_handle nf_def_dim_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_74_ nf__enddef_ nf90_open$netcdf_clone_115290_53_ addlayers_ nf90_open$netcdf_clone_115290_33_ _OPEN $data_init_1$getnetcdfaverage_ __cray2_EXP $data_init$getvarnetcdf_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_3_ nf90_put_var_1d_fourbytereal$netcdf_clone_116140_71_ $data_init$read_averages_ sub_temp_exp_ $data_init$save_out_1d_ nf90_open$netcdf_clone_115290_1_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_8_ memmove _FWF $data_init$mainprogram_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_9_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_91_ ini_dens_ nf90_create$netcdf_clone_116140_78_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_34_ nf90_put_var_1d_eightbytereal$netcdf_clone_116140_34_ nf90_get_var_1d_eightbytereal$netcdf_clone_115290_57_ nf90_inq_dimid$netcdf_ nf90_create$netcdf_clone_116140_15_ nf_put_var1_int_ nf90_open$netcdf_clone_115290_7_ _COS_W nf90_put_var_fourbyteint$netcdf_clone_116140_38_ to_out_1d_ to_out_2ddetail_ _END _F90_STRING_COMPARE nf90_get_var_2d_eightbytereal$netcdf_clone_115290_36_ _ATP_Text_Globals nf90_get_var_2d_eightbytereal$netcdf_clone_115290_22_ timeloop_var_ nf90_def_var_onedim$netcdf_ nf90_enddef$netcdf_clone_116140_85_ clearall_ _IO_stdin_used nf90_get_var_2d_eightbytereal$netcdf_clone_115290_30_ nf90_open$netcdf_clone_115290_40_ nf90_put_var_2d_fourbytereal$netcdf_clone_116140_90_ $data_init_3$ini_fromfile_ nf90_get_var_2d_eightbytereal$netcdf_clone_115290_10_ _ATP_Data_Globals sub_temp_imp_ $data_init_3$save_out_1d_ nf90_open$netcdf_clone_115290_39_ tridiag_ __cray2_EXP_W nf90_put_var_1d_fourbytereal$netcdf_clone_116140_70_  .symtab .strtab .shstrtab .interp .note.ABI-tag .note.SuSE .hash .dynsym .dynstr .gnu.version .gnu.version_r .rela.dyn .rela.plt .init .text .fini .rodata .eh_frame_hdr .eh_frame .preinit_array .init_array .ctors .dtors .jcr .dynamic .got .got.plt .data .bss .comment .debug_aranges .debug_pubnames .debug_info .debug_abbrev .debug_line .debug_frame .debug_str .debug_loc .debug_ranges                                                                                 @                                          #             @                                          1             <@     <                                    <             X@     X      �                           B             �@     �      �                          J             �	@     �	      �                             R   ���o       *@     *      x                            _   ���o       �@     �                                   n             �@     �                                 x      B       �@     �      8                          �             @                                         }              @            �                            �              @            �K                             �             �dD     �d                                   �              eD      e     �$                              �             ��D     ��     t                             �             ��D     ��     �                             �             �d     �                                   �             �d     �                                   �              �d      �                                   �             �d     �                                   �              �d      �                                   �             (�d     (�                                �             8�d     8�                                  �             @�d     @�     �                            �             ��d     ��     h             @               �             @�e     (�     
              @                    0               (�     �                                                  г     �                                                  ��     _                              ,                     ߵ     �                             8                     ��     �                             F                     ��     �E                             R                     (     X                              _     0               �     V                            j                     �     =                             u                          p                                                    ��     �                                                   �     �{      (   (                	                      0�     [A                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             