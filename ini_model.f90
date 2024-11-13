!------------------------------------
!  SUBROUTINES TO INITIALIZE THE FIRN MODEL
!------------------------------------

!----------------------------
	subroutine model_settings(dtSnow,nyears,nyearsSU,dtmodelImp,dtmodelExp,ImpExp,dtobs, &
		kkinit,NoR,startasice,beginT,writeinprof,writeinspeed,writeindetail, &
		proflayers,detlayers,detthick,dzmax,initdepth,th,Grid,settingsfile, &
		username,domain)
	
	implicit none
	
	integer :: dtSnow,nyears,nyearsSU,dtmodelImp,dtmodelExp,ImpExp,dtobs,kkinit
	integer :: NoR,startasice,beginT,proflayers,detlayers
	integer :: writeinprof,writeinspeed,writeindetail
	double precision :: dzmax,initdepth,th,Grid(2),detthick
	character*255 :: settingsfile,username,domain,pad

	pad = "/scratch/ms/nl/"//trim(username)//"/data/ms_files/"
	open(unit=11,file=trim(pad)//"model_settings_"//trim(domain)//"_"//trim(settingsfile)//".txt")

	read (11,*)
	read (11,*)
        ! read model settings
    	read (11,*) nyears       ! simulation time [yr]
	read (11,*) nyearsSU     ! simulation time during the spin up period [yr]
	read (11,*) dtmodelExp   ! timestep in model with explicit T-scheme [s]
    	read (11,*) dtmodelImp	 ! timestep in model with implicit T-scheme [s]
	read (11,*) ImpExp	 ! Impicit or Explicit scheme (1=Implicit/fast, 2= Explicit/slow)
    	read (11,*) dtobs        ! timestep in input data [s]
        read (11,*) dtSnow       ! duration of running average of variables used in snow parameterisation [s]
	read (11,*)
    	read (11,*) dzmax        ! vertical model resolution [m]
    	read (11,*) initdepth    ! initial depth of firn profile [m]
    	read (11,*) th           ! theta (if theta=0.5 , it is a Crank Nicolson scheme) 
	read (11,*) startasice	 ! indicates the initial rho-profile (1=linear, 2=ice)
    	read (11,*) beginT	 ! indicates the inital T-profile (0=winter, 1=summer, 2=linear)
	read (11,*) NoR          ! number of times the data series is repeated for the initial rho profile is constructed
	read (11,*) 
	read (11,*) writeinspeed ! frequency of writing speed components to file
    	read (11,*) writeinprof  ! frequency of writing of firn profiles to file
	read (11,*) proflayers    ! number of output layers in the prof file
	read (11,*) writeindetail ! frequency of writing to detailed firn profiles to file
	read (11,*) detlayers     ! number of output layers in the detailed prof file
	read (11,*) detthick	  ! thickness of output layer in the detailed prof file
      	read (11,*)
        read (11,*) Grid(1) 	 ! Lon; indicates the longitude gridpoint
	read (11,*) Grid(2)	 ! Lat; indicates the latitude gridpoint

        kkinit = initdepth/dzmax  ! initial amount of vertical layers

	close(11)
	
	print *, "testest"
	
	end subroutine


!---------------------------
	subroutine input_settings(numLons,numLats,numTimes,domain,username)
	
	implicit none
	
	integer :: numLons,numLats,numTimes
	character*255 :: domain,username

	print *, "/perm/ms/nl/"//trim(username)//"/code/DATA/input_settings_"//trim(domain)//".txt"
	open(unit=12,file="/perm/ms/nl/"//trim(username)//"/code/DATA/input_settings_"//trim(domain)//".txt")
	
	read(12,*)
	read(12,*)
	read(12,*) numLons
	read(12,*) numLats
	read(12,*) numTimes

        print *, "Read input settings"
	
	close(12)
	
	end subroutine

	
!----------------------------
	subroutine constants(maxpore,rho0max,rhoi,R,pi)
	
	implicit none
	
	double precision :: maxpore,rho0max,rhoi,R,pi
		
	maxpore = 0.02	! the maximum irreducable water content [%]
    	rho0max = 470.	! maximum density of surface snow [kg m-3]
    	rhoi = 917.	! density of ice [kg m-3]
    	R = 8.3145      ! gas constant [J mole-1 K-1]
    	pi=asin(1.)*2	! pi = 3.1415
	
	end subroutine
		

!---------------------------
	subroutine FindGrid(Points,Grids,Latitude,Longitude,LSM,numLons,numLats)
	
	implicit none
	
	integer :: i,j,numLons,numLats
	integer, dimension(2) :: Points
	
	double precision :: dmax,dist,lon,lat
	double precision, dimension(2) :: Grids
	double precision, dimension(numLons,numLats) :: Latitude,Longitude,LSM
	
	print *, Grids(1),Grids(2)
	
	dmax=999.
	do i=1,numLons
	 do j=1,numLats
	  dist = sqrt(((Grids(1)-Longitude(i,j))*cos(((Grids(2)+Latitude(i,j))/2.)/360.* &
      		asin(1.)*4.))**2. + (Grids(2)-Latitude(i,j))**2.)
	  if (dist.lt.dmax.and.LSM(i,j).gt.0.5) then
	    dmax = dist
            Points(1) = i
            Points(2) = j
            lon = Longitude(i,j)
            lat = Latitude(i,j)    
          endif
         enddo
	enddo

	write (*,*) " "
	write (*,*) " Lon: ",Grids(1)," and Lat: ",Grids(2)
	write (*,*) "------------------------------------"
	write (*,*) " Closest gridpoint: ",Points(1),", ",Points(2)
	write (*,*) "          with Lon: ",lon," and Lat: ",lat
 	write (*,*) "-----------------------------------------------------------------"
	write (*,*) " "

	end subroutine


!---------------------------------
	subroutine interpol(TempSurf,PreSol,PreLiq,Sublim,SnowMelt,SnowDrif,FF10m, &
		TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM,Rho0FM, &
        numTimes,numSteps,numPoints,dtSnow,dtmodel,domain)
	
	integer :: step,a,b
	double precision :: part1,part2
	
	integer :: numTimes,numSteps,numPoints,dtmodel,dtSnow,numSnow
	double precision, dimension(numTimes) :: TempSurf,PreSol,PreLiq,Sublim
	double precision, dimension(numTimes) :: SnowMelt,SnowDrif,FF10m
	double precision, dimension(numPoints) :: TempFM,PsolFM,PliqFM,SublFM,MeltFM,DrifFM
	double precision, dimension(numPoints) :: ff10FM,Rho0FM
    double precision ::  TempSnow,ff10Snow
    character*255 :: domain

	do a = 1,numTimes-1
	  do  b = 1,numSteps
	    step = (a-1)*numSteps+b
	    part1 = (b-1.)/numSteps
	    part2 = 1. - part1
	    TempFM(step) = part2*TempSurf(a) + part1*TempSurf(a+1)
	    PSolFM(step) = (part2*PreSol(a) + part1*PreSol(a+1))/numSteps
	    PliqFM(step) = (part2*PreLiq(a) + part1*PreLiq(a+1))/numSteps
	    SublFM(step) = (part2*Sublim(a) + part1*Sublim(a+1))/numSteps
	    MeltFM(step) = (part2*SnowMelt(a) + part1*SnowMelt(a+1))/numSteps
	    DrifFM(step) = (part2*SnowDrif(a) + part1*SnowDrif(a+1))/numSteps
	    ff10FM(step) = part2*FF10m(a) + part1*FF10m(a+1)
	  end do
	end do

	do b = 1,numSteps
	    step = (numTimes-1)*numSteps+b
	    TempFM(step) = TempSurf(numTimes)
	    PSolFM(step) = PreSol(numTimes)/numSteps
	    PliqFM(step) = PreLiq(numTimes)/numSteps
	    SublFM(step) = Sublim(numTimes)/numSteps
	    MeltFM(step) = SnowMelt(numTimes)/numSteps
	    DrifFM(step) = SnowDrif(numTimes)/numSteps
	    ff10FM(step) = FF10m(numTimes)
	end do

    numSnow = max( int(dtSnow/dtmodel), 1)

    if (numSnow == 1) then
        ! Use current temperature and wind speed for snow parameterisations
        if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
            ! Greenland
            do step = 1,numPoints
                Rho0FM(step) = 362.1 + 2.78*(TempFM(step) - 273.15)     ! Fasuto et al. (2018)
            end do
        else
            ! Not Greenland
            do step = 1,numPoints
                Rho0FM(step) = 83. + 0.77*TempFM(step) + 11.67*ff10FM(step) ! Veldhuijsen (2022) - adjusted from Lenearts(2012)
                Rho0FM(step) = min(470., Rho0FM(step))
            end do
        end if
    else
        ! Use mean temperature and wind speed, averaged over numSnow time steps
        if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
            ! Greenland
            TempSnow = sum( TempFM(1:numSnow) )/numSnow
            do step = 1,numSnow
                Rho0FM(step) = 362.1 + 2.78*(TempSnow - 273.15)         ! Fausto et al. (2018)
            end do
            do step = numSnow+1,numPoints
                TempSnow = sum( TempFM(step-numSnow:step) )/numSnow
                Rho0FM(step) = 362.1 + 2.78*(TempSnow - 273.15)         ! Fausto et al. (2018)
            end do
        else
            ! Not Greenland
            TempSnow = sum( TempFM(1:numSnow) )/numSnow
            ff10Snow = sum( ff10FM(1:numSnow) )/numSnow
            do step = 1,numSnow
                Rho0FM(step) = 83. + 0.77*TempSnow + 11.67*ff10Snow     ! Veldhuijsen (2022) - adjusted from Lenearts(2012)
                Rho0FM(step) = min(470., Rho0FM(step))
            end do
            do step = numSnow+1,numPoints
                TempSnow = sum( TempFM(step-numSnow:step) )/numSnow
                ff10Snow = sum( ff10FM(step-numSnow:step) )/numSnow
                Rho0FM(step) = 83. + 0.77*TempSnow + 11.67*ff10Snow     ! Veldhuijsen (2022) - adjusted from Lenearts(2012)
                Rho0FM(step) = min(470., Rho0FM(step))
            end do
        end if
    end if
	
	print *, dtmodel

	print *, numSteps

	print *, FF10m(1:10)
	
	print *, ff10FM(1:10)
	
	print *, Rho0FM(1:10)
	

	end subroutine


!----------------------------
	subroutine read_averages(AveTsurf,AveAcc,AveWind,AveMelt,tsav,acav,ffav, &
		numLons,numLats,lon,lat,nyears)
	
	implicit none
	
	integer :: numLons,numLats,lon,lat,nyears
	double precision :: tsav,acav,ffav
	double precision, dimension(numLons,numLats) :: AveTsurf,AveAcc,AveWind,AveMelt

	print *, AveTsurf(lon,lat)
	print *, AveMelt(lon,lat)
	print *, AveWind(lon,lat)

	acav = AveAcc(lon,lat)
	ffav = AveWind(lon,lat)
	tsav = AveTsurf(lon,lat)
	if (tsav .gt. 273.15) tsav = 273.15    

	print *, "  "
	print *, "------------------------------------"	
	write(*,'(A19,1X,F8.3,1X,A1)') " Average Tsurf:     ", tsav, "K"
	write(*,'(A19,1X,F8.3,1X,A13)') " Average Acc:       ", acav, "mm w.e. yr-1"
	write(*,'(A19,1X,F8.3,1X,A5)') " Average Wind:      ", ffav, "m s-1"
	write(*,'(A19,1X,F8.3,1X,A13)') " Total Melt:        ", AveMelt(lon,lat), "mm w.e. yr-1"
	print *, "------------------------------------"
	print *, "  "

	end subroutine


!----------------------------- 	
    subroutine ini_dens(kk,kUL,initdepth,dzmax,rho0,rhoi,R,acav,tsav, &
        zs,DZ,Depth,Rho,DenRho,Mlwc,Refreeze,domain)
    
    implicit none
    
    ! declare arguments
    integer, intent(in) :: kk, kUL, initdepth    
    double precision, intent(in) :: dzmax, rho0, rhoi, R, acav, tsav
    double precision, intent(inout) :: zs
    double precision, dimension(kk), intent(inout) :: DZ, Depth, Rho, DenRho, Mlwc, Refreeze
    character*255, intent(in) :: domain

    ! declare local variables
    integer :: k
    double precision :: drho, Ec, Eg, g, part1, cons
    
    Depth(kUL) = 0.5*dzmax
    DZ(kUL) = dzmax
    Rho(kUL) = rho0
    Mlwc(kUL) = 0.
    DenRho(kUL) = 0.
    
    Ec = 60000.
    Eg = 42400.
    g = 9.81
    
    do k=kUL-1,1,-1
        part1 = (rhoi-Rho(k+1))*exp((-Ec/(R*tsav))+(Eg/(R*tsav)))
        if (Rho(k+1) .le. 550.) then
            if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
                cons = 0.6688 + 0.0048*log(acav)
            else
                cons = 1.288 - 0.117*log(acav)
            endif
            drho = 0.07*dzmax*Rho(k+1)*g*part1
        else
            if (trim(domain) .eq. "FGRN11" .or. trim(domain) .eq. "FGRN055") then
                cons = 1.7465 - 0.2045*log(acav)
            else
                cons = 6.387 *(acav*(-0.477))+0.195
            endif
            if (cons .lt. 0.25) cons = 0.25
            drho = 0.03*dzmax*Rho(k+1)*g*part1
        endif
        Rho(k) = drho + Rho(k+1)
        if (Rho(k) .gt. rhoi) Rho(k) = rhoi
    enddo
    
    do k = kUL-1,1,-1
        DZ(k) = dzmax
        Mlwc(k) = 0.
        DenRho(k) = 0.
        Refreeze(k) = 0.
        Depth(k) = Depth(k+1) + 0.5*DZ(k+1) + 0.5*DZ(k)
    enddo
    
    zs = 0.
               
    end subroutine ini_dens


!----------------------------- 
    	subroutine ini_temp(kk,kUL,beginT,tsav,pi,DZ,M,T,Rho,Depth,Year,rhoi)
        
    	implicit none
    
    	double precision :: ki,ci,period,om,diff
    
    	integer :: k,kk,kUL,beginT
    	double precision :: tsav,ampts,pi,kice,rhoi,kice_ref,kcal,kf,kair,kair_ref,theta
    	double precision,dimension(kk) :: DZ,M,T,Depth,Rho,Year
    
        ampts = 10.
    
    	ci = 152.5 + 7.122 * tsav                  ! heat capacity, Paterson (1994)
    	period = 365.*24.*3600.                    ! seconds in 1 yr          
    	om = 2.*pi/period                          ! rads per second

    	do k=kUL,1,-1    	!temperature - depth loop
    	  M(k) = Rho(k) * DZ(k)                      ! mass in layer k

          !ki = 0.021+2.5*(Rho(k)/1000.)**2           ! Anderson (1976)
          kice = 9.828 * exp(-0.0057*T(k))              ! Paterson et al., 1994
          kice_ref = 9.828 * exp(-0.0057*270.15)              ! Paterson et al., 1994
          kcal = 0.024 - 1.23E-4*Rho(k) + 2.5E-6*Rho(k)**2.
          kf = 2.107 + 0.003618*(Rho(k)-rhoi)           ! Calonne (2019)
          kair = (2.334E-3*T(k)**(3/2))/(164.54 + T(k)) ! Reid (1966)
          kair_ref = (2.334E-3*270.15**(3/2))/(164.54 + 270.15)
          theta = 1./(1.+exp(-0.04*(Rho(k)-450.)))
          ki = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf   ! Calonne (2019)

    	  Diff = ki/(Rho(k)*ci)                      ! thermal diffusivity
    	  if (k.eq.1) then
    	    T(k) = tsav
	  else
      	    if (beginT .eq. 1) then ! winter
              T(k) = tsav - ampts*exp(-(om/(2.*Diff))**0.5*Depth(k))* &
      	 	cos(-(om/(2.*Diff))**0.5*Depth(k))
      	    elseif (beginT .eq. 2) then !summer
              T(k) = tsav + ampts*exp(-(om/(2.*Diff))**0.5*Depth(k))* &
      	 	cos(-(om/(2.*Diff))**0.5*Depth(k))
      	    else 
      	      T(k) = tsav
      	    endif
      	    if (T(k) .gt. 272.15) T(k) = 272.15
          endif
	  
	  Year(k) = -999.
     	enddo
        
    	end subroutine
    

!---------------------------------
	subroutine timeloop_var(kk,dtmodel,nyears,DenRho,Refreeze,zs,Totvice, &
		Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy,Mrunoff, &
		TotRunoff,Mrefreeze,Totrefreeze,Mrain,TotRain,Msurfmelt, &
		TotSurfmelt,Msolin,TotSolIn,writeinprof,writeinspeed,writeindetail, &
		numOutputProf,numOutputSpeed,numOutputDetail,outputProf, &
		outputSpeed,outputDetail,Rho0out)
	
	implicit none
	
	integer :: kk,dtmodel,nyears
	integer :: writeinprof,writeinspeed,writeindetail
	integer :: numOutputProf,numOutputSpeed,numOutputDetail
	integer :: outputProf,outputSpeed,outputDetail

	double precision :: zs,Rho0out
	double precision :: Totvice,Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy
  	double precision :: Mrunoff,TotRunoff,Mrefreeze,Totrefreeze
  	double precision :: Mrain,TotRain,Msurfmelt,TotSurfmelt,Msolin,TotSolIn
	
	double precision,dimension(kk) :: DenRho,Refreeze
	
	
	DenRho(:) = 0.
	Refreeze(:) = 0.
	
	zs = 0.
	Totvice = 0.
	Totvacc = 0.
	Totvsub = 0.
	Totvsnd = 0.
	Totvmelt = 0.
	Totvfc = 0.
	Totvbouy = 0.
	
	Mrunoff = 0.
	TotRunoff = 0.
	Mrefreeze = 0.
	TotRefreeze = 0.
	Mrain = 0.
	TotRain = 0.
	Msurfmelt = 0.
	TotSurfmelt = 0.
	Msolin = 0.
	TotSolIn = 0.
	Rho0out = 0.	
	
	numOutputProf = writeinprof / dtmodel
	numOutputSpeed = writeinspeed / dtmodel
	numOutputDetail = writeindetail / dtmodel
	outputProf = nyears*3600*24*365 / writeinprof
	outputSpeed = nyears*3600*24*365 / writeinspeed
	outputDetail = nyears*3600*24*365 / writeindetail
	
	print *, " "
	print *, "Output variables"
	print *, "Prof, Speed, Detail"	
	print *, "writein...",writeinprof, writeinspeed, writeindetail
	print *, "numOutput...",numOutputProf, numOutputSpeed, numOutputDetail
	print *, "output...",outputProf,outputSpeed, outputDetail
	print *, " "

	
	end subroutine
	
	
!---------------------------------
	subroutine spinup_var(kUL,kkinit,zs,Mrunoff,Mrefreeze,Mrain,Msurfmelt,Msolin)
	
	implicit none
	
	integer :: kUL,kkinit
	double precision :: zs,Mrunoff,Mrefreeze,Mrain,Msurfmelt,Msolin

	Mrunoff = 0.
	Mrefreeze = 0.
	Mrain = 0.
	Msurfmelt = 0.
	Msolin = 0.
	
	zs = 0.
	kUL = kkinit

	
	end subroutine
	
