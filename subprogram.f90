!------------------------------------
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
	
