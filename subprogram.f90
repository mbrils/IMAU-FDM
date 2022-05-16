!------------------------------------
!  SUBROUTINE WITH TIME LOOP OF FIRN MODEL
!------------------------------------

!----------------------------
    subroutine subprogram(beginT,dtmodel,dtobs,ImpExp,numTimes,numSteps, &
        numPoints,numPointsSU,dtSnow,nyears,nyearsSU,kk,kkinit,initdepth,numberofrepeat,th,R, &
        pi,dzmax,rhoi,maxpore,writeinprof,writeinspeed,writeindetail, &
        proflayers,detlayers,detthick,acav,tsav,TempSurf,PreSol,PreLiq, &
        Sublim,SnowMelt,SnowDrif,FF10m,IceShelf,settingsfile,fname_p1,username, &
        domain,ini_fname)
        
    IMPLICIT NONE
    
    integer :: dtmodel,dtobs,ImpExp,IceShelf,kk,kkinit,kUL,initdepth
    integer :: numTimes,numPoints,numPointsSU,numSteps,dtSnow,nyears,nyearsSU,beginT,i,time
    integer :: numberofrepeat,writeinprof,writeinspeed,writeindetail
    integer :: proflayers,detlayers
    double precision :: th,R,pi,dzmax,rho0,rhoi,maxpore,acav,tsav,detthick
    
    integer :: numOutputSpeed,numOutputProf,numOutputDetail
    integer :: outputSpeed,outputProf,outputDetail
    double precision :: zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy
    double precision :: Totvice,Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy
    double precision :: Mrunoff,TotRunoff,Mrefreeze,Totrefreeze,FirnAir,TotLwc,IceMass
    double precision :: Mrain,TotRain,Msurfmelt,TotSurfmelt,Msolin,TotSolIn,Rho0out
    double precision :: Ts,Psol,Pliq,Su,Me,Sd
    
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
    
    ! Spin up the model to a 'steady state'
    call spinup(numTimes,numPoints,numPointsSU,kk,kUL,dtmodel, &
        R,rhoi,acav,tsav,th,dzmax,rho0,maxpore,zs,Msurfmelt,Mrain,Msolin, &
        Mrunoff,Mrefreeze,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year,TempFM, &
        PSolFM,PLiqFM,SublFM,MeltFM,DrifFM,Rho0FM,IceShelf,ImpExp,nyears,nyearsSU,domain)
        
    ! Get variables for the time integration and open outputfile
    call timeloop_var(kk,dtmodel,nyears,DenRho,Refreeze,zs,Totvice, &
        Totvfc,Totvacc,Totvsub,Totvsnd,Totvmelt,Totvbouy,Mrunoff, &
        TotRunoff,Mrefreeze,Totrefreeze,Mrain,TotRain,Msurfmelt, &
        TotSurfmelt,Msolin,TotSolIn,writeinprof,writeinspeed,writeindetail, &
        numOutputProf,numOutputSpeed,numOutputDetail,outputProf, &
        outputSpeed,outputDetail,Rho0out)
    
    ! Write intitial profile to NetCDF-file and prepare output arrays
    call write_initial(kk,kUL,Rho,M,T,Depth,Mlwc,Year,settingsfile,fname_p1, &
        username)
    
    allocate(out_1D((outputSpeed+50),18))
    allocate(out_2D_dens((outputProf+50),proflayers))
    allocate(out_2D_temp((outputProf+50),proflayers))
    allocate(out_2D_lwc((outputProf+50),proflayers))
    allocate(out_2D_depth((outputProf+50),proflayers))
    allocate(out_2D_dRho((outputProf+50),proflayers))
    allocate(out_2D_year((outputProf+50),proflayers))
    allocate(out_2D_det_dens((outputDetail+50),detlayers))
    allocate(out_2D_det_temp((outputDetail+50),detlayers))
    allocate(out_2D_det_lwc((outputDetail+50),detlayers))
    allocate(out_2D_det_refreeze((outputDetail+50),detlayers))
    
    !set missing value
    out_1D(:,:) = 9.96921e+36
    out_2D_dens(:,:) = 9.96921e+36 
    out_2D_temp(:,:) = 9.96921e+36 
    out_2D_lwc(:,:) = 9.96921e+36 
    out_2D_depth(:,:) = 9.96921e+36 
    out_2D_dRho(:,:) = 9.96921e+36 
    out_2D_year(:,:) = 9.96921e+36 
    out_2D_det_dens(:,:) = 9.96921e+36
    out_2D_det_temp(:,:) = 9.96921e+36
    out_2D_det_lwc(:,:) = 9.96921e+36
    out_2D_det_refreeze(:,:) = 9.96921e+36
    
    ! Time integration
    print *, "Start of time loop"
    do time = 1, numPoints
        rho0 = Rho0FM(time)
        Ts = TempFM(time)
        Psol = PSolFM(time)
        Pliq = PliqFM(time)
        Su = SublFM(time)
        Me = MeltFM(time)
        Sd = DrifFM(time)
        
        ! Calculate the densification	  
        call densific(kk,kUL,dtmodel,R,rhoi,acav,tsav,Rho,T,domain)
        
        ! Re-calculate the Temp-profile (explicit or implicit)		  
        if (ImpExp .eq. 1) call sub_temp_imp(kk,kUL,dtmodel,th,Ts,T,Rho,DZ,rhoi)
        if (ImpExp .eq. 2) call sub_temp_exp(kk,kUL,dtmodel,Ts,T,Rho,DZ,rhoi)
        
        ! Re-caluclate DZ/M-values according to new Rho-/T-values
        ! And calculate all liquid water processes
        
        if (mod(time,200000) .eq. 0) print *, time, zs
        
        call vertgrid(time,kk,kUL,dtmodel,numPoints,nyears,dzmax,rho0,rhoi,acav, &
            maxpore,zs,dzd,vmelt,vacc,vsub,vsnd,vfc,vbouy,Ts,PSol,PLiq, &
            Su,Me,Sd,M,T,DZ,Rho,DenRho,Depth,Mlwc,Refreeze,Year, &
            ImpExp,IceShelf,Msurfmelt,Mrain,Msolin,Mrunoff,Mrefreeze)
        
        ! Calculate the firn air content and total liquid water content of the firn column
        FirnAir = 0.
        TotLwc = 0.
        IceMass = 0.
        do i = 1,kUL
            if (Rho(i).le.910.) FirnAir = FirnAir + DZ(i)*(rhoi-Rho(i))/(rhoi)
            TotLwc = TotLwc + Mlwc(i)
            IceMass = IceMass + M(i)
        enddo
        
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
    enddo
    
    Year(:) = Year(:) - nyears
    
    ! Finished time loop
    print *, "End of Time Loop"
    print *, " "
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

103 continue

    end subroutine

