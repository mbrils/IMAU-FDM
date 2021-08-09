	
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

