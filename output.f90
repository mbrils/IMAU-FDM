!------------------------------------
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
