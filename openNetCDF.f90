!------------------------------------
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

