#
#
#FC = ftn
#NetCDF = /usr/local/apps/netcdf4/4.1.3_libmfix/LP64/
#INCLUDE = -I/$(NetCDF)/include 
#LIB = -L/usr/local/apps/netcdf4/4.1.3_libmfix/LP64/lib -lnetcdff -L/usr/local/apps/hdf5/1.8.8_libmfix/LP64/lib -L/usr/local/apps/szip/2.1_libmfix/LP64/lib -L/usr/local/apps/zlib/1.2.6/LP64/lib -L/usr/local/lib/libm.aix7_rounding_issue.IZ87564s01 -lnetcdf -lm -lhdf5_hl -lhdf5 -lsz -lz
#FFLAGS = -O3 -qstrict -qtune=auto -qcache=auto -qhot -q64 -qarch=auto -qessl
#OFLAGS = -lessl
##########################
#NetCDF = /usr/local/apps/netcdf4/4.4.1/CRAY/82/
#INCLUDE = -I/$(NetCDF)/include
#LIB2 = -L/$(NetCDF)/lib -lnetcdf
#LIB = $(LIB2) -lhdf5_hl -lm -lz   ##$(LIB2) -lhdf5_hl -lm -lsz -lz  #-lhdf5 -lhdf5_hl -lm -lsz -lz
#FFLAGS = -O 3
#EXTRA_0901 = -Wl,"--as-needed"
##########################
##########################
#FC = gfortran
#NetCDF = /usr/local/apps/netcdf4/4.4.1/CRAY/82/
#INCLUDE = -I/$(NetCDF)/include
#LIB2 = -L/$(NetCDF)/lib -lnetcdf
#LIB = -lhdf5_hl -lm -lz   ##$(LIB2) -lhdf5_hl -lm -lsz -lz  #-lhdf5 -lhdf5_hl -lm -lsz -lz
#FFLAGS = -O 3
#EXTRA_0901 = -Wl,"--as-needed"
#FFLAGS = -O3 -qstrict -qtune=auto -qcache=auto -qhot -q64 -qarch=auto -qessl
#OFLAGS = -lessl
##########################
##########################
#FC = ifort
#NetCDF = /usr/local/apps/netcdf4/4.3.2/INTEL/150/
#INCLUDE = -I/$(NetCDF)/include
#LIB2 = -L/$(NetCDF)/lib -lnetcdf
#LIB = $(LIB2) -lhdf5 -lm -lz   ##$(LIB2) -lhdf5_hl -lm -lsz -lz  #-lhdf5 -lhdf5_hl -lm -lsz -lz
##FFLAGS = -O 3
##EXTRA_0901 = -Wl,"--as-needed"
##FFLAGS = -O3 -qstrict -qtune=auto -qcache=auto -qhot -q64 -qarch=auto -qessl
##OFLAGS = -lessl
##########################
FC = ftn
#NetCDF = /usr/local/apps/netcdf4/4.4.1/CRAY/82/
#INCLUDE = -I/$(NetCDF)/include
#LIB2 = -L/$(NetCDF)/lib -lnetcdf
#LIB = -lhdf5_hl -lm -lz  
FFLAGS = -O3
EXTRA_0901 = -Wl,"--as-needed"
############################

OBJECTS = \
  mainprogram_firnmodel.f90 \
  openNetCDF.f90 \
  ini_model.f90 \
  time_loop.f90 \
  output.f90 \
  subprogram.f90 

#IMAU-FDM_v0901.x: $(OBJECTS)
#	$(FC) $(INCLUDE) -o $@ $(OBJECTS) ${OFLAGS} $(FFLAGS) $(EXTRA_0901) $(LIB)
#	cp IMAU-FDM_v0901.x ../NP_Run/

#IMAU-FDM_v0901.x: $(OBJECTS)
#	$(FC) $(INCLUDE) -o $@ $(OBJECTS) $(LIB)
#	cp IMAU-FDM_v0901.x ../NP_Run/

IMAU-FDM_np.x: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(FFLAGS) $(EXTRA_0901)
	#cp IMAU-FDM_np.x ../NS_Run/


# .SUFFIXES:       .f90

.f90.o:
#	$(FC)  $*.f90 -o $@
	$(FC) -c $*.f90

clean:
	rm IMAU-FDM_v0901.x
	
  


