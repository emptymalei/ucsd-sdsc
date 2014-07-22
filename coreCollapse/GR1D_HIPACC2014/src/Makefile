include ../make.inc

SOURCES=GR1D_module.F90 \
	timers.F90 \
	GR1D.F90 \
	atmosphere.F90 \
	eos.F90 \
	ppm.F90 \
	ye_of_rho.F90 \
	allocate_vars.F90 \
	initialize_vars.F90 \
	input_parser.F90 \
	problem.F90 \
	outinfo.F90 \
	output.F90 \
	output_control.F90 \
	grid.F90 \
	map_profile.F90 \
	mass_interior.F90 \
	shocktube.F90 \
	sedov.F90 \
	collapse.F90 \
	OSC.F90 \
	boundaries.F90 \
	driver.F90 \
	tools.F90 \
	Step.F90 \
	prim2con.F90 \
	reconstruct.F90 \
	flux_differences_HLLE.F90 \
	gravity.F90 \
	TVD_reconstruction.F90 \
	con2prim.F90 \
	press_sources.F90 \
	restart_H5.F90 \
	analysis.F90 \
	neutrino_pressure.F90 \
	start.F90 

MODULES=GR1D_module.mod
OBJECTS=$(SOURCES:.F90=.o )
	
EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS)

ifeq ($(HAVE_NUC_EOS),1)
	DEFS += -DHAVE_NUC_EOS
	EXTRADEPS += nuc_eos/nuc_eos.a
	EXTRAOBJECTS += nuc_eos/nuc_eos.a
endif
ifeq ($(HAVE_LEAK_ROS),1)
	DEFS += -DHAVE_LEAK_ROS
	EXTRAINCS += -Ileakage_rosswog
	EXTRAOBJECTS += leakage_rosswog/leakage_rosswog.a 
	EXTRADEPS += leakage_rosswog/leakage_rosswog.a 
endif	

GR1D:  GR1D_module.o $(EXTRADEPS) $(OBJECTS)  
	$(F90) $(LDFLAGS) -o ../GR1D $(OBJECTS) $(EXTRAOBJECTS)


$(OBJECTS): %.o: %.F90 
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) $(EXTRAINCS) -c $< -o $@

GR1D.o: GR1D_module.o

nuc_eos/nuc_eos.a: nuc_eos/*.F90 nuc_eos/*.f
	$(MAKE) -C nuc_eos

leakage_rosswog/leakage_rosswog.a: leakage_rosswog/leak_rosswog.F90
	$(MAKE) -C leakage_rosswog

clean:
	rm -rf ../GR1D
	rm -rf *.o 	
	rm -rf *.mod
	rm -rf nuc_eos/*.o
	rm -rf nuc_eos/*.mod
	rm -rf nuc_eos/*.a
	rm -rf leakage_rosswog/*.a
	rm -rf leakage_rosswog/*.o
	rm -rf leakage_rosswog/*.mod 

