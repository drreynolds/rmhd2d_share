###################################################################
#  Makefile for unigrid resistive MHD code
#  Daniel R. Reynolds, UCSD Mathematics
###################################################################

# Automatic compilation of .o files from Fortran source
.SUFFIXES: .F .o

# Make.config is used for setting compile-time problem 
# configuration settings (boundary conditions, etc).
include Make.config

# Make.machine is used for setting machine-dependent compilation
# variables: e.g. compilers, libraries, F90-to-C conversion vars
include Make.machine

# append problem and interface preprocessor defs to compilation defs.
FFLAGS_ = $(FFLAGS) $(PROB_DEFS) $(SUND_DEFS)
CFLAGS_ = $(CFLAGS) $(PROB_DEFS) $(SUND_DEFS)



####### Executable source dependencies #######

# all executables depend on these files
ALLSRC = modules profiling main setupgrid setuplocaldomain   \
         setboundaryvalues constructlrstates output          \
         bdryexchange setroevariables setalphas              \
         seteigenvectors_prim seteigenvectors_cons           \
         seteigenvalues setdurl setsls inviscidfluxlf        \
         inviscidfluxroe inviscidfluxrp inviscidfluxcd2      \
         inviscidfluxcd4 inviscidfluxtcd inviscidfluxzip     \
         projection newdt xchange diags drivertools mhdtools

# executables including viscous fluxes depend on these
VISCSRC = viscousfluxcd2 viscousfluxcd4 viscousfluxtcd entsource 

# executables using explicit time-stepping depend on these
EXPSRC = mhdsolveRKTVD2 expRMHD_driver

# executables using KINSOL depend on these
KINSRC = impRMHD_interface impRMHD_driver nvector_mhd fnvector_mhd

# executables using CVODE depend on these
CVSRC = impCVODE_interface impCVODE_driver nvector_mhd fnvector_mhd

# header files (for dependency check when recompiling)
HEADS = nvector_mhd.h mesh.h mesh_common.h mesh_parms.h  \
        mesh_parms_static.h mesh_uparms.h properties.h   \
        mesh_uparms_static.h mgparams.h mpistuff.h boundaries.h
HEADERS = $(addprefix include/, $(HEADS))

# necessary lapack/blas files
SRC_LPK = dgemm dgemv dger dgetf2 dgetrf dgetrs dlaswp dscal \
          dswap dtrsm idamax ieeeck ilaenv xerbla lsame
OBJ_LPK = $(addprefix source/Lapack/, $(addsuffix .o, $(SRC_LPK)))

# preconditioning files
PRECSRC = ptri_parallel fastwave_prec viscous_prec-Du \
          viscous_prec-Db viscous_prec-De \
          vprec_solver vprec_mult combo_prec

# expRecon source files and corresp. object files
SRC1 = $(ALLSRC) $(VISCSRC) $(EXPSRC) init2
OBJ1 = $(addprefix source/, $(addsuffix .o, $(SRC1)))

# kinRecon source files and corresp object files
SRC2 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(KINSRC) init2
OBJ2 = $(addprefix source/, $(addsuffix .o, $(SRC2))) $(OBJ_LPK)

# cvRecon source files and corresp object files
SRC4 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(CVSRC) init2
OBJ4 = $(addprefix source/, $(addsuffix .o, $(SRC4))) $(OBJ_LPK)

# expLinWave source files and corresp. object files
SRC6 = $(ALLSRC) $(EXPSRC) initwave3
OBJ6 = $(addprefix source/, $(addsuffix .o, $(SRC6)))

# cvLinWave source files and corresp. object files
SRC7 = $(ALLSRC) $(PRECSRC) $(CVSRC) initwave3
OBJ7 = $(addprefix source/, $(addsuffix .o, $(SRC7))) $(OBJ_LPK)

# kinLinWave source files and corresp. object files
SRC8 = $(ALLSRC) $(PRECSRC) $(KINSRC) initwave3
OBJ8 = $(addprefix source/, $(addsuffix .o, $(SRC8))) $(OBJ_LPK)


# expKH source files and corresp. object files
SRC9 = $(ALLSRC) $(VISCSRC) $(EXPSRC) initKH
OBJ9 = $(addprefix source/, $(addsuffix .o, $(SRC9)))

# cvKH source files and corresp. object files
SRC10 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(CVSRC) initKH
OBJ10 = $(addprefix source/, $(addsuffix .o, $(SRC10))) $(OBJ_LPK)

# kinKH source files and corresp. object files
SRC11 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(KINSRC) initKH 
OBJ11 = $(addprefix source/, $(addsuffix .o, $(SRC11))) $(OBJ_LPK)


# kinTiltMode source files and corresp. object files
SRC12 = $(ALLSRC) $(PRECSRC) $(KINSRC) init_tiltmode bessel
OBJ12 = $(addprefix source/, $(addsuffix .o, $(SRC12))) $(OBJ_LPK)


# kinTearingMode source files and corresp. object files
SRC13 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(KINSRC) \
        init_tearingmode 
OBJ13 = $(addprefix source/, $(addsuffix .o, $(SRC13))) $(OBJ_LPK)


# expRT source files and corresp. object files
SRC14 = $(ALLSRC) $(VISCSRC) $(EXPSRC) initRT
OBJ14 = $(addprefix source/, $(addsuffix .o, $(SRC14)))

# kinRT source files and corresp. object files
SRC15 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(KINSRC) initRT
OBJ15 = $(addprefix source/, $(addsuffix .o, $(SRC15))) $(OBJ_LPK)

# cvRT source files and corresp. object files
SRC16 = $(ALLSRC) $(VISCSRC) $(PRECSRC) $(CVSRC) initRT
OBJ16 = $(addprefix source/, $(addsuffix .o, $(SRC16))) $(OBJ_LPK)


# test_psol_scaling source files and corresp. object files
PSOL_SCALING_TEST_SRC = $(addprefix source/, $(ALLSRC)) \
                        source/init2 source/ptri_parallel \
                        source/fastwave_prec tests/test_psol_scaling 
PSOL_SCALING_TEST_OBJ = $(addsuffix .o, $(PSOL_SCALING_TEST_SRC)) $(OBJ_LPK)


# test_psol_accuracy source files and corresp. object files
PSOL_ACCURACY_TEST_SRC = $(addprefix source/, $(ALLSRC)) \
                         source/init2 source/ptri_parallel \
                         source/ptri_product tests/test_psol_accuracy
PSOL_ACCURACY_TEST_OBJ = $(addsuffix .o, $(PSOL_ACCURACY_TEST_SRC)) $(OBJ_LPK)


# fnvec_mhd_test source files and corresp. object files
FNVEC_TEST_SRC = source/modules source/nvector_mhd source/fnvector_mhd \
                 tests/test_nvec_mhd tests/test_fnvec_mhd 
FNVEC_TEST_OBJ = $(addsuffix .o, $(FNVEC_TEST_SRC))


# test_driver source files and corresp. object files
TEST_DRIVER_SRC = $(ALLSRC) $(PRECSRC) init2 test_driver \
                  impCVODE_interface nvector_mhd fnvector_mhd
TEST_DRIVER_OBJ = $(addprefix source/, $(addsuffix .o, $(TEST_DRIVER_SRC))) $(OBJ_LPK)

# comparison_driver source files and corresp. object files
COMP_DRIVER_SRC = $(ALLSRC) $(PRECSRC) init2   \
                  mhdsolveRKTVD2 impCVODE_interface         \
                  comparison_driver nvector_mhd fnvector_mhd
COMP_DRIVER_OBJ = $(addprefix source/, $(addsuffix .o, $(COMP_DRIVER_SRC))) $(OBJ_LPK)


####### Remaining definitions (no adjustment necessary) #######

# SUNDIALS interface requirements
SUND_INCD   = -I${SUNDROOT}/include
SUND_LIBD   = -L${SUNDROOT}/lib
SUND_LIBS   = -lm
CVODE_LIBS  = -lsundials_fcvode -lsundials_cvode 
KINSOL_LIBS = -lsundials_fkinsol -lsundials_kinsol

# HYPRE interface requirements
HYPRE_INCD   = -I${HYPREROOT}/include
HYPRE_LIBD   = -L${HYPREROOT}/lib
HYPRE_LIBS   = -lHYPRE

# Final usable variables
INCS    = -Iinclude ${MPI_INCD} ${SUND_INCD} ${HYPRE_INCD} ${C_INCD}
LIBS    = ${MPI_LIBD} ${MPI_LIBS}
KLIBS   = ${SUND_LIBD} ${KINSOL_LIBS} ${SUND_LIBS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${LIBS}
CVLIBS  = ${SUND_LIBD} ${CVODE_LIBS} ${SUND_LIBS} ${HYPRE_LIBD} ${HYPRE_LIBS} ${LIBS}

# compilation cleanup files
TRASH = source/*.o source/Lapack/*.o tests/*.o *.mod source/*.il



####### Makefile targets #######

help :
	@echo " "
	@echo "usage: 'make <command>', where <command> is one of:"
	@echo "  expRecon             [explicit method, GEM reconnection]"
	@echo "  kinRecon             [KINSOL, GEM reconnection]"
	@echo "  cvRecon              [CVODE, GEM reconnection]"
	@echo "  kinTiltMode          [KINSOL, tearing-mode instability]"
	@echo "  kinTearingMode       [KINSOL, tearing-mode instability]"
	@echo "  expLinWave           [explicit, linear wave prop. test]"
	@echo "  cvLinWave            [CVODE linear wave prop. test]"
	@echo "  kinLinWave           [KINSOL linear wave prop. test]"
	@echo "  expKH                [explicit, Kelvin-Helmholtz test]"
	@echo "  cvKH                [CVODE Kelvin-Helmholtz test]"
	@echo "  kinKH                [KINSOL Kelvin-Helmholtz test]"
	@echo "  expRT                [explicit, Rayleigh-Taylor test]"
	@echo "  kinRT                [KINSOL Rayleigh-Taylor test]"
	@echo "  cvRT                 [CVODE Rayleigh-Taylor test]"
	@echo "  test_driver          [driver to run customized tests]"
	@echo "  comparison_driver    [driver to compare explicit vs implicit]"
	@echo "  fnvec_mhd_test       [checks F90/C vector kernel interface]"
	@echo "  clean                [cleans temporary build files]"
	@echo "  outclean             [cleans output from a run]"
	@echo "  tempclean            [cleans temporary build and ~ files]"
	@echo " "
	@echo "Note: controls over spatial dimension, parallel execution,"
	@echo "      boundary conditions, dynamic meshing, upwind vs CD"
	@echo "      inviscid fluxes, order of CD inviscid fluxes,"
	@echo "      primitive vs conservative upwind fluxing, use of"
	@echo "      Harten's entropy fix, and use of viscous fluxes"
	@echo "      are all decided through compiler definition flags"
	@echo "      in the file Make.config.  See that file for info."

expRecon : ${HEADERS} ${OBJ1}
	${F90} ${FFLAGS_} -o $@ ${OBJ1} ${INCS} ${LIBS}

kinRecon : ${HEADERS} ${OBJ2}
	${F90} ${FFLAGS_} -o $@ ${OBJ2} ${INCS} ${KLIBS} 

kinRecon2 : ${HEADERS} ${OBJ2N}
	${F90} ${FFLAGS_} -o $@ ${OBJ2N} ${INCS} ${KLIBS} 

cvRecon : ${HEADERS} ${OBJ4}
	${F90} ${FFLAGS_} -o $@ ${OBJ4} ${INCS} ${CVLIBS} 

cvRecon2 : ${HEADERS} ${OBJ4N}
	${F90} ${FFLAGS_} -o $@ ${OBJ4N} ${INCS} ${CVLIBS} 

kinTiltMode : ${HEADERS} ${OBJ12}
	${F90} ${FFLAGS_} -o $@ ${OBJ12} ${INCS} ${KLIBS} 

kinTiltMode2 : ${HEADERS} ${OBJ12N}
	${F90} ${FFLAGS_} -o $@ ${OBJ12N} ${INCS} ${KLIBS} 

kinTearingMode : ${HEADERS} ${OBJ13}
	${F90} ${FFLAGS_} -o $@ ${OBJ13} ${INCS} ${KLIBS} 

kinTearingMode2 : ${HEADERS} ${OBJ13N}
	${F90} ${FFLAGS_} -o $@ ${OBJ13N} ${INCS} ${KLIBS} 

expLinWave : ${HEADERS} ${OBJ6}
	${F90} ${FFLAGS_} -o $@ ${OBJ6} ${INCS} ${LIBS} 

cvLinWave : ${HEADERS} ${OBJ7}
	${F90} ${FFLAGS_} -o $@ ${OBJ7} ${INCS} ${CVLIBS} 

cvLinWave2 : ${HEADERS} ${OBJ7N}
	${F90} ${FFLAGS_} -o $@ ${OBJ7N} ${INCS} ${CVLIBS} 

kinLinWave : ${HEADERS} ${OBJ8}
	${F90} ${FFLAGS_} -o $@ ${OBJ8} ${INCS} ${KLIBS} 

kinLinWave2 : ${HEADERS} ${OBJ8N}
	${F90} ${FFLAGS_} -o $@ ${OBJ8N} ${INCS} ${KLIBS} 

expKH : ${HEADERS} ${OBJ9}
	${F90} ${FFLAGS_} -o $@ ${OBJ9} ${INCS} ${LIBS} 

cvKH : ${HEADERS} ${OBJ10}
	${F90} ${FFLAGS_} -o $@ ${OBJ10} ${INCS} ${CVLIBS} 

cvKH2 : ${HEADERS} ${OBJ10N}
	${F90} ${FFLAGS_} -o $@ ${OBJ10N} ${INCS} ${CVLIBS} 

kinKH : ${HEADERS} ${OBJ11}
	${F90} ${FFLAGS_} -o $@ ${OBJ11} ${INCS} ${KLIBS} 

kinKH2 : ${HEADERS} ${OBJ11N}
	${F90} ${FFLAGS_} -o $@ ${OBJ11N} ${INCS} ${KLIBS} 

expRT : ${HEADERS} ${OBJ14}
	${F90} ${FFLAGS_} -o $@ ${OBJ14} ${INCS} ${LIBS} 

kinRT : ${HEADERS} ${OBJ15}
	${F90} ${FFLAGS_} -o $@ ${OBJ15} ${INCS} ${KLIBS} 

kinRT2 : ${HEADERS} ${OBJ15N}
	${F90} ${FFLAGS_} -o $@ ${OBJ15N} ${INCS} ${KLIBS} 

cvRT : ${HEADERS} ${OBJ16}
	${F90} ${FFLAGS_} -o $@ ${OBJ16} ${INCS} ${CVLIBS} 

cvRT2 : ${HEADERS} ${OBJ16N}
	${F90} ${FFLAGS_} -o $@ ${OBJ16N} ${INCS} ${CVLIBS} 


test_psol_scaling : ${HEADERS} ${PSOL_SCALING_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${PSOL_SCALING_TEST_OBJ} ${INCS}

test_psol_accuracy : ${HEADERS} ${PSOL_ACCURACY_TEST_OBJ}
	${F90} ${FFLAGS_} -o $@ ${PSOL_ACCURACY_TEST_OBJ} ${INCS}


test_driver : ${HEADERS} ${TEST_DRIVER_OBJ}
	${F90} ${FFLAGS_} -o $@ ${TEST_DRIVER_OBJ} ${INCS} ${CVLIBS} 


comparison_driver : ${HEADERS} ${COMP_DRIVER_OBJ}
	${F90} ${FFLAGS_} -o $@ ${COMP_DRIVER_OBJ} ${INCS} ${CVLIBS} 


fnvec_mhd_test : ${HEADERS} ${FNVEC_TEST_OBJ}
	${F90} -o $@ ${FNVEC_TEST_OBJ} ${INCS} ${KLIBS}


clean:
	\rm -f ${TRASH}

outclean:
	\rm -f C_test* F_test* *.txt *.history output.0* outavs.00* dump.* gmon.out precmat.* rhsvec.*

tempclean: clean
	\rm -i *~ source/*~ include/*~



####### Makefile pattern rule redefinitions #######

.F.o :
	${F90} -c ${FFLAGS_} ${MPI_INCD} ${INCS} $< -o $@

.f.o :
	${F77} -c ${F77FLAGS} ${INCS} $< -o $@

.c.o :
	${CC} -c ${CFLAGS_} ${INCS} $< -o $@


####### End of Makefile #######
