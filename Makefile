OBJS =  functions_real.o functions_complex.o binary_decimal.o Base_to_Decimal.o reading_input.o Basis_KondoModel.o Basis_2_orb_Hubbard_chain.o Basis_2_orb_Hubbard_chain_KSector.o Basis_1_orb_Hubbard_2D_KSector.o Basis_3_orb_Hubbard_chain.o Basis_3_orb_Hubbard_chain_two_SzSectors.o Basis_3_orb_Hubbard_chain_GC.o Basis_multi_orb_Hubbard_chain_GC.o Basis_1_orb_Hubbard_GC.o  Basis_3_orb_Hubbard_chain_GC_restricted.o Basis_SpinlessFermionsFockSpace.o Model_SpinlessFermionsFockSpace.o Model_2_orb_Hubbard_chain.o Model_2_orb_Hubbard_chain_KSector.o Model_3_orb_Hubbard_chain.o Model_3_orb_Hubbard_chain_GC.o Model_multi_orb_Hubbard_chain_GC.o Model_1_orb_Hubbard_GC.o Model_3_orb_Hubbard_chain_two_SzSectors.o Model_3_orb_Hubbard_chain_two_SzSectors_complex.o Model_3_orb_Hubbard_chain_complex.o Basis_1_orb_Hubbard_chain.o Basis_1_orb_tJ.o Basis_Spins.o Model_Spins.o Basis_Spins_Target_Sz.o Model_Spins_Target_Sz.o  Model_1_orb_Hubbard_chain.o Model_1_orb_Hubbard_chain_complex.o Model_1_orb_Hubbard_2D_KSector.o Model_1_orb_tJ.o Model_KondoModel.o Model_1_orb_tJ_complex.o Lanczos_engine.o FTLM_Static.o LTLM_Static.o FTLM_Dynamics.o LTLM_Dynamics.o main.o 
DEBUG = -g3
#OPTFLAG = -O3
CC = g++ $(OPTFLAG) -std=c++11
CFLAGS = -c $(DEBUG) #-DUSE_COMPLEX
LFLAGS = $(DEBUG)
MKL_LIB = #/opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a /opt/intel/mkl/lib/libmkl_sequential.a
MKL_LIB += -llapack -lblas #-ldl -lpthread -lm
MKL_LIB2 = #-L/opt/intel/mkl/lib/ -lmkl_core -lmkl_intel_lp64 -lmkl_sequential /opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp65.a /opt/intel/mkl/lib/libmkl_sequential.a
MKL_LIB2 = #/opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a /opt/intel/mkl/lib/libmkl_sequential.a
MKL_include = #-I/opt/intel/mkl/include
OPENMP = #/opt/intel/compilers_and_libraries_2016.3.170/mac/compiler/lib/
LIBS_1 = -fopenmp #-L$(OPENMP)
LIBS_1 += #-liomp5 #-qopenmp
 

all : clean $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o lanczos $(MKL_include) $(MKL_LIB2) $(MKL_LIB) $(LIBS_1)

functions_real.o : functions_real.cpp
	$(CC) $(CFLAGS) functions_real.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

functions_complex.o : functions_complex.cpp
	$(CC) $(CFLAGS) functions_complex.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

reading_input.o : reading_input.cpp
	$(CC) $(CFLAGS) reading_input.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

binary_decimal.o : binary_decimal.cpp
	$(CC) $(CFLAGS) binary_decimal.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Base_to_Decimal.o : Base_to_Decimal.cpp
	$(CC) $(CFLAGS) Base_to_Decimal.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_SpinlessFermionsFockSpace.o : basis/Basis_SpinlessFermionsFockSpace.cpp
	$(CC) $(CFLAGS) basis/Basis_SpinlessFermionsFockSpace.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_2_orb_Hubbard_chain_KSector.o : basis/Basis_2_orb_Hubbard_chain_KSector.cpp
	$(CC) $(CFLAGS) basis/Basis_2_orb_Hubbard_chain_KSector.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_1_orb_Hubbard_2D_KSector.o : basis/Basis_1_orb_Hubbard_2D_KSector.cpp
	$(CC) $(CFLAGS) basis/Basis_1_orb_Hubbard_2D_KSector.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_2_orb_Hubbard_chain.o : basis/Basis_2_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) basis/Basis_2_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_3_orb_Hubbard_chain.o : basis/Basis_3_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) basis/Basis_3_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_3_orb_Hubbard_chain_two_SzSectors.o : basis/Basis_3_orb_Hubbard_chain_two_SzSectors.cpp
	$(CC) $(CFLAGS) basis/Basis_3_orb_Hubbard_chain_two_SzSectors.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_3_orb_Hubbard_chain_GC.o : basis/Basis_3_orb_Hubbard_chain_GC.cpp
	$(CC) $(CFLAGS) basis/Basis_3_orb_Hubbard_chain_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_multi_orb_Hubbard_chain_GC.o : basis/Basis_multi_orb_Hubbard_chain_GC.cpp
	$(CC) $(CFLAGS) basis/Basis_multi_orb_Hubbard_chain_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_1_orb_Hubbard_GC.o : basis/Basis_1_orb_Hubbard_GC.cpp
	$(CC) $(CFLAGS) basis/Basis_1_orb_Hubbard_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_3_orb_Hubbard_chain_GC_restricted.o : basis/Basis_3_orb_Hubbard_chain_GC_restricted.cpp
	$(CC) $(CFLAGS) basis/Basis_3_orb_Hubbard_chain_GC_restricted.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_KondoModel.o : basis/Basis_KondoModel.cpp
	$(CC) $(CFLAGS) basis/Basis_KondoModel.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_SpinlessFermionsFockSpace.o : models/Model_SpinlessFermionsFockSpace.cpp
	$(CC) $(CFLAGS) models/Model_SpinlessFermionsFockSpace.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_3_orb_Hubbard_chain_two_SzSectors.o : models/Model_3_orb_Hubbard_chain_two_SzSectors.cpp
	$(CC) $(CFLAGS) models/Model_3_orb_Hubbard_chain_two_SzSectors.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_3_orb_Hubbard_chain_two_SzSectors_complex.o : models/Model_3_orb_Hubbard_chain_two_SzSectors_complex.cpp
	$(CC) $(CFLAGS) models/Model_3_orb_Hubbard_chain_two_SzSectors_complex.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_2_orb_Hubbard_chain.o : models/Model_2_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) models/Model_2_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_2_orb_Hubbard_chain_KSector.o : models/Model_2_orb_Hubbard_chain_KSector.cpp
	$(CC) $(CFLAGS) models/Model_2_orb_Hubbard_chain_KSector.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_Hubbard_2D_KSector.o : models/Model_1_orb_Hubbard_2D_KSector.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_Hubbard_2D_KSector.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_3_orb_Hubbard_chain.o : models/Model_3_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) models/Model_3_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_3_orb_Hubbard_chain_complex.o : models/Model_3_orb_Hubbard_chain_complex.cpp
	$(CC) $(CFLAGS) models/Model_3_orb_Hubbard_chain_complex.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_3_orb_Hubbard_chain_GC.o : models/Model_3_orb_Hubbard_chain_GC.cpp
	$(CC) $(CFLAGS) models/Model_3_orb_Hubbard_chain_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_multi_orb_Hubbard_chain_GC.o : models/Model_multi_orb_Hubbard_chain_GC.cpp
	$(CC) $(CFLAGS) models/Model_multi_orb_Hubbard_chain_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_Hubbard_GC.o : models/Model_1_orb_Hubbard_GC.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_Hubbard_GC.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)


Model_KondoModel.o : models/Model_KondoModel.cpp
	$(CC) $(CFLAGS) models/Model_KondoModel.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)


Basis_1_orb_Hubbard_chain.o : basis/Basis_1_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) basis/Basis_1_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_Spins.o : basis/Basis_Spins.cpp
	$(CC) $(CFLAGS) basis/Basis_Spins.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_Spins_Target_Sz.o : basis/Basis_Spins_Target_Sz.cpp
	$(CC) $(CFLAGS) basis/Basis_Spins_Target_Sz.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Basis_1_orb_tJ.o : basis/Basis_1_orb_tJ.cpp
	$(CC) $(CFLAGS) basis/Basis_1_orb_tJ.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_Spins.o : models/Model_Spins.cpp
	$(CC) $(CFLAGS) models/Model_Spins.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_Spins_Target_Sz.o : models/Model_Spins_Target_Sz.cpp
	$(CC) $(CFLAGS) models/Model_Spins_Target_Sz.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_Hubbard_chain.o : models/Model_1_orb_Hubbard_chain.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_Hubbard_chain.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_tJ_complex.o : models/Model_1_orb_tJ_complex.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_tJ_complex.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_tJ.o : models/Model_1_orb_tJ.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_tJ.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Model_1_orb_Hubbard_chain_complex.o : models/Model_1_orb_Hubbard_chain_complex.cpp
	$(CC) $(CFLAGS) models/Model_1_orb_Hubbard_chain_complex.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

Lanczos_engine.o : Lanczos_engine.cpp
	$(CC) $(CFLAGS) Lanczos_engine.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

FTLM_Static.o : FTLM_Static.cpp
	$(CC) $(CFLAGS) FTLM_Static.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

LTLM_Static.o : LTLM_Static.cpp
	$(CC) $(CFLAGS) LTLM_Static.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

FTLM_Dynamics.o : FTLM_Dynamics.cpp
	$(CC) $(CFLAGS) FTLM_Dynamics.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

LTLM_Dynamics.o : LTLM_Dynamics.cpp
	$(CC) $(CFLAGS) LTLM_Dynamics.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp $(MKL_include) $(MKL_LIB) $(LIBS_1)

clean:
	rm -f *.o lanczos
