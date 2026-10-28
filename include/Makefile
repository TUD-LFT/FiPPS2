# Makefile
# erzeugt die Library für alle Module

include ../makeconfig.mk

LIBRARY := libmodules.a

OBJ_PATH = ./obj
MOD_PATH = ./mod

OBJECTS := \
 globale_variablen.o \
 eigenschaften_typen.o \
 element_typen.o \
 failurecriteria_typen.o \
 knoten_typen.o \
 koordinatensystem_typen.o \
 laminat_typen.o \
 last_typen.o \
 material_typen.o \
 randbedingung_typen.o \
 residual_typen.o \
 internals_typen.o \
 ausgabe_typen.o \
 ergebnis_typen.o \
 fesimulation_typen.o \
 pre_assemble_types.o \
 vtk_variablen.o \
 vec_func.o \
 dof_types.o \
 konstanten.o \
 mat_func.o \
 integration_schemes.o \
 failure_criteria.o

OBJECTS_ABS = $(patsubst %,$(OBJ_PATH)/%,$(OBJECTS))

###############################################################################
## rules follow
###############################################################################

.PHONY: mkobj default debug

default : FCFLAGS = $(FC_DEF_FLAGS) -fpp
debug : FCFLAGS = $(FC_DBG_FLAGS) -fpp
default debug : mkobj all

mkobj:
	mkdir -p $(OBJ_PATH)
	mkdir -p $(MOD_PATH)

all: $(LIBRARY)

$(OBJ_PATH)/%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -module $(MOD_PATH) $(SLEPC_INCLUDE) $(PETSC_FC_INCLUDES) -o $@

$(LIBRARY): $(OBJECTS_ABS)
	$(AR) $(AR_FLAGS) $(LIBRARY) $(OBJECTS_ABS) 

clean2: clean
	$(RM) -rf $(LIBRARY) $(OBJ_PATH) $(MOD_PATH)
