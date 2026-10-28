
include ../makeconfig.mk

PROG = PANEL2D

LIBRARY  = panel2d.a

SRC_PATH = ./src
OBJ_PATH = ./obj
INC_PATH = ./include
MOD_PATH = ./mod

#VPATH = $(SRC_PATH) $(OBJ_PATH) $(INC_PATH)

MAIN = panel2d.o

OBJECTS  = konst_var.o \
	netz_variablen.o \
	functions.o \
	bl_functions.o \
	runge.o \
	ppval.o \
	spline.o \
	drag_squire_young.o \
	solve_bl_moran.o \
	stagnation_point.o \
	write_vtk.o \
	boundary_layer.o \
	pressure_calculation.o \
	calculate_ABmat.o \
	panel_geometry.o \
	calculate_geometry.o \
	sort_nodes.o \
	add_displacements.o \
	unify_airfoil.o \
	check_airfoil.o \
	read_input.o \
	splineInterpol.o \
	panel2d_main.o \
	xfoil_deps.o \
	refine_mesh.o

OBJECTS_ABS = $(patsubst %,$(OBJ_PATH)/%,$(OBJECTS))
MAIN_ABS = $(patsubst %,$(OBJ_PATH)/%,$(MAIN))

.PHONY: mkobj default debug

default : FCFLAGS = $(FC_DEF_FLAGS)
debug : FCFLAGS = $(FC_DBG_FLAGS)
default debug : mkobj all

mkobj:
	mkdir -p $(OBJ_PATH)
	mkdir -p $(MOD_PATH)

all: $(PROG)

$(OBJ_PATH)/xfoil_deps.o: $(SRC_PATH)/xfoil_deps.f
	$(FC) -c $(FCFLAGS) -r8 $(SRC_PATH)/xfoil_deps.f -module $(MOD_PATH) -o $@

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.f90
	$(FC) $(FCFLAGS) -c $< -module $(MOD_PATH) -I$(INC_PATH) -o $@

$(PROG): $(LIBRARY) $(MAIN_ABS)
	$(FC) -o $@ $(MAIN_ABS) $(LIBRARY) $(BLASLAPACK_LIB) $(LIBS)

$(LIBRARY): $(OBJECTS_ABS)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

clean::
	$(RM) -rfv $(PROG) $(LIBRARY) $(OBJ_PATH) $(MOD_PATH)
