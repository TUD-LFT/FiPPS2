# Makefile
# erzeugt das Programm apame

include ../makeconfig.mk

SRC_PATH = ./src
OBJ_PATH = ./obj
INC_PATH = ./include
MOD_PATH = ./mod
BIN_PATH = ./bin

PROGRAM = apame

LIBRARY = apame.a

MAIN = \
	apame.o

OBJECTS =  \
	module_kind_and_konst.o \
	module_grid.o \
	module_induced_drag.o \
	module_influence.o \
	module_input.o \
	module_modify_wake.o \
	module_pressure.o \
	module_rhs.o \
	module_velocity.o \
	func_IDW.o \
	func_message.o \
	func_new_line.o \
	func_tic.o \
	func_title.o \
	func_toc.o \
	func_vect_ang.o \
	func_vect_inprod.o \
	func_vect_mod.o \
	func_vect_prod.o \
	subr_apame_main.o \
	subr_vtk_results.o \
	subr_results.o \
	subr_speeds.o \
	subr_write_case.o \
	subr_write_coef.o \
	subr_write_key.o \
	subr_write_ptot.o \
	subr_write_time.o

OBJECTS_ABS = $(patsubst %,$(OBJ_PATH)/%,$(OBJECTS))
MAIN_ABS = $(patsubst %,$(OBJ_PATH)/%,$(MAIN))

###############################################################################
## rules follow
###############################################################################

.PHONY: mkobj default debug

default : FCFLAGS = $(FC_DEF_FLAGS)
debug : FCFLAGS = $(FC_DBG_FLAGS)
default debug : mkobj all

mkobj:
	mkdir -p $(OBJ_PATH)
	mkdir -p $(MOD_PATH)

all: $(PROGRAM) move

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.f90
	$(FC) $(FCFLAGS) -c $< -module $(MOD_PATH) -I$(INC_PATH) -o $@
	
$(PROGRAM): $(LIBRARY) $(MAIN_ABS)
	$(FC) -o $(PROGRAM) $(MAIN_ABS) $(BLASLAPACK_LIB) $(LIBRARY) $(LIB)

$(LIBRARY): $(OBJECTS_ABS)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

move:
	mkdir -p $(BIN_PATH)
	$(MV) $(PROGRAM) $(BIN_PATH)

clean::
	$(RM) -rf $(LIBRARY) $(BIN_PATH) $(OBJ_PATH) $(MOD_PATH)
	$(RM) -f $(SRC_PATH)/*__genmod.f90
