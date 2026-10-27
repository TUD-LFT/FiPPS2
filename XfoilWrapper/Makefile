
include ../makeconfig.mk

PROG = XfoilWrapper

LIBRARY  = xfoilwrapper.a

SRC_PATH = ./src
OBJ_PATH = ./obj
INC_PATH = ./include
MOD_PATH = ./mod

#VPATH = $(SRC_PATH) $(OBJ_PATH) $(INC_PATH)

OBJECTS  = xfw_konst_var.o \
        xfw_netz_variablen.o \
        xfw_functions.o \
        xfw_read_input.o \
        xfw_write_xfoil_input.o \
        xfw_write_af_coords.o \
        xfw_start_xfoil.o \
        xfw_read_xfoil.o \
        xfw_ma_crit.o \
        xfoilwrapper_main.o

OBJECTS_ABS = $(patsubst %,$(OBJ_PATH)/%,$(OBJECTS))


.PHONY: mkobj default debug

default : FCFLAGS = $(FC_DEF_FLAGS)
debug : FCFLAGS = $(FC_DBG_FLAGS)
default debug : mkobj all

mkobj:
	mkdir -p $(OBJ_PATH)
	mkdir -p $(MOD_PATH)

all: $(LIBRARY)

$(OBJ_PATH)/%.o: $(SRC_PATH)/%.f90
	$(FC) $(FCFLAGS) -c $< -module $(MOD_PATH) -I$(INC_PATH) -o $@

$(LIBRARY): $(OBJECTS_ABS)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

clean::
	$(RM) -rfv $(LIBRARY) $(OBJ_PATH) $(MOD_PATH)
