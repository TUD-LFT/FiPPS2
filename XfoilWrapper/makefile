
include ../makeconfig.mk

PROG = XfoilWrapper

LIBRARY  = xfoilwrapper.a

SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include
MOD_PATH = ./mod/

VPATH = $(SRC_PATH) $(OBJ_PATH) $(INC_PATH)

OBJS  = xfw_konst_var.o \
        xfw_netz_variablen.o \
        xfw_functions.o \
        xfw_read_input.o \
        xfw_write_xfoil_input.o \
        xfw_write_af_coords.o \
        xfw_start_xfoil.o \
        xfw_read_xfoil.o \
        xfw_ma_crit.o \
        xfoilwrapper_main.o

OBJ = $(patsubst %,$(OBJ_PATH)%,$(OBJS))

MAIN = $(patsubst %,$(OBJ_PATH)%,$(MAINS))

all: $(LIBRARY)

$(OBJ_PATH)%.o: $(SRC_PATH)%.f90
	$(FC) $(FC_FLAGS) -c $< $(INC_PATH) -o $@

$(LIBRARY): $(OBJ)
	$(AR) $(AR_FLAGS) $(LIBRARY) $^

clean::
	$(RM) -rfv $(LIBRARY) $(OBJ) $(MAIN) $(MOD_PATH)*.mod
