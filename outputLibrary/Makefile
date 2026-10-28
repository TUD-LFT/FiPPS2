# Makefile
# erzeugt die dynamische Bibliothek der Benutzerdefinierten Ergebnisausgabe

include ../makeconfig.mk

LIB = \
	-L../include/ -lmodules \
	-L../source/ -lfipps

LIBRARY := libuseroutput.so

OBJ_PATH = ./obj
MOD_PATH = ./mod
INCL_PATH = ../include/mod

OBJECTS := \
 sol1_output_vtk.o \
 sol1_output_user.o \
 sol2_output_vtk.o \
 sol2_output_user.o

OBJECTS_ABS = $(patsubst %,$(OBJ_PATH)/%,$(OBJECTS))

###############################################################################
## rules follow
###############################################################################

.PHONY: default debug

default : FCFLAGS = $(FC_DEF_FLAGS) -fpp
debug : FCFLAGS = $(FC_DBG_FLAGS) -fpp
default debug : mkobj all

mkobj:
	mkdir -p $(OBJ_PATH)
	mkdir -p $(MOD_PATH)
	
all: $(LIBRARY) move

$(LIBRARY): $(OBJECTS_ABS)
	$(FC) $(FCFLAGS) -shared -o $(LIBRARY) $(OBJECTS_ABS) $(LIB)

move:
	$(MV) $(LIBRARY) ../
	
# rule for objects
$(OBJ_PATH)/%.o : %.f90
	$(FC) $(FCFLAGS) -c $< -module $(MOD_PATH) -I$(INCL_PATH) -o $@

clean2: clean
	$(RM) -rf ../$(LIBRARY) $(OBJ_PATH) $(MOD_PATH)
