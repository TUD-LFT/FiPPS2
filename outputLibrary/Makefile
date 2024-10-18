# Makefile
# erzeugt die dynamische Bibliothek der Benutzerdefinierten Ergebnisausgabe

include ../makeconfig.mk

INCL   = -I../include/

LIB = \
	-L../include/ -lmodule \
	-L../source/ -lfipps

LIBRARY := libuseroutput

OBJECTS := \
 sol1_output_vtk.o \
 sol1_output_user.o \
 sol2_output_vtk.o \
 sol2_output_user.o

###############################################################################
## rules follow
###############################################################################
.PHONY: all clean

all: $(LIBRARY) move

$(LIBRARY): $(OBJECTS)
	$(FC) $(FOPTS) -shared -o $(LIBRARY).so $(OBJECTS) $(LIB)

move:
	$(MV) $(LIBRARY).so ../
	
# rule for objects
%.o : %.f90
	$(FC) $(FOPTS) -fpp -fPIC $(INCL) -c $<

clean2: clean
	$(RM) ../$(LIBRARY).so
