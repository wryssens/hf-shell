OBJDIR :=   obj
SRCDIR :=   src
LIBDIR :=   lib
MODDIR :=   mod

#-------------------------------------------------------------------------------
# Compilers and some recommended options
#Default compiler is gfortran
CXX :=   gfortran

ifeq ($(CXX),gfortran)
	CXXFLAGS := -O3 -J$(MODDIR) -Wall 
else ifeq ($(CXX),ifort)
	CXXFLAGS := -O3 -no-wrap-margin -module $(MODDIR)
endif
#-------------------------------------------------------------------------------
# Source code specification
TARGET  :=  hf_shell.exe
SRC     :=  clebsches.f90 generic.f90   smbasis.f90     interaction.f90        
SRC     +=  HF.f90        BCS.f90         HFB.f90
SRC     +=  evolution.f90 observables.f90 thermal_projection.f90 
SRC     +=  printing.f90  IO.f90          hf_shell.f90 
#Make the lists of objects
OBJ     :=  $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))
LIBS   := -llapack -lblas

#-------------------------------------------------------------------------------
# Recipes
.PHONY: all clean

$(TARGET):  $(OBJ_LIB) $(OBJ) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_LIB) $(OBJ) $(LIBS)    

clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(LIBDIR)/*.o
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90 | $(OBJDIR)/ $(MODDIR)/  
	$(CXX) $(CXXFLAGS) -c  $< -o $@

$(OBJDIR)/:
	mkdir -p  obj/

$(MODDIR)/:
	mkdir -p  mod/

#-------------------------------------------------------------------------------

