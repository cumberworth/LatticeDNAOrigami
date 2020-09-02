OPTLEVEL = -O3
BUILDDIR = ../../build/
TARGET = calc_num_walks
TARGETDIR = ./
SRCDIR = ../../src
INCLUDEDIR = ../../include

vpath %.h $(INCLUDEDIR)/
vpath %.cpp $(SRCDIR)/

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(subst .cpp,.o,$(SOURCES))
OBJECTS := $(subst $(SRCDIR),$(BUILDDIR),$(OBJECTS))

CPP = mpicxx
CPPFLAGS = -std=c++14 -I ../../include $(OPTLEVEL)
LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -lboost_system -lboost_filesystem $(OPTLEVEL)

all: $(TARGET)

$(TARGET): calc_num_walks.o $(BUILDDIR)origami_system.o $(BUILDDIR)hash.o $(BUILDDIR)order_params.o $(BUILDDIR)bias_functions.o $(BUILDDIR)origami_potential.o $(BUILDDIR)utility.o $(BUILDDIR)random_gens.o $(BUILDDIR)nearest_neighbour.o $(BUILDDIR)jsoncpp.o $(BUILDDIR)files.o $(BUILDDIR)domain.o $(BUILDDIR)simulation.o $(BUILDDIR)top_constraint_points.o $(BUILDDIR)parser.o $(BUILDDIR)random_gens.o $(BUILDDIR)movetypes.o $(BUILDDIR)ideal_random_walk.o $(BUILDDIR)met_movetypes.o $(BUILDDIR)orientation_movetype.o $(BUILDDIR)rg_movetypes.o $(BUILDDIR)transform_movetypes.o $(BUILDDIR)cb_movetypes.o $(BUILDDIR)us_simulation.o $(BUILDDIR)ptmc_simulation.o
	$(CPP) $(LDFLAGS) -o $@ $^

$calc_num_walks.cpp: calc_num_walks.cpp
	$(CCP) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean
clean:
	rm *.o
