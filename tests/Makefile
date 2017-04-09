CC = mpicxx

# Without ceres
CPPFLAGS = -I../include -g
LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -g
#CPPFLAGS = -I ../include -O3
#LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -O3

# With ceres
#CPPFLAGS = -I ../include -I/usr/include/eigen3 -g
#LDFLAGS = -lblas -llapack -lcholmod -lcxsparse -lgomp -lceres -lglog -lboost_program_options -lboost_mpi -lboost_serialization -g
#CPPFLAGS = -I ../include -I/usr/include/eigen3 -O3
#LDFLAGS = -lblas -llapack -lcholmod -lcxsparse -lgomp -lceres -lglog -lboost_program_options -lboost_mpi -lboost_serialization -O3

BUILD = ../build/
BUILDTEST = build/
vpath %.h ../include
vpath %.cpp src

TARGET = bin/testLatticeDNAOrigami

all: $(TARGET)

$(TARGET): $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)movetypes.o $(BUILD)random_gens.o $(BUILD)ideal_random_walk.o $(BUILD)parser.o $(BUILD)order_params.o $(BUILDTEST)test_main.o $(BUILDTEST)test_ideal_random_walk.o $(BUILDTEST)test_movetypes.o $(BUILDTEST)test_nearest_neighbour.o $(BUILDTEST)test_order_params.o $(BUILDTEST)test_origami_system.o $(BUILDTEST)test_random_gens.o $(BUILDTEST)test_simulation.o $(BUILD)constant_temp_simulation.o $(BUILD)annealing_simulation.o $(BUILD)ptmc_simulation.o $(BUILD)us_simulation.o #$(BUILD)mezei_us_simulation.o
	$(CC) -o $@ $^ $(LDFLAGS)
	$(TESTTARGET)

$(BUILDTEST)test_main.o: test_main.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_ideal_random_walk.o: test_ideal_random_walk.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_movetypes.o: test_movetypes.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_nearest_neighbour.o: test_nearest_neighbour.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_order_params.o: test_order_params.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_origami_system.o: test_origami_system.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_random_gens.o: test_random_gens.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILDTEST)test_simulation.o: test_simulation.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean install
clean:
	rm $(BUILDTEST)*.o
