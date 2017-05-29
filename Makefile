#CC = g++
CC = mpicxx

# Without ceres
CPPFLAGS = -Iinclude -g
LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -g
#CPPFLAGS = -Iinclude -O3
#LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -O3

# With ceres
#CPPFLAGS = -Iinclude -I/usr/include/eigen3 -g
#LDFLAGS = -lblas -llapack -lcholmod -lcxsparse -lgomp -lceres -lglog -lboost_program_options -lboost_mpi -lboost_serialization -g
#CPPFLAGS = -Iinclude -I/usr/include/eigen3 -O3
#LDFLAGS = -lblas -llapack -lcholmod -lcxsparse -lgomp -lceres -lglog -lboost_program_options -lboost_mpi -lboost_serialization -O3

BUILD = build/
vpath %.h include
vpath %.cpp src

INSTALL_LOC = ../../bin/latticeDNAOrigami
TARGET = bin/latticeDNAOrigami

all: $(TARGET)

$(TARGET): $(BUILD)main.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)movetypes.o $(BUILD)random_gens.o $(BUILD)ideal_random_walk.o $(BUILD)parser.o $(BUILD)order_params.o $(BUILD)enumerate.o $(BUILD)constant_temp_simulation.o $(BUILD)annealing_simulation.o $(BUILD)ptmc_simulation.o $(BUILD)us_simulation.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(TESTTARGET): $(BUILD)test.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)movetypes.o $(BUILD)random_gens.o $(BUILD)ideal_random_walk.o $(BUILD)parser.o $(BUID)order_params.o
	$(CC) -o $@ $^ $(LDFLAGS)
	$(TESTTARGET)

$(BUILD)main.o: main.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)origami_system.o: origami_system.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)utility.o: utility.cpp utility.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)nearest_neighbour.o: nearest_neighbour.cpp nearest_neighbour.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)files.o: files.cpp files.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)json.o: jsoncpp.cpp json/json.h json/json-forwards.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)domain.o: domain.cpp domain.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)simulation.o: simulation.cpp simulation.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)movetypes.o: movetypes.cpp movetypes.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)random_gens.o: random_gens.cpp random_gens.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)ideal_random_walk.o: ideal_random_walk.cpp ideal_random_walk.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)parser.o: parser.cpp parser.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)order_params.o: order_params.cpp order_params.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)enumerate.o: enumerate.cpp enumerate.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)constant_temp_simulation.o: constant_temp_simulation.cpp constant_temp_simulation.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)annealing_simulation.o: annealing_simulation.cpp annealing_simulation.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)ptmc_simulation.o: ptmc_simulation.cpp ptmc_simulation.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)us_simulation.o: us_simulation.cpp us_simulation.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean install
clean:
	rm $(BUILD)*.o

install:
	cp $(TARGET) $(INSTALL_LOC)
