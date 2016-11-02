#CC = g++
CC = mpicxx
#CPPFLAGS = -I ../../include -g
#LDFLAGS =  -lboost_program_options -g
CPPFLAGS = -I ../../include -O3
LDFLAGS =  -lboost_program_options -lboost_mpi -lboost_serialization -O3
BUILD = ../../build/
vpath %.h ../../include
#vpath %.cpp src

TARGET = calc_num_walks

all: $(TARGET)

testing: $(TESTTARGET)

$(TARGET): calc_num_walks.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)parser.o $(BUILD)random_gens.o $(BUILD)movetypes.o $(BUILD)ideal_random_walk.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(BUILD)calc_num_walks.o: calc_num_walks.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean
clean:
	rm *.o
