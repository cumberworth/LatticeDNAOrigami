CC = g++
CPPFLAGS = -O3 -I include
LDFLAGS =

BUILD = build/
vpath %.h include
vpath %.cpp src

TARGET = bin/latticeDNAOrigami
TESTTARGET = tests/testLatticeDNAOrigami

all: $(TARGET)

testing: $(TESTTARGET)

$(TARGET): $(BUILD)main.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)movetypes.o $(BUILD)random_gens.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(TESTTARGET): $(BUILD)test.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o $(BUILD)files.o $(BUILD)json.o $(BUILD)domain.o $(BUILD)simulation.o $(BUILD)movetypes.o
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

$(BUILD)test.o: test.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean
clean:
	rm $(BUILD)*.o
