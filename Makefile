CC = g++
CPPFLAGS = -I include -O3
LDFLAGS = -O3

BUILD = build/
vpath %.h include
vpath %.cpp src

TARGET = bin/latticeDNAOrigami

all: $(TARGET)

$(TARGET): $(BUILD)main.o $(BUILD)origami_system.o $(BUILD)utility.o $(BUILD)nearest_neighbour.o
	$(CC) -o $@ $^ $(LDFLAGS)

$(BUILD)main.o: main.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)origami_system.o: origami_system.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)utility.o: utility.cpp utility.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

$(BUILD)nearest_neighbour.o: nearest_neighbour.cpp nearest_neighbour.h
	$(CC) -o $@ -c $(CPPFLAGS) $<

.PHONY: clean
clean:
	rm $(BUILD)*.o
