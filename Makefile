# AUTODEPENDENCY GENERATION DOES NOT WORK
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# It is surprising because the guy who wrote the above post in the make maintainer

OPTLEVEL = -g
BUILDDIR = build
PREFIX = ~/bin/latticeDNAOrigami
TARGET = latticeDNAOrigami
TARGETDIR = bin
SRCDIR = src
INCLUDEDIR = include

vpath %.h $(INCLUDEDIR)/
vpath %.cpp $(SRCDIR)/

SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(subst .cpp,.o,$(SOURCES))
OBJECTS := $(subst $(SRCDIR),$(BUILDDIR),$(OBJECTS))

CPP = mpicxx
CPPFLAGS = -std=c++14 -I include $(OPTLEVEL)
LDFLAGS = -lboost_program_options -lboost_mpi -lboost_serialization -lboost_system -lboost_filesystem $(OPTLEVEL)

# For compiling on clusters with local Boost installation
#CPPFLAGS = -std=c++14 -I/home/amc226/include -Iinclude $(OPTLEVEL)
#LDFLAGS = -L/home/amc226/lib -lboost_program_options -lboost_mpi -lboost_serialization -lboost_system -lboost_filesystem $(OPTLEVEL)

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

COMPILE.cpp = $(CPP) $(DEPFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: $(TARGETDIR)/$(TARGET)

$(TARGETDIR)/$(TARGET): $(OBJECTS)
	$(CPP) $(LDFLAGS) -o $@ $^

$(BUILDDIR)/%.o: %.cpp
$(BUILDDIR)/%.o: %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

.PHONY: clean install
clean:
	rm $(BUILDDIR)/*.o
	rm .d/*.d

install:
	cp $(TARGETDIR)/$(TARGET) $(PREFIX)

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SOURCES))))
