# AUTODEPENDENCY GENERATION DOES NOT WORK
# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/
# It is surprising because the guy who wrote the above post in the make maintainer

OPTLEVEL = -O3
BUILDDIR = build
PREFIX = ~/bin/latticeDNAOrigami
TARGET = latticeDNAOrigami
TARGETDIR = bin
SRCDIR = src
INCLUDEDIR = include

vpath %.h $(INCLUDEDIR)/
vpath %.cpp $(SRCDIR)/

# Put git commit hash in code
$(shell echo -e "#include \"version.h\"\n\nchar const *const GIT_COMMIT = \"$$(git rev-parse HEAD)\";" > src/version.cpp.tmp; if diff -q src/version.cpp.tmp src/version.cpp >/dev/null 2>&1; then rm src/version.cpp.tmp; else mv src/version.cpp.tmp src/version.cpp; fi)

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
	$(CPP) -o $@ $^ $(LDFLAGS)

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
