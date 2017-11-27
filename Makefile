CC := g++

SRCDIR := src
TARGETDIR := bin
BUILDDIR := build
TARGET := bin/ASElux

# Normal CPP files
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

$(info $(SOURCES))
$(info $(OBJECTS))

CFLAGS := --std=c++11 -pthread -O3

$(TARGET): $(OBJECTS)
	mkdir -p $(TARGETDIR)
	echo " Linking..."
	echo " $(CC) $^ -o $(TARGET)"; $(CC) $(CFLAGS) $^ -o $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	mkdir -p $(BUILDDIR)
	echo "Compiling $<..."; $(CC) $(CFLAGS) -c -o $@ $<

