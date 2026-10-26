# GNU Make; GCC with system LAPACK/BLAS.
# Run from this directory (containing src/ and include/).
#   make -f Makefile.gnu -j8
#   make -f Makefile.gnu debug-parallel
#   make -f Makefile.gnu release-serial LTO=full NATIVE=1
#   make -f Makefile.gnu CXX=g++-16      # Homebrew GCC on macOS
# Defaults: -O2, OpenMP, no LTO, portable CPU target.
# On macOS, /usr/bin/g++ is Apple Clang; select an installed GCC explicitly.
# Compiler upgrades / loaded module changes require a clean rebuild.
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

# Capture the main makefile before dependency files are included.
THIS_MAKEFILE := $(abspath $(firstword $(MAKEFILE_LIST)))

# Replace Make's built-in CXX default while preserving user overrides.
ifeq ($(origin CXX),default)
CXX := g++
endif
CXX ?= g++

# Each configuration has separate objects and an executable. GCC supports
# full LTO, so the Clang/Intel LTO=thin setting is not accepted here.
MODE   ?= release
OPENMP ?= 1
LTO    ?= off
NATIVE ?= 0

ifneq ($(words $(MODE)),1)
$(error MODE must be release or debug)
endif
ifeq ($(filter release debug,$(MODE)),)
$(error MODE must be release or debug)
endif
ifneq ($(words $(OPENMP)),1)
$(error OPENMP must be 0 or 1)
endif
ifeq ($(filter 0 1,$(OPENMP)),)
$(error OPENMP must be 0 or 1)
endif
ifneq ($(words $(LTO)),1)
$(error LTO must be off or full; GCC does not support ThinLTO)
endif
ifeq ($(filter off full,$(LTO)),)
$(error LTO must be off or full; GCC does not support ThinLTO)
endif
ifneq ($(words $(NATIVE)),1)
$(error NATIVE must be 0 or 1)
endif
ifeq ($(filter 0 1,$(NATIVE)),)
$(error NATIVE must be 0 or 1)
endif

# Keep GCC objects and binaries separate from Clang and Intel outputs.
SRC_DIR    := src
BIN_DIR    := bin/gnu
BUILD_ROOT := build/gnu
BUILD_TAG  := $(MODE)-omp$(OPENMP)-lto$(LTO)-native$(NATIVE)
BUILD_DIR  := $(BUILD_ROOT)/$(BUILD_TAG)
OBJ_DIR    := $(BUILD_DIR)/obj
DEP_DIR    := $(BUILD_DIR)/dep
EXE        := $(BIN_DIR)/$(BUILD_TAG)/AS3
CONFIG     := $(BUILD_DIR)/flags.txt

SRCS := $(sort $(wildcard $(SRC_DIR)/*.cpp))
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))
DEPS := $(patsubst $(SRC_DIR)/%.cpp,$(DEP_DIR)/%.d,$(SRCS))

# User flags are separate from required build flags. Examples:
# make CXXFLAGS=-Wshadow
# make EXTRA_CXXFLAGS=-pg EXTRA_LDFLAGS=-pg
# LAPACK/BLAS must use the LP64 interface (32-bit Fortran integers).
CPPFLAGS ?=
CXXFLAGS ?=
LDFLAGS  ?=
LDLIBS   ?= -llapack -lblas
EXTRA_CPPFLAGS ?=
EXTRA_CXXFLAGS ?=
EXTRA_LDFLAGS  ?=

PROJECT_CPPFLAGS := -Iinclude
PROJECT_CXXFLAGS := -std=c++20 -Wall -Wextra -Wconversion
PROJECT_LDFLAGS  :=

ifeq ($(MODE),debug)
# Enable assertions and solver NaN checks without dropping OpenMP settings.
PROJECT_CPPFLAGS += -UNDEBUG -DDEBUG -DENABLE_NAN_CHECK
PROJECT_CXXFLAGS += -O0 -g
else
PROJECT_CPPFLAGS += -DNDEBUG
PROJECT_CXXFLAGS += -O2
endif

ifeq ($(NATIVE),1)
# Targets the build host's instruction set: use only on compatible CPUs.
# AArch64 selects architecture and tuning together with -mcpu; this also
# avoids GCC's conflicting native architecture/tuning expansion on macOS.
GCC_TARGET := $(shell $(CXX) -dumpmachine)
ifneq ($(filter aarch64% arm64%,$(GCC_TARGET)),)
PROJECT_CXXFLAGS += -mcpu=native
else
PROJECT_CXXFLAGS += -march=native -mtune=native
endif
endif

ifeq ($(OPENMP),1)
# GCC supplies its own OpenMP headers and libgomp runtime. Enable both the
# solver's guarded code and compiler support, including at link time.
PROJECT_CPPFLAGS += -DHAVE_OPENMP
PROJECT_CXXFLAGS += -fopenmp
PROJECT_LDFLAGS  += -fopenmp
else
PROJECT_CPPFLAGS += -UHAVE_OPENMP
endif

ifeq ($(LTO),full)
# GCC uses -flto, not Clang's -flto=full. Generate LTO data while compiling
# and enable interprocedural optimization when linking through g++.
PROJECT_CXXFLAGS += -flto
PROJECT_LDFLAGS  += -flto
endif

COMPILE_CPPFLAGS = $(CPPFLAGS) $(EXTRA_CPPFLAGS) $(PROJECT_CPPFLAGS)
COMPILE_CXXFLAGS = $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(PROJECT_CXXFLAGS)
LINK_FLAGS       = $(LDFLAGS) $(EXTRA_LDFLAGS) $(PROJECT_LDFLAGS)

.PHONY: all debug release debug-serial debug-parallel release-serial \
        release-parallel print clean FORCE

all: $(EXE)

# Select the configuration before parsing its output paths. Recursive make
# inherits command-line overrides and the parallel jobserver.
debug:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=debug all
release:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=release all
debug-serial:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=debug OPENMP=0 all
debug-parallel:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=debug OPENMP=1 all
release-serial:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=release OPENMP=0 all
release-parallel:
	+$(MAKE) --file="$(THIS_MAKEFILE)" MODE=release OPENMP=1 all

# Track effective flags and source membership so incremental builds cannot
# reuse objects after a configuration change. Keep the timestamp unchanged
# when settings match. Flag values should not contain literal apostrophes.
$(CONFIG): FORCE
	@mkdir -p $(BUILD_DIR)
	@printf '%s\n' 'CXX=$(CXX)' \
	    'CPPFLAGS=$(COMPILE_CPPFLAGS)' \
	    'CXXFLAGS=$(COMPILE_CXXFLAGS)' \
	    'LDFLAGS=$(LINK_FLAGS)' 'LDLIBS=$(LDLIBS)' \
	    'SRCS=$(SRCS)' > $@.tmp
	@cmp -s $@.tmp $@ && rm -f $@.tmp || mv -f $@.tmp $@

# Dependency files include project headers and .inl files; -MP tolerates
# headers that have been removed or renamed since the previous build.
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(CONFIG)
	@mkdir -p $(OBJ_DIR) $(DEP_DIR)
	$(CXX) $(COMPILE_CPPFLAGS) $(COMPILE_CXXFLAGS) \
	    -MMD -MP -MF $(DEP_DIR)/$*.d -MT $@ -c $< -o $@

# Keep libraries after objects for linkers that resolve symbols left to right.
$(EXE): $(OBJS) $(CONFIG)
	@mkdir -p $(dir $@)
	@printf '%s\n' '$(strip $(CXX) $(COMPILE_CXXFLAGS) $(LINK_FLAGS)) \' \
	    $(foreach obj,$(OBJS),'    $(obj) \') \
	    '    $(LDLIBS) -o $@'
	@$(CXX) $(COMPILE_CXXFLAGS) $(LINK_FLAGS) $(OBJS) $(LDLIBS) -o $@

print:
	@printf '%s\n' 'Executable: $(EXE)' 'Compiler: $(CXX)' \
	    'CPPFLAGS: $(COMPILE_CPPFLAGS)' 'CXXFLAGS: $(COMPILE_CXXFLAGS)' \
	    'LDFLAGS: $(LINK_FLAGS)' 'LDLIBS: $(LDLIBS)' \
	    'Sources: $(SRCS)' 'Objects: $(OBJS)' 'Dependencies: $(DEPS)'

# Remove only GCC outputs. Legacy obj/, dep/ and bin/AS3 are left alone.
clean:
	$(RM) -r $(BUILD_ROOT) $(BIN_DIR)

FORCE:

-include $(DEPS)
