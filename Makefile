# the desired compiler
CXX = g++

# check operating system
ifeq ($(OS),Windows_NT)
	$(error Windows is not supported by libFWTransform at this time)
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
		EXT = so
		LIB_TYPE = shared
    endif
    ifeq ($(UNAME_S),Darwin)
		EXT = dylib
		LIB_TYPE = dynamiclib
    endif
endif

# name of the shared library
LIB_NAME = libWTools.$(EXT)

# the directory to store the final library
LIB_DIR = lib

# the directory to store any executables
BIN_DIR = bin

# the name of the test bin
TEST_BIN = $(BIN_DIR)/testFWTransform

# the directory to look for library source code
SRC_DIR = src

# the directory to look for test source code
TEST_DIR = test

# the directory to store objects
OBJ_DIR = obj

# the directory to search for includes
INCLUDE_DIR = include

# make the appropriate directories
DIRS = $(shell mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR))

# find all library source code files
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# and the corresponding object files
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# find all test source code files
TEST_SRC = $(wildcard $(TEST_DIR)/*.cpp)

# and the object files for test sources
TEST_OBJ = $(TEST_SRC:$(TEST_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# include the relevant directories
CPPFLAGS = -I$(INCLUDE_DIR)
CPPFLAGS += -isystemtest

# some compiler flags
CFLAGS = -Wall -Wextra -Wshadow -Wold-style-cast -Wcast-align
CFLAGS += -Wpointer-arith -Wundef -Wno-unused -Warray-bounds
CFLAGS += -march=native -O3 -flto -fPIC -Wmissing-braces -Wmaybe-uninitialized
CFLAGS += -Wstrict-overflow=5 -Wno-sign-conversion -Wmisleading-indentation
CFLAGS += -fmax-errors=5

# and some dependency generation flags
DEP_FLAGS = -MMD -MP

# the libraries we want to link against
LIBS = -lfftw3 #/usr/lib/x86_64-linux-gnu/libfftw3.so

# the primary target builds both the library and the test executables
all: $(DIRS) $(LIB_DIR)/$(LIB_NAME) $(TEST_BIN)

# the shared library target
$(LIB_DIR)/$(LIB_NAME): $(OBJ)
	$(CXX) $(LDFLAGS) -$(LIB_TYPE) $(LIBS) $^ -o $@

# the test binary target
$(TEST_BIN): $(TEST_OBJ)
	$(CXX) $(LDFLAGS) -Llib/ $^ -o $@ -lWTools $(LIBS) -lWTools

# build library object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $< $(DEP_FLAGS) -o $@

# build test binary object files
$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $< $(DEP_FLAGS) -o $@

# and cleanup all targets and dependencies
clean:
	rm -rf $(OBJ) $(TEST_OBJ) $(LIB_DIR)/* $(TEST_BIN) $(OBJ_DIR)/*.d
