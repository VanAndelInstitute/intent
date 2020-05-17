
CXX  := c++ 
CXXFLAGS := -std=c++17
APP_DIR := ./bin
LDFLAGS  := -lz
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
TARGET   := intent
INCLUDE  := -Iinclude
SRC      := main.cpp \
   $(wildcard src/*.cpp) 

# on macOS, sanitize=leak not supported out of the box by xcode. 
# consider:
# brew install llvm@8make CXX=/usr/local/opt/llvm/bin/clang++ debug
# make CXX=/usr/local/opt/llvm@8/bin/clang++ debug
SANITIZE_LEAK_SUPPORT := $(shell touch _.c && $(CXX) -fsanitize=leak -c _.c -o _.o &> /dev/null && echo 1; rm -f _.c _.o)
ifeq ($(SANITIZE_LEAK_SUPPORT), 1)
	DEBUGFLAGS=-g -fsanitize=leak
else
	DEBUGFLAGS=-g -fsanitize=address
endif

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean all test

clean:
	-@rm -rvf $(OBJ_DIR)/*

install:
	cp $(APP_DIR)/$(TARGET) /usr/local/bin

test:
	$(APP_DIR)/$(TARGET) ./test/test_R1.fastq.gz ./test/test_R2.fastq.gz	
