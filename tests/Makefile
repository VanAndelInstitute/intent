CXX      := c++
CXXFLAGS := -std=c++17 -lz
APP_DIR := ./
LDFLAGS  := 
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
TARGET   := intent_tests
INCLUDE  := -I../include/
SRC      := $(wildcard ../src/*.cpp) \
			$(wildcard ./*.cpp) 
OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET) run

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

clean:
	-@rm -rvf $(OBJ_DIR)/*

run:
	./intent_tests --verbose
