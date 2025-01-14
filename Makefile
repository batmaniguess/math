CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

SRC_DIR = src
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
BINS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%)

all: $(BUILD_DIR) $(BINS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR)/*

rebuild: clean all

.PHONY: all clean rebuild
