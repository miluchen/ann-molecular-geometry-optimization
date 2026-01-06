CXX = g++
CXXFLAGS = -std=c++11 -O3 -DADD_
INCLUDES = -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
LDFLAGS = -framework Accelerate
TARGET = chemnet_main
SOURCES = chemnet_main.cpp chemnet.cpp

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) $(LDFLAGS) -o $(TARGET)

.PHONY: clean run
clean:
	rm -f $(TARGET) *.o

run: $(TARGET)
	./$(TARGET)
