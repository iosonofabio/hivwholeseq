CPPFLAGS=-I/ebio/ag-neher/share/programs/include

all: seqan align_global

seqan: main.cpp
	g++ $(CPPFLAGS) main.cpp -o seqan

align_global: align_global.cpp
	g++ $(CPPFLAGS) $^ -o $@

clean:
	rm -rf seqan
