cpp_source=sequence_batch.cc index.cc ksw.cc chromap.cc
src_dir=src
objs_dir=objs
objs+=$(patsubst %.cc,$(objs_dir)/%.o,$(cpp_source))

cxx=g++
cxxflags=-std=c++11 -Wall -O3 -fopenmp -march=native
ldflags=-lm -lz

exec=chromap

all: dir $(exec) 
	
dir:
	mkdir -p $(objs_dir)

$(exec): $(objs)
	$(cxx) $(cxxflags) $(objs) -o $(exec) $(ldflags)
	
$(objs_dir)/%.o: $(src_dir)/%.cc
	$(cxx) $(cxxflags) -c $< -o $@ $(ldflags)

.PHONY: clean
clean:
	-rm -r $(exec) $(objs_dir)
