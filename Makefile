CXX=g++
CXXFLAGS=-std=c++11 -Wall -O3 -fopenmp -msse4.1
LDFLAGS=-lm -lz

cpp_source=sequence_batch.cc index.cc candidate_processor.cc feature_barcode_matrix.cc ksw.cc output_tools.cc chromap.cc
src_dir=src
objs_dir=objs
objs+=$(patsubst %.cc,$(objs_dir)/%.o,$(cpp_source))

exec=chromap

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

all: dir $(exec) 
	
dir:
	mkdir -p $(objs_dir)

$(exec): $(objs)
	$(CXX) $(CXXFLAGS) $(objs) -o $(exec) $(LDFLAGS)
	
$(objs_dir)/%.o: $(src_dir)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm -rf $(exec) $(objs_dir)
