# type make. needs gtk3 (same as a1) + umfpack.
# umfpack: run  bash build_umfpack.sh  once first (builds it from the tarball into ./suitesparse, no blas needed)
CXX = g++
FLAGS = -O2 -Wall -std=c++17 -I. -Iezgl -Isuitesparse/include
GTK = $(shell pkg-config --cflags --libs gtk+-3.0 x11)
UMF = suitesparse/lib/libumfpack.a suitesparse/lib/libamd.a suitesparse/lib/libsuitesparseconfig.a -lm

placer: placer.cpp ezgl/*.cpp suitesparse/lib/libumfpack.a
	$(CXX) $(FLAGS) placer.cpp ezgl/*.cpp $(GTK) $(UMF) -o placer

# no graphics version (doesnt need gtk), for quick runs
placer_nogfx: placer.cpp suitesparse/lib/libumfpack.a
	$(CXX) $(FLAGS) -DNOGFX placer.cpp $(UMF) -o placer_nogfx

suitesparse/lib/libumfpack.a:
	bash build_umfpack.sh

clean:
	rm -f placer placer_nogfx
