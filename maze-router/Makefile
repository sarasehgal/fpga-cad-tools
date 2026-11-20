# just type make. needs gtk3 (same as ezgl example)
CXX = g++
FLAGS = -O2 -Wall -std=c++14 -I. -Iezgl
GTK = $(shell pkg-config --cflags --libs gtk+-3.0 x11)

router: router.cpp ezgl/*.cpp
	$(CXX) $(FLAGS) router.cpp ezgl/*.cpp $(GTK) -o router

clean:
	rm -f router
