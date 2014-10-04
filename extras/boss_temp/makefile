all: minimal_1_0 minimal_1_1 minimal_1_2

minimal_1_0: libminimal_1_0.so
libminimal_1_0.so:
	g++ -shared minimal_1_0/*.cpp -I minimal_1_0 -fPIC -o libminimal_1_0.so

minimal_1_1: libminimal_1_1.so
libminimal_1_1.so:
	g++ -shared minimal_1_1/*.cpp -I minimal_1_1 -fPIC -o libminimal_1_1.so

minimal_1_2: libminimal_1_2.so
libminimal_1_2.so:
	g++ -shared minimal_1_2/*.cpp -I minimal_1_2 -fPIC -o libminimal_1_2.so

clean:
	rm -rf libminimal_*.so
