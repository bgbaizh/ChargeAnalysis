objs= ChargeAna.o ChargeAna_binding.o
soft=g++
libdir=-L/home/bzhao/.conda/envs/bzh/lib/
lib=-lpython3
includedir= -I /home/bzhao/soft/pyscal/lib/pybind11/include/ -I /home/bzhao/.conda/envs/bzh/include/python3.9/

cChargeAnalysis.so:${objs}
	${soft} -o $@  -shared  $^ ${libdir} ${lib} -O3 -Wall -fopenmp -std=c++2a -shared -fPIC
ChargeAna.o:ChargeAna.cpp ChargeAna.h
	${soft} -O3 -Wall -fopenmp -std=c++2a -c $< ${includedir} -o $@  -shared -fPIC
ChargeAna_binding.o:ChargeAna_binding.cpp ChargeAna.h
	${soft} -O3 -Wall  -fopenmp -std=c++2a -c $< ${includedir} -o $@ -shared -fPIC 
clean:
	rm cChargeAna.so ${objs}

