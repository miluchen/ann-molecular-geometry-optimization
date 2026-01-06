ann:
	g++ -c -O3 -DADD_ -I../CBLAS/src main.cpp chemnet.cpp -std=c++11
	g++ -o test main.o chemnet.o ../CBLAS/lib/cblas_LINUX.a -lblas -std=c++11
