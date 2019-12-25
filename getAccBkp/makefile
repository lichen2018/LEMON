objects= get_acc_bkp.o BamHeader.o BamReader.o BamRecord.o BamWalker.o GenomicRegion.o GetAccBkp.o RefGenome.o ssw_cpp.o ssw.o
get_acc_bkp:$(objects)
	g++ -g -Wall -o get_acc_bkp $(objects) -I/home/lichen/htslib/include -L/home/lichen/htslib/lib -lhts
ssw_cpp.o : ssw_cpp.h
	g++ -g -Wall -c ssw_cpp.cpp -std=c++11 -o ssw_cpp.o
ssw.o : ssw.h
	gcc -g -Wall -c ssw.c -std=c++11 -o ssw.o
get_acc_bkp.o : GetAccBkp.h
	g++ -g -Wall -c get_acc_bkp.cpp -o get_acc_bkp.o -std=c++11  
RefGenome.o : RefGenome.h
	g++ -g -Wall -c RefGenome.cpp -std=c++11 -o RefGenome.o 
BamHeader.o :  BamHeader.h
	g++ -g -Wall -c BamHeader.cpp -std=c++11 -o BamHeader.o 
BamReader.o : ThreadPool.h GenomicRegionCollection.h
	g++ -g -Wall -c BamReader.cpp -std=c++11 -o BamReader.o 
BamRecord.o : BamHeader.h
	g++ -g -Wall -c BamRecord.cpp -std=c++11 -o BamRecord.o 
BamWalker.o : BamReader.h
	g++ -g -Wall -c BamWalker.cpp -std=c++11 -o BamWalker.o 
GenomicRegion.o : BamHeader.h
	g++ -g -Wall -c GenomicRegion.cpp -std=c++11 -o GenomicRegion.o 
GetAccBkp.o : BamWalker.h RefGenome.h
	g++ -g -Wall -c GetAccBkp.cpp -std=c++11 -pthread -o GetAccBkp.o


clean :
	rm -rf get_acc_bkp $(objects)
