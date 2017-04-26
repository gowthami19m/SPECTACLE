CC=g++
CFLAGS=-Wall -O3 -std=c++11 -I ./boost/include
LDFLAGS=-std=c++11 ./boost/lib/libboost_iostreams.a ./zlib/install/lib/libz.a
SRC_DIR=src
LIB_DIR=lib
BIN_DIR=bin
ZLIB=ZLIB

all: $(ZLIB) generate-a-single generate-q-single reconstruct q-to-q-paired q-to-q-single q-to-a-paired q-to-a-single remove-postfix-lsc remove-postfix-proovread sam-paired evaluate

generate-a-single: $(SRC_DIR)/generate-map.from-fasta.single.common.o
	$(CC) $(SRC_DIR)/generate-map.from-fasta.single.common.o $(LDFLAGS) -o $(BIN_DIR)/generate-map.from-fasta.single.common

generate-q-single: $(SRC_DIR)/generate-map.from-fastq.single.common.o
	$(CC) $(SRC_DIR)/generate-map.from-fastq.single.common.o $(LDFLAGS) -o $(BIN_DIR)/generate-map.from-fastq.single.common

reconstruct: $(SRC_DIR)/reconstruct-genome.rna.o
	$(CC) $(SRC_DIR)/reconstruct-genome.rna.o $(LDFLAGS) -o $(BIN_DIR)/reconstruct-genome.rna

q-to-q-paired: $(SRC_DIR)/write-order-file.from-fastq.to-fastq.paired.common.o
	$(CC) $(SRC_DIR)/write-order-file.from-fastq.to-fastq.paired.common.o $(LDFLAGS) -o $(BIN_DIR)/write-order-file.from-fastq.to-fastq.paired.common

q-to-q-single: $(SRC_DIR)/write-order-file.from-fastq.to-fastq.single.common.o
	$(CC) $(SRC_DIR)/write-order-file.from-fastq.to-fastq.single.common.o $(LDFLAGS) -o $(BIN_DIR)/write-order-file.from-fastq.to-fastq.single.common

q-to-a-paired: $(SRC_DIR)/write-order-file.from-fastq.to-fasta.paired.common.o
	$(CC) $(SRC_DIR)/write-order-file.from-fastq.to-fasta.paired.common.o $(LDFLAGS) -o $(BIN_DIR)/write-order-file.from-fastq.to-fasta.paired.common

q-to-a-single: $(SRC_DIR)/write-order-file.from-fastq.to-fasta.single.common.o
	$(CC) $(SRC_DIR)/write-order-file.from-fastq.to-fasta.single.common.o $(LDFLAGS) -o $(BIN_DIR)/write-order-file.from-fastq.to-fasta.single.common

remove-postfix-lsc: $(SRC_DIR)/remove-postfix.fasta.single.lsc.o
	$(CC) $(SRC_DIR)/remove-postfix.fasta.single.lsc.cpp $(LDFLAGS) -o $(BIN_DIR)/remove-postfix.fasta.single.lsc

remove-postfix-proovread: $(SRC_DIR)/remove-postfix.fasta.single.proovread.o
	$(CC) $(SRC_DIR)/remove-postfix.fasta.single.proovread.cpp $(LDFLAGS) -o $(BIN_DIR)/remove-postfix.fasta.single.proovread

sam-paired: $(SRC_DIR)/write-order-file.sam.paired.common.o
	$(CC) $(SRC_DIR)/write-order-file.sam.paired.common.o $(LDFLAGS) -o $(BIN_DIR)/write-order-file.sam.paired.common

$(SRC_DIR)/generate-map.from-fasta.single.common.o: $(SRC_DIR)/generate-map.from-fasta.single.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/generate-map.from-fastq.single.common.o: $(SRC_DIR)/generate-map.from-fastq.single.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/reconstruct-genome.rna.o: $(SRC_DIR)/reconstruct-genome.rna.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/remove-postfix.fasta.single.lsc.o: $(SRC_DIR)/remove-postfix.fasta.single.lsc.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/remove-postfix.fasta.single.proovread.o: $(SRC_DIR)/remove-postfix.fasta.single.proovread.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/write-order-file.from-fastq.to-fastq.paired.common.o: $(SRC_DIR)/write-order-file.from-fastq.to-fastq.paired.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/write-order-file.from-fastq.to-fastq.single.common.o: $(SRC_DIR)/write-order-file.from-fastq.to-fastq.single.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/write-order-file.from-fastq.to-fasta.paired.common.o: $(SRC_DIR)/write-order-file.from-fastq.to-fasta.paired.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/write-order-file.from-fastq.to-fasta.single.common.o: $(SRC_DIR)/write-order-file.from-fastq.to-fasta.single.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/write-order-file.sam.paired.common.o: $(SRC_DIR)/write-order-file.sam.paired.common.cpp
	$(CC) $(CFLAGS) -c -o $@ $?

$(SRC_DIR)/evaluate.o: $(SRC_DIR)/evaluate.cpp
	$(CC) -O3 -std=c++11 -c `perl -MConfig -e 'print join(" ", @Config{qw(ccflags optimize cccdlflags)}, "-I$$Config{archlib}/CORE")'` -o $@ $?

$(SRC_DIR)/evaluate-wrap.o: swig
	$(CC) -O3 -std=c++11 -c `perl -MConfig -e 'print join(" ", @Config{qw(ccflags optimize cccdlflags)}, "-I$$Config{archlib}/CORE")'` -o $@ $(SRC_DIR)/evaluate-wrap.cpp

swig:
	swig -perl5 -c++ -o $(SRC_DIR)/evaluate-wrap.cpp $(LIB_DIR)/evaluate.i
	mv $(SRC_DIR)/evaluate.pm $(LIB_DIR)

$(ZLIB):
	cd zlib; ./compile

evaluate: $(SRC_DIR)/evaluate.o $(SRC_DIR)/evaluate-wrap.o
	$(CC) -O3 -std=c++11 `perl -MConfig -e 'print $$Config{lddlflags}'` -o $(LIB_DIR)/evaluate.so $?
	chmod 644 $(LIB_DIR)/evaluate.so
	cd ncurses; ./compile; cd ..
	cd samtools; ./compile

clean:
	rm -f $(BIN_DIR)/generate-map.from-fasta.single.common
	rm -f $(BIN_DIR)/generate-map.from-fastq.single.common
	rm -f $(BIN_DIR)/reconstruct-genome.rna
	rm -f $(BIN_DIR)/remove-postfix.fasta.single.lsc
	rm -f $(BIN_DIR)/remove-postfix.fasta.single.proovread
	rm -f $(BIN_DIR)/write-order-file.from-fastq.to-fastq.paired.common
	rm -f $(BIN_DIR)/write-order-file.from-fastq.to-fastq.single.common
	rm -f $(BIN_DIR)/write-order-file.from-fastq.to-fasta.paired.common
	rm -f $(BIN_DIR)/write-order-file.from-fastq.to-fasta.single.common
	rm -f $(BIN_DIR)/write-order-file.sam.paired.common
	rm -f $(SRC_DIR)/*.o
	rm -f $(LIB_DIR)/evaluate.so
	rm -f $(SRC_DIR)/evaluate-wrap.cpp
	rm -f $(LIB_DIR)/evaluate.pm
	rm -r samtools/install;cd samtools/samtools-1.2; make clean
	cd zlib; rm -rf install; cd zlib-1.2.8; make clean
