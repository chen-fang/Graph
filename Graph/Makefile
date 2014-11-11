GCC=g++
GFLAGS=-Wall -g -std=c++11 -O3 -DNDEBUG -w


#SUPERLUINC= -I../../SuperLU_MT_2.1
ADETLINC = -I /home/blackicefc/ADETL_DISTRO
#ADETLINC = -I /home/chenfang/Library/AD_Library/ADETL_DISTRO

SRC = ./*.cpp

#LIB= ../../SuperLU_MT_2.1/lib/libsuperlu_mt_PTHREAD.a ../../SuperLU_MT_2.1/lib/libblas_PTHREAD.a -lm -lpthread

#LDFLAGS   :=../../SuperLU_MT_2.1/lib/libsuperlu_mt_PTHREAD.a ../../SuperLU_MT_2.1/lib/libblas_PTHREAD.a -lm -lpthread

OUT= run_ex1.out


g:
	$(GCC) $(GFLAGS) $(SUPERLUINC) $(ADETLINC) $(LIB) $(SRCDIR) $(SRC) $(LDFLAGS)  -o $(OUT)

clean:
	rm -rf $(OUT) ./*~ ./*.o RECORD_RESULT.txt

