#include"Kseq_Cpp.h"
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

Kseq_Cpp::Kseq_Cpp(std::string file_path,unsigned int buffer_size){
    _fp=gzopen(file_path.c_str(), "r");
    seq = kseq_init(_fp);
}
Kseq_Cpp::~Kseq_Cpp(){
    kseq_destroy(seq);
    gzclose(_fp);
}
bool Kseq_Cpp::get_seq(){
    return kseq_read(seq)>=0?true:false;
}