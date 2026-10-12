#include"Kseq_Cpp.h"
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <stdexcept>

Kseq_Cpp::Kseq_Cpp(std::string file_path,unsigned int buffer_size){
    _fp=gzopen(file_path.c_str(), "r");
    if (!_fp)
        throw std::runtime_error{"Could not open input file: " + file_path};
    seq = kseq_init(_fp);
}
Kseq_Cpp::~Kseq_Cpp(){
    kseq_destroy(seq);
    gzclose(_fp);
}
bool Kseq_Cpp::get_seq(){
    int const result = kseq_read(seq);
    if (result < -1)
        throw std::runtime_error{"Failed to read sequence: malformed record or stream error."};
    return result >= 0;
}
