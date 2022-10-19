#ifndef KSEQ_CPP_H
#define KSEQ_CPP_H
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <string>
KSEQ_INIT(gzFile, gzread)

using record_type=kseq_t;
class Kseq_Cpp
{
public:
    Kseq_Cpp()=delete;
    explicit Kseq_Cpp(std::string  file_path,unsigned int buffer_size=16384);
    Kseq_Cpp(Kseq_Cpp&)=delete;
    ~Kseq_Cpp();

public:
    bool get_seq();
    record_type * seq;

private:
    gzFile _fp;
    
};
#endif // KSEQ_CPP_H
