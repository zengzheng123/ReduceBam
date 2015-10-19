//
//  VariantFile.h
//  ReduceBam
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __ReduceBam__VariantFile__
#define __ReduceBam__VariantFile__

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include "gzstream.h"
using namespace std;

class VcfFile
{
public:
    
    VcfFile(string input_vcf_file);
    
    ~VcfFile();
    
    bool get_next(string &line);
    
    void roll_back(string &line);
    
    void close();
    
    bool eof();
    
    ifstream vcf_txt_fs;
    igzstream vcf_gz_fs;
    istream* vcf_fs;
    string cur_line;
    string file_type;
};


#endif /* defined(__ReduceBam__VariantFile__) */
