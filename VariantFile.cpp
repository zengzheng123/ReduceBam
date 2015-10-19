//
//  VariantFile.cpp
//  ReduceBam
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "VariantFile.h"


VcfFile::VcfFile(string input_vcf_file)
{
    if(input_vcf_file.length() > 7 && input_vcf_file.substr(input_vcf_file.length() - 7) == ".vcf.gz")
    {
        cout << "[INFO] Input VCF is gzip file" << endl;
        vcf_gz_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_gz_fs;
        file_type = "vcf";
    }
    else if(input_vcf_file.length() > 4 && input_vcf_file.substr(input_vcf_file.length() - 4) == ".vcf")
    {
        cout << "[INFO] Input VCF is plain text file" << endl;
        vcf_txt_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_txt_fs;
        file_type = "vcf";
    }
    else if(input_vcf_file.length() > 7 && input_vcf_file.substr(input_vcf_file.length() - 7) == ".maf.gz")
    {
        cout << "[INFO] Input MAF is gzip file" << endl;
        vcf_gz_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_gz_fs;
        file_type = "maf";
    }
    else if(input_vcf_file.length() > 4 && input_vcf_file.substr(input_vcf_file.length() - 4) == ".maf")
    {
        cout << "[INFO] Input MAF is plain text file" << endl;
        vcf_txt_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_txt_fs;
        file_type = "maf";
    }
    else
    {
        cerr << "[ERROR] Input variant is in unrecognized format: " << input_vcf_file << endl;
        exit(1);
    }
    if(!(*vcf_fs))
    {
        cerr << "[ERROR] Fail to open input target file: " << input_vcf_file << endl;
        exit(1);
    }
}

VcfFile::~VcfFile()
{
    close();
}

bool VcfFile::get_next(string &line)
{
    if(cur_line.empty())
    {
        if(!getline((*vcf_fs), line))
            return false;
    }
    else
    {
        line = cur_line;
        cur_line.clear();
    }
    return true;
}

void VcfFile::roll_back(string &line)
{
    cur_line = line;
}

void VcfFile::close()
{
    if(vcf_txt_fs)
        vcf_txt_fs.close();
    if(vcf_gz_fs)
        vcf_gz_fs.close();
}

bool VcfFile::eof()
{
    return vcf_fs->eof();
}

