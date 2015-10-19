//
//  main.cpp
//  ReduceBam
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 10/8/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//


#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "gzstream.h"
#include "VariantFile.h"


using namespace std;
using namespace BamTools;

int opt_buffer_size = 0;
vector<string> opt_input_variant_file;
string opt_input_bam_file;
string opt_output_bam_file;
string opt_output_bed_file;

const string VERSION = "ReduceBam 1.0.0";


bool DEBUG_FLAG = true;


void printUsage(string msg = "")
{
    cout << endl;
    cout << VERSION << endl;
    cout << "Usage: " << endl;
    cout << "\t--input_variant               <string>                     Input variant file, support .vcf/.vcf.gz/.maf/.maf.gz format, can be specified multiple times." << endl;
    cout << "\t--input_bam                   <string>                     Input BAM file" << endl;
    cout << "\t--output_bam                  <string>                     Output BAM file" << endl;
    cout << "\t--output_bed                  <string>                     Optional output BED3 file computed from the input VCF/MAF file and the buffer" << endl;
    cout << "\t--buffer                      <int>                        Buffer size to be added to each side of the VCF/MAF entry. Default is 0" << endl;
    cout << "\t--help                                                     Print command line usage" << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"input_variant",               required_argument,      0,     'v'},
    {"input_bam",                   required_argument,      0,     'i'},
    {"output_bam",                  required_argument,      0,     'o'},
    {"output_bed",                  required_argument,      0,     'd'},
    {"buffer",                      required_argument,      0,     'b'},
    {"help",                        no_argument,            0,     'h'},
    {0, 0, 0, 0}
};


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "v:i:o:d:b:h", long_options, &option_index);
        switch(next_option)
        {
            case 'v':
                opt_input_variant_file.push_back(optarg);
                break;
            case 'i':
                opt_input_bam_file = optarg;
                break;
            case 'o':
                opt_output_bam_file = optarg;
                break;
            case 'd':
                opt_output_bed_file = optarg;
                break;
            case 'b':
                opt_buffer_size = atoi(optarg);
                break;
            case 'h':
                printUsage();
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR] Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);
    
    if(opt_input_variant_file.empty())
        printUsage("[ERROR] Please specify input VCF file or MAF file");
    if(opt_input_bam_file.empty())
        printUsage("[ERROR] Please specify input BAM file");
    if(opt_output_bam_file.empty())
        printUsage("[ERROR] Please specify output BAM file");
    if(opt_buffer_size < 0)
        printUsage("[ERROR] --buffer need to >= 0");
}


void split(const string& line, char delim, vector<std::string>& parsed_item, bool ignore_empty_item = false)
{
    stringstream my_ss(line);
    string item;
    while (getline(my_ss, item, delim))
    {
#ifdef _PARSING_DEBUG
        cerr << "[DEBUG] Parsed item: " << item << endl;
#endif
        if (ignore_empty_item && item.empty())
            continue;    // use to skip empty item
        parsed_item.push_back(item);
    }
}


bool sort_interval_by_pos(const pair<int, int> &lhs, const pair<int, int> &rhs)
{
    if(lhs.first != rhs.first)
        return lhs.first < rhs.first;
    else
        return lhs.second < rhs.second;
}

int check_interval_overlap(int lstart, int lend, int rstart, int rend)
{
    if(lend < rstart)
        return -1; //lhs < rhs
    else if(lstart > rend)
        return 1;  //lhs > rhs
    else
        return 0;  //lhs and rhs overlap

}

int check_interval_overlap(const pair<int, int> &lhs, const pair<int, int> &rhs)
{
    return check_interval_overlap(lhs.first, lhs.second, rhs.first, rhs.second);
}

int check_interval_connected(const pair<int, int> &lhs, const pair<int, int> &rhs)
{
    return check_interval_overlap(lhs.first, lhs.second + 1, rhs.first, rhs.second);
}


void get_vcf_interval(VcfFile &my_vcf, map<string, vector< pair<int, int> > > &unsorted_interval_vec, map<string, int> &chrom_len_map, int buffer_size)
{

    string line;
    while(my_vcf.get_next(line))
    {
        if(line[0] == '#')
            continue;
        int start_pos = 0, end_pos = 0;
        char chrom_tmp[5000], ref_tmp[5000], alt_tmp[5000];
        sscanf(line.c_str(), "%s\t%d\t%*s\t%s\t%s", chrom_tmp, &(start_pos), ref_tmp, alt_tmp);
        string chrom = chrom_tmp;
        string ref = ref_tmp;
        string alt = alt_tmp;
        if(chrom_len_map.find(chrom) == chrom_len_map.end())
        {
            cerr << "[ERROR] Could not find VCF chrom: " << chrom << " in BAM header" << endl;
            exit(1);
        }
        size_t max_interest_length = ref.length();
        vector<string> alt_vec;
        split(alt, ',', alt_vec, true);
        for(size_t alt_index = 0; alt_index < alt_vec.size(); alt_index++)
        {
            if(alt_vec[alt_index].length() > max_interest_length)
                max_interest_length = alt_vec[alt_index].length();
        }
        end_pos = start_pos + (int)max_interest_length - 1;
        int buffered_start_pos = max(1, start_pos - buffer_size) - 1;                      //add buffer, convert to 0-indexed to match BAM entry
        int buffered_end_pos = min(chrom_len_map[chrom], end_pos + buffer_size) - 1;       //add buffer, convert to 0-indexed to match BAM entry
        if(unsorted_interval_vec.find(chrom) == unsorted_interval_vec.end())
        {
            vector< pair<int, int> > new_vec_unsorted;
            unsorted_interval_vec.insert(make_pair(chrom, new_vec_unsorted));
        }
        unsorted_interval_vec[chrom].push_back(make_pair(buffered_start_pos, buffered_end_pos));
    }
}



void get_maf_interval(VcfFile &my_maf, map<string, vector< pair<int, int> > > &unsorted_interval_vec, map<string, int> &chrom_len_map, int buffer_size)
{
    string line;
    while(my_maf.get_next(line))
    {
        if(line[0] != '#') // first non-comment line is the header line
            break;
    }
    while(my_maf.get_next(line))
    {
        
        vector<string> maf_items;
        split(line, '\t', maf_items);
        string chrom = maf_items[4];
        string variant_type = maf_items[9];
        string ref = maf_items[10];
        string alt = maf_items[12];
        if(alt.empty())
            alt = maf_items[11];
        int start_pos = atoi(maf_items[5].c_str());
        int end_pos = atoi(maf_items[6].c_str());
        
        if(chrom_len_map.find(chrom) == chrom_len_map.end())
        {
            cerr << "[ERROR] Could not find MAF chrom: " << chrom << " in BAM header" << endl;
            exit(1);
        }
        if(variant_type == "SNP" || variant_type == "DNP" || variant_type == "TNP" || variant_type == "ONP")
        {
            // do nothing
        }
        else if(variant_type == "INS")
        {
            end_pos += ((int)alt.length() - 1);
        }
        else if(variant_type == "DEL")
        {
            start_pos--;
        }
        else
        {
            cout << "[WARNING] Ignoring variant with unrecognized type: " << variant_type << endl;
            continue;
        }

        int buffered_start_pos = max(1, start_pos - buffer_size) - 1;                      //add buffer, convert to 0-indexed to match BAM entry
        int buffered_end_pos = min(chrom_len_map[chrom], end_pos + buffer_size) - 1;       //add buffer, convert to 0-indexed to match BAM entry
        
        if(unsorted_interval_vec.find(chrom) == unsorted_interval_vec.end())
        {
            vector< pair<int, int> > new_vec_unsorted;
            unsorted_interval_vec.insert(make_pair(chrom, new_vec_unsorted));
        }
        unsorted_interval_vec[chrom].push_back(make_pair(buffered_start_pos, buffered_end_pos));
    }
}


void merge_interval(map<string, vector< pair<int, int> > > &unsorted_interval_vec, map<string, vector< pair<int, int> > > &sorted_interval_vec)
{
    for(map<string, vector< pair<int, int> > >::iterator it = unsorted_interval_vec.begin(); it != unsorted_interval_vec.end(); it++)
    {
        sort(it->second.begin(), it->second.end(), sort_interval_by_pos);
        vector<pair<int, int> > new_vec_sorted;
        sorted_interval_vec.insert(make_pair(it->first, new_vec_sorted));
        for(size_t i = 0; i < it->second.size(); i++)
        {
            if(sorted_interval_vec[it->first].size() == 0 || check_interval_connected(sorted_interval_vec[it->first][sorted_interval_vec[it->first].size() - 1], it->second[i]) != 0)
            {
                sorted_interval_vec[it->first].push_back(it->second[i]);
            }
            else
            {
                sorted_interval_vec[it->first][sorted_interval_vec[it->first].size() - 1].second = it->second[i].second;
            }
        }
    }
}



void output_interval_to_bed(map<string, vector< pair<int, int> > > &sorted_interval_vec, string output_bed_file)
{
    ofstream output_fs(output_bed_file.c_str());
    if(!output_fs)
    {
        cerr << "[ERROR] Fail to write to output bed file: " << output_bed_file << endl;
        exit(1);
    }
    
    for(map<string, vector< pair<int, int> > >::iterator it = sorted_interval_vec.begin(); it != sorted_interval_vec.end(); it++)
    {
        for(size_t i = 0; i < it->second.size(); i++)
        {
            output_fs << it->first << "\t" << it->second[i].first << "\t" << (it->second[i].second + 1) << endl;   //bed file is half-open
        }
    }
    output_fs.close();
    cout << "[INFO] Output BED file complete: " << output_bed_file << endl;
}


bool binary_search_overlap(pair<int, int> &pattern, vector< pair<int, int> > &context_vec)
{
    if(context_vec.size() == 0)
        return false;
    int min_index = 0;
    int max_index = (int)context_vec.size() - 1;
    while (min_index <= max_index)
    {
        int mid_index = (min_index + max_index) / 2;
        int compare_result = check_interval_overlap(pattern, context_vec[mid_index]);
        
        //if(DEBUG_FLAG)
        //{
        //cerr << "[DEBUG] min_index: " << min_index << " max_index: " << max_index << " mid_index: " << mid_index << endl;
        //cerr << "[DEBUG] pattern: (" << pattern.first << "," << pattern.second << "), context: (" << context_vec[mid_index].first << "," << context_vec[mid_index].second << ")" << endl;
        //cerr << "[DEBUG] compare result: " << compare_result << endl;
        //}
        
        if(compare_result == 0)
            return true;
        else if (compare_result == 1)
            min_index = mid_index + 1;
        else
            max_index = mid_index - 1;
    }
    return false;
}


void create_index(string input_bam_file)
{
    cout << "[INFO] Creating index for BAM file: " << input_bam_file << endl;
    BamReader my_bam_reader;
    if(!my_bam_reader.Open(input_bam_file))
    {
        cerr << "[WARNING] Fail to open input bam file for indexing: " << input_bam_file << endl;
        return;
    }
    if(!my_bam_reader.CreateIndex())
    {
        cerr << "[WARNING] Fail to indexing bam: " << input_bam_file << endl;
        return;
    }
    my_bam_reader.Close();
}

void reduce_bam(string input_bam_file, vector<string> &input_variant_file, string output_bam_file, string output_bed_file, int buffer_size)
{
    BamReader my_bam_reader;
    if(!my_bam_reader.Open(input_bam_file))
    {
        cerr << "[ERROR] Fail to open input bam file: " << input_bam_file << endl;
        exit(1);
    }
    string input_bam_index_file1 = input_bam_file + ".bai";
    string input_bam_index_file2 = input_bam_file.substr(0, input_bam_file.length() - 3) + "bai";
    if(!my_bam_reader.OpenIndex(input_bam_index_file1))
    {
        if(!my_bam_reader.OpenIndex(input_bam_index_file2))
        {
            cerr << "[ERROR] Fail to open input bam index file: " << input_bam_index_file1 << " or " << input_bam_index_file2 << endl;
            exit(1);
        }
    }
    
    map<string, int> chrom_len_map;
    const RefVector &reference_data =  my_bam_reader.GetReferenceData();
    for(size_t i = 0; i < reference_data.size(); i++)
    {
        chrom_len_map[reference_data[i].RefName] = reference_data[i].RefLength;
    }
    
    map<string, vector< pair<int, int> > > unsorted_interval_vec;
    map<string, vector< pair<int, int> > > sorted_interval_vec;
    for(size_t i = 0; i < input_variant_file.size(); i++)
    {
        VcfFile my_vf(input_variant_file[i]);
        if(my_vf.file_type == "vcf")
            get_vcf_interval(my_vf, unsorted_interval_vec, chrom_len_map, buffer_size);
        else if(my_vf.file_type == "maf")
            get_maf_interval(my_vf, unsorted_interval_vec, chrom_len_map, buffer_size);
        else
        {
            cerr << "[ERROR]: Input variant is in unrecognized format" << endl;
            exit(1);
        }
        my_vf.close();
        cout << "[INFO] Load VCF file complete: " << input_variant_file[i] << endl;
    }
    merge_interval(unsorted_interval_vec, sorted_interval_vec);
    unsorted_interval_vec.clear();
    
    if(!output_bed_file.empty())
    {
        output_interval_to_bed(sorted_interval_vec, output_bed_file);
    }
    
    SamHeader original_header = my_bam_reader.GetHeader();
    BamWriter my_bam_writer;
    if(!my_bam_writer.Open(output_bam_file, original_header, reference_data))
    {
        cerr << "[ERROR] Fail to open output bam file: " << output_bam_file << endl;
        exit(1);
    }
    
    BamAlignment my_bam_alignment;
    while(my_bam_reader.GetNextAlignmentCore(my_bam_alignment))
    {
        if(!my_bam_alignment.IsMapped())
            continue;
        string bam_alignment_chrom = reference_data[my_bam_alignment.RefID].RefName;
        pair<int, int> bam_alignment_range = make_pair(my_bam_alignment.Position, my_bam_alignment.GetEndPosition() - 1); //GetEndPosition is half-opened

        //if(DEBUG_FLAG)
        //    cerr << "[DEBUG] searching overlap for BAM entry: " << bam_alignment_chrom << " " << bam_alignment_range.first << " " << bam_alignment_range.second << endl;
        
        if(sorted_interval_vec.find(bam_alignment_chrom) == sorted_interval_vec.end())
            continue;
        if(!binary_search_overlap(bam_alignment_range, sorted_interval_vec[bam_alignment_chrom]))
            continue;
        my_bam_alignment.BuildCharData();
        my_bam_writer.SaveAlignment(my_bam_alignment);
    }

    my_bam_reader.Close();
    my_bam_writer.Close();
    cout << "[INFO] Intersect BAM file with buffered VCF complete" << endl;
}



int main(int argc, const char * argv[]) {

    parseOption(argc, argv);
    reduce_bam(opt_input_bam_file, opt_input_variant_file, opt_output_bam_file, opt_output_bed_file, opt_buffer_size);
    create_index(opt_output_bam_file);
    return 0;
}
