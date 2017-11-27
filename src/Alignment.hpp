//
//  Alignment.hpp
//  ASElux
//
//  Created by MiaoZong on 1/11/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#ifndef Alignment_hpp
#define Alignment_hpp

#include <stdio.h>
#include "DynamicIndex.hpp"
#include "ThreadPool.h"
#include "BinaryFile.h"
#include <assert.h>

class Alignment{
public:
    Alignment(DynamicIndex* index, int n);
    ~Alignment(){
        delete m_dynamic_index;
        for(auto target:m_result)
            delete target.second;
    };
    
    void align_paired(string file_read1, string file_read2, bool is_fq);
    void align_single(string file, bool is_fq);
    void output_result(string file_name);
    
    vector<string>* fetch_reads(BinaryFile& b_file, BinaryFile::Offset& index, string& residual, unsigned long& read_n, unsigned long& line_n, bool is_fq);
    
//    void fetch_paired_reads(BinaryFile& b_file1, BinaryFile& b_file2, BinaryFile::Offset& index1, BinaryFile::Offset& index2, string& residual1, string& residual2, unsigned long& read_n, unsigned long& line_n, bool is_fq);
    
    DynamicIndex* m_dynamic_index;
    mutex m_mtx;
    unordered_map<string, int*> m_result;
    
    queue<string> m_read_names;
private:
//    DynamicIndex* m_dynamic_index;
    int m_threadn;
    
    
};



#endif /* Alignment_hpp */
