//
//  SuffixArray.hpp
//  ASElux
//
//  Created by MiaoZong on 12/20/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#ifndef SuffixArray_hpp
#define SuffixArray_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "Annotation.hpp"
using namespace std;


class SuffixArray{
public:
    enum{
        LESS,
        EQUAL,
        MORE
    };
    SuffixArray(const string& input);
    SuffixArray(const char* input, const int& n);
    SuffixArray(const char* sequence, unsigned int length, unsigned int* suffix_array);
    ~SuffixArray();
    int binary_search(const char* read, unsigned int read_l, vector<unsigned int>* posi);
//    int debug_binary_search(const char* read, unsigned int read_l, vector<unsigned int>* posi);
    bool local_search(const char* read, unsigned int read_l);
//    char* get_string(){return m_string;}
//    int binary_search(const char* read, unsigned int read_l, vector<unsigned int>* posi,
//                    vector<int>* geneID, Annotation* genome);
    void show_sa();
    void verify_sa();
    
    unsigned int m_length;
    char* m_string;
    unsigned int* m_sarray;
    const int m_anchor_l = 10;
protected:
    
    bool less_than(const char* str, const unsigned int& len, const unsigned int& posi1, const unsigned int& posi2);
};

#endif /* SuffixArray_hpp */
