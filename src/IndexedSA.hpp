//
//  IndexedSA.hpp
//  ASElux
//
//  Created by MiaoZong on 1/6/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#ifndef IndexedSA_hpp
#define IndexedSA_hpp

#include <stdio.h>
#include "SuffixArray.hpp"
#include <unordered_map>
#include "Funtions.hpp"


class IndexedSA : public SuffixArray{
public:
    IndexedSA(const string& input, int length);
    IndexedSA(const char* input, const int& n, const int length);
    ~IndexedSA();
    int quick_search(const char *read, unsigned int read_l, vector<unsigned int> *posi);
    void set_mis(string mis){m_mis_allowed = stoi(mis) + 1;};
    
protected:
    int m_len;
    int m_mis_allowed;
    unordered_map<string, unsigned int *> m_index;
    
    unsigned int * binary_search_region(const string& seed);
    int binary_search(const char *read, unsigned int read_l, vector<unsigned int> *posi, int l, int r);
    
    bool count_mismatch(const char* read, const unsigned int& posi, const int& read_l);
};

#endif /* IndexedSA_hpp */
