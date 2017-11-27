//
//  StaticIndex.hpp
//  ASElux
//
//  Created by MiaoZong on 1/3/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#ifndef StaticIndex_hpp
#define StaticIndex_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <set>
#include "SuffixArray.hpp"
#include "Annotation.hpp"
#include "BinaryFile.h"
using namespace std;

class StaticIndex{
public:
    StaticIndex(string gtf_file, string genome_file);
    StaticIndex(string file_prefix);
    ~StaticIndex();
    void save_to_disk(string file_prefix);
    void output_txt(string file_name);
    bool is_exonic(string chr, int posi)
    {return m_annotation->is_exonic(chr, posi);}
    
    size_t get_gene_num()
    {return m_gene_index.size();}
    
    char* get_gene_str(int i)
    {return m_gene_index[i]->m_string;}
    
    int masked_search(const char* read, unsigned int read_l, vector<unsigned int>* posi, const set<uint16_t>& genes);
    
    bool align_to_gene(const char* read, unsigned int read_l, set<uint16_t>& genes);
    
    void global_align(const char *read, unsigned int read_l, set<uint16_t>& forward_genes, set<uint16_t>& backward_genes);
    
    bool locate_in_gene(const unsigned int& locus, const set<uint16_t>& genes){
        if (genes.find(m_annotation->get_id(locus)) != genes.end())
            return true;
        else
            return false;
    }
    
    bool is_in_gene(const char* read, const int& len, uint16_t gene_id);
    
    Annotation* m_annotation;
protected:
    vector<SuffixArray*> m_gene_index;
    SuffixArray* m_genome_index;
    const int m_iteration_n = 2;
};

#endif /* StaticIndex_hpp */
