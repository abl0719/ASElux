//
//  DynamicIndex.hpp
//  ASElux
//
//  Created by MiaoZong on 1/9/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#ifndef DynamicIndex_hpp
#define DynamicIndex_hpp

#include <stdio.h>
#include "StaticIndex.hpp"
#include "IndexedSA.hpp"
#include "Parameters.hpp"
#include <vector>
#include <sstream>
#include "Funtions.hpp"
#include <set>
#include <mutex>
#include <assert.h>

class DynamicIndex{
public:
    class SNP{
    public:
        string m_chr;
        string m_id;
        int m_posi;
        char m_ref;
        char m_alt;
        SNP(string chr,string name, int posi, char ref, char alt)
        :m_chr(chr), m_posi(posi), m_ref(ref), m_alt(alt), m_id(name){}
        SNP(const string& input, int& genotype);
        SNP(const SNP& snp):m_chr(snp.m_chr), m_posi(snp.m_posi), m_ref(snp.m_ref), m_alt(snp.m_alt), m_id(snp.m_id){}
    };
    class Genotype{
    public:
        string m_id;
        int m_geno;
        uint16_t m_gene_id;
        Genotype(const SNP& in, char geno, uint16_t gene_id)
        :m_geno(geno), m_id(in.m_id), m_gene_id(gene_id){
//            cout << geno << endl;
            m_geno = geno - '0';
        }
        Genotype(){}
    };
    
    DynamicIndex(string vcf_file, StaticIndex* static_index, int len, string mis_allowed, Parameters* p);
    ~DynamicIndex();
    void output_txt(string file_name);
    
    bool align_paired(const string& read1, const string& read2, unordered_map<string, int*>& result, mutex& m_lock);
    void align_single(const string& read, unordered_map<string, int*>& result, mutex& m_lock);
    
private:
    int m_read_len;
    StaticIndex* m_static_index;
    IndexedSA* m_ASE_reads;
    unordered_map<string, vector<SNP>*> m_hetro_SNPs;
    unordered_map<string, vector<SNP>*> m_homo_SNPs;
    unordered_map<string, vector<SNP>*> m_INDELs;
    vector<Genotype> m_SNP_index;
    // The backward read always follows the forward read
    // with the same genotype. So one genotype for two read.
    vector<vector<int>*> m_same_read;
    // The read point to the same vector which contains all
    // related reads (including itself). The forward reads
    // have even id, and the backward reads have odd id.
    
    void add_SNP(const SNP* snp,
                 unordered_map<string, vector<SNP>*>& snp_map);
    void all_hap(int n, vector<string>& hap);
//    string reverse_read(string read);
    int get_snp_n(const vector<SNP>::iterator& it,
                  vector<SNP>::iterator& it_start,
                  vector<SNP>::iterator& it_end,
                  const vector<SNP>* snps);
    int locus_info(const unsigned int& posi, uint16_t& gene_id, int& dist, int& read_i);
    void snp_info(const unsigned int& posi, Genotype& snp){
        int read_i = posi/(2*m_read_len);
        snp = m_SNP_index[read_i/2];
    }
    
    //add for transcript base dynamic index
    bool locate_ASE_region(const SNP& snp, Annotation::transcript trans, pair<string, vector<pair<int, int>>>& ASE_regions);
    Parameters* m_p;
};

#endif /* DynamicIndex_hpp */
