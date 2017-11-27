//
//  Annotation.hpp
//  ASElux
//
//  Created by MiaoZong on 12/20/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#ifndef Annotation_hpp
#define Annotation_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include "BinaryFile.h"

using namespace std;
class Annotation{
public:
    Annotation(string gtf_file);
    Annotation(string binary_file, int i);
    ~Annotation();                  // to be done...
    uint16_t get_id(const unsigned int& posi)
    {return m_ref_genome[posi];}
    
    void get_region(const int id, int& start, int& end);
    bool is_exonic(string chr, int posi){
        if (m_exon.find(chr) == m_exon.end()){return false;}
        return (*((m_exon.find(chr))->second))[posi];
    }
    void save_to_disk(string file_name);
    
    
    void output_anno(string file_name);
    void output_exon(string file_name);
    
    class transcript{
    public:
        string m_chr;
        int m_start;
        int m_end;
        vector<pair<int, int>> m_exons;
        
        transcript(string& chr, int& start, int& end):m_chr(chr),m_start(start),m_end(end){}
        void add_exon(int start, int end){
            m_exons.push_back({start, end});
        }
    };
    
    class gene{
    public:
        string m_chr;
        int m_start;
        int m_end;
        uint16_t m_id;
        vector<transcript> m_transcripts;
        
        gene(string& chr, int& start, int& end, uint16_t& id):m_chr(chr),m_start(start),m_end(end),m_id(id){}
        gene(char* chr, int& start, int& end, uint16_t& id):m_chr(chr),m_start(start),m_end(end),m_id(id){}
        void add_transcript(const transcript& data){
            m_transcripts.push_back(data);
        }
    };
    vector<gene> m_gene_list;
    void check_anno();
private:
    vector<uint16_t> m_ref_genome;
    vector<int> m_start;
    vector<int> m_end;
    unordered_map<string, int> m_chr_l;
    unordered_map<string, vector<bool>*> m_exon;
    
    void add_gene(unordered_map<string, vector<int*>*>& annotation,
                  const string& chr, const int& start, const int& end);
    
};

#endif /* Annotation_hpp */
