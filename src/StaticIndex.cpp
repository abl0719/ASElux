//
//  StaticIndex.cpp
//  ASElux
//
//  Created by MiaoZong on 1/3/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#include "StaticIndex.hpp"
#include <unordered_map>
#include <assert.h>
#include <fstream>
#include "Funtions.hpp"
#include <cstring>

StaticIndex::StaticIndex(string gtf_file, string genome_file){
    //create genome annotation
    m_annotation = new Annotation(gtf_file);
    cout << "annotation created...\n";
    //load genome sequence
    ifstream in_file(genome_file);
    string line;
    unordered_map<string, string> genome_temp;
    if (in_file.is_open()){
        string chr;
        while ( getline (in_file,line)){
            if (line.substr(0,1) == ">"){
                chr = line.substr(1,(line.length()-1));
                chr = chr.substr(0, chr.find(' '));
                continue;
            }
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            auto sequence = genome_temp.find(chr);
            if (sequence == genome_temp.end()){
                genome_temp.insert({chr, line});
            }else{
                (sequence->second).append(line);
            }
        }
        in_file.close();
    }
    //create suffix array and new genome
    string new_genome;
    for(auto gene_temp:(m_annotation->m_gene_list)){
        auto seq = genome_temp.find(gene_temp.m_chr);
	if (seq == genome_temp.end())
	    continue;
        string seq_temp = seq->second.substr(gene_temp.m_start-1, gene_temp.m_end-gene_temp.m_start+1);
        new_genome.append(seq_temp);
        new_genome.append(" ");
        SuffixArray* temp_sa = new SuffixArray(seq_temp);
        m_gene_index.push_back(temp_sa);
    }
    m_genome_index = new SuffixArray(new_genome);
}

StaticIndex::~StaticIndex(){
    delete m_annotation;
    delete m_genome_index;
    for (auto sa:m_gene_index){
        delete sa;
    }
}

void append_bfile_by_sa(SuffixArray* input, BinaryFile& bfile,
                          BinaryFile::Offset& index){
    unsigned int length = input->m_length;
    bfile.write(length, index);
    index += sizeof(unsigned int);
    bfile.write(input->m_string, length, index);
    index += length;
    bfile.write_array(input->m_sarray, length, index);
    index += length*sizeof(unsigned int);
}

SuffixArray* load_from_file(BinaryFile& bfile, BinaryFile::Offset& index){
    unsigned int length;
    bfile.read(length, index);
    index += sizeof(unsigned int);
    char* sequence = new char[length];
    bfile.read(sequence, length, index);
    index += length;
    unsigned int* sa = new unsigned int[length];
    bfile.read_array(sa, length, index);
    index += length*sizeof(unsigned int);
    SuffixArray* result = new SuffixArray(sequence, length, sa);
    delete[] sequence;
    delete[] sa;
    return result;
}

void StaticIndex::save_to_disk(string file_prefix){
//    cout << "saving to disk..." << endl;
    // save genome annotation to one file
    string annotation_file = file_prefix + ".annotation";
    m_annotation->save_to_disk(annotation_file);
//    cout << "annotaiton finished..." << endl;
    // save local sa to another file
    BinaryFile bfile;
    string sa_local_file = file_prefix + "_gene.sa";
    if (! bfile.createNew(sa_local_file)){
        cout << "can not create" << sa_local_file << endl;
        return;
    }
    BinaryFile::Offset index = 0;
    for (auto it = m_gene_index.begin(); it != m_gene_index.end(); it ++){
        append_bfile_by_sa(*it, bfile, index);
    }
    bfile.close();
//    cout << "genes finished..." << endl;
    // save genome sa to the last file
    string sa_genome_file = file_prefix + "_genome.sa";
    if (! bfile.createNew(sa_genome_file)){
        cout << "can not create" << sa_genome_file << endl;
        return;
    }
    index = 0;
    append_bfile_by_sa(m_genome_index, bfile, index);
    bfile.close();
//    cout << "genome finished..." << endl;
}

StaticIndex::StaticIndex(string file_prefix){
    //load annotaiton from file
    string annotation_file = file_prefix + ".annotation";
    m_annotation = new Annotation(annotation_file, 0);
    //load genome sa from file
    BinaryFile bfile;
    string sa_genome_file = file_prefix + "_genome.sa";
    if (! bfile.openExisting(sa_genome_file)){
        cout << "can not open file " << sa_genome_file << endl;
        return;
    }
    BinaryFile::Offset index(0);
    m_genome_index = load_from_file(bfile, index);
    bfile.close();
    //load gene sa from file
    string sa_local_file = file_prefix + "_gene.sa";
    if (! bfile.openExisting(sa_local_file)){
        cout << "can not open file " << sa_local_file << endl;
        return;
    }
    index = 0;
    while(index < bfile.fileLength()){
        SuffixArray* temp = load_from_file(bfile, index);
        m_gene_index.push_back(temp);
    }
    bfile.close();
}

void StaticIndex::output_txt(string file_name){
    ofstream out_file;
    out_file.open(file_name);
    for (auto sa:m_gene_index){
        unsigned int n = sa->m_length;
        for (unsigned int i = 0; i < n; i ++){
            out_file << *(sa->m_string + i);
        }
        out_file << endl;
        for (unsigned int i = 0; i < n; i ++){
            out_file << *(sa->m_sarray + i) << " ";
        }
        out_file << endl;
    }
    out_file.close();
}

inline int my_strncmp(const char *s1, const char *s2, int& n){
    for (size_t i=n ; i > 0; s1 ++, s2 ++, --i){
        if (*s1 != *s2){
            return ((*(unsigned char *)s1 < *(unsigned char *)s2) ? -1 : +1);
        }
    }
    for (; *s1 != '\0' && *s1 == *s2; s1 ++, s2 ++){
        ++n;
    }
    return 0;
}

int StaticIndex::masked_search(const char *read, unsigned int read_l, vector<unsigned int> *posi, const set<uint16_t> &genes){
    if (!posi->empty())
        return 0;
    long r = m_genome_index->m_length -1;
    long l = 0;
    // if l or r is in the masked genes, push the bundry
    for (;l <= r && locate_in_gene(m_genome_index->m_sarray[l], genes);++l){}
    for (;l <= r && locate_in_gene(m_genome_index->m_sarray[r], genes);--r){}
    
    long mid = (r+l)/2;
    int matched_l = 1;
    unsigned int result[2] = {0, 0};
    while(l <= r && matched_l <= read_l){
        for (;locate_in_gene(m_genome_index->m_sarray[mid], genes);
             ++mid){}
        int c = my_strncmp(m_genome_index->m_string+m_genome_index->m_sarray[mid], read, matched_l);
        if (c == 0){
            result[0] = (unsigned int) mid;
            result[1] = matched_l;
            matched_l ++;
        }else if(c < 0){
            l = mid + 1;
            for (;l <= r && locate_in_gene(m_genome_index->m_sarray[l], genes);
                 ++l){}
            mid = (l+r)/2;
        }else{
            r = mid-1;
            for (;l <= r && locate_in_gene(m_genome_index->m_sarray[r], genes);
                 --r){}
            mid = (l+r)/2;
        }
    }
    if (result[1] < m_genome_index->m_anchor_l){return result[1];}
    posi->push_back(m_genome_index->m_sarray[result[0]]);
    if (result[0] + 1 < m_genome_index->m_length){
        for (int i = 1; result[0] + i < m_genome_index->m_length && strncmp(read, m_genome_index->m_string+m_genome_index->m_sarray[result[0] + i], result[1]) == 0 ;i ++){
            posi->push_back(m_genome_index->m_sarray[result[0] + i]);
        }
    }
    if (result[0] >= 1){
        for (int i = 1; result[0] >= i && strncmp(read, m_genome_index->m_string+m_genome_index->m_sarray[result[0] - i], result[1]) == 0 ;i ++){
            posi->push_back(m_genome_index->m_sarray[result[0] - i]);
        }
    }
    return result[1];
}

bool StaticIndex::align_to_gene(const char *read, unsigned int read_l, set<uint16_t> &genes){
    int aligned_l = 0;
    bool result = false;
//    cout << "start:" << endl;
    while (aligned_l + 20 < read_l){
        vector<unsigned int>* posis = new vector<unsigned int>;
        int mmp_l = masked_search(read+aligned_l, read_l-aligned_l, posis, genes);
//        cout << mmp_l << endl;
        aligned_l += mmp_l;
        if (mmp_l < 20){
            delete posis;
            continue;
        }
        if (posis->size() > 10){
            delete posis;
            return false;
        }
        result = true;
        for (unsigned int locus: *posis){
            uint16_t gene_id = m_annotation->get_id(locus);
//            cout << "global align to gene " << gene_id << endl;
            if (m_gene_index[gene_id]->local_search(read, read_l))
                genes.insert(gene_id);
        }
        delete posis;
    }
    return result;
}

void StaticIndex::global_align(const char *read, unsigned int read_l, set<uint16_t>& forward_genes, set<uint16_t>& backward_genes){
    for (int i = 0; i < m_iteration_n; ++i){
        if (!align_to_gene(read, read_l, forward_genes))
            break;
    }
    char* rev_read = reverse_read(read, read_l);
    for (int i = 0; i < m_iteration_n; ++i){
        if (!align_to_gene(rev_read, read_l, backward_genes))
            break;
    }
    delete rev_read;
}

bool StaticIndex::is_in_gene(const char *read, const int &len, uint16_t gene_id){
//	cout << "local align to gene " << gene_id << endl;
    return m_gene_index[gene_id]->local_search(read, len);
}









