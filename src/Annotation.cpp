//
//  Annotation.cpp
//  ASElux
//
//  Created by MiaoZong on 12/20/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#include "Annotation.hpp"
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <assert.h>
#include "Funtions.hpp"
#include <cstring>

//void split_str(const std::string &s, char delim,
//               std::vector<std::string> &elems) {
//    std::stringstream ss;
//    ss.str(s);
//    std::string item;
//    while (std::getline(ss, item, delim)) {
//        elems.push_back(item);
//    }
//}

void Annotation::add_gene(unordered_map<string, vector<int*>*>& annotation,
              const string& chr, const int& start, const int& end){
    // add gene to the gene list
    int* region = new int[2] {start, end};
//    cout << region[0] << " " << region[1] << endl;
    auto search = annotation.find(chr);
    if (search == annotation.end()){
        vector<int*>* temp = new vector<int*>;
        temp->push_back(region);
        annotation.insert({chr, temp});
    }else{
        for (vector<int*>::reverse_iterator it = (search->second)->rbegin(); it != (search->second)->rend(); ++it){
            if ((*it)[1] < region[0]){
//                cout << "insert" << endl;
                (search->second)->insert(it.base(), region);
                break;
            }else if ((*it)[0] > region[1]){
//                cout << "pass" << endl;
                continue;
            }else{
//                cout << "merge" << endl;
                region[0] = min(region[0], (*it)[0]);
                region[1] = max(region[1], (*it)[1]);
                it = vector<int*>::reverse_iterator((search->second)->erase(next(it).base()));
                if (search->second->empty() || it == search->second->rbegin()){
                    search->second->push_back(region);
                    break;
                }
            }
        }
    }
//    cout << chr << ": " << end << endl;
    // get the end of each chromosome
    auto search_2 = m_chr_l.find(chr);
    if (search_2 == m_chr_l.end()){
        m_chr_l.insert({chr, end});
    }else{
        search_2->second = max(search_2->second, end);
    }
}

Annotation::Annotation(string gtf_file){
    ifstream in_file(gtf_file);
    string line;
    unordered_map<string, vector<int*>*> annotation;
    //  load gtf file
    uint16_t gene_id = 0;
    gene* temp_gene = nullptr;
    transcript* temp_transcript = nullptr;
    if (in_file.is_open()){
        while ( getline (in_file,line)){
	    if (line[0] == '#')
		continue;
            vector<string> data;
            split_str(line, '\t', data);
            int start = stoi(data[3]);
            int end = stoi(data[4]);
            if (data[2] == "gene"){
                //record the gene
                if (temp_gene == nullptr){
                    temp_gene = new gene(data[0], start, end, gene_id);
                    ++ gene_id;
                }else if(data[0] != temp_gene->m_chr || start > temp_gene->m_end){
                    temp_gene->m_transcripts.push_back(*temp_transcript);
                    m_gene_list.push_back(*temp_gene);
                    delete temp_transcript;
                    temp_transcript = nullptr;
                    delete temp_gene;
                    temp_gene = new gene(data[0], start, end, gene_id);
                    ++ gene_id;
                }else{
                    temp_gene->m_end = max(temp_gene->m_end, end);
                }
                //record the chromosome length
                auto search_it = m_chr_l.find(data[0]);
                if (search_it == m_chr_l.end()){
                    m_chr_l.insert({data[0], end});
                }else{
                    search_it->second = max(search_it->second, end);
                }
            }else if (data[2] == "transcript"){
                if (temp_transcript == nullptr){
                    temp_transcript = new transcript(data[0], start, end);
                }else{
                    temp_gene->add_transcript(*temp_transcript);
                    delete temp_transcript;
                    temp_transcript = new transcript(data[0], start, end);
                }
            }else if (data[2] == "exon"){
                temp_transcript->add_exon(start, end);
            }
        }
        in_file.close();
        assert(temp_gene != nullptr);
        assert(temp_transcript != nullptr);
        temp_gene->m_transcripts.push_back(*temp_transcript);
        m_gene_list.push_back(*temp_gene);
        delete temp_transcript;
        delete temp_gene;
    }
    //reverse exons if they are in the - strand
    for (auto&& g:m_gene_list){
        for (auto&& tr:g.m_transcripts){
            if (tr.m_exons.empty())
                continue;
            if (tr.m_exons.front().first > tr.m_exons.back().first){
                reverse(tr.m_exons.begin(), tr.m_exons.end());
            }
        }
    }
//    //test output
//    output_anno(gtf_file + ".out.txt");
//    abort();
//    uint16_t gene_id = 0;
    //  initialize the data members
    gene_id = 0;
    int start = 0;
    int end = 0;
    for (auto&& temp_gene:m_gene_list){
        end = start + temp_gene.m_end - temp_gene.m_start;
        m_start.push_back(start);
        m_end.push_back(end);
        for (int i = 0; i + temp_gene.m_start <= temp_gene.m_end + 1; i ++){
            m_ref_genome.push_back(gene_id);
        }
        gene_id++;
        start = end + 2;
    }
//    for (auto chromosome:annotation){
//        string chr = chromosome.first;
//        vector<int*>* regions = chromosome.second;
//        for(int* region:*regions){
//            m_gene_list.push_back(gene(chr, region[0], region[1], gene_id));
//            end = start + region[1]-region[0];
//            m_start.push_back(start);
//            m_end.push_back(end);
//            for (int i = 0; i + region[0] <= region[1] + 1; i ++){
//                m_ref_genome.push_back(gene_id);
//            }
//            gene_id++;
//            start = end + 2;
//            delete region;
//        }
//        delete regions;
//    }
    // initialize each chromosome
    for (auto chromosome:m_chr_l){
        string chr = chromosome.first;
        int length = chromosome.second;
//        cout << chr << ": " << length << endl;
        vector<bool>* temp = new vector<bool>((length+1), false);
        m_exon.insert({chr, temp});
    }
    // add exonic region mark
    for (auto&& temp_gene:m_gene_list){
        string chr = temp_gene.m_chr;
        auto search = m_exon.find(chr);
        for (auto&& temp_transcript:temp_gene.m_transcripts){
            for (auto&& exon:temp_transcript.m_exons){
                for (int i = exon.first; i <= exon.second ; i ++){
                    (*(search->second))[i] = true;
                }
            }
        }
    }
//    ifstream in_file2(gtf_file);
////    in_file2.open(gtf_file);
//    if (in_file2.is_open()){
//        while ( getline (in_file2,line)){
//            vector<string> data_temp;
//            split_str(line, '\t', data_temp);
//            if (data_temp[2] != "exon"){continue;}
//            string chr = data_temp[0];
//            int start = stoi(data_temp[3]);
//            int end = stoi(data_temp[4]);
////            cout << chr << ": " << start << "-" << end << endl;
//            for (int i = start; i <= end ; i ++){
//                auto search = m_exon.find(chr);
//                (*(search->second))[i] = true;
//            }
//        }
//        in_file2.close();
//    }
    cout << "Anootation finished..." << endl;
    
    //output_anno(gtf_file + ".out.txt");
    //output_exon(gtf_file + ".exon.txt");
}

Annotation::~Annotation(){
    for (auto chr:m_exon){
        delete chr.second;
    }
    //need to delete the newed vectors
}

void Annotation::save_to_disk(string file_name){
    BinaryFile bfile;
    if (!bfile.createNew(file_name)){
        cout << "can not open file " << file_name << endl;
        return;
    }
    // the file contains No. of genes, m_start, m_end,
    //   gene list, region of exons (in this order).
    // write gene No.
    bfile.write((int) m_gene_list.size(), 0);
    // write m_start & m_end
    BinaryFile::Offset index_end = sizeof(int);
    bfile.write_vector(m_start, index_end);
    index_end += sizeof(int)*m_start.size();
    bfile.write_vector(m_end, index_end);
    index_end += sizeof(int)*m_end.size();
    // wirte gene list
    char chr_name_p[20];
    for (auto it = m_gene_list.begin(); it != m_gene_list.end(); it ++){
        strcpy(chr_name_p, it->m_chr.c_str());
        bfile.write(chr_name_p, 20, index_end);
        index_end += 20;
        bfile.write(it->m_start, index_end);
        index_end += sizeof(int);
        bfile.write(it->m_end, index_end);
        index_end += sizeof(int);
        bfile.write(it->m_id, index_end);
        index_end += 2;
        bfile.write((int) it->m_transcripts.size(), index_end);
        index_end += sizeof(int);
        for (auto&& temp_transcript:it->m_transcripts){
            bfile.write(temp_transcript.m_start, index_end);
            index_end += sizeof(int);
            bfile.write(temp_transcript.m_end, index_end);
            index_end += sizeof(int);
            bfile.write((int) temp_transcript.m_exons.size(), index_end);
            index_end += sizeof(int);
            for (auto&& temp_exon:temp_transcript.m_exons){
                bfile.write(temp_exon.first, index_end);
                index_end += sizeof(int);
                bfile.write(temp_exon.second, index_end);
                index_end += sizeof(int);
            }
        }
    }
    // save the space for the end of file
    BinaryFile::Offset file_end_index = index_end;
    index_end += sizeof(BinaryFile::Offset);
    // write region of exons
    for (auto chr_l:m_chr_l){
        int length = chr_l.second;
        string chr_name = chr_l.first;
        // name of chromosome
        strcpy(chr_name_p, chr_name.c_str());
        bfile.write(chr_name_p, 20, index_end);
        index_end += 20;
        // length of chromosome
        bfile.write(length, index_end);
        index_end += sizeof(int);
        // pre set the space for No. of exons
        BinaryFile::Offset index_exon_n = index_end;
        index_end += sizeof(int);
        // exon info
        auto exon_temp = m_exon.find(chr_name);
        bool x = false;
        int exon_n = 0;
        for (int i = 0; i < (exon_temp->second)->size(); i ++){
            if ((*(exon_temp->second))[i] != x){
                //this way, the file record the start and (end+1),
                //    need to minus 1 when reading the file
                bfile.write(i, index_end);
                index_end += sizeof(int);
                x = (*(exon_temp->second))[i];
                exon_n ++;
            }
        }
        if ((exon_temp->second)->back() == true){
            bfile.write((int)((exon_temp->second)->size()-1), index_end);
            index_end += sizeof(int);
            exon_n ++;
        }
        assert(exon_n%2 == 0);
        bfile.write((int) (exon_n/2), index_exon_n);
    }
    //finish
    bfile.write(index_end, file_end_index);
    assert(index_end == bfile.fileLength());
    bfile.close();
}

Annotation::Annotation(string binary_file, int nouse){
    //load annotation from binary file
    BinaryFile bfile;
    if (! bfile.openExisting(binary_file)){
        cout << "can not open file " << binary_file << endl;
        return;
    }
    // read gene No.
    int gene_n;
    bfile.read(gene_n, 0);
    BinaryFile::Offset index = sizeof(int);
    //read m_start
    m_start = vector<int>(gene_n);
    bfile.read_vector(m_start, index);
    //read m_end
    index += gene_n*sizeof(int);
    m_end = vector<int>(gene_n);
    bfile.read_vector(m_end, index);
    index += gene_n*sizeof(int);
    //read gene list
    int gene_left = gene_n;
    while (gene_left > 0){
        --gene_left;
        char chr_name[20];
        int start;
        int end;
        uint16_t id;
        bfile.read(chr_name, 20, index);
        bfile.read(start, index + 20);
        bfile.read(end, index + sizeof(int) + 20);
        bfile.read(id, index + 2*sizeof(int) + 20);
        string chr(chr_name);
        gene temp_gene(chr, start, end, id);
        for (int i = 0; i + start <= end + 1; i ++){
            m_ref_genome.push_back(id);
        }
        index += sizeof(int)*2 + 22;
        int transcript_left;
        bfile.read(transcript_left, index);
        index += sizeof(int);
        while (transcript_left > 0){
            --transcript_left;
            int start;
            bfile.read(start, index);
            index += sizeof(int);
            int end;
            bfile.read(end, index);
            index += sizeof(int);
            transcript temp_transcript(chr, start, end);
            int exon_left;
            bfile.read(exon_left, index);
            index += sizeof(int);
            while (exon_left > 0){
                --exon_left;
                int start;
                bfile.read(start, index);
                index += sizeof(int);
                int end;
                bfile.read(end, index);
                index += sizeof(int);
                temp_transcript.m_exons.push_back({start, end});
            }
            temp_gene.m_transcripts.push_back(temp_transcript);
        }
        m_gene_list.push_back(temp_gene);
    }
//    for (int i = 0; i < gene_n; i ++){
//        char chr_name[20];
//        int start;
//        int end;
//        uint16_t id;
//        bfile.read(chr_name, 20, index);
//        bfile.read(start, index + 20);
//        bfile.read(end, index + sizeof(int) + 20);
//        bfile.read(id, index + 2*sizeof(int) + 20);
//        m_gene_list.push_back(Annotation::gene(chr_name, start, end, id));
//        for (int i = 0; i + start <= end + 1; i ++){
//            m_ref_genome.push_back(id);
//        }
//        index += sizeof(int)*2 + 22;
//    }
    //read exons
    BinaryFile::Offset file_end;
    bfile.read(file_end, index);
    index += sizeof(BinaryFile::Offset);
    while(index < file_end){
        char chr_name[20];
        bfile.read(chr_name, 20, index);
        index += 20;
        int chr_l;
        bfile.read(chr_l, index);
        index += sizeof(int);
        int exon_n;
        bfile.read(exon_n, index);
        index += sizeof(int);
        string chr_name_s(chr_name);
//        cout << chr_name_s << '\t' << chr_l << endl;
        m_chr_l.insert({chr_name_s, chr_l});
        vector<bool>* exon_temp = new vector<bool>(chr_l+1, false);
        int start_temp;
        int end_temp;
        for (int i = 0; i < exon_n; i ++){
            bfile.read(start_temp, index);
            bfile.read(end_temp, index + sizeof(int));
            for (int j = start_temp; j < end_temp; j ++){
                (*exon_temp)[j] = true;
            }
            index += 2*sizeof(int);
        }
        m_exon.insert({chr_name_s, exon_temp});
    }
}

void Annotation::output_anno(string file_name){
    ofstream out_file;
    out_file.open(file_name);
    for (auto g:m_gene_list){
        int id = (int) g.m_id;
        out_file << g.m_chr << ":\t" << g.m_start << "-" << g.m_end << "\t" << id << "\n";
        for (auto tr:g.m_transcripts){
            assert(tr.m_start == tr.m_exons.front().first);
            assert(tr.m_end == tr.m_exons.back().second);
            out_file << "\t" << tr.m_start << "-" << tr.m_end << " :";
            for (auto exon:tr.m_exons){
                out_file << " " << exon.first << '-' << exon.second;
            }
            out_file << endl;
        }
//        out_file << m_start[id] << "\t" << m_end[id] << endl;
//        for (int i = m_start[id]; i <= m_end[id]; i ++){
//            assert(m_ref_genome[i] == g.m_id);
//        }
    }
    out_file.close();
}

void Annotation::output_exon(string file_name){
    ofstream out_file;
    out_file.open(file_name);
    for (auto chr_l:m_chr_l){
        string chr_name = chr_l.first;
        out_file << chr_name << endl;
        // exon info
        auto exon_temp = m_exon.find(chr_name);
        bool x = false;
        int exon_n = 0;
//        cout << chr_name << endl << (exon_temp->second)->size() << endl;
        for (int i = 0; i < (exon_temp->second)->size(); i ++){
            if ((*(exon_temp->second))[i] != x){
                out_file << i << " ";
                //this way, the file record the start and (end+1),
                //    need to minus 1 when reading the file
                x = (*(exon_temp->second))[i];
                exon_n ++;
            }
        }
        if ((exon_temp->second)->back() == true){
            out_file << (exon_temp->second)->size()-1 << " ";
            exon_n ++;
        }
        assert(exon_n%2 == 0);
        out_file << endl;
    }
    out_file.close();
}

void Annotation::check_anno(){
    int last_end(0);
    for (auto g:m_gene_list){
        assert(g.m_end > g.m_start);
        assert(g.m_start > last_end);
        last_end = g.m_end;
    }
}














