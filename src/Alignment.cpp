//
//  Alignment.cpp
//  ASElux
//
//  Created by MiaoZong on 1/11/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#include "Alignment.hpp"

void Alignment::output_result(string file_name){
    ofstream outfile;
    outfile.open(file_name);
    for (auto target:m_result) {
        double ratio = (double) *(target.second)/ ((double) *(target.second+1) + (double) *(target.second));
        outfile << target.first << '\t';
        outfile << *(target.second) << '\t';
        outfile << *(target.second+1) << '\t' << ratio << endl;
    }
    outfile.close();
}

Alignment::Alignment(DynamicIndex* index, int n, Parameters* p)
:m_dynamic_index(index), m_threadn(n), m_p(p){}

void wrapper_single(Alignment* this_class, const string& read){
    this_class->m_dynamic_index->align_single(read, this_class->m_result, this_class->m_mtx);
}

void wrapper_single_vector(Alignment* this_class, const vector<string>* reads){
    for (string read: *reads){
        transform(read.begin(), read.end(), read.begin(), ::toupper);
        this_class->m_dynamic_index->align_single(read, this_class->m_result, this_class->m_mtx);
    }
    delete reads;
}

void wrapper_paired_vector(Alignment* this_class, const vector<string>* read_box1, const vector<string>* read_box2){
    if (read_box1->size() != read_box2->size()){
        cout << "read reading went wrong.." << endl;
        return;
    }
    for (int i = 0; i < read_box1->size(); ++i){
//		cout << i << "th read:\n";
        string r_name = this_class->m_read_names.front();
        this_class->m_read_names.pop();
        string read1 = (*read_box1)[i];
        string read2 = (*read_box2)[i];
        transform(read1.begin(), read1.end(), read1.begin(), ::toupper);
        transform(read2.begin(), read2.end(), read2.begin(), ::toupper);
        if (this_class->m_dynamic_index->align_paired(read1, read2, this_class->m_result, this_class->m_mtx)){
			//~ cout << r_name << '\t' << read1 << '\t' << read2 << endl;
            //cout << r_name << endl;
        }
    }
    delete read_box1;
    delete read_box2;
}

void wrapper_paired(Alignment* this_class, const string& read1,
                    const string& read2){
    this_class->m_dynamic_index->align_paired(read1, read2, this_class->m_result, this_class->m_mtx);
}

vector<string>* Alignment::fetch_reads(BinaryFile& b_file, BinaryFile::Offset& index, string& residual, unsigned long& read_n, unsigned long& line_n, bool is_fq){
    vector<string>* result = new vector<string>;
    char* data;
    size_t reading_size = 50000000;
    if (index + reading_size >= b_file.fileLength())
        reading_size = b_file.fileLength() - index;
    data = new char[reading_size];
    b_file.read(data, reading_size, index);
    int start_posi = 0;
    
    bool r1_done = false;
    if (m_read_names.size() > 1)
        r1_done = true;
    
    for (int i = 0; i < reading_size; ++i){
        if (data[i] != '\n')
            continue;
        residual.append(data+start_posi, i-start_posi);
        int temp = (line_n++)%4;
        if (temp == 1 || (!is_fq && temp == 3)){
            result->push_back(residual);
            ++ read_n;
        }else if (temp == 0 || (!is_fq && temp == 2)){
//            cout << residual << endl;
            if (!r1_done)
                m_read_names.push(residual);
        }
        start_posi = i+1;
        residual.clear();
    }
    if (start_posi < reading_size)
        residual = string(data+start_posi, reading_size-start_posi);
    delete data;
    index += reading_size;
    return result;
}

void Alignment::align_single(string file, bool is_fq){
    BinaryFile b_file;
    if (!b_file.openExisting(file)){
        cout << "can not open " << file << endl;
        return;
    }
//    ifstream read_input;string line;
//    read_input.open(file);
    time_t start_t;
    time(&start_t);
//    auto start_t = clock();
    unsigned long read_n = 0;
    unsigned long line_i = 0;
    BinaryFile::Offset index = 0;
    string residual;
    if (m_threadn > 1){
        ThreadPool m_pool(m_threadn);
        while (index < b_file.fileLength()){
            vector<string>* reads_box = fetch_reads(b_file, index, residual, read_n, line_i, is_fq);
            m_pool.enqueue(wrapper_single_vector, this, reads_box);
//            wrapper_single_vector(this, reads_box);
        }
    }else{
        while (index < b_file.fileLength()){
            vector<string>* reads_box = fetch_reads(b_file, index, residual, read_n, line_i, is_fq);
            wrapper_single_vector(this, reads_box);
        }
    }
    time_t end_t;
    time(&end_t);
    double seconds = difftime(end_t, start_t);
    cout << "Processed " << read_n << " reads" << endl;
    double speed = ((double) read_n)/seconds;
    cout << "average speed is: " << speed << "reads/s" << endl;
}

void Alignment::align_paired(string file_read1, string file_read2, bool is_fq){
    BinaryFile b_file1;
    BinaryFile b_file2;
    if (!b_file1.openExisting(file_read1) || !b_file2.openExisting(file_read2)){
        cout << "can not open reads file." << endl;
        return;
    }
    time_t start_t;
    time(&start_t);
    unsigned long read_n = 0;
    unsigned long line_i = 0;
    BinaryFile::Offset index = 0;
    string residual1;
    string residual2;
    if (m_threadn > 1){
        ThreadPool m_pool(m_threadn);
        while (index < b_file1.fileLength()){
            unsigned long read_n2(read_n);
            unsigned long line_i2(line_i);
            BinaryFile::Offset index2(index);
            vector<string>* reads_box1 = fetch_reads(b_file1, index, residual1, read_n, line_i, is_fq);
            vector<string>* reads_box2 = fetch_reads(b_file2, index2, residual2, read_n2, line_i2, is_fq);
            m_pool.enqueue(wrapper_paired_vector, this, reads_box1, reads_box2);
        }
    }else{
        while (index < b_file1.fileLength()){
            unsigned long read_n2(read_n);
            unsigned long line_i2(line_i);
            BinaryFile::Offset index2(index);
            vector<string>* reads_box1 = fetch_reads(b_file1, index, residual1, read_n, line_i, is_fq);
            vector<string>* reads_box2 = fetch_reads(b_file2, index2, residual2, read_n2, line_i2, is_fq);
            wrapper_paired_vector(this, reads_box1, reads_box2);
        }
    }
    time_t end_t;
    time(&end_t);
    double seconds = difftime(end_t, start_t);
    cout << "Processed " << read_n << " reads" << endl;
    double speed = ((double) read_n)/seconds;
    cout << "average speed is: " << speed << "reads/s" << endl;
}
