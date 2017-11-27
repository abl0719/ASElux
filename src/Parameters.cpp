//
//  Parameters.cpp
//  ASElux
//
//  Created by MiaoZong on 1/17/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#include "Parameters.hpp"
#include <fstream>
#include <iostream>

inline bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

Parameters::Parameters(int argc,  char** argv){
    m_parameter.insert({"--readLen",""});
    m_parameter.insert({"--index",""});
    m_parameter.insert({"--vcf",""});
    m_parameter.insert({"--gtf",""});
    m_parameter.insert({"--ref",""});
    m_parameter.insert({"--out",""});
    m_parameter.insert({"--seqFiles",""});
    m_parameter.insert({"--nthread","1"});
    m_parameter.insert({"--mis","2"});
    m_b_paramter.insert({"--fq",false});
    m_b_paramter.insert({"--fa",false});
    m_b_paramter.insert({"--se",false});
    m_b_paramter.insert({"--pe",false});
    if (argc < 2)
        return;
    m_runMode = string(argv[1]);
    string par_str;
    for (int i = 2; i < argc; ++i){
        string this_str = string(argv[i]);
        if (m_parameter.find(this_str) != m_parameter.end()){
            par_str = this_str;
        }else if (m_b_paramter.find(this_str) != m_b_paramter.end()){
            m_b_paramter[this_str] = true;
        }else if (par_str == "--seqFiles"){
            m_read_file.push_back(this_str);
        }else{
            m_parameter[par_str] = this_str;
        }
    }
}

bool Parameters::check_parameters(){
    if (build_index()){
        if (!file_exists(m_parameter["--gtf"])){
            cout << "Can not open the annotation file." << endl;
            return false;
        }else if(!file_exists(m_parameter["--ref"])){
            cout << "Can not open the reference genome." << endl;
            return false;
        }else if(m_parameter["--out"].size() == 0){
            cout << "Please indicate the output file." << endl;
            return false;
        }
        return true;
    }else if (run_alignment()){
        if (m_b_paramter["--fq"] == m_b_paramter["--fa"]){
            cout << "Please indicate the read type (fa/fq)." << endl;
            return false;
        }else if (m_b_paramter["--pe"] == m_b_paramter["--se"]){
            cout << "Please indicate the read type (pe/se)." << endl;
            return false;
        }else if(!file_exists(m_parameter["--vcf"])){
            cout << "Can not open the SNP file." << endl;
            return false;
        }else if(!file_exists(m_parameter["--index"]+".annotation")){
            cout << "Can not open the index file." << endl;
            return false;
        }else if(!file_exists(m_parameter["--index"]+"_genome.sa")){
            cout << "Can not open the index file." << endl;
            return false;
        }else if(!file_exists(m_parameter["--index"]+"_gene.sa")){
            cout << "Can not open the index file." << endl;
            return false;
        }else if(m_parameter["--readLen"].size() == 0){
            cout << "Please indicate the read length." << endl;
            return false;
        }else if(m_read_file.size() == 0){
            cout << "Please indicate the read sequences." << endl;
            return false;
        }else if(is_pe() && m_read_file.size()%2 != 0){
            cout << "The read files are not paired end reads." << endl;
            return false;
        }else if(m_parameter["--out"].size() == 0){
            cout << "Please indicate the output file." << endl;
            return false;
        }
        return true;
    }
    cout << "The running mode is unkown." << endl;
    return false;
}
