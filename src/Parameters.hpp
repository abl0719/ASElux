//
//  Parameters.hpp
//  ASElux
//
//  Created by MiaoZong on 1/17/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#ifndef Parameters_hpp
#define Parameters_hpp

#include <map>
#include <vector>
#include <string>
#include <stdio.h>
using namespace std;

class Parameters{
public:
    Parameters(int argc,  char** argv);
    ~Parameters(){};
    
    bool check_parameters();
    string get_par(string key){
        return m_parameter[key];
    }
    bool is_fq(){
        return m_b_paramter["--fq"];
    }
    bool is_pe(){
        return m_b_paramter["--pe"];
    }
    vector<string>& get_read_files(){
        return m_read_file;
    }
    bool build_index(){
        if (m_runMode == "build")
            return true;
        else
            return false;
    }
    bool run_alignment(){
        if (m_runMode == "align")
            return true;
        else
            return false;
    }
private:
    string m_runMode;
    map<string, string> m_parameter;
    map<string, bool> m_b_paramter;
    vector<string> m_read_file;
};

#endif /* Parameters_hpp */
