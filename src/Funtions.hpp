//
//  Funtions.hpp
//  ASElux
//
//  Created by MiaoZong on 12/27/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#ifndef Funtions_hpp
#define Funtions_hpp

#include <stdio.h>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <sstream>

using namespace std;

class Functions{
public:
    Functions(){
        m_alphabet_rev = vector<string>(256);
        unsigned char value = 0;
        for (int i = 0; i <= 3333; i ++){
            string key = to_string(i);
            key.insert(0, "0000");
            key = key.substr(key.length() - 4);
            replace(key.begin(), key.end(), '0', 'A');
            replace(key.begin(), key.end(), '1', 'T');
            replace(key.begin(), key.end(), '2', 'G');
            replace(key.begin(), key.end(), '3', 'C');
            if (has_number(key)) {continue;}
            m_alphabet.insert({key, value});
            m_alphabet_rev[value] = key;
            while (key.back() == 'A'){
                key.pop_back();
                m_alphabet.insert({key, value});
            }
            value ++;
        }
    }
    
    ~Functions(){}
    
    bool has_number(const std::string& s){
        return !s.empty() && find_if(s.begin(), s.end(), [](char c) { return isdigit(c); }) != s.end();
    }
    
    string compress_read(const string& read){
        size_t read_len = read.length();
        string result((read_len+3)/4+1, '\0');
        result.back() = (unsigned char) (4*result.length() - read_len);
        for (int i = 0; i*4 < read_len; i ++){
            auto searcher = m_alphabet.find(read.substr(i*4, 4));
            result[i] = searcher->second;
        }
        return result;
    }
    
    string* compress_read(const string* read){
        size_t read_len = read->length();
        string* result = new string((read_len+3)/4+1, '\0');
        result->back() = (unsigned char) (4*result->length() - read_len);
        for (int i = 0; i*4 < read_len; i ++){
            auto searcher = m_alphabet.find(read->substr(i*4, 4));
            (*result)[i] = searcher->second;
        }
        return result;
    }
    
    string decompress_read(const string& compressed_read){
        size_t new_len = 4*compressed_read.length() - (unsigned char) compressed_read.back();
        string result(4*compressed_read.length(), '\0');
        for (int i = 0; i < (compressed_read.length()-1); i ++){
            result.replace(i*4, i*4+4, m_alphabet_rev[(unsigned char) compressed_read[i]]);
        }
        result.resize(new_len);
        return result;
    }
    
    string* decompress_read(const string* compressed_read){
        size_t new_len = 4*compressed_read->length() - (unsigned char) compressed_read->back();
        string* result = new string(4*compressed_read->length(), '\0');
        for (int i = 0; i < (compressed_read->length()-1); i ++){
            result->replace(i*4, 4, m_alphabet_rev[(unsigned char) compressed_read->at(i)]);
        }
        result->resize(new_len);
        return result;
    }
private:
    unordered_map<string, unsigned char> m_alphabet;
    vector<string> m_alphabet_rev;
};

inline void split_str(const std::string &s, char delim,
               std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

inline string reverse_read(string read){
    string result(read);
    for(int i = 0; i < read.length(); i ++){
        switch (read[read.length()-i-1]) {
            case 'A':
                result[i] = 'T';
                break;
            case 'T':
                result[i] = 'A';
                break;
            case 'G':
                result[i] = 'C';
                break;
            case 'C':
                result[i] = 'G';
                break;
            default:
                break;
        }
    }
    return result;
}

inline char* reverse_read(const char* read, const int& len){
    char* result = new char[len];
    for(int i = 0; i < len; i ++){
        switch (read[len-i-1]) {
            case 'A':
                result[i] = 'T';
                break;
            case 'T':
                result[i] = 'A';
                break;
            case 'G':
                result[i] = 'C';
                break;
            case 'C':
                result[i] = 'G';
                break;
            default:
                break;
        }
    }
    return result;
}


#endif /* Funtions_hpp */
