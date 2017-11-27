//
//  IndexedSA.cpp
//  ASElux
//
//  Created by MiaoZong on 1/6/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#include "IndexedSA.hpp"
#include <cstring>

IndexedSA::IndexedSA(const char* input, const int& n, const int length):SuffixArray(input, n), m_len(length){
    for (int i = 0; i + length <= m_length; i ++){
        string seed = string(input+i, length);
        auto point = m_index.find(seed);
        if (point != m_index.end()){continue;}
        unsigned int* current = binary_search_region(seed);
        m_index.insert({seed, current});
    }
}

IndexedSA::IndexedSA(const string& input, int length):SuffixArray(input), m_len(length){
    for (int i = 0; i + length <= m_length; i ++){
        string seed = input.substr(i, length);
        if (seed.find(' ') != string::npos){continue;}
        auto point = m_index.find(seed);
        if (point != m_index.end()){continue;}
        unsigned int* current = binary_search_region(seed);
//        cout << *(current) << " : " << *(current + 1) << endl;
        m_index.insert({seed, current});
    }
    m_mis_allowed = 3;
}

IndexedSA::~IndexedSA(){
    for (auto x:m_index)
        delete x.second;
}

int IndexedSA::quick_search(const char *read, unsigned int read_l, vector<unsigned int> *posi){
    if (read_l < m_len){
        fputs("The read is too short.", stderr);
        abort();
//        return SuffixArray::binary_search(read, read_l, posi);
    }
    int mis(0), l(0);
    vector<unsigned int> temp_posi;
    set<unsigned int> candidates;
    for (int i=0; i < read_l; ++i){
	if (read[i] == 'N')
	    return 0;
    }
    while (mis < 3 && l + m_anchor_l < read_l) {
//        cout << "start: " << mis << ' ' << l << endl;
        string seed = string(read+l, m_len);
        auto point = m_index.find(seed);
        if (point == m_index.end()){
            ++ mis;
            l += m_len;
        }else{
            int temp_l = binary_search(read+l, read_l-l, &temp_posi, *(point->second), *(point->second + 1));
            if (temp_l == read_l){
                posi->swap(temp_posi);
                return read_l;
            }else if (temp_l >= m_anchor_l && temp_posi.size() < 6){
                for (unsigned int x: temp_posi)
                    candidates.insert(x);
            }
            l += (temp_l+1);
            ++ mis;
        }
//        cout << "end: " << mis << ' ' << l << endl;
    }
    for (unsigned int posi_i: candidates){
//        cout << "check mismatch\n";
        if (count_mismatch(read, posi_i, read_l))
            posi->push_back(posi_i);
    }
//    cout << "quichk_search ends\n";
    if (posi->empty())
        return 0;
    else
        return read_l;
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

int IndexedSA::binary_search(const char *read, unsigned int read_l, vector<unsigned int> *posi, int l, int r){
    if (!posi->empty())
        return 0;
//    cout << l << " : " << r << endl;
    int mid = (r+l)/2;
    int matched_l = 1;
    unsigned int result[2] = {0, 0};
    //    cout << "stage 0" << endl;
    while(l <= r && matched_l <= read_l){
        int c = my_strncmp(m_string+m_sarray[mid], read, matched_l);
        if (c == 0){
            result[0] = mid;
            result[1] = matched_l;
            matched_l ++;
        }else if(c < 0){
            l = mid + 1;
            mid = (l+r)/2;
        }else{
            r = mid-1;
            mid = (l+r)/2;
        }
    }
    //    cout << "stage 1" << endl;
    if (result[1] < m_anchor_l){return result[1];}
    posi->push_back(m_sarray[result[0]]);
    if (result[0] + 1 < m_length){
        for (int i = 1; result[0] + i < m_length && strncmp(read, m_string+m_sarray[result[0] + i], result[1]) == 0 ;++ i){
            posi->push_back(m_sarray[result[0] + i]);
        }
    }
    if (result[0] >= 1){
        for (int i = 1; result[0] >= i && strncmp(read, m_string+m_sarray[result[0] - i], result[1]) == 0 ;++ i){
            posi->push_back(m_sarray[result[0] - i]);
        }
    }
    //    cout << "stage 2" << endl;
    return result[1];
}

unsigned int* IndexedSA::binary_search_region(const string& seed){
    unsigned int* result = new unsigned int[2];
    int r = m_length -1;
    int l = 0;
    int mid = (r+l)/2;
    int matched_l = 1;
    //    cout << "stage 0" << endl;
    while(matched_l <= m_len){
        int c = strncmp(m_string+m_sarray[mid], seed.c_str(), matched_l);
        if (c == 0){
            matched_l ++;
        }else if(c < 0){
            l = mid + 1;
            mid = (l+r)/2;
        }else{
            r = mid-1;
            mid = (l+r)/2;
        }
    }
    r = mid;
    l = mid;
    -- matched_l;
//    cout << matched_l << endl;
    while (l >= 0 && strncmp(m_string+m_sarray[l], seed.c_str(), matched_l) == 0){
        *result = l--;
//        cout << *result << " : ";
    }
    while (r < m_length && strncmp(m_string+m_sarray[r], seed.c_str(), matched_l) == 0){
        *(result+1) = r++;
//        cout << *(result+1) << endl;
    }
    return result;
}

bool IndexedSA::count_mismatch(const char* read, const unsigned int& posi, const int& read_l){
    if ((posi/read_l)%2 == 1)
        return false;
    int dist = read_l - 1 - (posi%read_l);
    if (m_string[posi + dist] != read[dist])
        return false;
    int mis(0);
    for (int i = 0; i < read_l && mis < m_mis_allowed; ++i){
        if (m_string[posi + i] != read[i]){
            ++ mis;
            if (m_string[posi + i] == ' '){
                return false;
            }
        }
    }
    return mis < m_mis_allowed;
}

