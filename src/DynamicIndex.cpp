//
//  DynamicIndex.cpp
//  ASElux
//
//  Created by MiaoZong on 1/9/17.
//  Copyright © 2017 UCLA. All rights reserved.
//

#include "DynamicIndex.hpp"
#include "ThreadPool.h"

DynamicIndex::~DynamicIndex(){
    delete m_static_index;
    delete m_ASE_reads;
    for (auto target:m_hetro_SNPs)
        delete target.second;
    for (auto target:m_homo_SNPs)
        delete target.second;
    for (auto target:m_INDELs)
        delete target.second;
    unordered_map<vector<int>*, int> garbage;
    for (auto target:m_same_read){
        garbage.insert({target, 1});
    }
    for (auto target:garbage)
        delete target.first;
}

DynamicIndex::DynamicIndex(string vcf_file, StaticIndex* static_index, int len, string mis_allowed):
m_static_index(static_index), m_read_len(len){
    // load vcf file
    ifstream in_file(vcf_file);
    string line;
    if (in_file.is_open()){
        while ( getline (in_file,line)){
            if (line[0] == '#'){continue;}
            int geno;
            SNP* snp = new SNP(line, geno);
            if (!m_static_index->is_exonic(snp->m_chr, snp->m_posi))
                continue;
            if (geno == 1){
                add_SNP(snp, m_hetro_SNPs);
            }else if (geno == 2){
                add_SNP(snp, m_homo_SNPs);
            }else if (geno == 3){
                add_SNP(snp, m_INDELs);
            }
        }
        in_file.close();
    }
    // clean heterozygous SNPs near indel
    for (auto snps:m_hetro_SNPs){
        if (m_INDELs.find(snps.first) == m_INDELs.end())
            continue;
        auto hetro = snps.second->begin();
        auto indel = m_INDELs[snps.first]->begin();
        while (hetro != snps.second->end() &&
               indel != m_INDELs[snps.first]->end()) {
            if (abs(indel->m_posi - hetro->m_posi) < m_read_len){
                hetro = snps.second->erase(hetro);
            }else if (indel->m_posi < hetro->m_posi){
                indel ++;
            }else{
                hetro ++;
            }
        }
    }
    // construct ASE_reads
    if (m_hetro_SNPs.size() == 0){
	    cout << "Can not find any SNPs." << endl;
	    exit(1);
    }
    string ASE_reads;
    vector<SNP>::iterator hetro = m_hetro_SNPs.begin()->second->begin();
    vector<SNP>::iterator homo = m_homo_SNPs.begin()->second->begin();
    vector<SNP>::iterator it_start;
    vector<SNP>::iterator it_end;
    int read_i = 0;
    //new added
    for (auto gene:m_static_index->m_annotation->m_gene_list){
        if (m_hetro_SNPs.find(gene.m_chr) == m_hetro_SNPs.end())
            continue;
        //re-locate iterators in new chromosome
        if (hetro->m_chr != gene.m_chr){
            hetro = m_hetro_SNPs[gene.m_chr]->begin();
            homo = m_homo_SNPs[gene.m_chr]->begin();
        }
        //move to next gene if no SNPs
        if (hetro->m_posi > gene.m_end)
            continue;
        //get gene str
        char* gene_cstr = m_static_index->get_gene_str(gene.m_id);
        string gene_str(gene_cstr, gene.m_end-gene.m_start+1);
        for(;homo !=  m_homo_SNPs[gene.m_chr]->end() &&
            homo->m_posi <= gene.m_end; homo ++){
            if (homo->m_posi >= gene.m_start)
                gene_str[homo->m_posi - gene.m_start] = homo->m_alt;
        }
        //get the hetro SNPs withing the gene
        vector<SNP> SNPs_in_gene;
        for (;hetro != m_hetro_SNPs[gene.m_chr]->end() &&
             hetro->m_posi <= gene.m_end; hetro ++){
            SNPs_in_gene.push_back(*hetro);
        }
        //record ASE regions from all transcripts
        unordered_map<string, vector<pair<int, int>>> ASE_regions;
        for (auto&& transcript:gene.m_transcripts){
            for (int i = 0; i < SNPs_in_gene.size();++i){
//                cout << "SNP at " << SNPs_in_gene[i].m_posi;
                pair<string, vector<pair<int, int>>> x;
                if (locate_ASE_region(SNPs_in_gene[i], transcript, x)){
                    ASE_regions.insert(x);
//                    cout << " found in this transcript\n";
                }else{
//                    cout << " not found in this transcript\n";
                }
            }
        }
        //generate ASE reads for each SNP
        for (auto&& x:ASE_regions){
            pair<int, int> blank_region = x.second.back();
            x.second.pop_back();
            int ASE_l = blank_region.first;
            int snp_i = 0;
            //temp data for the next step
            string ASE_read(ASE_l, ' ');
            int snp_n = 0;
            vector<pair<int, SNP>> insert_snps;
//            cout << "start counting SNPs\n";
            for (auto&& exon:x.second){
                while (snp_i < SNPs_in_gene.size() && SNPs_in_gene[snp_i].m_posi < exon.first) {
                    ++snp_i;
                }
                while (snp_i < SNPs_in_gene.size() && SNPs_in_gene[snp_i].m_posi <= exon.second){
                    ++snp_n;
                    auto&& temp_snp = SNPs_in_gene[snp_i];
                    insert_snps.push_back({(temp_snp.m_posi-exon.first+ASE_l), temp_snp});
                    ++snp_i;
                }
                ASE_read.append(string(gene_str.c_str() + (exon.first-gene.m_start), exon.second - exon.first + 1));
                ASE_l += exon.second - exon.first + 1;
            }
            ASE_read.append(string(blank_region.second,' '));
//            assert(ASE_l == ASE_read.size());
            assert(ASE_read.size() == 2*m_read_len-1);
            assert(snp_n == insert_snps.size());
//            cout << ASE_read << endl;
            if (snp_n > 5)
                continue;
            vector<string> haps;
            all_hap(snp_n, haps);
            vector<int>* forward = new vector<int>;
            vector<int>* backward = new vector<int>;
            for (auto hap:haps){
                auto snp_i = hap.begin();
                // for each hetrozygous SNP
                for (auto&& hetro_snp:insert_snps){
                    if (*snp_i == '1')
                        ASE_read[hetro_snp.first] = hetro_snp.second.m_alt;
                    else if (*snp_i == '0')
                        ASE_read[hetro_snp.first] = hetro_snp.second.m_ref;
                    if (hetro_snp.first == m_read_len - 1)
                        m_SNP_index.push_back(Genotype(hetro_snp.second, *snp_i, gene.m_id));
                }
                ASE_reads.append(ASE_read);
                ASE_reads.append(" ");
                forward->push_back(read_i ++);
                m_same_read.push_back(forward);
                ASE_reads.append(reverse_read(ASE_read));
                ASE_reads.append(" ");
                backward->push_back(read_i ++);
                m_same_read.push_back(backward);
            }
        }
    }
    m_ASE_reads = new IndexedSA(ASE_reads, 20);
    m_ASE_reads->set_mis(mis_allowed);
}

int DynamicIndex::get_snp_n(const vector<SNP>::iterator& it,
                            vector<SNP>::iterator& it_start,
                            vector<SNP>::iterator& it_end,
                            const vector<SNP>* snps){
    for (it_start = it; it_start != snps->begin() &&
         it->m_posi - it_start->m_posi < m_read_len; it_start --){}
    if (it->m_posi - it_start->m_posi >= m_read_len)
        it_start ++;
    for (it_end = it; it_end != snps->end() &&
         it_end->m_posi - it->m_posi < m_read_len; it_end ++){}
    return (int) (it_end - it_start);
}

void DynamicIndex::all_hap(int n, vector<string> &hap){
    if (n <= 1){
        hap = {"0","1"};
    }else{
        all_hap(n-1, hap);
        vector<string> temp = hap;
        for (int i = 0; i < hap.size(); ++i){
            temp[i].append("1");
            hap[i].append("0");
        }
        hap.insert(hap.end(), temp.begin(), temp.end());
    }
}

DynamicIndex::SNP::SNP(const string& input, int& genotype){
    //genotype = 0 : homozygous ref
    //genotype = 1 : heterozygous
    //genotype = 2 : homozygous alt
    //genotype = 3 : indel
    vector<string> data;
    split_str(input, '\t', data);
    m_chr = data[0];
    m_posi = stoi(data[1]);
    m_id = data[2];
    m_ref = data[3][0];
    m_alt = data[4][0];
    if (data[3].length() > 1 || data[4].length() > 1){
        genotype = 3;
        return;
    }
    vector<string> info;
    vector<string> info_content;
    split_str(data[8], ':', info);
    split_str(data[9], ':', info_content);
    string geno;
    while (!info.empty()) {
        if (info.back() == "GT"){
            geno = info_content.back();
            break;
        }
        info.pop_back();
        info_content.pop_back();
    }
//    cout << geno[0] << "/" << geno[2] << endl ;
    genotype = geno[0] - '0';
    genotype += geno[2] - '0';
}


void DynamicIndex::add_SNP(const DynamicIndex::SNP* snp, unordered_map<string, vector<SNP>*>& snp_map){
    auto search = snp_map.find(snp->m_chr);
    if (search == snp_map.end()){
        vector<SNP>* temp = new vector<SNP>;
        temp->push_back(*snp);
        snp_map.insert({snp->m_chr, temp});
    }else{
        for (auto it = (search->second)->rbegin(); it != (search->second)->rend(); ++it){
            if (it->m_posi < snp->m_posi){
                (search->second)->insert(it.base(), *snp);
                break;
            }else if (it->m_posi > snp->m_posi){
                if (it+1 == search->second->rend()){
                    search->second->insert(search->second->begin(), *snp);
                    break;
                }
                continue;
            }
        }
    }
}

void DynamicIndex::output_txt(string file_name){
    ofstream out_file;
    out_file.open(file_name);
    out_file << "read length: " << m_read_len << endl;
    auto same_reads = m_same_read.begin();
    char* read = m_ASE_reads->m_string;
    for (auto geno:m_SNP_index){
        out_file << geno.m_id << ':' << geno.m_geno << '\t' << geno.m_gene_id << endl;
        for (int id:*(*same_reads))
            out_file << id << '\t';
        out_file << endl;
        same_reads ++;
        out_file << string(read, m_read_len*2) << endl;
        read += m_read_len*2;
        
        for (int id:*(*same_reads))
            out_file << id << '\t';
        out_file << endl;
        same_reads ++;
        out_file << string(read, m_read_len*2) << endl;
        read += m_read_len*2;
    }
    out_file.close();
}

int DynamicIndex::locus_info(const unsigned int &posi, uint16_t& gene_id, int& dist, int& read_i){
    read_i = posi/(2*m_read_len);
    dist = posi%(2*m_read_len);
    gene_id = m_SNP_index[read_i/2].m_gene_id;
    return read_i%2;
}

bool DynamicIndex::align_paired(const string &read1, const string &read2, unordered_map<string, int *> &result, mutex& m_lock){
    bool debug = false;
//    return;
    //treat read1 as the main read
    vector<unsigned int>* m_main_posi = new vector<unsigned int>;
    int aligned_l = m_ASE_reads->quick_search(read1.c_str(), m_read_len, m_main_posi);
    if (aligned_l == m_read_len){
//        cout << "main read found!\n";
		//~ cout << "r1: " << read1 << endl;
        set<uint16_t> f_gene;
        set<uint16_t> b_gene;
        set<const char*> same_reads;
//        set<char*> b_reads;
        set<uint16_t> all_gene;
        set<unsigned int> unique_posi(m_main_posi->begin(), m_main_posi->end());
        
        //get all reads with different genotypes
        for (unsigned int posi:unique_posi) {
            int dist; int read_i;uint16_t gene_id;
            int direction = locus_info(posi, gene_id, dist, read_i);
            for (auto target:*m_same_read[read_i]){
				if (read_i == target){
                    same_reads.insert(read1.c_str());
				}else{
					same_reads.insert(m_ASE_reads->m_string + dist + target*m_read_len*2);
				}
            }
            if (direction == 0){
                f_gene.insert(gene_id);
            }else{
                b_gene.insert(gene_id);
            }
        }
        //align all reads against the genome
        m_static_index->global_align(read1.c_str(), m_read_len, f_gene, b_gene);
        for (auto f_read:same_reads){
            m_static_index->global_align(f_read, m_read_len, f_gene, b_gene);
        }
        // align the paired read
        if (!f_gene.empty()){
            for (auto gene:f_gene) {
                char* read2_rev = reverse_read(read2.c_str(), m_read_len);
//                cout << "forward1:" << gene << endl;
                if (m_static_index->is_in_gene(read2_rev, m_read_len, gene))
                    all_gene.insert(gene);
                delete read2_rev;
            }
        }
        if (!b_gene.empty()){
            for (auto gene:b_gene) {
//                cout << "reverse1:" << gene << endl;
                if (m_static_index->is_in_gene(read2.c_str(), m_read_len, gene))
                    all_gene.insert(gene);
            }
        }
//        cout << "gene read 1:" << all_gene.size() << endl;
        // output result
        if (all_gene.size() ==  1) {
            Genotype snp;
            set<string> done_snps;
            for (auto posi:unique_posi){
                //cout << posi << "\t";
                snp_info(posi, snp);
//                cout << "aligned to gene: " << *all_gene.begin() << endl;
                if (snp.m_gene_id != *all_gene.begin())
                    continue;
                if (done_snps.find(snp.m_id) != done_snps.end())
                    continue;
                done_snps.insert(snp.m_id);
                m_lock.lock();
                //~ cout << posi << "\t";
                //~ cout << "r1: " << snp.m_id << endl;
                debug = true;
                auto target = result.find(snp.m_id);
                if (target == result.end()){
                    int* temp = new int[2];
                    temp[snp.m_geno] = 1;
                    temp[1-snp.m_geno] = 0;
                    result.insert({snp.m_id, temp});
                }else{
                    result[snp.m_id][snp.m_geno] ++;
                }
                m_lock.unlock();
            }
        }
    }
    m_main_posi->clear();
    // treat read2 as the main read
//    cout << "start alignment2\n";
    aligned_l = m_ASE_reads->quick_search(read2.c_str(), m_read_len, m_main_posi);
    if (aligned_l == m_read_len){
        set<uint16_t> f_gene;
        set<uint16_t> b_gene;
        set<string> snp_id;
        set<const char*> same_reads;
//        set<char*> b_reads;
        set<uint16_t> all_gene;
        set<unsigned int> unique_posi(m_main_posi->begin(), m_main_posi->end());
        //get all reads with different genotypes
        for (unsigned int posi:unique_posi) {
            //            cout << "posi: " << posi << endl;
            int dist; int read_i;uint16_t gene_id;
            int direction = locus_info(posi, gene_id, dist, read_i);
            for (auto target:*m_same_read[read_i]){
				if (read_i == target){
                    same_reads.insert(read2.c_str());
				}else{
					same_reads.insert(m_ASE_reads->m_string + dist + target*m_read_len*2);
				}
                //~ same_reads.insert(m_ASE_reads->m_string
                                  //~ + dist + target*m_read_len*2);
            }
            if (direction == 0){
                f_gene.insert(gene_id);
            }else{
                b_gene.insert(gene_id);
            }
        }
        //align all reads against the genome
        m_static_index->global_align(read2.c_str(), m_read_len, f_gene, b_gene);
        for (auto f_read:same_reads){
            m_static_index->global_align(f_read, m_read_len, f_gene, b_gene);
        }
        // align the paired read
        if (!f_gene.empty()){
            for (auto gene:f_gene) {
//                cout << "forward2:" << gene << endl;
                char* read1_rev = reverse_read(read1.c_str(), m_read_len);
                if (m_static_index->is_in_gene(read1_rev, m_read_len, gene))
                    all_gene.insert(gene);
                delete read1_rev;
            }
        }
        if (!b_gene.empty()){
            for (auto gene:b_gene) {
//                cout << "reverse2:" << gene << endl;
                if (m_static_index->is_in_gene(read1.c_str(), m_read_len, gene))
                    all_gene.insert(gene);
            }
        }
//        cout << "gene read 2:" << all_gene.size() << endl;
        if (all_gene.size() ==  1) {
            Genotype snp;
            set<string> done_snps;
            for (auto posi:unique_posi){
                //cout << posi << "\t";
                snp_info(posi, snp);
//		cout << "aligned to gene: " << *all_gene.begin() << endl;
                if (snp.m_gene_id != *all_gene.begin())
                    continue;
                if (done_snps.find(snp.m_id) != done_snps.end())
                        continue;
                done_snps.insert(snp.m_id);
                m_lock.lock();
                debug = true;
                //~ cout << posi << "\t";
                //~ cout << "r2: " << snp.m_id << endl;
                auto target = result.find(snp.m_id);
                if (target == result.end()){
                    int* temp = new int[2];
                    temp[snp.m_geno] = 1;
                    temp[1-snp.m_geno] = 0;
                    result.insert({snp.m_id, temp});
                }else{
                    result[snp.m_id][snp.m_geno] ++;
                }
                m_lock.unlock();
            }
        }
    }
    delete m_main_posi;
    return debug;
//    cout << "alignment ends\n";
}

void DynamicIndex::align_single(const string &read, unordered_map<string, int *> &result, mutex& m_lock){
    vector<unsigned int>* m_main_posi = new vector<unsigned int>;
    int aligned_l = m_ASE_reads->quick_search(read.c_str(), m_read_len, m_main_posi);
    if (aligned_l == m_read_len){
        set<uint16_t> f_gene;
        set<uint16_t> b_gene;
        set<const char*> same_reads;
        set<uint16_t> all_gene;
        set<unsigned int> unique_posi(m_main_posi->begin(), m_main_posi->end());
        //get all reads with different genotypes
        for (unsigned int posi:unique_posi) {
            int dist; int read_i;uint16_t gene_id;
            int direction = locus_info(posi, gene_id, dist, read_i);
            for (auto target:*m_same_read[read_i]){
                if (read_i == target){
                    same_reads.insert(read.c_str());
                }else{
                    same_reads.insert(m_ASE_reads->m_string + dist + target*m_read_len*2);
                }
            }
            if (direction == 0){
                f_gene.insert(gene_id);
            }else{
                b_gene.insert(gene_id);
            }
        }
        //align all reads against the genome
        m_static_index->global_align(read.c_str(), m_read_len, f_gene, b_gene);
        for (auto f_read:same_reads){
            m_static_index->global_align(f_read, m_read_len, f_gene, b_gene);
        }
        // align the paired read
        for (auto gene:f_gene) {
            all_gene.insert(gene);
        }
        for (auto gene:b_gene) {
            all_gene.insert(gene);
        }
        // output result
        if (all_gene.size() ==  1) {
            Genotype snp;
            set<string> done_snps;
            for (auto posi:unique_posi){
                snp_info(posi, snp);
                if (snp.m_gene_id != *all_gene.begin())
                    continue;
                if (done_snps.find(snp.m_id) != done_snps.end())
                    continue;
                done_snps.insert(snp.m_id);
                m_lock.lock();
                auto target = result.find(snp.m_id);
                if (target == result.end()){
                    int* temp = new int[2];
                    temp[snp.m_geno] = 1;
                    temp[1-snp.m_geno] = 0;
                    result.insert({snp.m_id, temp});
                }else{
                    result[snp.m_id][snp.m_geno] ++;
                }
                m_lock.unlock();
            }
        }
    }
    delete m_main_posi;
}

bool DynamicIndex::locate_ASE_region(const SNP& snp, Annotation::transcript trans, pair<string, vector<pair<int, int>>>& ASE_regions){
    //identify the ASE exon
    auto exon = trans.m_exons.begin();
    for (;exon != trans.m_exons.end() && exon->second <= snp.m_posi; exon ++){}
    if (exon == trans.m_exons.end() || snp.m_posi < exon->first){return false;}
    //get the regions around the SNP
    int start(snp.m_posi), end(snp.m_posi);
    pair<int, int> blank_region = {0,0};
    ASE_regions.first = "";
    ASE_regions.second.clear();
    auto temp_exon = exon;
    for(int i = m_read_len-1; i > 0;){
        if (start - i >= temp_exon->first){
            start -= i;
            i = 0;
        }else if (temp_exon != trans.m_exons.begin()){
            i -= (start - temp_exon->first + 1);
            temp_exon --;
            start = temp_exon->second;
        }else{
            blank_region.first = temp_exon->first - (start - i);
            i = 0;
            start = temp_exon->first;
        }
    }
    auto start_exon = temp_exon;
    temp_exon = exon;
    for(int i = m_read_len-1; i > 0;){
        if (end + i <= temp_exon->second){
            end += i;
            i = 0;
        }else if (temp_exon+1 != trans.m_exons.end()){
            i -= (temp_exon->second - end + 1);
            temp_exon ++;
            end = temp_exon->first;
        }else{
            blank_region.second = end + i - temp_exon->second;
            i = 0;
            end = temp_exon->second;
        }
    }
    auto end_exon = temp_exon;
    for (auto exon = start_exon; exon != end_exon + 1; exon ++){
        int my_start(start);
        int my_end(end);
        if (my_start < exon->first)
            my_start = exon->first;
        if (my_end > exon->second)
            my_end = exon->second;
        string id = to_string(my_start) + "-" + to_string(my_end) + ";";
        ASE_regions.first.append(id);
        ASE_regions.second.push_back({my_start, my_end});
    }
    ASE_regions.second.push_back(blank_region);
    //not tested yet
    return true;
}

