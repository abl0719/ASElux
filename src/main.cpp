//
//  main.cpp
//  ASElux
//
//  Created by MiaoZong on 12/16/16.
//  Copyright © 2016 UCLA. All rights reserved.
//


#include <iostream>
#include "Parameters.hpp"
#include "StaticIndex.hpp"
#include "DynamicIndex.hpp"
#include "Alignment.hpp"


int main(int argc,  char** argv){
    Parameters p(argc, argv);
    if (!p.check_parameters()){
        return 1;
    }
    if (p.build_index()){
        cout << "Building static index..." << endl;
        StaticIndex * static_index = new StaticIndex(p.get_par("--gtf"), p.get_par("--ref"));
        static_index->save_to_disk(p.get_par("--out"));
    }else if (p.run_alignment()){
        cout << "loading index: " << p.get_par("--index") << endl;
        StaticIndex * static_index = new StaticIndex(p.get_par("--index"));
        int read_len = stoi(p.get_par("--readLen"));
        cout << "Building dynamic index... " << endl;
        DynamicIndex * dynamic_index = new DynamicIndex(p.get_par("--vcf"), static_index, read_len, p.get_par("--mis"), &p);
        int thread_n = stoi(p.get_par("--nthread"));
        Alignment* align = new Alignment(dynamic_index, thread_n, &p);
        cout << "Start alignment... " << endl;
        auto read_files = p.get_read_files();
        if (p.is_pe()){
            for (int i = 0; i*2 < read_files.size(); ++i){
                align->align_paired(read_files[2*i], read_files[2*i+1], p.is_fq());
            }
        }else{
            for (int i = 0; i < read_files.size(); ++i){
                align->align_single(read_files[i], p.is_fq());
            }
        }
        align->output_result(p.get_par("--out"));
        delete align;
    }
    return 0;
}

//#include "Annotation.hpp"
//int main(){
////    StaticIndex * static_index = new StaticIndex("/Users/zmiao/code/ASElux/data/shrinked_annotation.gtf", "/Users/zmiao/code/ASElux/data/shrinked_genome.fa");
////    static_index->save_to_disk("/Users/zmiao/code/ASElux/data/junc_static");
//    StaticIndex * static_index = new StaticIndex("/Users/zmiao/code/ASElux/data/junc_static");
//    DynamicIndex * dynamic_index = new DynamicIndex("/Users/zmiao/code/ASElux/data/new_vcf.vcf", static_index, 50);
////    Annotation test("/Users/zmiao/code/ASElux/data/shrinked_annotation.gtf");
////    test.save_to_disk("/Users/zmiao/code/ASElux/data/new_anno_test.anno");
////    Annotation test2("/Users/zmiao/code/ASElux/data/new_anno_test.anno", 0);
////    test2.output_anno("/Users/zmiao/code/ASElux/data/test2.anno.txt");
////    test2.output_exon("/Users/zmiao/code/ASElux/data/test2.exon.txt");
//}
