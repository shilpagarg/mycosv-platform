// reannotate_elements.cpp - re-classify element_class of existing MycoSV calls
// with the fixed classify_repeat_element (structural-first, Pezizomycotina
// Starship gate, filamentous-only RIP gate, adaptive HGT threshold).
// Reads VARSEQ from a calls.vcf; emits query-coord element labels.
// Build: g++ -O2 -std=c++17 -pthread -I<repo_root> reannotate_elements.cpp -o reannotate_elements
// Usage: ./reannotate_elements <calls.vcf> <cladeGc> <phylum> [class]

#include "layer1_clade_graph.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
using namespace tol;
int main(int argc, char** argv){
    // args: <vcf> <cladeGc> <phylum>
    std::ifstream in(argv[1]);
    double cg = std::stod(argv[2]);
    std::string phylum = argv[3];
    std::string cls = (argc>4)? argv[4] : ".";
    std::string line;
    std::regex vs("VARSEQ=([ACGTNacgtn]+)");
    std::regex rend("END=([0-9]+)"), sv("SVTYPE=([^;]+)");
    while(std::getline(in,line)){
        if(line.empty()||line[0]=='#') continue;
        std::smatch m;
        if(!std::regex_search(line,m,vs)) continue;
        std::string seq=m[1];
        ElementClass ec=classify_repeat_element(seq,cg,phylum,cls);
        // query coordinates = VCF CHROM (field 0) + POS (field 1), END from INFO
        std::istringstream ss(line); std::string chrom,pos,rest;
        std::getline(ss,chrom,'\t'); std::getline(ss,pos,'\t');
        std::smatch me,ms;
        std::string end=std::regex_search(line,me,rend)?me[1].str():pos;
        std::string svt=std::regex_search(line,ms,sv)?ms[1].str():".";
        std::cout<<chrom<<"\t"<<pos<<"\t"<<end<<"\t"<<svt<<"\t"<<element_class_name(ec)<<"\t"<<seq.size()<<"\n";
    }
    return 0;
}
