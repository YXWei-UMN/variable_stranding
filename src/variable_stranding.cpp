//
// Created by eason on 3/8/21.
//

#include "../include/variable_stranding.h"
#include <sstream>
#include <set>

variable_stranding::variable_stranding(string blastfile) {

    nts_=vector<int8_t>(g_total_strand_number,0);


    fstream result_file(blastfile,ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << blastfile << "!\n";
    }
    string line;
    while(getline(result_file,line)){
        if (line.size()<=1 || line[0]== '#')
            continue;

        string delimiter = "\t";
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete payload#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete % identity
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete alignment length
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete mismatches
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete gap opens
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer start
        line.erase(0, line.find(delimiter) + delimiter.length());
        //delete primer end
        line.erase(0, line.find(delimiter) + delimiter.length());

        string start_pos, end_pos;
        start_pos = line.substr(0, line.find(delimiter));
        line.erase(0, line.find(delimiter) + delimiter.length());
        end_pos = line.substr(0, line.find(delimiter));
        collision_num_++;
        long start = stol(start_pos)>stol(end_pos)?stol(end_pos):stol(start_pos);
        long end = stol(start_pos)<stol(end_pos)?stol(end_pos):stol(start_pos);
        for (long i = start; i <= end; ++i) {
            nts_[i-1]++;
        }
    }
}


void variable_stranding::greedy() {
    int reduced_collision=0;
    int strand_point=0;
    while (strand_point<g_total_nt_number-220){
        int cuts[5];
        cuts[0]=(nts_[strand_point+180]);
        cuts[1]=(nts_[strand_point+190]);
        cuts[2]=(nts_[strand_point+200]);
        cuts[3]=(nts_[strand_point+210]);
        cuts[4]=(nts_[strand_point+220]);
        int largest = *max_element(cuts,cuts+5);
        reduced_collision += largest;
        if (largest==cuts[0])
            strand_point+=180;
        if (largest==cuts[1])
            strand_point+=190;
        if (largest==cuts[2])
            strand_point+=200;
        if (largest==cuts[3])
            strand_point+=210;
        if (largest==cuts[4])
            strand_point+=220;
    }
    cout<<"reduced_collision"<<reduced_collision<<" "<<(reduced_collision/(collision_num_*1.0))<<endl;
}