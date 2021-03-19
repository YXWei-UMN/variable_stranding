//
// Created by eason on 3/8/21.
//

#include "../include/variable_stranding.h"
#include <sstream>
#include <set>

variable_stranding::variable_stranding(string blastfile) {

    nts_.resize(g_total_nt_number,0);


    fstream result_file(blastfile,ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << blastfile << "!\n";
    }
    string line;
    while(getline(result_file,line)){
        if (line.size()<=1 || line[0]== '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID =line.substr(0, line.find(delimiter));
        if (primers_.find(primerID)==primers_.end()){
            primers_.emplace(primerID,make_pair(false,1));
        } else{
            primers_.find(primerID)->second.second++;
        }
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

        // the following code is used to tell whether the cut can reduce the collision
        // if the collision length is > 19, no matter how to cut, there is always a collision >= 10
        if (end-start>=19) {
            primers_.find(primerID)->second.first= true;
            continue;
        }
        // if the coolision length is <=10, it's ok to be cut arbitrily
        // thus we only have to process situation with 20> length >10
        else if(end-start>10){
            // the real cut range is (end-10,start+10)
            // since end-start>10, end-10>start && start+10<end
            // it's impossible that end-10>start+10, since end-start>20 violate our primer length
            long tmp = start;
            start = end-9;
            end = tmp+9;

        }


        //TODO for primers which has collision longer than 20, all its collision shouldn't be in the nts_
        //TODO for primer which has too many collisions, all its collision shouldn't be in the nts_
        /*for (long i = start; i <= end; ++i) {
            nts_[i-1]++;
        }*/
    }
}

void variable_stranding::data_analysis() {
    // check the portion of primers that has collision longer than 20
    int primer_with_long_collision_num=0;
    for(auto n:primers_){
        if (n.second.first) primer_with_long_collision_num++;
    }

    cout<<"number of primers that has long collision is: "<<primer_with_long_collision_num<<"    portion:"<<(primer_with_long_collision_num/(primers_.size()*1.0))<<endl;

    // check the collision number of each primer
    vector<int> primer_distribution;
    for(auto i:primers_){
        // primers with different collision number
        int x_axis = i.second.second;
        /*if (if_count_intra_redundant_collision){
            x_axis+=i.second->redundant_collision;
        }*/
        if(primer_distribution.size()<x_axis+1){
            for (int j = primer_distribution.size(); j <= x_axis+1; ++j) {
                primer_distribution.push_back(0);
            }
        }
        primer_distribution[x_axis]++;
    }

    ofstream myfile;
    double portion_primer=0;
    double portion_collision=0;
    myfile.open ("primer_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < primer_distribution.size(); i++){
        // Cumulative distribution function (cdf) of primer
        portion_primer+=(primer_distribution[i]/(primers_.size()*1.0));

        // Cumulative distribution function (cdf) of collision
        portion_collision+=(primer_distribution[i]*i/(collision_num_*1.0));

        // write into file
        myfile<<portion_primer<<","<<portion_collision<<endl;
    }
    myfile.close();


}

void variable_stranding::fixed_length() {
    int reduced_collision=0;
    int strand_point=0;
        while(strand_point<g_total_nt_number-g_strand_len_1) {
            strand_num_++;
            reduced_collision += nts_[strand_point + g_strand_len_1];
            strand_point += g_strand_len_1;
        }
    cout << "strand num:" << strand_num_ << endl;
    cout << "reduced collision:" << reduced_collision << " " << reduced_collision << " "
         << 100 * (reduced_collision / (collision_num_ * 1.0)) << "%" << endl;
}
void variable_stranding::greedy() {
    int reduced_collision=0;
    int strand_point=0;
    long g_strand_len_1_num=0;
    long g_strand_len_2_num=0;
    long g_strand_len_3_num=0;
    long g_strand_len_4_num=0;
    while (strand_point<g_total_nt_number-g_strand_len_1){
        strand_num_++;
        int cuts[4];
        cuts[0]=(nts_[strand_point+g_strand_len_1]);
        cuts[1]=(nts_[strand_point+g_strand_len_2]);
        cuts[2]=(nts_[strand_point+g_strand_len_3]);
        cuts[3]=(nts_[strand_point+g_strand_len_4]);
        //cuts[4]=(nts_[strand_point+220]);
        int largest = *max_element(cuts,cuts+4);
        reduced_collision += largest;
        if (largest==cuts[0]){
            strand_point+=g_strand_len_1;
            g_strand_len_1_num++;
        }
        else if (largest==cuts[1]){
            strand_point+=g_strand_len_2;
            g_strand_len_2_num++;
        }
        else if (largest==cuts[2]){
            strand_point+=g_strand_len_3;
            g_strand_len_3_num++;
        }
        else if (largest==cuts[3]){
            strand_point+=g_strand_len_4;
            g_strand_len_4_num++;
        }
        /*if (largest==cuts[4])
            strand_point+=220;*/
    }
    cout<<"strand num:"<<strand_num_<<endl;
    cout<<"len-"<<g_strand_len_1<<":"<<g_strand_len_1_num<<endl;
    cout<<"len-"<<g_strand_len_2<<":"<<g_strand_len_2_num<<endl;
    cout<<"len-"<<g_strand_len_3<<":"<<g_strand_len_3_num<<endl;
    cout<<"len-"<<g_strand_len_4<<":"<<g_strand_len_4_num<<endl;
    cout<<"reduced collision:"<<reduced_collision<<" "<<100*(reduced_collision/(collision_num_*1.0))<<"%"<<endl;
}