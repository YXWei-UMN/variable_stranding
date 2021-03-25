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
        total_collision_num_++;
        long start = stol(start_pos)>stol(end_pos)?stol(end_pos):stol(start_pos);
        long end = stol(start_pos)<stol(end_pos)?stol(end_pos):stol(start_pos);

        // code for detect overlong collision but not screen out
        if (end-start>=11) {
            primers_.find(primerID)->second.first= true;
        }

        /* code for overlong screen out, not used currently
         * // the following code is used to tell whether the cut can reduce the collision
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

        }*/
        //TODO for primers which has collision longer than 20, all its collision shouldn't be in the nts_
        //TODO for primer which has too many collisions, all its collision shouldn't be in the nts_


        for (long i = start; i <= end; ++i) {
            nts_[i-1]++;
        }
    }
}

void variable_stranding::primer_analysis() {
    // check the portion of primers that has collision longer than 20
    cout<<"collided primers:"<<primers_.size()<<endl;
    cout<<"collision number:"<<total_collision_num_<<endl;
    cout<<"average collision/primer:"<<total_collision_num_/(primers_.size()*1.0)<<endl;

    vector<int> primer_distribution;
    vector<int> primer_distribution_overlong;
    long overlong_collision=0;
    for(auto i:primers_){
        // primers with different collision number
        int x_axis = i.second.second;
        /*if (if_count_intra_redundant_collision){
            x_axis+=i.second->redundant_collision;
        }*/
        if(primer_distribution.size()<x_axis+1){
            for (int j = primer_distribution.size(); j <= x_axis+1; ++j) {
                primer_distribution.push_back(0);
                primer_distribution_overlong.push_back(0);
            }
        }
        primer_distribution[x_axis]++;

        if (i.second.first){
            primer_distribution_overlong[x_axis]++;
            overlong_collision++;
        }
    }

    ofstream myfile;
    double portion_primer=0;
    double portion_collision=0;
    myfile.open ("primer_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < primer_distribution.size(); i++){
        // Cumulative distribution function (cdf) of primer
        if (primer_distribution[i]==0) continue;

        portion_primer+=(primer_distribution[i]/(primers_.size()*1.0));

        // Cumulative distribution function (cdf) of collision
        portion_collision+=(primer_distribution[i]*i/(total_collision_num_*1.0));


        // write into file
        myfile<<portion_primer<<","<<portion_collision<<","<<primer_distribution_overlong[i]<<","<<primer_distribution[i]-primer_distribution_overlong[i]<<endl;

    }
    myfile.close();
    cout<<"primers with overlong collision:"<<overlong_collision<<endl;
    cout<<"portion of good primer:"<<(primers_.size()-overlong_collision)/(primers_.size()*1.0)<<endl;
}

void variable_stranding::strand_analysis() {
    cout<<"strand number:"<<nts_.size()/200<<endl;
    vector<long int> collision_number_distribution;
    vector<long int> collision_length_distribution(22,0);
    vector<long int> collision_coverage_distribution(200,0);

    int nt=0;
    int strand=0;
    int collided_strand=0;
    long int total_length=0;
    long int total_collision=0;
    int total_coverage=0;
    while (strand < nts_.size()/200){

        int collision_number=0;
        int total_collision_length_inside_strand=0;
        int collision_coverage=0;
        for (int i = 0; i < 200; ++i) {
            //every 200 nt is a strand, current nt is nt+i
            // length
            total_collision_length_inside_strand+=nts_[nt+i];
            // coverage
            if (nts_[nt+1]!=0) collision_coverage++;
            // number = the increase on every element in nt array
            if (nt==0 && i==0) continue;
            if(nts_[nt+i]>nts_[nt+i-1]) collision_number+=nts_[nt+i]-nts_[nt+i-1];
        }
        total_length+= total_collision_length_inside_strand;
        total_collision+=collision_number;
        total_coverage+=collision_coverage;

        int average_length = 0;
        if (collision_number!=0) {
            collided_strand++;
            average_length = total_collision_length_inside_strand/collision_number;
        }


        // it's possible average length is 0 (free strand)
        collision_length_distribution[average_length]++;
        collision_coverage_distribution[collision_coverage]++;
        // padding number_distribution vector
        if(collision_number_distribution.size()<collision_number+1){
            for (int j = collision_number_distribution.size(); j <= collision_number+1; ++j) {
                collision_number_distribution.push_back(0);
            }
        }
        collision_number_distribution[collision_number]++;

        // adjust strand and nt for next iteration
        strand++; // finish one strand
        nt+=200; // finish 200 nts
    }

    ofstream myfile;
    myfile.open ("strand_collision_number_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collision_number_distribution.size(); i++){
        if (collision_number_distribution[i]==0) continue;
        // write into file
        myfile<<i<<","<<collision_number_distribution[i]<<endl;
    }
    myfile.close();

    myfile.open ("strand_collision_length_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collision_length_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collision_length_distribution[i]<<endl;
    }
    myfile.close();

    myfile.open ("strand_collision_coverage_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collision_coverage_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collision_coverage_distribution[i]<<endl;
    }
    myfile.close();

    cout<<"collision num in strand_analysis function:"<<total_collision<<" in real data:"<<total_collision_num_<<endl;

    cout<<"overall average collision length (total length/collision number):"<<total_length/(total_collision*1.0)<<endl;
    cout<<"overall average collision coverage (total coverage/total collided strand):"<<(total_coverage)/(collided_strand*1.0)<<endl;
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
         << 100 * (reduced_collision / (total_collision_num_ * 1.0)) << "%" << endl;
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
    cout<<"reduced collision:"<<reduced_collision<<" "<<100*(reduced_collision/(total_collision_num_*1.0))<<"%"<<endl;
}