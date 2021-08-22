//
// Created by eason on 3/8/21.
//

#include "../include/variable_stranding.h"
#include <sstream>
#include <set>

bool isDir(string dir)
{
    struct stat fileInfo;
    stat(dir.c_str(), &fileInfo);
    if (S_ISDIR(fileInfo.st_mode)) {
        return true;
    } else {
        return false;
    }
}

void variable_stranding::listFiles(string baseDir, bool recursive)
{
    DIR *dp;
    struct dirent *dirp;
    if ((dp = opendir(baseDir.c_str())) == NULL) {
        cout << "[ERROR: " << errno << " ] Couldn't open " << baseDir << "." << endl;
        return;
    } else {
        while ((dirp = readdir(dp)) != NULL) {
            if (dirp->d_name != string(".") && dirp->d_name != string("..")) {
                if (isDir(baseDir + dirp->d_name) == true && recursive == true) {
                    //all_files_.push_back(baseDir + dirp->d_name);
                    listFiles(baseDir + dirp->d_name + "/", true);
                } else {
                    all_files_.push_back(baseDir + dirp->d_name);
                }
            }
        }
        closedir(dp);
    }
}


variable_stranding::variable_stranding(string blast_result_path) {
    listFiles(blast_result_path, true);
    fstream result_file;

    int x=0;
    for(auto n:all_files_){
        cout<<n<<endl;
        result_file.open(n,ios::in);
        if (result_file.fail()) {
            cerr << "fail to open blastfile:" << n << "!\n";
        }
        string line;
        bool next_primer= true;
        while(getline(result_file,line)){
            if (line.size()<=1 || line[0]== '#'){
                next_primer= true;
                continue;
            }

            if (!next_primer) continue;
            next_primer= false;
            string delimiter = "\t";
            //record primerID and how many collisions this primer has
            string primerID =line.substr(0, line.find(delimiter));
            //delete primer#
/*            line.erase(0, line.find(delimiter) + delimiter.length());
            //record payloadID to indicate corresponding collision on nt sequences
            string payloadID = line.substr(0, line.find(delimiter));
            long strand_ID = stol(payloadID.substr(7));

            total_collision_num_++;


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


        long start = stol(start_pos)>stol(end_pos)?stol(end_pos):stol(start_pos);
        long end = stol(start_pos)<stol(end_pos)?stol(end_pos):stol(start_pos);
        start+=g_strand_len_1*strand_ID;
        end+=g_strand_len_1*strand_ID;
*/


            if (primers_.find(primerID)==primers_.end()){
                primer p;
                //p.collided_file_.emplace(strand_ID/100);
                primers_.emplace(primerID,p);
            } /*else{
                primers_.find(primerID)->second.collided_file_.emplace(strand_ID/100);
            }*/


            // one chunk 4KB  one strand 40 bytes. One chunk = 100 strand
            /*if (chunks_.find(strand_ID/100)==chunks_.end()){
                chunk f;
                f.collided_primer_.emplace(primerID);
                chunks_.emplace(strand_ID/100,f);
            } else{
                chunks_.find(strand_ID/100)->second.collided_primer_.emplace(primerID);
            }*/
        }
        cout<<x++<<"th payload, accumulated collided primers:"<<primers_.size()<<endl;
        result_file.close();
    }

    //different_primers("/home/eason/CLionProjects/variable_stranding/PDF_GF47/blast_PDF_GF47");
    //different_primers("/home/eason/CLionProjects/variable_stranding/trans_3/blast_PDF_trans_3");

    /*ofstream myfile;
    myfile.open (g_blast_result_path+"primer_distribution.csv",ios::out | ios::trunc);
    vector<int> primer_distribution(140,0);
    for(auto n:primers_){
        // every row has 200 primers, totally has 140 rows
        long primer_ID = stol(n.first.substr(6));
        int x = primer_ID/200;
        primer_distribution[x]++;
    }
    for(auto n:primer_distribution){
        myfile<<n<<endl;
    }
    myfile.close();*/
}


void variable_stranding::different_primers(string file){

    fstream result_file;

    result_file.open(file,ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << file << "!\n";
    }
    int different=0;
    string line;
    while(getline(result_file,line)){
        if (line.size()<=1 || line[0]== '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID =line.substr(0, line.find(delimiter));

        if (primers_.find(primerID)==primers_.end()){
            different++;
            primer p;
            //p.collided_file_.emplace(strand_ID/100);
            primers_.emplace(primerID,p);
        }
    }

    result_file.close();
    cout<<different<<endl;
}



void variable_stranding::collisions_among_primer() {
    cout<<"collided primer number: "<<primers_.size()<<endl;
    ofstream myfile;
    // collision number of each chunks
    /*vector<int> collision_distribution_among_primers(chunks_.size(),0);
    for(auto n:primers_){
        collision_distribution_among_primers[n.second.collided_file_.size()]++;
    }



    myfile.open ("collision distribution among primers "+g_blast_Result+".csv",ios::out | ios::trunc);
    for(int i=0; i < collision_distribution_among_primers.size(); i++){
        // write into file
        myfile<<i<<","<<collision_distribution_among_primers[i]<<endl;
    }
    myfile.close();*/


    /*vector<int> common_collision_primer_degree(28010,0);
    int i=0;
    for(auto n:primers_){
        i++;
        if (i%100==0) {
            cout<<i<<endl;
        }
        set<string> common_primers;
        for(auto m:n.second.collided_file_){
            auto p = chunks_.find(m);
            for(auto f:p->second.collided_primer_){
                common_primers.emplace(f);
            }
        }
        common_collision_primer_degree[common_primers.size()]++;
        common_primers.clear();
    }
    cout<<"Start to write file"<<endl;
    myfile.open ("common_collision_primer_degree_linux_200strand.csv",ios::out | ios::trunc);
    for(int i=0; i < common_collision_primer_degree.size(); i++){
        // write into file
        myfile<<i<<","<<common_collision_primer_degree[i]<<endl;
    }
    myfile.close();*/
}





