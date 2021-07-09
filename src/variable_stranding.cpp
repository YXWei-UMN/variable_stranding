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

    for(auto n:all_files_){
        cout<<n<<endl;
        result_file.open(n,ios::in);
        if (result_file.fail()) {
            cerr << "fail to open blastfile:" << n << "!\n";
        }
        string line;
        while(getline(result_file,line)){
            if (line.size()<=1 || line[0]== '#')
                continue;

            string delimiter = "\t";
            //record primerID and how many collisions this primer has
            string primerID =line.substr(0, line.find(delimiter));
            //delete primer#
            line.erase(0, line.find(delimiter) + delimiter.length());
            //record payloadID to indicate corresponding collision on nt sequences
            string payloadID = line.substr(0, line.find(delimiter));
            long strand_ID = stol(payloadID.substr(7));

            total_collision_num_++;
/*

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
                p.collided_file_.emplace(strand_ID/100);
                primers_.emplace(primerID,p);
            } else{
                primers_.find(primerID)->second.collided_file_.emplace(strand_ID/100);
            }


            // one chunk 4KB  one strand 40 bytes. One chunk = 100 strand
            if (chunks_.find(strand_ID/100)==chunks_.end()){
                chunk f;
                f.collided_primer_.emplace(primerID);
                chunks_.emplace(strand_ID/100,f);
            } else{
                chunks_.find(strand_ID/100)->second.collided_primer_.emplace(primerID);
            }
        }
        result_file.close();
    }
}

void variable_stranding::collisions_among_chunks() {
    cout<<"chunk number: "<<chunks_.size()<<endl;
    ofstream myfile;

    /*int all_collide_chunks = 0;
    for(auto n:primers_){
        all_collide_chunks+=n.second.collided_file_.size();
    }
    cout<<"average collide chunk per primer: "<<all_collide_chunks/(primers_.size()*1.0)<<endl;*/

    // collision number of each chunks
    vector<int> collision_distribution_among_chunks(28002,0);
    for(auto n:chunks_){
        collision_distribution_among_chunks[n.second.collided_primer_.size()]++;
    }



    myfile.open ("collide_primer_per_chunk_video.csv",ios::out | ios::trunc);
    for(int i=0; i < collision_distribution_among_chunks.size(); i++){
        // write into file
        myfile<<i<<","<<collision_distribution_among_chunks[i]<<endl;
    }
    myfile.close();




    // for each chunk/file, go over others to see whether they have common collision
    /*vector<int> common_collision_chunk_degree(chunks_.size(),0);
    int i=0;
    for(auto n:chunks_){
        i++;
        cout<<i<<endl;
        unordered_set<int> commons;
        for(auto m:n.second.collided_primer_){
            for(auto f:primers_[m].collided_file_){
                if (commons.find(f)!=commons.end()) {
                    cout<<"find duplicate one "<<f<<endl;
                    cout<<commons.size()<<endl;
                    commons.emplace(f);
                    cout<<commons.size()<<endl;
                }
                commons.emplace(f);
            }
        }
        common_collision_chunk_degree[commons.size()]++;
        cout<<commons.size()<<" "<<common_collision_chunk_degree[commons.size()]<<endl;
    }
    cout<<"Start to write file"<<endl;
    myfile.open ("common_collision_chunk_degree_video1.3G_200strand.csv",ios::out | ios::trunc);
    for(int i=0; i < common_collision_chunk_degree.size(); i++){
        // write into file
        myfile<<i<<","<<common_collision_chunk_degree[i]<<endl;
    }
    myfile.close();*/
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


    vector<int> common_collision_primer_degree(28010,0);
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
    myfile.close();
}




void variable_stranding::compare_nostrand_with_fixed200() {
    cout<<primers_.size()<<endl;

    unordered_set<int> original_strand;
    unordered_set<int> double_strand;
    unordered_set<int> triple_strand;
    unordered_set<int> overall_strand;
    unordered_set<string> primer_upperbound;

    for(auto n:collisions_){
        if(n.first/200 != n.second/200) continue;
        int id = n.first/200;
        original_strand.emplace(id);
    }

    //fstream result_file("./blast_test_map_A2T_evalue50_200strand",ios::in);
    fstream result_file("./blast_random_C_200strand",ios::in);
    //fstream result_file("./blast_map_A2T_evalue50_200strand",ios::in);
    //fstream result_file("./blast_swap_1_200strand",ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    string line;
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
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
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (original_strand.find(id)!=original_strand.end()) double_strand.emplace(id);
        }
    }
    result_file.close();


    //result_file.open("./blast_test_map_C2G_evalue50_200strand",ios::in);
    result_file.open("./blast_random_2_C_200strand",ios::in);
    //result_file.open("./blast_map_C2G_evalue50_200strand",ios::in);
    //result_file.open("./blast_swap_3_200strand",ios::in);
    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
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
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (double_strand.find(id)!=double_strand.end()) triple_strand.emplace(id);
        }
    }
    result_file.close();

    //result_file.open("./blast_test_map_C2G_A2T_evalue50_200strand",ios::in);
    result_file.open("./blast_random_3_C_200strand",ios::in);
    //result_file.open("./blast_map_C2G_A2T_evalue50_200strand",ios::in);
    //result_file.open("./blast_swap_5_200strand",ios::in);

    if (result_file.fail()) {
        cerr << "fail to open blastfile:" << "blast_result_evalue50_no_strand" << "!\n";
    }
    while(getline(result_file,line)) {
        if (line.size() <= 1 || line[0] == '#')
            continue;

        string delimiter = "\t";
        //record primerID and how many collisions this primer has
        string primerID = line.substr(0, line.find(delimiter));
        //delete primer#
        line.erase(0, line.find(delimiter) + delimiter.length());
        //record payloadID
        string payloadID = line.substr(0, line.find(delimiter));
        long strand_ID = stol(payloadID.substr(7));
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
        long start = stol(start_pos) > stol(end_pos) ? stol(end_pos) : stol(start_pos);
        long end = stol(start_pos) < stol(end_pos) ? stol(end_pos) : stol(start_pos);

        start += strand_ID * g_strand_len;
        end += strand_ID * g_strand_len;

        // even if 10k long strand, still divide 200 to see whether crossed by 200-fixed cutoff point
        int end_strand = end / 200;
        int start_strand = start / 200;
        if (end_strand != start_strand) continue;
        else {
            int id = start_strand;
            if (triple_strand.find(id)!=triple_strand.end()) {
                overall_strand.emplace(id);
                primer_upperbound.emplace(primerID);
            }
        }
    }
    result_file.close();

    cout<<"original strand num: "<<original_strand.size()<<endl;
    cout<<"double strand num: "<<double_strand.size()<<endl;
    cout<<"triple strand num: "<<triple_strand.size()<<endl;
    cout<<"overall strand num: "<<overall_strand.size()<<"  "<<overall_strand.size()/(3293065*1.0)<<endl;
    cout<<"primer upperbound: "<<primer_upperbound.size()<<endl;
}


void variable_stranding::collision_analysis() {
    ofstream myfile;
    // collision distance distribution
    sort(collisions_.begin(), collisions_.end());

    int last = 0;
    int portion_cut=0;
    int portion_trans=0;
    int cut_primers=0;
    int trans_primers=0;
    for(auto n:collisions_){
        if ((n.first-last)>450) {
            portion_cut+=(n.first-last);
            cut_primers++;
        }
        else {
            portion_trans+=(n.first-last);
            trans_primers++;
        }
        last = n.first;
    }
    cout<<portion_cut<<" "<<portion_trans<<endl;
    cout<<cut_primers<<" "<<trans_primers<<endl;
    //collision length distribution
    /*vector<int> collsion_length_distribution(22,0);
    for(auto n:collisions_){
        int len = n.second-n.first;
        collsion_length_distribution[len]++;
    }


    myfile.open ("collsion_length_distribution.csv",ios::out | ios::trunc);
    for(int i=1; i < collsion_length_distribution.size(); i++){
        // write into file
        myfile<<i<<","<<collsion_length_distribution[i]<<endl;
    }
    myfile.close();*/




    // distance distribution
    /*sort(collisions_.begin(),collisions_.end());
    unordered_map<int,int> collision_distance_distribution;

    // calculate distance of collisions inside one primer if it has more than one collision
    int next_start,last_end;
    auto n = collisions_.begin();

    int close=0;
    int far=0;
    while(n!=collisions_.end()){
        last_end=n->second;
        n++;
        next_start=n->first;
        int dis = next_start-last_end;
        dis=dis/100*100;

        if (collision_distance_distribution.find(dis)!=collision_distance_distribution.end()){
            collision_distance_distribution.find(dis)->second++;
        } else{
            collision_distance_distribution.emplace(dis,1);
        }
    }

    myfile.open ("collision_distance.csv",ios::out | ios::trunc);
    for(auto n:collision_distance_distribution){
        // write into file
        myfile<<n.first<<","<<n.second<<endl;
    }
    myfile.close();*/


    // collision distribution on nt sequences
    /*myfile.open ("collision_distribution.csv",ios::out | ios::trunc);
    for(int i=0; i<nts_.size()/1024000; i++){
        // every 200nt print collision distribution
        int collision_nt=0;
        for (int j = 0; j < 1024000; ++j) {
            if (nts_[j+i*1024000]>0) collision_nt++;
        }
        myfile<<i<<"    "<<collision_nt<<endl;
    }
    myfile.close();*/

}




/*
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
       */
/*if (largest==cuts[4])
            strand_point+=220;*//*

    }
    cout<<"strand num:"<<strand_num_<<endl;
    cout<<"len-"<<g_strand_len_1<<":"<<g_strand_len_1_num<<endl;
    cout<<"len-"<<g_strand_len_2<<":"<<g_strand_len_2_num<<endl;
    cout<<"len-"<<g_strand_len_3<<":"<<g_strand_len_3_num<<endl;
    cout<<"len-"<<g_strand_len_4<<":"<<g_strand_len_4_num<<endl;
    cout<<"reduced collision:"<<reduced_collision<<" "<<100*(reduced_collision/(total_collision_num_*1.0))<<"%"<<endl;
}*/
