#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <iostream>
#include <cstring>

#include "protein.h"

using namespace std;

int getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (string(dirp->d_name)!="." and string(dirp->d_name)!=".."){
            files.push_back(string(dirp->d_name));
        }
    }
    closedir(dp);
    return 0;
}

#define KINK myoVI.size/2+1

int main(){

    string traj_directory = "../Struct_data/1/";
    vector<string> files = vector<string>();

    getdir(traj_directory,files);
    
    for (int i=0;i<files.size();i++){
        Protein myoVI((traj_directory + "myoVI1.structure_"+ to_string(i) +".pdb").c_str());
        cout << myoVI.chain[KINK].z << endl;
    }
    

return 0;
}
