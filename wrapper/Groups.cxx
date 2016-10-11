/*************************************************
 *	main.cpp
 *
 *	Release: Sep 2016
 *	Update: Sep 2016
 *
 *	University of North Carolina at Chapel Hill
 *	Department of Computer Science
 *
 *	Ilwoo Lyu, ilwoolyu@cs.unc.edu
 *************************************************/

#include <cstdlib>
#include <vector>
#include <string>
#include <dirent.h>
#include "GroupsCLP.h"
#include "GroupwiseRegistration.h"

void getListFile(string path, vector<string> &list, const string &suffix)
{
    DIR *dir = opendir(path.c_str());
    if (dir != NULL)
    {
        while (dirent *entry = readdir(dir))
        {
            string filename = entry->d_name;
            if (filename.size() >= suffix.size() && equal(suffix.rbegin(), suffix.rend(), filename.rbegin()))
            list.push_back(path + "/" + filename);
        }
    }
    closedir(dir);
    sort(list.begin(), list.begin() + list.size());
}

void getTrimmedList(vector<string> &list, const vector<string> &name)
{
    int i = 0;
    while (i < list.size())
    {
        size_t found;
        for (int j = 0; j < name.size(); j++)
        {
            found = list[i].find(name[j].substr(name[j].rfind('/') + 1));
            if (string::npos != found) break;
        }
        if (string::npos == found) list.erase(list.begin() + i);
        else i++;
    }
    sort(list.begin(), list.begin() + list.size());
}

int main(int argc, char *argv[])
{
    PARSE_ARGS;

    if (argc == 1)
    {
        cout << "Usage: " << argv[0] << " --help" << endl;
        return EXIT_SUCCESS;
    }
    
    // update list files from the directory information
    if (!dirSphere.empty() && listSphere.empty()) getListFile(dirSphere, listSphere, "vtk");// listSphere.erase(listSphere.begin() + 30, listSphere.begin() + listSphere.size());
    if (!dirProperty.empty() && listProperty.empty()) getListFile(dirProperty, listProperty, "txt");
    if (!dirSurf.empty() && listSurf.empty()) getListFile(dirSurf, listSurf, "vtk");
    if (!dirLandmark.empty() && listLandmark.empty()) getListFile(dirLandmark, listLandmark, "txt");
    if (!dirCoeff.empty() && listCoeff.empty()) getListFile(dirCoeff, listCoeff, "coeff");
    
    // subject names
    int nSubj = listSphere.size();
    vector<string> subjName;
    if (nSubj > 0)
    {
        for (int i = 0; i < nSubj; i++)
        {
            int pivot = listSphere[i].rfind('/') + 1;
            string name = listSphere[i].substr(pivot, listSphere[i].length() - 4 - pivot);
            pivot = name.find('.');
            if (pivot == string::npos) pivot = name.length();
            subjName.push_back(name.substr(0, pivot));
        }
    }
    //for (int i = 0; i < nSubj; i++) cout << subjName[i] << endl;
    if (listOutput.empty())
    for (int i = 0; i < nSubj; i++) listOutput.push_back(dirOutput + "/" + subjName[i] + ".coeff");
    
    // trim all irrelevant files to the sphere files
    if (!dirProperty.empty()) getTrimmedList(listProperty, subjName);
    if (!dirLandmark.empty()) getTrimmedList(listLandmark, subjName);
    if (!dirCoeff.empty()) getTrimmedList(listCoeff, subjName);
    if (!dirProperty.empty()) getTrimmedList(listProperty, listFilter);
    if (listWeight.empty())
    for (int i = 0; i < listProperty.size() / nSubj; i++)
    listWeight.push_back(1);
    
    int nProperties = listProperty.size();
    int nOutput = listOutput.size();
    int nWeight = listWeight.size();
    int nLandmark = listLandmark.size();
    int nCoeff = listCoeff.size();
    int nSurf = listSurf.size();
    
    const char **property = NULL;
    const char **sphere = NULL;
    const char **output = NULL;
    const char **landmark = NULL;
    const char **coeff = NULL;
    const char **surf = NULL;
    
    if (nProperties > 0) property = new const char*[nProperties];
    if (nSubj > 0) sphere = new const char*[nSubj];
    if (nOutput > 0) output = new const char*[nOutput];
    if (nLandmark > 0) landmark = new const char*[nLandmark];
    if (nCoeff > 0) coeff = new const char*[nCoeff];
    if (nSurf > 0) surf = new const char*[nSurf];
    if (surf == NULL) weightLoc = 0;
    float *weight = new float[nWeight];
    
    // exception handling
    if (nSubj == 0)
    {
        cout << "Fatal error: no subject is provided!" << endl;
        return EXIT_FAILURE;
    }
    else if (nSubj != nOutput)
    {
        cout << "Fatal error: # of subjects is incosistent with # of outputs!" << endl;
        return EXIT_FAILURE;
    }
    else if (nLandmark == 0 && nProperties == 0)
    {
        cout << "Fatal error: neither landmarks nor properties are provided!" << endl;
        return EXIT_FAILURE;
    }
    else if (nProperties / nSubj != nWeight)
    {
        cout << "Fatal error: # of properties is incosistent with # of weighting factors!" << endl;
        return EXIT_FAILURE;
    }
    
    for (int i = 0; i < nSubj; i++) sphere[i] = listSphere[i].c_str();
    for (int i = 0; i < nProperties; i++) property[i] = listProperty[i].c_str();
    for (int i = 0; i < nOutput; i++) output[i] = listOutput[i].c_str();
    for (int i = 0; i < nLandmark; i++) landmark[i] = listLandmark[i].c_str();
    for (int i = 0; i < nCoeff; i++) coeff[i] = listCoeff[i].c_str();
    for (int i = 0; i < nSurf; i++) surf[i] = listSurf[i].c_str();
    for (int i = 0; i < nWeight; i++) weight[i] = listWeight[i];
    if (nWeight == 0) for (int i = 0; i < nProperties / nSubj; i++) weight[i] = 1;
    
    // display for lists of files
    cout << "Property: " << nProperties / nSubj << endl;	for (int i = 0; i < nProperties; i++) cout << property[i] << endl;
    cout << "Sphere: " << nSubj << endl;					for (int i = 0; i < nSubj; i++) cout << sphere[i] << endl;
    cout << "Output: " << nOutput << endl;					for (int i = 0; i < nOutput; i++) cout << output[i] << endl;
    cout << "Landmark: " << nLandmark << endl;				for (int i = 0; i < nLandmark; i++) cout << landmark[i] << endl;
    cout << "Coefficient: " << nCoeff << endl;				for (int i = 0; i < nCoeff; i++) cout << coeff[i] << endl;
    cout << "Surface: " << nSurf << endl;					for (int i = 0; i < nSurf; i++) cout << surf[i] << endl;
    
    try{
        GroupwiseRegistration groups(sphere, nSubj, property, nProperties / nSubj, output, weight, degree, landmark, weightLoc, coeff, surf, maxIter);
        groups.run();
        
        // delete memory allocation
        delete [] property;
        delete [] sphere;
        delete [] output;
        delete [] landmark;
        delete [] coeff;
        delete [] surf;
        delete [] weight;    
    }catch(exception e){

        cerr<<e.what()<<endl;
    }
    
    
    return EXIT_SUCCESS;
}
