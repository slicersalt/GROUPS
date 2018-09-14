/*************************************************
 *  main.cpp
 *
 *  Release: Sep 2016
 *  Update: Sep 2016
 *
 *  University of North Carolina at Chapel Hill
 *  Department of Computer Science
 *
 *  Ilwoo Lyu, ilwoolyu@cs.unc.edu
 *************************************************/

#include <cstdlib>
#include <vector>
#include <string>

#ifdef WIN32
#include "dirent.h"
#else
#include <dirent.h>
#endif
#include "GroupsCLP.h"
#include "GroupwiseRegistration.h"
#include <iostream>

bool getListFile(string path, vector<string> &list, const string &suffix)
{

    DIR *dir = opendir(path.c_str());
    if (dir != NULL)
    {
        while (dirent *entry = readdir(dir))
        {
            string filename = entry->d_name;
            if(filename.find(suffix) != string::npos && filename.find_last_of(suffix) == filename.size() - 1)
            {
                list.push_back(path + "/" + filename);
            }
        }
        closedir(dir);
        sort(list.begin(), list.begin() + list.size());
        return true;
    }else{
        
        cerr<<"The directory does not exist! "<<path<<endl;
        return false;
    }
}

void getTrimmedList(vector<string> &list, const vector<string> &name)
{
    size_t i = 0;
    while (i < list.size())
    {
        size_t found;
        for (size_t j = 0; j < name.size(); j++)
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
    if (!dirSphere.empty() && listSphere.empty()){
        if(! getListFile(dirSphere, listSphere, ".vtk")){
            return EXIT_FAILURE;
        }
    }
    
    if (!dirSurf.empty() && listSurf.empty()){
        if(! getListFile(dirSurf, listSurf, ".vtk")){
            return EXIT_FAILURE;
        }
    }
    
    if (!dirCoeff.empty() && listCoeff.empty()){
        if(! getListFile(dirCoeff, listCoeff, "coeff")){
            return EXIT_FAILURE;
        }
    }

    if(dirOutput.compare("") == 0){
        cout<<"Setting dirOutput to ./"<<endl;
        dirOutput = "./";
    }else{
        if (opendir(dirOutput.c_str()) == NULL)
        {
            cout<<"The output directory does not exist! "<<dirOutput<<endl;
            return EXIT_FAILURE;
        }
    }

    // subject names
    size_t nSubj = listSphere.size();
    vector<string> subjName;
    if (nSubj > 0)
    {
        for (size_t i = 0; i < nSubj; i++)
        {
            size_t pivot = listSphere[i].rfind('/') + 1;

            std::string filename = listSphere[i].substr(pivot);

            std::string name;
            std::string suffixe_surf_para = "_pp_surf_para.vtk";
            std::string suffixe_para = "_pp_para.vtk";
            size_t suffixe_size;
            
            if ( filename.substr(filename.length() - suffixe_para.length()) == suffixe_para )        
                suffixe_size = suffixe_para.length();

            else if ( filename.substr(filename.length() - suffixe_surf_para.length()) == suffixe_surf_para )
                suffixe_size = suffixe_surf_para.length();

            name = filename.substr(0, filename.length() - suffixe_size);
            subjName.push_back(name);
            std::cout << "Name du subject " << i << " :: " << name << std::endl;
        }
    }
    for (size_t i = 0; i < nSubj; i++) cout << subjName[i] << endl;
    if (listOutput.empty())

    for (size_t i = 0; i < nSubj; i++){
        listOutput.push_back(dirOutput + "/" + subjName[i].substr(0, subjName[i].find_last_of(".")) + ".coeff");
        cout<<listOutput[i]<<endl;
    } 
    
    // trim all irrelevant files to the sphere files
    if (!dirCoeff.empty()) getTrimmedList(listCoeff, subjName);

    std::map<std::string, float> mapProperty;

    try{
        cout<<"ReadingnNumber of properties: "<<modelProperty.size()/2<<endl;
        for(size_t i = 0; i < modelProperty.size(); i+=2){
            string propertyname = modelProperty[i];
            float weight = atof(modelProperty[i+1].c_str());
            cout<<"Reading property: "<<propertyname<<endl;
            cout<<"Reading weight: "<<weight<<endl;
            
            mapProperty[propertyname] = weight;
        }
    }catch(exception e){
        cerr<<e.what()<<endl;
        return EXIT_FAILURE;
    }
    

    // if ( phiOn )
    //     mapProperty["Color_Map_Phi"] = phiWeight;
    // if ( thetaOn )
    //     mapProperty["Color_Map_Theta"] = thetaWeight;
    // if ( CurvednessOn )
    //     mapProperty["Curvedness"] = CurvednessWeight;
    // if ( ShapeIndexOn )
    //     mapProperty["Shape_Index"] = ShapeIndexWeight;
    // if ( MaxCurvOn )
    //     mapProperty["Maximum_Curvature"] = MaxCurvWeight;
    // if ( MinCurvOn )
    //     mapProperty["Minimum_Curvature"] = MinCurvWeight;
    // if ( MeanCurvOn )
    //     mapProperty["Meanimum_Curvature"] = MeanCurvWeight;
    // if ( GaussCurvOn )
    //     mapProperty["Gaussian_Curvature"] = GaussCurvWeight;

    size_t nOutput = listOutput.size();
    // int nCoeff = listCoeff.size();
    size_t nSurf = listSurf.size();
    weightLoc = 0;

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
    // else if (nLandmark == 0 && nProperties == 0)
    // {
    //     cout << "Fatal error: neither landmarks nor properties are provided!" << endl;
    //     return EXIT_FAILURE;
    // }
    // else if (nProperties / nSubj != nWeight)
    // {
    //     cout << "Fatal error: # of properties is incosistent with # of weighting factors!" << endl;
    //     return EXIT_FAILURE;
    // }
    
    for (size_t i = 0; i < nSubj; i++) listSphere[i] = listSphere[i].c_str();
    for (size_t i = 0; i < nOutput; i++) listOutput[i] = listOutput[i].c_str();
    // for (size_t i = 0; i < nCoeff; i++) listCoeff[i] = listCoeff[i].c_str();
    for (size_t i = 0; i < nSurf; i++) listSurf[i] = listSurf[i].c_str();
    // for (size_t i = 0; i < nWeight; i++) listWeight[i] = listWeight[i];
    // if (nWeight == 0) for (size_t i = 0; i < nProperties / nSubj; i++) listWeight[i] = 1;
    

    // display for lists of files
    cout << "Sphere: " << nSubj << endl;                    for (size_t i = 0; i < nSubj; i++) cout << listSphere[i] << endl;
    cout << "Output: " << nOutput << endl;                  for (size_t i = 0; i < nOutput; i++) cout << listOutput[i] << endl;
    cout << "Surface: " << nSurf << endl;                   for (size_t i = 0; i < nSurf; i++) cout << listSurf[i] << endl;
    
// GroupwiseRegistration(vector<string> sphere, vector<string> surf, vector<string> propertiesnames, vector<string> outputcoeff, vector<double> weight, double weightLoc, int deg, vector<string> inputcoeff, int maxIter);

    try{
        GroupwiseRegistration groups(listSphere, listSurf, mapProperty, listOutput, landmarksOn, weightLoc, degree, listCoeff, maxIter);
        // GroupwiseRegistration groups(listSphere, nSubj, listProperty, nProperties / nSubj, listOutput, listWeight, degree, listLandmark, weightLoc, listCoeff, listSurf, maxIter);
        groups.run();
        
        // delete memory allocation
        // delete [] property;
        // delete [] sphere;
        // delete [] output;
        // delete [] landmark;
        // delete [] coeff;
        // delete [] surf;
        // delete [] weight;    
    }catch(exception e){
        cerr<<e.what()<<endl;
        return EXIT_FAILURE;
    }
    
    
    return EXIT_SUCCESS;
}
