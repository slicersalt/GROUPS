#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <dirent.h>
#include "GroupsCLP.h"
#include <sstream>
#include "GroupwiseRegistration.h"

const char * VTKSUFFIX = ".vtk";
const char * TXTSUFFIX = ".txt";
const char * COEFFSUFFIX = ".coef";

std::vector<std::string> openDirectory (std::string path)
{
	DIR* dir;
	dirent* pdir;
	std::vector<std::string> files;

	dir = opendir(path.c_str());

	while(pdir = readdir(dir)) {
	  files.push_back(std::string(pdir->d_name));
	}
	return files;
}

std::string removePathAndSuffix(std::string fileName){
	std::string fileWithoutPath = fileName.substr(fileName.find_last_of("/")+1, fileName.size());
	fileWithoutPath = fileWithoutPath.substr(0, fileWithoutPath.find("."));
	return fileWithoutPath;
}

bool checkPropList(std::string propDirPath, std::string fileName){
	std::vector<std::string> propDirList;
	propDirList = openDirectory(propDirPath);
	bool hasProp = false;
	for(int i = 0; i < propDirList.size(); i++){
		if(propDirList[i].find(fileName) != std::string::npos)
			hasProp = true;
	}
	return hasProp;
}

std::vector<std::string> sortSurfFiles (std::string surfPath, std::string propPath, std::string ID){
	
	std::vector<std::string> directoryFileList;
	directoryFileList = openDirectory(surfPath);
	std::vector<std::string> sortedFileList;
	std::string absPathFileName;
	int substring;
//	std::vector<std::string>::iterator it;
	if (directoryFileList.size() > 0) {
		bool doesHaveProp = false;
		for (int listIndex = 0; listIndex < directoryFileList.size(); listIndex++) {
			std::string currentString = directoryFileList[listIndex];
//			it = directoryFileList.begin() + listIndex;
			absPathFileName = surfPath;
			substring = currentString.find(ID);
		  	if (substring != std::string::npos) {
				absPathFileName = absPathFileName + "/" + currentString;
		    	sortedFileList.push_back(absPathFileName);
				doesHaveProp = checkPropList(propPath, currentString);
				if(!doesHaveProp){
					std::cout << "\nWARNING: PROPERTIES FOR " <<currentString<< " DO NOT EXIST. IT HAS BEEN REMOVED\n" << std::endl ;
					sortedFileList.erase(sortedFileList.begin() + listIndex);
				}
		  	}
		}
	}
	std::sort(sortedFileList.begin(), sortedFileList.end());
	return sortedFileList;
}
std::vector<std::string> sortFiles (std::string path, std::string ID){
	
	std::vector<std::string> directoryFileList;
	directoryFileList = openDirectory(path);
	std::vector<std::string> sortedFileList;
	std::string absPathFileName;
	int substring;
//	std::vector<std::string>::iterator it;
	if (directoryFileList.size() > 0) {
		for (int listIndex = 0; listIndex < directoryFileList.size(); listIndex++) {
			std::string currentString = directoryFileList[listIndex];
//			it = directoryFileList.begin() + listIndex;
			absPathFileName = path;
			substring = currentString.find(ID);
		  	if (substring != std::string::npos) {
				absPathFileName = absPathFileName + "/" + currentString;
		    		sortedFileList.push_back(absPathFileName);
		  	}
		}
	}
	std::sort(sortedFileList.begin(), sortedFileList.end());
	return sortedFileList;
}

std::vector<std::string> setOutputFiles (std::vector<std::string> surfFiles, std::string outPath, std::string ext, std::string coeff){
	std::vector<std::string> outFiles;
	std::string surfFileName;
	std::string outFileName;
	for(int i = 0; i < surfFiles.size(); i++){
		surfFileName = removePathAndSuffix(surfFiles[i]);
		outFileName = outPath + "/" + surfFileName + coeff + ext;
		outFiles.push_back(outFileName);
	}
	return outFiles;
}

std::vector<std::string> sortFiles (std::vector<std::string> list, std::string ID){
	
	std::vector<std::string> sortedFileList;

	if (list.size() > 0) {
		for (int listIndex = 0; listIndex < list.size(); listIndex++) {
			std::string currentString = list[listIndex];
			int substring = currentString.find(ID);
		  	if (substring != std::string::npos) {
		    		sortedFileList.push_back(currentString);
		  	}
		}
	}
	std::sort(sortedFileList.begin(), sortedFileList.end());
	return sortedFileList;
}

char **putFilesInCharList (std::vector<std::string> files, int size){
	char **fileList = NULL;
	fileList = (char **) malloc(sizeof(char *) * size);
	for (int i = 0; i < size; i++)
		fileList[i] = (char *)malloc(sizeof(char) * 1024);
	for (int i = 0; i < size; i++)
		strcpy(fileList[i], files[i].c_str());
	return fileList;
}

int main(int argc , char* argv[])
{
	PARSE_ARGS;
	
	std::vector<std::string> tempPropFiles;
	std::vector<std::string> extList;
	std::vector<std::string> extSpecifiedTempPropList ;
	std::vector<std::string> extSpecifiedTempPropFiles ;
	std::string ext;
	
	int tempSize = tempPropFiles.size();
	extList = extensions;
	int extSize = extList.size();
	std::ostringstream ss0 ;
	ss0 << tempDir ;
	std::string tempDirString = ss0.str() ;
	char **tempProp;
	if(!tempDirString.empty()){
		for(int i = 0; i < extSize; i++){				//filters based on desired ext
			extSpecifiedTempPropList = sortFiles(tempDir, extList[i]);
			for(int j = 0; j < extSpecifiedTempPropList.size(); j++){
				extSpecifiedTempPropFiles.push_back(extSpecifiedTempPropList[j]);
			}
		}
		int extTempPropSize = extSpecifiedTempPropFiles.size();
		tempProp = putFilesInCharList(extSpecifiedTempPropFiles, extTempPropSize);
		std::cout<<"Template model properties: "<< tempDir << std::endl ;
		for(int i = 0; i < extTempPropSize; i++){
			std::cout<< "	" << tempProp[i] << ' ' << std::endl ;
		}
	}	
	else{
		std::cout << "Warning: No input for '-t'." << std::endl ;
		tempProp = NULL;
	}
	
	char *sph = (char *) sphere.c_str();
	char *log = (char *) logfile.c_str();
	if(logfile.empty()) log = NULL;

 	std::cout<<"Sphere: " << sph << std::endl ;
	std::cout<<"Degree: " << degree << std::endl ;
	std::cout<<"Logfile: " << logfile << std::endl ;
 	std::cout<<"Max Iterations: " << maxIter << std::endl ;
 	std::cout<<"addLoc: " << addLoc << std::endl ;
 	std::cout<<"Temp Surf: " << tempSurface << std::endl ;
	std::vector<std::string> surfaceFiles, coefficientFiles, allPropertyFiles, extSpecifiedPropFiles, landmarkFiles, asphereFiles, outputFiles;
	
//Handling of Surface directory**********************************************************

	std::ostringstream ss ;
	ss << surfDir ;
	std::string surfDirString = ss.str() ;
	std::ostringstream ss1 ;
	ss1 << propDir ;
	std::string propDirString = ss1.str() ;
	char **surfFileList;
	int surfSize;
	if(!surfDirString.empty()){
		if(!propDirString.empty()){
			surfaceFiles = sortSurfFiles(surfDir, propDir, VTKSUFFIX);
			surfSize = surfaceFiles.size();
			surfFileList = putFilesInCharList(surfaceFiles, surfSize);
			std::cout<<"Surface Directory: "<< surfDir << std::endl ;
			for(int i = 0; i < surfSize; i++)
			std::cout<< "	" << surfFileList[i] << ' ' << std::endl ;
		}
		else
		{
			surfSize = 0;
			surfFileList = NULL;
			std::cout << "Warning: No input for '-p'." << std::endl ;
		}
	}
	else
		surfFileList = NULL;

//Handling of Property directory**********************************************************

	std::string surfFileWithoutPath;
	std::vector<std::string> foundPropFiles ;
	char **propFileList;
	if(!propDirString.empty()){
		for(int i = 0; i < surfSize; i++){				//puts all props w/corresponding surf in surfFileList
			surfFileWithoutPath = removePathAndSuffix(surfFileList[i]);
			foundPropFiles = sortFiles(propDir, surfFileWithoutPath);
			for(int j = 0; j < foundPropFiles.size(); j++)
				allPropertyFiles.push_back(foundPropFiles[j]);
		}
		std::vector<std::string> extSpecifiedPropList ;
		for(int i = 0; i < extSize; i++){				//filters props based on desired ext
			extSpecifiedPropList = sortFiles(allPropertyFiles, extList[i]);
			for(int j = 0; j < extSpecifiedPropList.size(); j++){
			extSpecifiedPropFiles.push_back(extSpecifiedPropList[j]);
			}
		}
		int extPropSize = extSpecifiedPropFiles.size();
		int newExtPropSize = extSpecifiedPropFiles.size();
		propFileList = putFilesInCharList(extSpecifiedPropFiles, newExtPropSize);
		std::cout<<"Property Directory: "<< propDir << std::endl ;
		for(int i = 0; i < newExtPropSize; i++){
			std::cout<< "	" << propFileList[i] << ' ' << std::endl ;
		}
	}
	else
		propFileList = NULL;

//Handling of Coefficient directory********************************************************
	
	std::ostringstream ss2 ;
	ss2 << spharmDir ;
	std::string spharmDirString = ss2.str() ;
	char **coeffFileList;
	if(!spharmDirString.empty()){
		coefficientFiles = sortFiles(spharmDir, COEFFSUFFIX);
		int coeffSize = coefficientFiles.size();
		coeffFileList = putFilesInCharList(coefficientFiles, coeffSize);
		std::cout<<"Coefficient Directory: "<< spharmDir << std::endl ;
		for(int i = 0; i < coeffSize; i++)
			std::cout<< "	" << coeffFileList[i] << ' ' << std::endl ;
	}
	else
		coeffFileList = NULL;

//Output Directory*************************************************************************
	
	outputFiles = setOutputFiles(surfaceFiles, outputDir, TXTSUFFIX, ".coeff");
	int outSize = outputFiles.size();
	char **outFileList = putFilesInCharList(outputFiles, outSize);
	std::cout << "Output Directory: " << outputDir << std::endl ;
	for(int i = 0; i < outSize; i++){
		std::cout << "	" << outFileList[i] << ' ' << std::endl ;
	}
	
	char *tempSurf = NULL;
	if (tempSurface.length() > 0)
	{
		tempSurf = new char[tempSurface.length() + 1];
		std::strcpy(tempSurf, tempSurface.c_str());
	}
	
	// landmark directory
	char **landmark = NULL;
	if(!lmDir.empty()){
		landmarkFiles = sortFiles(lmDir, TXTSUFFIX);
		int lmSize = landmarkFiles.size();
		landmark = putFilesInCharList(landmarkFiles, lmSize);
		std::cout<<"Landmark Directory: "<< lmDir << std::endl ;
		for(int i = 0; i < lmSize; i++)
			std::cout<< "	" << landmark[i] << ' ' << std::endl ;
	}
	
	// aligned sphere directory
	char **alignedSphere = NULL;
	if(!asphereDir.empty()){
		asphereFiles = sortFiles(asphereDir, TXTSUFFIX);
		int asphereSize = asphereFiles.size();
		alignedSphere = putFilesInCharList(asphereFiles, asphereSize);
		std::cout<<"Alinged Sphere Directory: "<< asphereDir << std::endl ;
		for(int i = 0; i < asphereSize; i++)
			std::cout<< "	" << alignedSphere[i] << ' ' << std::endl ;
	}
	
	float *w = new float[extSize];
	if (!weight.empty())
	{
		if (weight.size() != extSize)
		{
			cout << "# of the weighting factors is not matched to # of the properties" << endl;
			return 1;
		}
		for (int i = 0; i < weight.size(); i++)
			w[i] = weight[i];
	}
	else
	{
		for (int i = 0; i < extSize; i++)
			w[i] = 1;
	}
	cout << "Weighting Factors\n";
	for (int i = 0; i < weight.size(); i++)
	{
		cout << extensions[i] << " - " << w[i] << endl;
	}
	if (addLoc != 0)
		cout << "Spatial weighting factor: " << addLoc << endl;

// Call main procedure
	// char *sphere, char **tmpDepth, char **subjDepth, int nSubj, int deg, int nProperties, char *coeffLog, char **coeff
	GroupwiseRegistration *r = new GroupwiseRegistration(sph, tempProp, propFileList, surfSize, landmark, alignedSphere, w, degree, extSize, addLoc, tempSurf, surfFileList, log, coeffFileList, maxIter, outFileList);
	//GroupwiseRegistration *r = new GroupwiseRegistration(sph, NULL, NULL, surfSize, degree, 0, addLoc, tempSurf, surfFileList, log, coeffFileList, maxIter);	// location only
	
	for (int i = 0; i < outSize; i++)
	{
		r->saveLCoeff(outFileList[i], i);
	}
	
	return 0 ;
	
}//end of main
