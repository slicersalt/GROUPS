#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <dirent.h>
#include "GRTestProjectCLP.h"
#include "GroupwiseRegistration.h"

const char * VTKSUFFIX = ".vtk";
const char * TXTSUFFIX = ".txt";
const char * CORRPREFIX = "corr";
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

std::vector<std::string> sortFiles (std::string path, std::string ID){
	
	std::vector<std::string> directoryFileList;
	directoryFileList = openDirectory(path);
	std::vector<std::string> sortedFileList;
	std::string absPathFileName;
	if (directoryFileList.size() > 0) {
		bool doesHaveProp = false;
		for (int listIndex = 0; listIndex < directoryFileList.size(); listIndex++) {
			std::string currentString = directoryFileList[listIndex];
			absPathFileName = path;
			int substring = currentString.find(ID);
		  	if (substring != std::string::npos) {
				absPathFileName = absPathFileName + "/" + currentString;
		    		sortedFileList.push_back(absPathFileName);
				doesHaveProp = true;
		  	}
		}
		if(!doesHaveProp)
			std::cout << "\nWARNING: PROPERTIES FOR " << ID << " DO NOT EXIST.\n" << std::endl ;
	}

	return sortedFileList;
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

	return sortedFileList;
}

char **putFilesInCharList (std::vector<std::string> files, int size){
	
//	int size = files.size();
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
	
	for(int i = 0; i < extSize; i++){						//filters based on desired ext
		extSpecifiedTempPropList = sortFiles(tempDir, extList[i]);
		for(int j = 0; j < extSpecifiedTempPropList.size(); j++){
			extSpecifiedTempPropFiles.push_back(extSpecifiedTempPropList[j]);
		}
	}
	int extTempPropSize = extSpecifiedTempPropFiles.size();
	char **tempProp = NULL;
	tempProp = putFilesInCharList(extSpecifiedTempPropFiles, extTempPropSize);
	std::cout<<"Template model properties: "<< tempDir << std::endl ;
	for(int i = 0; i < extTempPropSize; i++){
		std::cout<< "	" << tempProp[i] << ' ' << std::endl ;
	}

	char *sph = (char *) sphere.c_str();
	char *log = (char *) logfile.c_str();

 	std::cout<<"Sphere: " << sph << std::endl ;
//	std::cout<<"tmpDepth: " << tmpDepth << std::endl ;
	std::cout<<"Degree: " << degree << std::endl ;
	std::cout<<"Logfile: " << logfile << std::endl ;
	
	std::vector<std::string> surfaceFiles, coefficientFiles, allPropertyFiles, extSpecifiedPropFiles, mapFiles;
	
//Handling of Surface directory**********************************************************
	
	surfaceFiles = sortFiles(surfDir, VTKSUFFIX);
	int surfSize = surfaceFiles.size();
	char **surfFileList = NULL;
	surfFileList = putFilesInCharList(surfaceFiles, surfSize);
	std::cout<<"Surface Directory: "<< surfDir << " " << surfSize << std::endl ;
	for(int i = 1; i < surfSize; i++)
		std::cout<< "	" << surfFileList[i] << ' ' << std::endl ;

//Handling of Property directory**********************************************************

	std::string surfFileWithoutPath;
	std::vector<std::string> foundPropFiles ;
	for(int i = 0; i < surfSize; i++){					//puts all props w/corresponding surf in list
		surfFileWithoutPath = removePathAndSuffix(surfFileList[i]);
		foundPropFiles = sortFiles(propDir, surfFileWithoutPath);
		for(int j = 0; j < foundPropFiles.size(); j++)
			allPropertyFiles.push_back(foundPropFiles[j]);
	}
	std::vector<std::string> extSpecifiedPropList ;
	for(int i = 0; i < extSize; i++){						//filters props based on desired ext
		extSpecifiedPropList = sortFiles(allPropertyFiles, extList[i]);
		for(int j = 0; j < extSpecifiedPropList.size(); j++){
			extSpecifiedPropFiles.push_back(extSpecifiedPropList[j]);
		}
	}
	int extPropSize = extSpecifiedPropFiles.size();
	char **propFileList = NULL;
	propFileList = putFilesInCharList(extSpecifiedPropFiles, extPropSize);
	std::cout<<"Property Directory: "<< propDir << std::endl ;
	for(int i = 1; i < extPropSize; i++){
		std::cout<< "	" << propFileList[i] << ' ' << std::endl ;
	}
	
//Handling of Coefficient directory********************************************************

	coefficientFiles = sortFiles(spharmDir, COEFFSUFFIX);
	int coeffSize = coefficientFiles.size();
	char **coeffFileList = NULL;
	coeffFileList = putFilesInCharList(coefficientFiles, coeffSize);
	std::cout<<"Coefficient Directory: "<< spharmDir << std::endl ;
	for(int i = 1; i < coeffSize; i++)
		std::cout<< "	" << coeffFileList[i] << ' ' << std::endl ;

/*Handling of Map directory (delete) *****************************************************************

	mapFiles = sortFiles(mapDir, CORRPREFIX);
	int mapSize = mapFiles.size();
	char **mapFileList = putFilesInCharList(mapFiles, mapSize);
	std::cout<<"Map Directory: "<< mapDir << std::endl ;
	for(int i = 1; i < mapSize; i++)
		std::cout<< "	" << mapFileList[i] << ' ' << std::endl ;
*/
	
	// char *sphere, char **tmpDepth, char **subjDepth, int nSubj, int deg, int nProperties, char *coeffLog, char **coeff
	GroupwiseRegistration *r = new GroupwiseRegistration(sph, tempProp, propFileList, surfSize, degree, extSize, log);
	
	return 0 ;	
} 													//end of main

