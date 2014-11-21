#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <dirent.h>
#include "GRTestProjectCLP.h"
#include <sstream>

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
	std::vector<std::string>::iterator it;
	if (directoryFileList.size() > 0) {
		bool doesHaveProp = false;
		for (int listIndex = 0; listIndex < directoryFileList.size(); listIndex++) {
			std::string currentString = directoryFileList[listIndex];
			it = directoryFileList.begin() + listIndex;
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
/*			if(!doesHaveProp){
				std::cout << "\nWARNING: PROPERTIES FOR " << ID << " DO NOT EXIST.\n" << std::endl ;
				sortedFileList.erase(it);
//				std::cout<< sortedFileList[] << std::endl ;
			}*/
		}

/*		std::cout << "size: " << sortedFileList.size() << std::endl;
		if(!doesHaveProp){
			std::cout << "\nWARNING: PROPERTIES FOR " << ID << " DO NOT EXIST.\n" << std::endl ;
			sortedFileList.erase(it);
//			std::cout<< sortedFileList[] << std::endl ;
		}
		std::cout << "new size: " << sortedFileList.size() << std::endl;*/
	}
	return sortedFileList;
}
std::vector<std::string> sortFiles (std::string path, std::string ID){
	
	std::vector<std::string> directoryFileList;
	directoryFileList = openDirectory(path);
	std::vector<std::string> sortedFileList;
	std::string absPathFileName;
	int substring;
	std::vector<std::string>::iterator it;
	if (directoryFileList.size() > 0) {
		bool doesHaveProp = false;
		for (int listIndex = 0; listIndex < directoryFileList.size(); listIndex++) {
			std::string currentString = directoryFileList[listIndex];
			it = directoryFileList.begin() + listIndex;
			absPathFileName = path;
			substring = currentString.find(ID);
		  	if (substring != std::string::npos) {
				absPathFileName = absPathFileName + "/" + currentString;
		    		sortedFileList.push_back(absPathFileName);
		  	}
		}
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
	std::ostringstream ss0 ;
	ss0 << tempDir ;
	std::string tempDirString = ss0.str() ;
	if(!tempDirString.empty()){
		for(int i = 0; i < extSize; i++){						//filters based on desired ext
			extSpecifiedTempPropList = sortFiles(tempDir, extList[i]);
			for(int j = 0; j < extSpecifiedTempPropList.size(); j++){
				extSpecifiedTempPropFiles.push_back(extSpecifiedTempPropList[j]);
			}
		}
	}	
	else
		std::cout << "Warning: No input for '-t'." << std::endl ;

	int extTempPropSize = extSpecifiedTempPropFiles.size();
	char **tempProp = putFilesInCharList(extSpecifiedTempPropFiles, extTempPropSize);
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

	std::ostringstream ss ;
	ss << surfDir ;
	std::string surfDirString = ss.str() ;
	std::ostringstream ss1 ;
	ss1 << propDir ;
	std::string propDirString = ss1.str() ;
	if(!surfDirString.empty()){
		if(!propDirString.empty())
			surfaceFiles = sortSurfFiles(surfDir, propDir, VTKSUFFIX);
		else
			std::cout << "Warning: No input for '-p'." << std::endl ;
	}
	int surfSize = surfaceFiles.size();
	std::cout << "Surf size: " << surfSize << std::endl;
	char **surfFileList = putFilesInCharList(surfaceFiles, surfSize);
	std::cout<<"Surface Directory: "<< surfDir << std::endl ;
	for(int i = 1; i < surfSize; i++)
		std::cout<< "	" << surfFileList[i] << ' ' << std::endl ;

//Handling of Property directory**********************************************************

	std::string surfFileWithoutPath;
	std::vector<std::string> foundPropFiles ;
	if(!propDirString.empty()){
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
		for(int i = 0; i < extPropSize; i++){					//removes files w/ "974644"
			std::vector<std::string>::iterator it = extSpecifiedPropFiles.begin() +i;
			std::string theCurrentFile = extSpecifiedPropFiles[i];
			if(theCurrentFile.find("974644") != std::string::npos){
				extSpecifiedPropFiles.erase(it);
//				std::cout<< theCurrentFile << std::endl;
			}
		}
		int newExtPropSize = extSpecifiedPropFiles.size();
		char **propFileList = putFilesInCharList(extSpecifiedPropFiles, newExtPropSize);
//		std::cout<< "HERE" << std::endl;
		std::cout<<"Property Directory: "<< propDir << std::endl ;
		for(int i = 0; i < newExtPropSize; i++){
			std::cout<< "	" << propFileList[i] << ' ' << std::endl ;
		}
	}
//Handling of Coefficient directory********************************************************
	
	std::ostringstream ss2 ;
	ss2 << spharmDir ;
	std::string spharmDirString = ss2.str() ;
	if(!spharmDirString.empty()){
		coefficientFiles = sortFiles(spharmDir, COEFFSUFFIX);
		int coeffSize = coefficientFiles.size();
		char **coeffFileList = putFilesInCharList(coefficientFiles, coeffSize);
		std::cout<<"Coefficient Directory: "<< spharmDir << std::endl ;
		for(int i = 1; i < coeffSize; i++)
			std::cout<< "	" << coeffFileList[i] << ' ' << std::endl ;
	}
/*Handling of Map directory (delete) *****************************************************************

	mapFiles = sortFiles(mapDir, CORRPREFIX);
	int mapSize = mapFiles.size();
	char **mapFileList = putFilesInCharList(mapFiles, mapSize);
	std::cout<<"Map Directory: "<< mapDir << std::endl ;
	for(int i = 1; i < mapSize; i++)
		std::cout<< "	" << mapFileList[i] << ' ' << std::endl ;
*/
	return 0 ;	
} 													//end of main

