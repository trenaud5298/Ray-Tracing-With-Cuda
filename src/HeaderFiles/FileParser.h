#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <filesystem>
#include <vector>
#include <string>
#include <map>
namespace fs = std::filesystem;

class FileParser {

public:
    FileParser(fs::path sourcePath, bool useRelativePaths);
    ~FileParser();
    fs::path sourcePath;
    std::map<std::string, fs::path> worldFiles;
    std::map<std::string, fs::path> fontFiles;
    std::map<std::string, fs::path> shapeFiles;
    std::map<std::string, fs::path> iniFiles;
    std::map<std::string, fs::path> folders;
    fs::path configFile;
    void printAllFiles();
    fs::path getDirectoryPath(std::string folderName);
    std::string getDirectoryPathString(std::string folderName);

    fs::path getFontFile(std::string fontName);
    fs::path getWorldFile(std::string worldName);
    fs::path getShapeFile(std::string shapeName);
    fs::path getIniFile(std::string iniFileName);



protected:



private:
    bool useRelativePaths;
    void searchFiles(fs::path searchDirectory);
    std::string getRelativePath(fs::path file);
};



#endif