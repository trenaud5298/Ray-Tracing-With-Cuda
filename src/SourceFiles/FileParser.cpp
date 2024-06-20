#include "../HeaderFiles/FileParser.h"
#include <iostream>

FileParser::FileParser(fs::path sourcePath, bool useRelativePaths) : sourcePath(sourcePath),useRelativePaths(useRelativePaths) {
    std::cout<<"Searching Files From: "<<this->sourcePath.string()<<std::endl;
    this->searchFiles(this->sourcePath);
    this->printAllFiles();
}

FileParser::~FileParser() {
    std::cout<<"File Parser Deconstructing"<<std::endl;
    
}


void FileParser::searchFiles(fs::path searchDirectory) {
    this->folders[searchDirectory.filename().string()] = searchDirectory;
    for (const auto& entry : fs::directory_iterator(searchDirectory)) {
        // Check if the entry is a regular file
        if (fs::is_regular_file(entry.path())) {
            std::string extension = entry.path().extension().string();

            if(extension == ".world") {
                /*ADDS RELATIVE/ABSOLUTE WORLD PATH TO WORLD FILES MAP*/
                this->worldFiles[entry.path().stem().string()] = this->useRelativePaths ? fs::relative(entry.path(),this->sourcePath) : entry.path();

            }else if(extension == ".ini") {
                /*CHECKS IF IT IS THE SETTINGS.INI FILE AND IF SO SETS THAT FOR THE CONFIG FILE*/
                this->configFile = this->useRelativePaths ? fs::relative(entry.path(),this->sourcePath) : entry.path();
                /*ADDS RELATIVE/ABSOLUTE INI PATH TO INI FILES MAP*/
                this->iniFiles[entry.path().stem().string()] = this->useRelativePaths ? fs::relative(entry.path(),this->sourcePath) : entry.path();


            }else if(extension == ".shape") {
                /*ADDS RELATIVE/ABSOLUTE SHAPE PATH TO SHAPE FILES MAP*/
                this->shapeFiles[entry.path().stem().string()] = this->useRelativePaths ? fs::relative(entry.path(),this->sourcePath) : entry.path();

            }else if(extension == ".ttf") {
                /*ADDS RELATIVE/ABSOLUTE FONT PATH TO FONT FILES MAP*/
                this->fontFiles[entry.path().stem().string()] = this->useRelativePaths ? fs::relative(entry.path(),this->sourcePath) : entry.path();
            }


        } else if (fs::is_directory(entry.path())) {
            this->searchFiles(entry.path());

        } else {
            std::cout<<"Error Unable to Recognize File:"<<entry.path()<<std::endl;
        };
    };
}

std::string FileParser::getRelativePath(fs::path file) {
    std::string relPath = fs::relative(file,this->sourcePath).string();
    std::replace(relPath.begin(),relPath.end(),'\\','/');
    if (relPath.front() != '.') {
        relPath = "./" + relPath;
    }

    return relPath;
}

void FileParser::printAllFiles() {
    /*PRINT ALL FOLDERS*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Folders Found:"<<std::endl;
    for(std::pair<std::string,fs::path> folder : this->folders) {
        std::cout<<"--"<<folder.first<<": "<<folder.second.string()<<std::endl;
    }
    std::cout<<"----------------------------------"<<std::endl<<std::endl;


    /*PRINT WORLD FILES*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"World Files:"<<std::endl;
    for(std::pair<std::string,fs::path> world : this->worldFiles) {
        std::cout<<"--"<<world.first<<": "<<world.second.string()<<std::endl;
    }
    std::cout<<"----------------------------------"<<std::endl<<std::endl;


    /*PRINT FONT FILES*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Font Files:"<<std::endl;
    for(std::pair<std::string,fs::path> font : this->fontFiles) {
        std::cout<<"--"<<font.first<<": "<<font.second.string()<<std::endl;
    }
    std::cout<<"----------------------------------"<<std::endl<<std::endl;


    /*PRINT SHAPE FILES*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Shape Files:"<<std::endl;
    for(std::pair<std::string,fs::path> shape : this->shapeFiles) {
        std::cout<<"--"<<shape.first<<": "<<shape.second.string()<<std::endl;
    }
    std::cout<<"----------------------------------"<<std::endl<<std::endl;


    /*PRINT INI FILES*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Ini Files:"<<std::endl;
    for(std::pair<std::string,fs::path> iniFile : this->iniFiles) {
        std::cout<<"--"<<iniFile.first<<": "<<iniFile.second.string()<<std::endl;
    }
    std::cout<<"----------------------------------"<<std::endl<<std::endl;


    /*PRINT CONFIG FILE*/
    std::cout<<"----------------------------------"<<std::endl;
    std::cout<<"Config File:"<<std::endl;
    std::cout<<"--"<<this->configFile.stem().string()<<": "<<this->configFile.string()<<std::endl;
    std::cout<<"----------------------------------"<<std::endl<<std::endl;
}

fs::path FileParser::getDirectoryPath(std::string folderName) {
    auto folderExists = this->folders.find(folderName);
    if(folderExists != this->folders.end()) {
        /*FOLDER WAS FOUND*/
        return this->folders[folderName];
    } else {
        /*FOLDER NOT FOUND*/
        return "";
    }
}

std::string FileParser::getDirectoryPathString(std::string folderName) {
    auto folderExists = this->folders.find(folderName);
    if(folderExists != this->folders.end()) {
        /*FOLDER WAS FOUND*/
        return this->folders[folderName].string();
    } else {
        /*FOLDER NOT FOUND*/
        return "";
    }
}

fs::path FileParser::getFontFile(std::string fontName) {
    auto fontExists = this->fontFiles.find(fontName);
    if(fontExists != this->fontFiles.end()) {
        /*FONT WAS FOUND*/
        return this->fontFiles[fontName];
    } else {
        /*FONT NOT FOUND*/
        return "";
    }
}

fs::path FileParser::getWorldFile(std::string worldName) {
    auto worldExists = this->worldFiles.find(worldName);
    if(worldExists != this->worldFiles.end()) {
        /*WORLD WAS FOUND*/
        return this->worldFiles[worldName];
    } else {
        /*WORLD NOT FOUND*/
        return "";
    }
}

fs::path FileParser::getShapeFile(std::string shapeName) {
    auto shapeExists = this->shapeFiles.find(shapeName);
    if(shapeExists != this->shapeFiles.end()) {
        /*SHAPE WAS FOUND*/
        return this->shapeFiles[shapeName];
    } else {
        /*SHAPE NOT FOUND*/
        return "";
    }
}

fs::path FileParser::getIniFile(std::string iniFileName) {
    auto iniFileExists = this->iniFiles.find(iniFileName);
    if(iniFileExists != this->iniFiles.end()) {
        /*INI FILE WAS FOUND*/
        return this->iniFiles[iniFileName];
    } else {
        /*INI FILE NOT FOUND*/
        return "";
    }
}

