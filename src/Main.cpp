#include <iostream>
#include "./HeaderFiles/Display.cuh"
#include "./HeaderFiles/Camera.h"
#include "./HeaderFiles/Settings.h"
#include "./HeaderFiles/FileParser.h"
#include "./HeaderFiles/EventManager.h"
#include <filesystem>

#include <thread>
#include <chrono>

#undef main

int main(int argc, char* argv[]) {
    std::cout<<"Hello World"<<std::endl;
    //Creates File Parse to Get all Important Files
    FileParser* mainFileReader = new FileParser(std::filesystem::current_path(), false);
    //Creates Settings Menu To Read/Write Ini File
    Settings* mainSettings = new Settings(mainFileReader->configFile.string(), mainFileReader->getDirectoryPathString("Settings"));
    mainSettings->printAllSettings();
    
    EventManager* programManager = new EventManager(mainSettings,mainFileReader);
    programManager->initializeProgram();
    programManager->runProgram();
    std::cout<<"GoodBye World"<<std::endl;
    SDL_Quit();
    TTF_Quit();
    return 0; 
}
