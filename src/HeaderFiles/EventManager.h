#ifndef EVENTMANAGER_H
#define EVENTMANAGER_H

#include <thread>
#include "./Camera.h"
#include "./Display.cuh"
#include "./Settings.h"
#include "./FileParser.h"
#include "./Menu.h"
#include "./WorldBuilding/World.h"
#include "./WorldBuilding/Material.h"
#include "./WorldBuilding/Sphere.h"
#include "./WorldBuilding/Triangle.h"
#include "./WorldBuilding/Vertex.h"
#include <ctime>
#include <sstream>
#include <iomanip>
#include <chrono>




class EventManager {

public:
    EventManager(Settings* mainSettings, FileParser* mainFileParser);
    void handleEvents();
    void initializeProgram();
    void runProgram();
    
protected:
    /*WILL RENDER SCENE BASED OFF CURRENT PROGRAM STATE*/
    void renderScene();

private:
    int maxRayBounces;
    int numOfFramesPerRender;
    const Uint8* keyBoardState;
    bool programAlive;
    int previouseProgramState;
    int programState;
    int expectedImageWidth;
    int expectedImageHeight;
    int mouseSensitivityX;
    int mouseSensitivityY;
    FileParser* mainFileParser;
    Settings* mainSettings;
    Camera* viewCamera;
    Display* mainDisplay;
    Menu* mainMenu;
    Menu* pauseMenu;
    SDL_Renderer* mainProgramRenderer;
    SDL_Event programEvent;
    void handleMainMenuEvents();
    void handleRunningEvents();
    void handlePausedEvents();

    /*Functions Involving Selected World*/
    World* selectedWorld;
    void buildWorld();

    /*General Functions To Be Called*/
    void updateProgramState(int newProgramState);
    void endProgram();


    /* Functions involving the Main Menu*/
    void buildMainMenu();
    void mainMenuClicked(int mousePosX, int mousePosY, bool isLeftClick);



    /* Functions involving the Paused Menu*/
    void buildPauseMenu();
    void pauseMenuClicked(int mousePosX, int mousePosY, bool isLeftClick);





    /* Functions involving the Running State*/
    void handleRunningKeyPresses(SDL_Keycode keyPressed);
    void updateViewCamera(float deltaX, float deltaY); /*Uses A Decimal Between 0 and 1 to represent change*/
    void handleKeyboardState(); /*DEALS WITH STUFF LIKE MOVEMENT AND MORE*/
    

};





#endif