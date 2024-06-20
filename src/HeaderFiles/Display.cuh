#ifndef DISPLAY_CUH
#define DISPLAY_CUH

#include <iostream>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include "../../lib/SDL/include/SDL.h"
#include "./Camera.h"
#include "./Button.h"
#include "./Settings.h"
#include "./WorldBuilding/World.h"
#include "./WorldBuilding/Sphere.h"
#include <random>


class Display {

public:
    Display(int displayWidth,int displayHeight,Camera* viewCamera);
    ~Display();
    void updateDisplayParameters(Settings* settings);

    void getPixelInfo(int x,int y);
    void renderImage();
    void seriousRenderOfImage(int numOfFrames,int maxRayBounces);
    void customResolutionRender(std::string fileName, int numOfFrames,int maxRayBounces,int resolutionX, int resolutionY);
    void renderMainMenu();
    void renderPauseMenu();
    void renderScene();
    void saveImage(const std::string& fileName);
    void updateDisplay();
    void frameUpdate();
    Camera* viewCamera;
    int displayState;
    /*
    0=MainMenu
    1=Running
    2=Paused
    */



protected:
    void allocateWorldData(World* world);
    int numOfWorldObjects;
    void allocateGPUMemory();
    void reAllocateGPUMemory();
    void updateGPUData();
    void copyImageData();
    void handleEvents();
    void resetCursor();
    void handleDisplaySizeChange(int newDisplayWidth, int newDisplayHeight);
    unsigned char* imageData; //RGB (0-255) VERSION OF IMAGE THAT IS COPIED BACK TO CPU
    unsigned char* gpuImage; 
    float* rawGpuImageData;//HOW THE ACTUAL FLOATING POITN VERSION OF IMAGE IS STORED AND EDITED ON GPU
    float* gpuCameraData;
    int displayWidth;
    int displayHeight;
    bool paused;
    SDL_Window* displayWindow;
    SDL_Renderer* imageRenderer;
    SDL_Texture* imageTexture;
    float** gpuWorldObjectData;

    float* gpuTroubleShootData;
    float* recieveGpuTroubleShootData;

private:
    int currentFrame;
    int previousMousePosX;
    int previousMousePosY;
    int mousePosX;
    int mousePosY;
    bool sizeChanged;

    /*Random Number Generator Stuff:*/
    std::mt19937 randomSeedGenerator;
    std::uniform_real_distribution<float> seedDistribution;
    
    

    /*
    Grants the EventManager Class Access to Display Class Protected and Private Members
    */
    friend class EventManager;



};



#endif