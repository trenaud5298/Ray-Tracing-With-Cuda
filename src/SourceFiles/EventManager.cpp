
#include "../HeaderFiles/EventManager.h"


#define MAIN_MENU 0
#define RUNNING 1
#define PAUSED 2
#define SETTINGS 3



EventManager::EventManager(Settings* mainSettings, FileParser* mainFileParser) : mainSettings(mainSettings),mainFileParser(mainFileParser),programState(MAIN_MENU),programAlive(true) {

    /* EVENT MANAGER SHOULD BE THE ONE STOP SHOP SO THAT THE MAIN FUNCTION IN A LOOP
     * SHOULD JUST HAVE TO CALL THE RUN PROGRAM (OR SOMETHING SIMILAR) METHOD TO HAVE EVERYTHING WORK
    */
    this->mouseSensitivityX = this->mainSettings->getMainSettingsValue("HorizontalCameraSensitivity");
    this->mouseSensitivityY = this->mainSettings->getMainSettingsValue("VerticalCameraSensitivity");
    this->numOfFramesPerRender = this->mainSettings->getGraphicsSettingsValue("FramesPerRender");
    this->maxRayBounces = this->mainSettings->getGraphicsSettingsValue("MaxRayBounces");
    this->keyBoardState = SDL_GetKeyboardState(NULL);

}




#pragma region Program Init

void EventManager::initializeProgram() {
    /*CREATES OBJECTS*/
    this->viewCamera = new Camera(this->mainSettings);
    this->viewCamera->setFocalPoint(0,-3,2.5);
    this->mainDisplay = new Display(this->mainSettings->getGraphicsSettingsValue("ScreenResolutionX"),this->mainSettings->getGraphicsSettingsValue("ScreenResolutionY"),this->viewCamera);
    this->mainProgramRenderer = this->mainDisplay->imageRenderer;
    this->expectedImageWidth = this->mainDisplay->displayWidth;
    this->expectedImageHeight = this->mainDisplay->displayHeight;
    this->buildMainMenu();
    this->buildPauseMenu();
    this->buildWorld();
    //NOTE
    /*
     * CONSIDER GAMMA CORRECTION
     * EVEN ON NORMAL 
     * 
     * 
     * 
    
    // // */
    this->mainDisplay->customResolutionRender("test",50,15,1920,1080);
    this->endProgram();
    /*NOT ONLY DO BUTTONS NEED TO BE CREATED BUT SO DO COLORS AND STUFF*/
}



void EventManager::buildMainMenu() {
    std::string fontPath = this->mainFileParser->getFontFile("Default").string();
    if(fontPath == "") {
        if(!this->mainFileParser->fontFiles.empty()) {
            auto font = this->mainFileParser->fontFiles.begin();
            fontPath = font->second.string();
        } else {
            std::cout<<"ERROR: NO FONTS FOUND: QUITTING PROGRAM"<<std::endl;
            this->endProgram();
        }
    }

    this->mainMenu = new Menu("Main Menu",this->mainDisplay->displayWidth,this->mainDisplay->displayHeight,this->mainProgramRenderer,fontPath);
    this->mainMenu->setupMenuColors(
        this->mainSettings->getMainMenuSettingsValue("BackgroundColorR"),this->mainSettings->getMainMenuSettingsValue("BackgroundColorG"),this->mainSettings->getMainMenuSettingsValue("BackgroundColorB"),this->mainSettings->getMainMenuSettingsValue("BackgroundOpacity"),
        this->mainSettings->getMainMenuSettingsValue("DefaultButtonColorR"),this->mainSettings->getMainMenuSettingsValue("DefaultButtonColorG"),this->mainSettings->getMainMenuSettingsValue("DefaultButtonColorB"),this->mainSettings->getMainMenuSettingsValue("DefaultButtonOpacity"),
        this->mainSettings->getMainMenuSettingsValue("SelectedButtonColorR"),this->mainSettings->getMainMenuSettingsValue("SelectedButtonColorG"),this->mainSettings->getMainMenuSettingsValue("SelectedButtonColorB"),this->mainSettings->getMainMenuSettingsValue("SelectedButtonOpacity"),
        this->mainSettings->getMainMenuSettingsValue("ButtonTextColorR"),this->mainSettings->getMainMenuSettingsValue("ButtonTextColorG"),this->mainSettings->getMainMenuSettingsValue("ButtonTextColorB"),this->mainSettings->getMainMenuSettingsValue("ButtonTextOpacity")
    );
    Button* startButton = new Button(0.5,0.35,0.2,0.15,"Start");
    Button* settingsButton = new Button(0.5,0.5,0.2,0.15,"Settings");
    Button* quitButton = new Button(0.5,0.65,0.2,0.15,"Quit");

    this->mainMenu->addButton(startButton,"Start Button");
    this->mainMenu->addButton(settingsButton,"Settings Button");
    this->mainMenu->addButton(quitButton,"Quit Button");
    std::cout<<"MAIN MENU BUILT"<<std::endl;
}




void EventManager::buildPauseMenu() {
    std::string fontPath = this->mainFileParser->getFontFile("Default").string();
    if(fontPath == "") {
        if(!this->mainFileParser->fontFiles.empty()) {
            auto font = this->mainFileParser->fontFiles.begin();
            fontPath = font->second.string();
        } else {
            std::cout<<"ERROR: NO FONTS FOUND: QUITTING PROGRAM"<<std::endl;
            this->endProgram();
        }
    }

    this->pauseMenu = new Menu("Pause Menu",this->mainDisplay->displayWidth,this->mainDisplay->displayHeight,this->mainProgramRenderer,fontPath);
    this->pauseMenu->setupMenuColors(
        this->mainSettings->getPauseMenuSettingsValue("BackgroundColorR"),this->mainSettings->getPauseMenuSettingsValue("BackgroundColorG"),this->mainSettings->getPauseMenuSettingsValue("BackgroundColorB"),this->mainSettings->getPauseMenuSettingsValue("BackgroundOpacity"),
        this->mainSettings->getPauseMenuSettingsValue("DefaultButtonColorR"),this->mainSettings->getPauseMenuSettingsValue("DefaultButtonColorG"),this->mainSettings->getPauseMenuSettingsValue("DefaultButtonColorB"),this->mainSettings->getPauseMenuSettingsValue("DefaultButtonOpacity"),
        this->mainSettings->getPauseMenuSettingsValue("SelectedButtonColorR"),this->mainSettings->getPauseMenuSettingsValue("SelectedButtonColorG"),this->mainSettings->getPauseMenuSettingsValue("SelectedButtonColorB"),this->mainSettings->getPauseMenuSettingsValue("SelectedButtonOpacity"),
        this->mainSettings->getPauseMenuSettingsValue("ButtonTextColorR"),this->mainSettings->getPauseMenuSettingsValue("ButtonTextColorG"),this->mainSettings->getPauseMenuSettingsValue("ButtonTextColorB"),this->mainSettings->getPauseMenuSettingsValue("ButtonTextOpacity")
    );
    Button* resumeButton = new Button(0.5,0.35,0.25,0.125,"Resume");
    Button* settingsButton = new Button(0.5,0.5,0.25,0.125,"Settings");
    Button* quitButton = new Button(0.5,0.65,0.25,0.125,"Main Menu");

    this->pauseMenu->addButton(resumeButton,"Resume Button");
    this->pauseMenu->addButton(settingsButton,"Settings Button");
    this->pauseMenu->addButton(quitButton,"Main Menu Button");
    std::cout<<"PAUSE MENU BUILT"<<std::endl;
}

void EventManager::buildWorld() {
    this->selectedWorld = new World("Demo World");
    /*
     * Need to Move Different Materials to the Material Class Instead
     * Also A Live World Editor Is Deffinitely Needed!
    */
    Material* perfectReflector = new Material(1.0,1.0,1.0,1.0,false);
    Material* basicBlack = new Material(0.0,0.0,0.0,0.0,false);
    Material* basicRed = new Material(1.0,0.0,0.0,0.0,false);
    Material* basicGreen = new Material(0.0,1.0,0.0,0.0,false);
    Material* basicBlue = new Material(0.0,0.0,1.0,0.0,false);
    Material* basicYellow = new Material(1.0,1.0,0.0,0.0,false);
    Material* basicCyan = new Material(0.0,0.6,0.6,0.0,false);


    /*NEED TO CHANGE HOW LIGHT WORKS LIGHT RECIENVED  +=  LIGHT STRENGTH * LIGHT COLOR * COLOR(OF RAY SO FAR)*/
    /*SHOULD JUST REMAKE THE MATERAIL SYSTEM AND MAYBE EVEN THE OBJECT AND WORLD SYSTEM IT BETTER AND MORE ROBUST*/
    Material* basicBlackLight = new Material(0.0,0.0,0.0,0.0,true);
    Material* basicRedLight = new Material(1.0,0.0,0.0,0.0,true);
    Material* basicGreenLight = new Material(0.0,1.0,0.0,0.0,true);
    Material* basicBlueLight = new Material(0.0,0.0,1.0,0.0,true);
    Material* basicYellowLight = new Material(1.0,1.0,0.0,0.0,true);
    Material* basicCyanLight = new Material(0.0,1.0,1.0,0.0,true);

    

    Material* metallicRed = new Material(1.0,0.0,0.0,0.5,false);
    Material* metallicGreen = new Material(0.0,1.0,0.0,0.5,false);
    Material* metallicBlue = new Material(0.0,0.0,1.0,1.0,false);
    Material* metallicWhite = new Material(1.0,1.0,1.0,0.5,false);
    Material* metallicCyan = new Material(0.0,1.0,1.0,0.5,false);

    Material* reflectiveRed = new Material(1.0,0.0,0.0,1.0,false);
    Material* reflectiveWhite = new Material(1.0,1.0,1.0,1.0,false);
    Material* reflectiveCyan = new Material(0.0,1.0,1.0,1.0,false);
    Material* reflectiveGreen = new Material(0.0,1.0,0.0,1.0,false);
    Material* reflectiveBlue = new Material(0.0,0.0,1.0,1.0,false);

    Material* mostlyReflectiveRed = new Material(1.0,0.0,0.0,0.98,false);
    Material* mostlyReflectiveWhite = new Material(1.0,1.0,1.0,0.98,false);
    Material* mostlyReflectiveCyan = new Material(0.0,1.0,1.0,0.98,false);
    Material* mostlyReflectiveGreen = new Material(0.0,1.0,0.0,0.98,false);
    Material* mostlyReflectiveBlue = new Material(0.0,0.0,1.0,0.98,false);


    Material* perfectReflectorRed = new Material(1.0,0.0,0.0,1.0,false);
    Material* perfectReflectorWhite = new Material(1.0,1.0,1.0,1.0,false);
    Material* perfectReflectorCyan = new Material(0.0,1.0,1.0,1.0,false);
    Material* perfectReflectorGreen = new Material(0.0,1.0,0.0,1.0,false);
    Material* perfectReflectorBlue = new Material(0.0,0.0,1.0,1.0,false);

    Material* reflectiveGrey = new Material(0.1603,0.1603,0.1603,1.0,false);

    //Need to Change Material So that there can be a seperate lgith color and lumeninance
    //So that the final light color is correct for the desired brightness of the light(use luminance mapping
    //from gpu code) to make sure that if you want a light with a 2 lumen brightness it is accurate

    float neonLightStrength = 1.5f;
    Material* basicNeonPinkLight = new Material(neonLightStrength*255/255.0f,neonLightStrength*59/255.0f,neonLightStrength*148/255.0f,0.0,true);
    Material* basicNeonPurpleLight = new Material(neonLightStrength*175/255.0f,neonLightStrength*61/255.0f,neonLightStrength*255/255.0f,0.0,true);
    Material* basicNeonCyanLight = new Material(neonLightStrength*85/255.0f,neonLightStrength*255/255.0f,neonLightStrength*225/255.0f,0.0,true);
    Material* basicNeonGreenLight = new Material(neonLightStrength*166/255.0f,neonLightStrength*253/255.0f,neonLightStrength*41/255.0f,0.0,true);

    Material* shinyNeonPinkLight = new Material(neonLightStrength*255/255.0f,neonLightStrength*59/255.0f,neonLightStrength*148/255.0f,0.5,true);
    Material* shinyNeonPurpleLight = new Material(neonLightStrength*175/255.0f,neonLightStrength*61/255.0f,neonLightStrength*255/255.0f,0.5,true);
    Material* shinyNeonCyanLight = new Material(neonLightStrength*85/255.0f,neonLightStrength*255/255.0f,neonLightStrength*225/255.0f,0.5,true);
    Material* shinyNeonGreenLight = new Material(neonLightStrength*166/255.0f,neonLightStrength*253/255.0f,neonLightStrength*41/255.0f,0.5,true);

    Material* reflectiveNeonPinkLight = new Material(neonLightStrength*255/255.0f,neonLightStrength*59/255.0f,neonLightStrength*148/255.0f,1.0,true);
    Material* reflectiveNeonPurpleLight = new Material(neonLightStrength*175/255.0f,neonLightStrength*61/255.0f,neonLightStrength*255/255.0f,1.0,true);
    Material* reflectiveNeonCyanLight = new Material(neonLightStrength*85/255.0f,neonLightStrength*255/255.0f,neonLightStrength*225/255.0f,1.0,true);
    Material* reflectiveNeonGreenLight = new Material(neonLightStrength*166/255.0f,neonLightStrength*253/255.0f,neonLightStrength*41/255.0f,1.0,true);
    Material* reflectiveWhiteLight = new Material(1.0,1.0,1.0,1.0,true);


    Material* basicWhite = new Material(1.0,1.0,1.0,0.0,false);
    Material* basicGrey = new Material(0.1603,0.1603,0.1603,0.0,false);

    float lightStrength = 4.0;
    float lightColorR = 255.0f/255.0f;
    float lightColorG = 244.0f/255.0f;
    float lightColorB = 229.0f/255.0f;

    Material* basicWhiteLight = new Material(lightStrength*lightColorR,lightStrength*lightColorG,lightStrength*lightColorB,0.0,true);

    Material* testWhite1 = new Material(1.0,1.0,1.0,0.02,false);
    Material* testWhite2 = new Material(1.0,1.0,1.0,0.04,false);
    Material* testWhite3 = new Material(1.0,1.0,1.0,0.06,false);
    Material* testWhite4 = new Material(1.0,1.0,1.0,0.9,false);


    Sphere* leftSphere = new Sphere(-3.0,-3.5,5.0,1.5,testWhite1);
    Sphere* rightSphere = new Sphere(3.0,-3.5,5.0,1.5,testWhite2);
    Sphere* centerSphere = new Sphere(0.0,-3.5,2.5,1.5,testWhite3);

    
    /*
     *NEED TO FIX BRIGHTNESS ISSUE FROM ITS ROOT FOR REFLECTIONS TO WORK NICELY
     *
    */

    //8 Verticies of Main Room:
    Vertex* bottomLeftCorner1 = new Vertex(-5.0,-5.0,-7.0);
    Vertex* bottomRightCorner1 = new Vertex(5.0,-5.0,-7.0);
    Vertex* topLeftCorner1 = new Vertex(-5.0,5.0,-7.0);
    Vertex* topRightCorner1 = new Vertex(5.0,5.0,-7.0);
    Vertex* bottomLeftCorner2 = new Vertex(-5.0,-5.0,7.0);
    Vertex* bottomRightCorner2 = new Vertex(5.0,-5.0,7.0);
    Vertex* topLeftCorner2 = new Vertex(-5.0,5.0,7.0);
    Vertex* topRightCorner2 = new Vertex(5.0,5.0,7.0);
    //Triangle Constructed From Verticies:
    Triangle* leftWall1 = new Triangle(bottomLeftCorner1,bottomLeftCorner2,topLeftCorner1,basicRed);
    Triangle* leftWall2 = new Triangle(bottomLeftCorner2,topLeftCorner2,topLeftCorner1,basicRed);
    Triangle* rightWall1 = new Triangle(bottomRightCorner1,bottomRightCorner2,topRightCorner1,basicCyan);
    Triangle* rightWall2 = new Triangle(topRightCorner1,bottomRightCorner2,topRightCorner2,basicCyan);
    Triangle* backWall1 = new Triangle(bottomLeftCorner2,topRightCorner2,bottomRightCorner2,basicGreen);
    Triangle* backWall2 = new Triangle(bottomLeftCorner2,topRightCorner2,topLeftCorner2,basicGreen);
    Triangle* frontWall1 = new Triangle(bottomLeftCorner1,topRightCorner1,bottomRightCorner1,basicWhite);
    Triangle* frontWall2 = new Triangle(bottomLeftCorner1,topRightCorner1,topLeftCorner1,basicWhite);
    Triangle* bottomWall1 = new Triangle(bottomLeftCorner1,bottomRightCorner2,bottomRightCorner1,basicWhite);
    Triangle* bottomWall2 = new Triangle(bottomLeftCorner1,bottomRightCorner2,bottomLeftCorner2,basicWhite);
    Triangle* topWall1 = new Triangle(topLeftCorner1,topRightCorner2,topRightCorner1,basicWhite);
    Triangle* topWall2 = new Triangle(topLeftCorner1,topRightCorner2,topLeftCorner2,basicWhite);
    //Create A simple light in the center of the room:
    Vertex* leftLightCorner1 = new Vertex(-2.0,4.99,-0.5);
    Vertex* leftLightCorner2 = new Vertex(-2.0,4.99,3.5);
    Vertex* rightLightCorner1 = new Vertex(2.0,4.99,-0.5);
    Vertex* rightLightCorner2 = new Vertex(2.0,4.99,3.5);

    Triangle* roomLight1 = new Triangle(leftLightCorner1,rightLightCorner2,rightLightCorner1,basicWhiteLight);
    Triangle* roomLight2 = new Triangle(leftLightCorner1,rightLightCorner2,leftLightCorner2,basicWhiteLight);

    Sphere* lightSphere = new Sphere(0,5.0,2.0,2.0,basicWhiteLight);

    this->selectedWorld->addTriangle(leftWall1);
    this->selectedWorld->addTriangle(leftWall2);
    this->selectedWorld->addTriangle(rightWall1);
    this->selectedWorld->addTriangle(rightWall2);
    this->selectedWorld->addTriangle(backWall1);
    this->selectedWorld->addTriangle(backWall2);
    this->selectedWorld->addTriangle(frontWall1);
    this->selectedWorld->addTriangle(frontWall2);
    this->selectedWorld->addTriangle(bottomWall1);
    this->selectedWorld->addTriangle(bottomWall2);
    this->selectedWorld->addTriangle(topWall1);
    this->selectedWorld->addTriangle(topWall2);


    this->selectedWorld->addSphere(centerSphere);
    this->selectedWorld->addSphere(rightSphere);
    this->selectedWorld->addSphere(leftSphere);

    //Adds a Sphere That Is emmitting Light in the Ceiling
    // this->selectedWorld->addSphere(lightSphere);

    this->selectedWorld->addTriangle(roomLight1);
    this->selectedWorld->addTriangle(roomLight2);


    
    this->mainDisplay->allocateWorldData(this->selectedWorld);
}


#pragma endregion




#pragma region Main Render Calls

void EventManager::runProgram() {
    while(this->programAlive){
        this->handleEvents();
        this->renderScene();
    }
}


# pragma endregion








#pragma region Event Displatcher
void EventManager::handleEvents() {
    /*NEEDS TO CALUCLATE DELTA TIME SINCE LAST CALL FOR USE IN
      STUFF LIKE CAMERA MOVEMENT AND STUFF
    */
    switch(this->programState) {
        case MAIN_MENU:
            this->handleMainMenuEvents();
            break;
        case RUNNING:
            this->handleRunningEvents();
            this->handleKeyboardState();
            break;

        case PAUSED:
            this->handlePausedEvents();
            break;

        default:
            std::cout<<"EVENTS NOT HANDLED FOR CURRENT DISPLAY STATE"<<std::endl;
            break;
    }
}
#pragma endregion




#pragma region  Main Menu Events
void EventManager::handleMainMenuEvents() {
    while(SDL_PollEvent(&this->programEvent)) {
        switch(this->programEvent.type) {
            case SDL_QUIT:
                this->endProgram();
                break;
            
            case SDL_MOUSEBUTTONDOWN:   
                if(this->programEvent.button.button== SDL_BUTTON_RIGHT) {
                    this->mainMenuClicked(this->programEvent.button.x,this->programEvent.button.y,false);
                } else {
                    this->mainMenuClicked(this->programEvent.button.x,this->programEvent.button.y,true);
                }
                break;

            case SDL_WINDOWEVENT:
                if(this->programEvent.window.event != SDL_WINDOWEVENT_RESIZED) {
                    break;
                }
                this->mainDisplay->handleDisplaySizeChange(this->programEvent.window.data1,this->programEvent.window.data2);
                break;
        }
    }
}

/* CURRENTLY NEEDS UPDATING:
 * USES OLD MAIN MENU METHOD
 * TO GET CLICKED BUTTON
 * NEEDS TO USE NEW MENU
 * 
*/


void EventManager::mainMenuClicked(int mousePosX, int mousePosY, bool isLeftClick) {
    float posX = static_cast<float>(mousePosX)/this->mainDisplay->displayWidth;
    float posY = static_cast<float>(mousePosY)/this->mainDisplay->displayHeight;
    std::string buttonClicked = this->mainMenu->getClickedButton(posX,posY);

    switch(isLeftClick) {
        case true:      /*LEFT CLICK*/
            /*CODE IF THE MOSUE CLICK WAS A LEFT CLICK*/
            std::cout<<"LEFT CLICKED: "<<buttonClicked<<std::endl;
            if(buttonClicked == "Start Button") {
                this->updateProgramState(RUNNING);

            } else if(buttonClicked == "Settings Button") {
                std::cout<<"Settings Not Made Yet!"<<std::endl;

            } else if(buttonClicked == "Quit Button") {
                this->endProgram();

            } else if(buttonClicked == "") {
                std::cout<<"No Button Clicked!"<<std::endl;

            } else {
                std::cerr<<"Unknown Button Clicked!"<<std::endl;

            }
            break;


        case false:     /*RIGHT CLICK*/ 
            /*CODE IF THE MOSUE CLICK WAS A RIGHT CLICK*/
            std::cout<<"LEFT CLICKED: "<<buttonClicked<<std::endl;

            break;  

    }

}


#pragma endregion



#pragma region Running Events
void EventManager::handleRunningEvents() {
    while(SDL_PollEvent(&this->programEvent)) {
        switch(this->programEvent.type) {
            case SDL_QUIT:
                exit(true);
                break;
            case SDL_KEYDOWN:
                this->handleRunningKeyPresses(this->programEvent.key.keysym.sym);
                break;

            case SDL_MOUSEMOTION:
                this->mainDisplay->currentFrame = 0;
                if(static_cast<float>(this->programEvent.motion.x)/this->mainDisplay->displayWidth > 0.75 || static_cast<float>(this->programEvent.motion.x)/this->mainDisplay->displayWidth < 0.25 || static_cast<float>(this->programEvent.motion.y)/this->mainDisplay->displayHeight > 0.75 || static_cast<float>(this->programEvent.motion.y)/this->mainDisplay->displayHeight < 0.25) {
                    SDL_WarpMouseInWindow(this->mainDisplay->displayWindow, this->mainDisplay->displayWidth*0.5, this->mainDisplay->displayHeight*0.5);
                }
                this->updateViewCamera(static_cast<float>(this->programEvent.motion.xrel)/this->mainDisplay->displayWidth , static_cast<float>(this->programEvent.motion.yrel)/this->mainDisplay->displayHeight);
                
                break;
            
            case SDL_WINDOWEVENT:
                if(this->programEvent.window.event != SDL_WINDOWEVENT_RESIZED) {
                    break;
                }
                this->mainDisplay->handleDisplaySizeChange(this->programEvent.window.data1,this->programEvent.window.data2);
                break;

            default:
                break;
        }
    }
}


void EventManager::handleRunningKeyPresses(SDL_Keycode keyPressed) {
    switch(keyPressed) {
        case SDLK_ESCAPE:
            /*CODE TO PAUSE PROGRAM*/
            this->updateProgramState(PAUSED);
            SDL_ShowCursor(SDL_ENABLE);
            break;
        case SDLK_r:
            if(this->keyBoardState[SDL_SCANCODE_LCTRL]) {
                std::cout<<"Starting Serious Render of Image!"<<std::endl;
                this->mainDisplay->seriousRenderOfImage(this->numOfFramesPerRender,this->maxRayBounces);
                auto now = std::chrono::system_clock::now();
                std::time_t time = std::chrono::system_clock::to_time_t(now);
                std::tm* tm = std::localtime(&time);
                std::ostringstream oss;
                oss << std::put_time(tm, "%Y.%m.%d_%H.%M.%S");
                this->mainDisplay->saveImage(oss.str());
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
            break;

        case SDLK_UP:
            if(this->keyBoardState[SDL_SCANCODE_LCTRL]) {
                this->viewCamera->incrimentFocalDistance(0.01);
                this->mainDisplay->currentFrame = 0;
            }
            break;
        
        case SDLK_DOWN:
            if(this->keyBoardState[SDL_SCANCODE_LCTRL]) {
                this->viewCamera->decreaseFocalDistance(0.01);
                this->mainDisplay->currentFrame = 0;
            }
            break;
    }
    
}


void EventManager::handleKeyboardState() {
    if(this->keyBoardState[SDL_SCANCODE_W])        {this->viewCamera->moveForward(0.03);    this->mainDisplay->currentFrame = 0;}
    if(this->keyBoardState[SDL_SCANCODE_A])        {this->viewCamera->moveLeft(0.03);       this->mainDisplay->currentFrame = 0;}
    if(this->keyBoardState[SDL_SCANCODE_S])        {this->viewCamera->moveBackward(0.03);   this->mainDisplay->currentFrame = 0;}
    if(this->keyBoardState[SDL_SCANCODE_D])        {this->viewCamera->moveRight(0.03);      this->mainDisplay->currentFrame = 0;}
    if(this->keyBoardState[SDL_SCANCODE_SPACE])    {this->viewCamera->moveUp(0.03);         this->mainDisplay->currentFrame = 0;}
    if(this->keyBoardState[SDL_SCANCODE_LSHIFT])   {this->viewCamera->moveDown(0.03);       this->mainDisplay->currentFrame = 0;}

    
}

void EventManager::updateViewCamera(float deltaX, float deltaY) {
    /* NEEDS TO WORK ON SETTINGS FIRST AND MAYBE CONSIDER GETTING EVERYTHING WORKING FIRST
     * BUT EITHER WAY THIS IS CURRENTLY INCOMPLETE
    */
    this->viewCamera->rotateCamera(deltaX * this->mouseSensitivityX, deltaY * this->mouseSensitivityY);
}


#pragma endregion








#pragma region Paused Menu Events


void EventManager::handlePausedEvents() {
    while(SDL_PollEvent(&this->programEvent)) {
        switch(this->programEvent.type) {
            case SDL_QUIT:
                this->endProgram();
                break;

            case SDL_MOUSEBUTTONDOWN:   
                if(this->programEvent.button.button== SDL_BUTTON_RIGHT) {
                    this->pauseMenuClicked(this->programEvent.button.x,this->programEvent.button.y,false);
                } else {
                    this->pauseMenuClicked(this->programEvent.button.x,this->programEvent.button.y,true);
                }
                break;

            case SDL_WINDOWEVENT:
                if(this->programEvent.window.event != SDL_WINDOWEVENT_RESIZED) {
                    break;
                }
                this->mainDisplay->handleDisplaySizeChange(this->programEvent.window.data1,this->programEvent.window.data2);
                break;

            default:
                break;
        }
    }
}


void EventManager::pauseMenuClicked(int mousePosX, int mousePosY, bool isLeftClick) {
    float posX = static_cast<float>(mousePosX)/this->mainDisplay->displayWidth;
    float posY = static_cast<float>(mousePosY)/this->mainDisplay->displayHeight;
    std::string buttonClicked = this->pauseMenu->getClickedButton(posX,posY);

    switch(isLeftClick) {
        case true:      /*LEFT CLICK*/
            /*CODE IF THE MOSUE CLICK WAS A LEFT CLICK*/
            std::cout<<"LEFT CLICKED: "<<buttonClicked<<std::endl;
            if(buttonClicked == "Resume Button") {
                this->updateProgramState(RUNNING);

            } else if(buttonClicked == "Settings Button") {
                std::cout<<"Settings Not Made!"<<std::endl;

            } else if(buttonClicked == "Main Menu Button") {
                this->updateProgramState(MAIN_MENU);

            } else if(buttonClicked == "") {
                std::cout<<"No Button Clicked!"<<std::endl;

            } else {
                std::cerr<<"Unknown Button Clicked!"<<std::endl;

            }
            break;


        case false:     /*RIGHT CLICK*/ 
            /*CODE IF THE MOSUE CLICK WAS A RIGHT CLICK*/
            std::cout<<"LEFT CLICKED: "<<buttonClicked<<std::endl;

            break;  

    }

}



#pragma endregion







#pragma region Render Calls


void EventManager::renderScene() {
    switch(this->programState) {
        /*CASE TO RENDER MAIN MENU*/
        case MAIN_MENU:
            SDL_RenderClear(this->mainProgramRenderer);
            this->mainMenu->renderMenu(this->mainDisplay->displayWidth,this->mainDisplay->displayHeight);
            SDL_RenderPresent(this->mainProgramRenderer);
            break;

        /*CASE TO RENDER RAYTRACING PART*/
        case RUNNING: 
            this->mainDisplay->renderScene();
            break;

        /*CASE TO RENDER PAUSED MENU*/
        case PAUSED:
            SDL_RenderClear(this->mainProgramRenderer);
            if(this->mainDisplay->displayWidth!=this->expectedImageWidth || this->mainDisplay->displayHeight!=this->expectedImageHeight) {
                std::cout<<"Updating Image!"<<std::endl;
                this->mainDisplay->renderImage();
                expectedImageWidth = this->mainDisplay->displayWidth;
                expectedImageHeight = this->mainDisplay->displayHeight;
            }
            /*NEED CODE HERE TO RENDER BACKGROUND IMAGE NOW*/
            SDL_UpdateTexture(this->mainDisplay->imageTexture, NULL, this->mainDisplay->imageData, this->mainDisplay->displayWidth*3);
            SDL_RenderCopy(this->mainProgramRenderer, this->mainDisplay->imageTexture, NULL, NULL);
            /* THIS SHOUDL WORK ^^^ HOWEVER IT CALLS THROUGH THE CLASS A LOT SO BE CAREFUL
             * AND CONSIDER MAKING A DISPLAY FUNCTION*/
            this->pauseMenu->renderMenu(this->mainDisplay->displayWidth,this->mainDisplay->displayHeight);
            SDL_RenderPresent(this->mainProgramRenderer);
            break;

        /*CASE TO RENDER SETTINGS MENU*/
        case SETTINGS:

            break;
    }
}



#pragma endregion
/*NEXT STEP IS TO HANDLE THE SWITCHING BETWEEN STATES*/









#pragma region General

void EventManager::updateProgramState(int newProgramState) {
    switch(newProgramState) {
        case MAIN_MENU:
            SDL_ShowCursor(SDL_ENABLE);
            this->previouseProgramState = this->programState;
            this->programState = MAIN_MENU;
            break;
        case RUNNING:
            SDL_ShowCursor(SDL_DISABLE);
            this->previouseProgramState = this->programState;
            this->programState = RUNNING;
            break;
        case PAUSED:
            SDL_ShowCursor(SDL_ENABLE);
            this->previouseProgramState = this->programState;
            this->programState = PAUSED;
            break;
            
        case -1:
            if(this->previouseProgramState >=0) {
                this->updateProgramState(this->previouseProgramState);
            } else {
                this->updateProgramState(MAIN_MENU);
            }

    }
}


void EventManager::endProgram() {
    std::cout<<"Program Terminating!"<<std::endl;
    delete this->mainDisplay;
    delete this->viewCamera;
    delete this->mainMenu;
    delete this->pauseMenu;
    delete this->mainSettings;
    delete this->mainFileParser;

    if(this->mainProgramRenderer)
        SDL_DestroyRenderer(this->mainProgramRenderer);
    this->programAlive = false;
}



#pragma endregion





