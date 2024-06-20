#include "../HeaderFiles/MainMenu.h"

#define NUMBER_OF_BUTTONS 3

MainMenu::MainMenu(int backgroundColorR,int backgroundColorG,int backgroundColorB) :
backgroundColorR(backgroundColorR),backgroundColorG(backgroundColorG),backgroundColorB(backgroundColorB) {
    this->menuButtons = new Button*[NUMBER_OF_BUTTONS];
    this->menuButtons[0] = new Button(0.5,0.3,0.2,0.15,"Start",255,255,255,150,255,0,0);
    this->menuButtons[1] = new Button(0.5,0.5,0.2,0.15,"Edit",255,255,255,150,255,0,0);
    this->menuButtons[2] = new Button(0.5,0.7,0.2,0.15,"Quit",255,255,255,150,255,0,0);
}

void MainMenu::renderMenu(SDL_Renderer* buttonRenderer, int imageWidth, int imageHeight, float posX, float posY) {
    this->backgroundRectangle = {0,0,imageWidth,imageHeight};
    SDL_SetRenderDrawColor(buttonRenderer,this->backgroundColorR,this->backgroundColorG,this->backgroundColorB,255);
    SDL_RenderFillRect(buttonRenderer, &this->backgroundRectangle);
    for(int i = 0; i < NUMBER_OF_BUTTONS; i ++) {
        /* RATHER THAN BEING PASSED IN THE POSX IT SHOULD BE CALCULATED
         * FROM DISPLAY SIZE AS WELL AS USING SDL TO GET CURRENT MOUSE POS
         * SO MOUSE POS DOES NOT ALWAYS NEED TO BE UPDATED AND STUFF
         * 
         * COULD ALSO BE RETRIEVED IN THE EVENT MANAGER AREA AND THEN PASSED LIKE NOW
         * JUST NO LONGER STORE THE MOUSE POS IN AN INSTANCE VARIABLE
        */
        this->menuButtons[i]->renderButton(buttonRenderer,imageWidth,imageHeight, posX, posY);
    }
}

int MainMenu::checkClickedButton(float posX, float posY) {
    for(int i = 0; i < NUMBER_OF_BUTTONS; i++) {
        if(this->menuButtons[i]->isClicked(posX,posY)) {
            return i;
        }
    }
    return -1;
}