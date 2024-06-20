#include "../HeaderFiles/Menu.h"






Menu::Menu(std::string menuName, int imageWidth, int imageHeight, SDL_Renderer* menuRenderer, std::string fontPath) : 
menuName(menuName),expectedImageWidth(0), expectedImageHeight(0),menuRenderer(menuRenderer),fontPath(fontPath) {
    this->menuFont = TTF_OpenFont(this->fontPath.c_str(), 24);
    std::cout<<"FONT FOUND AND MADE"<<std::endl;
    this->menuFontSize = 24;
    // this->updateMenuSize(this->expectedImageWidth,this->expectedImageHeight);
}

Menu::~Menu() {
    std::cout<<"Menu Deconstructing"<<std::endl;

    for(std::pair<std::string,Button*> menuButton: this->menuButtons) {
        delete menuButton.second;
    }
    if(this->menuFont) 
        TTF_CloseFont(this->menuFont);
    if(this->menuNameTextTexture)
        SDL_DestroyTexture(this->menuNameTextTexture);
}




void Menu::setupMenuColors(
uint8_t backgroundColorR,uint8_t backgroundColorG,uint8_t backgroundColorB,uint8_t backgroundColorA,
uint8_t defaultButtonColorR,uint8_t defaultButtonColorG,uint8_t defaultButtonColorB,uint8_t defaultButtonColorA,
uint8_t selectedButtonColorR,uint8_t selectedButtonColorG,uint8_t selectedButtonColorB,uint8_t selectedButtonColorA,
uint8_t buttonTextColorR,uint8_t buttonTextColorG,uint8_t buttonTextColorB,uint8_t buttonTextColorA) {

    this->backgroundColor = {backgroundColorR,backgroundColorG,backgroundColorB,backgroundColorA};
    this->defaultButtonColor = {defaultButtonColorR,defaultButtonColorG,defaultButtonColorB,defaultButtonColorA};
    this->selectedButtonColor = {selectedButtonColorR,selectedButtonColorG,selectedButtonColorB,selectedButtonColorA};
    this->buttonTextColor = {buttonTextColorR,buttonTextColorG,buttonTextColorB,buttonTextColorA};
}


void Menu::addButton(Button* newButton, std::string buttonName) {
    newButton->linkButtonToMenu(this);
    std::pair<std::string,Button*> newButtonPair = {buttonName,newButton};
    this->menuButtons.push_back(newButtonPair);
}




void Menu::updateMenuSize(int newImageWidth, int newImageHeight) {
    for(std::pair<std::string,Button*>& menuButton : this->menuButtons) {
        menuButton.second->updateButtonSize(newImageWidth,newImageHeight);
    }
    this->backgroundRectangle = {0,0,newImageWidth,newImageHeight};
    this->updateFontSize(static_cast<int>(newImageWidth*this->menuNameTextWidth),static_cast<int>(newImageWidth*this->menuNameTextHeight));
    
    SDL_Surface* menuNameTextSurface = TTF_RenderText_Solid(this->menuFont, this->menuName.c_str(), this->buttonTextColor);
    
    this->menuNameTextTexture = SDL_CreateTextureFromSurface(this->menuRenderer,menuNameTextSurface);
    int textWidth = menuNameTextSurface->w;
    int textHeight = menuNameTextSurface->h;
    int textCornerX = 0.5*newImageWidth - 0.5*textWidth;
    int textCornerY = 0.15*newImageHeight - 0.5*textHeight;
    this->menuNameTextRectangle = {textCornerX,textCornerY,textWidth,textHeight};
}




void Menu::updateFontSize(int mainTextWidth, int mainTextHeight) {\
    SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->menuFont, this->menuName.c_str(),{255,255,255});
    
    int textWidth = tempTextSurface->w;
    int textHeight = tempTextSurface->h;
    bool switchCase = textWidth < mainTextWidth && textHeight < mainTextHeight; 
    switch(switchCase) {
        case true:
            do {
                //INCRIMENT FONT SIZE UNTIL MAXIMUM
                this->incrimentFontSize();
                SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->menuFont, this->menuName.c_str(),{255,255,255});
                textWidth = tempTextSurface->w;
                textHeight = tempTextSurface->h;    
            } while(textWidth < mainTextWidth && textHeight < mainTextHeight);
                this->menuFontSize -=3;
            break;

        default:
            do {
                //DECREMENT FONT SIZE UNTIL SIZE FITS
                this->decrementFontSize();
                SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->menuFont, this->menuName.c_str(),{255,255,255});
                textWidth = tempTextSurface->w;
                textHeight = tempTextSurface->h;  
            } while(textWidth > mainTextWidth || textHeight > mainTextHeight);
                this->menuFontSize -= 2;
            break;
    }
    TTF_CloseFont(this->menuFont);
    this->menuFont = TTF_OpenFont(this->fontPath.c_str(), this->menuFontSize);
}



void Menu::incrimentFontSize() {
    this->menuFontSize += 1;
    TTF_CloseFont(this->menuFont);
    this->menuFont = TTF_OpenFont(this->fontPath.c_str(), this->menuFontSize);
}

void Menu::decrementFontSize() {
    this->menuFontSize -= 1;
    TTF_CloseFont(this->menuFont);
    this->menuFont = TTF_OpenFont(this->fontPath.c_str(), this->menuFontSize);
}




void Menu::renderMenu(int imageWidth, int imageHeight) {
    if(imageWidth != this->expectedImageWidth || imageHeight != this->expectedImageHeight) {
        this->updateMenuSize(imageWidth,imageHeight);
    }
    int tempMouseX, tempMouseY;
    float mousePosX, mousePosY;
    SDL_GetMouseState(&tempMouseX,&tempMouseY);
    mousePosX = static_cast<float>(tempMouseX)/imageWidth;
    mousePosY = static_cast<float>(tempMouseY)/imageHeight;

    /*RENDERS MENU BACKGROUND*/
    SDL_SetRenderDrawColor(this->menuRenderer,this->backgroundColor.r,this->backgroundColor.g,this->backgroundColor.b,this->backgroundColor.a);
    SDL_RenderFillRect(this->menuRenderer,&this->backgroundRectangle);
    /*RENDERS MENU NAME TEXT*/
    SDL_RenderCopy(this->menuRenderer,this->menuNameTextTexture, NULL, &this->menuNameTextRectangle);

    /*RENDER BUTTONS*/
    for(std::pair<std::string,Button*>& menuButton : this->menuButtons) {
        menuButton.second->renderButton(mousePosX, mousePosY);
    }
}

std::string Menu::getClickedButton(float posX, float posY) {
    for(std::pair<std::string,Button*> menuButton : this->menuButtons) {
        if(menuButton.second->isClicked(posX,posY)) {
            std::cout<<"A BUTTON WAS CLICKED2!"<<std::endl;
            return menuButton.first;
        }
    }
    return "";
}