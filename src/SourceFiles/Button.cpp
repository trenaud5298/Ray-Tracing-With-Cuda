#include "../HeaderFiles/Button.h"
#include <filesystem>

Button::Button(float posX, float posY, float widthScale, float heightScale, std::string buttonText) :
posX(posX),posY(posY),widthScale(widthScale),heightScale(heightScale),buttonText(buttonText) {

};


void Button::linkButtonToMenu(Menu* linkedMenu) {
    this->buttonRenderer = linkedMenu->menuRenderer;
    this->buttonTextColor = linkedMenu->buttonTextColor;
    this->defaultButtonColor = linkedMenu->defaultButtonColor;
    this->selectedButtonColor = linkedMenu->selectedButtonColor;
    this->font = TTF_OpenFont(linkedMenu->fontPath.c_str(),24);
    this->fontSize = 24;
}

Button::~Button() {
    std::cout<<"Button Deconstructing"<<std::endl;
    if(this->font)
        TTF_CloseFont(this->font);
    if(this->buttonTextTexture)
        SDL_DestroyTexture(this->buttonTextTexture);
}


void Button::renderButton(float mousePosX, float mousePosY) {
    /*
    Render Button Background
    */
    if(this->isMouseOver(mousePosX, mousePosY)) {
        SDL_SetRenderDrawColor(this->buttonRenderer,this->selectedButtonColor.r,this->selectedButtonColor.g,this->selectedButtonColor.b,this->selectedButtonColor.a);
    } else {
        SDL_SetRenderDrawColor(this->buttonRenderer,this->defaultButtonColor.r,this->defaultButtonColor.g,this->defaultButtonColor.b,this->defaultButtonColor.a);
    }

    SDL_RenderFillRect(this->buttonRenderer,&this->buttonBackgroundRectangle);

    /*
    Render Button Text
    */
    SDL_RenderCopy(this->buttonRenderer,this->buttonTextTexture, NULL, &this->buttonTextRectangle);

    return;
}



void Button::updateButtonSize(int imageWidth, int imageHeight) {
    int buttonWidth, buttonHeight, buttonCornerX, buttonCornerY;
    buttonWidth = this->widthScale * imageWidth;
    buttonHeight = this->heightScale * imageHeight;
    buttonCornerX = this->posX*imageWidth-0.5*buttonWidth;
    buttonCornerY = this->posY*imageHeight-0.5*buttonHeight;
    this->buttonBackgroundRectangle =  {buttonCornerX,buttonCornerY,buttonWidth,buttonHeight};

    this->setFontSize(buttonWidth,buttonHeight);
    SDL_Surface* buttonTextSurface = TTF_RenderText_Solid(this->font, this->buttonText.c_str(), this->buttonTextColor);
    this->buttonTextTexture = SDL_CreateTextureFromSurface(this->buttonRenderer,buttonTextSurface);
    int textWidth = buttonTextSurface->w;
    int textHeight = buttonTextSurface->h;
    int textCornerX = this->posX*imageWidth - 0.5*textWidth;
    int textCornerY = this->posY*imageHeight - 0.5*textHeight;
    this->buttonTextRectangle = {textCornerX,textCornerY,textWidth,textHeight};
}


void Button::setFontSize(int buttonWidth, int buttonHeight) {
    SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->font, this->buttonText.c_str(),{255,255,255});
    int textWidth = tempTextSurface->w;
    int textHeight = tempTextSurface->h;
    bool switchCase = textWidth < buttonWidth && textHeight < buttonHeight; 
    switch(switchCase) {
        case true:
            do {
                //INCRIMENT FONT SIZE UNTIL MAXIMUM
                this->incrimentFontSize();
                SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->font, this->buttonText.c_str(),{255,255,255});
                textWidth = tempTextSurface->w;
                textHeight = tempTextSurface->h;    
            } while(textWidth < buttonWidth && textHeight < buttonHeight);
                this->fontSize -=3;
            break;

        default:
            do {
                //DECREMENT FONT SIZE UNTIL SIZE FITS
                this->decrementFontSize();
                SDL_Surface* tempTextSurface = TTF_RenderText_Solid(this->font, this->buttonText.c_str(),{255,255,255});
                textWidth = tempTextSurface->w;
                textHeight = tempTextSurface->h;  
            } while(textWidth > buttonWidth || textHeight > buttonHeight);
                this->fontSize -= 2;
            break;
    }
    TTF_CloseFont(this->font);
    this->font = TTF_OpenFont("./Fonts/Default.ttf", this->fontSize);
}

void Button::incrimentFontSize() {
    this->fontSize += 1;
    TTF_CloseFont(this->font);
    this->font = TTF_OpenFont("./Fonts/Default.ttf", this->fontSize);

}

void Button::decrementFontSize() {
    this->fontSize -= 1;
    TTF_CloseFont(this->font);
    this->font = TTF_OpenFont("./Fonts/Default.ttf", this->fontSize);
}

bool Button::isClicked(float mousePosX, float mousePosY) {
    float halfWidth = this->widthScale/2;
    float halfHeight = this->heightScale/2;
   
    bool inBoundsX = mousePosX < (this->posX + halfWidth) && mousePosX > (this->posX - halfWidth);
    bool inBoundsY = mousePosY < (this->posY + halfHeight) && mousePosY > (this->posY - halfHeight);
    if(inBoundsX && inBoundsY) {
        return true;
    }

    return false;
}

bool Button::isMouseOver(float mousePosX, float mousePosY) {
    float halfWidth = this->widthScale/2;
    float halfHeight = this->heightScale/2;

    bool inBoundsX = mousePosX < (this->posX + halfWidth) && mousePosX > (this->posX - halfWidth);
    bool inBoundsY = mousePosY < (this->posY + halfHeight) && mousePosY > (this->posY - halfHeight);
    if(inBoundsX && inBoundsY) {
        return true;
    }

    return false;
}