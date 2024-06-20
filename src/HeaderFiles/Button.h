#ifndef Button_H
#define Button_H

#include <string.h>
#include <iostream>
#include "../../lib/SDL/include/SDL.h"
#include "../../lib/SDL/include/SDL_ttf.h"
#include "./Menu.h"

class Menu;

class Button {

public:
    /*
    Values are relatives compared to a screen that goes from 0 to 1 in both x and y
    For the widht and height a value of 1 would be the entire screen
    */
    Button(float posX, float posY, float widthScale, float heightScale, std::string buttonText);
    ~Button();
    void linkButtonToMenu(Menu* linkedMenu);
    void renderButton(float mousePosX, float mousePosY);
    void updateButtonSize(int imageWidth, int imageHeight);
    bool isClicked(float mousePosX, float mousePosY);

protected:
    bool isMouseOver(float mousePosX, float mousePosY);

private:
    /*
    Values are relatives compared to a screen that goes from 0 to 1 in both x and y
    For the widht and height a value of 1 would be the entire screen
    */
    float posX;
    float posY;
    float widthScale;
    float heightScale;
    std::string buttonText;
    /*
    Values are out of 255
    */
    int colorR;
    int colorG;
    int colorB;
    int opacity;
    /*
    Basic Font that just uses whatever is avalilable on
    the system and is initialized in the constructor:
    */
    TTF_Font* font;
    uint8_t fontColorR;
    uint8_t fontColorG;
    uint8_t fontColorB;
    /*
    Function to calculate proper font size
    */
    void setFontSize(int buttonWidth, int buttonHeight);
    void incrimentFontSize();
    void decrementFontSize();
    int fontSize;
    /*
    SDL COMPONENTS USED TO RENDER BUTTON
    */
    int expectedImageWidth;
    int expectedImageHeight;
    SDL_Rect buttonBackgroundRectangle;
    SDL_Texture* buttonTextTexture;
    SDL_Rect buttonTextRectangle;

    SDL_Renderer* buttonRenderer;
    SDL_Color defaultButtonColor;
    SDL_Color selectedButtonColor;
    SDL_Color buttonTextColor;

    friend class Menu;



};





#endif