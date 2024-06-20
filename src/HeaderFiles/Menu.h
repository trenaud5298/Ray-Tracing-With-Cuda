#ifndef MENU_H
#define MENU_H

#include "./Button.h"
#include <vector>
#include "../../lib/SDL/include/SDL.h"
#include "../../lib/SDL/include/SDL_ttf.h"

class Button;

class Menu {

public:
    Menu(std::string menuName, int imageWidth, int imageHeight, SDL_Renderer* menuRenderer, std::string fontPath);
    ~Menu();
    std::vector< std::pair<std::string,Button*> > menuButtons;
    void setupMenuColors(
        uint8_t backgroundColorR,uint8_t backgroundColorG,uint8_t backgroundColorB,uint8_t backgroundColorA,
        uint8_t defaultButtonColorR,uint8_t defaultButtonColorG,uint8_t defaultButtonColorB,uint8_t defaultButtonColorA,
        uint8_t selectedButtonColorR,uint8_t selectedButtonColorG,uint8_t selectedButtonColorB,uint8_t selectedButtonColorA,
        uint8_t buttonTextColorR,uint8_t buttonTextColorG,uint8_t buttonTextColorB,uint8_t buttonTextColorA);
    void addButton(Button* newButton, std::string buttonName);
    void renderMenu(int imageWidth, int imageHeight);
    void updateMenuSize(int newImageWidth, int newImageHeight);
    std::string getClickedButton(float posX, float posY);

protected:


private:
    std::string menuName;
    float menuNameTextWidth = 0.2;
    float menuNameTextHeight = 0.1;
    SDL_Texture* menuNameTextTexture;
    SDL_Rect menuNameTextRectangle;

    int expectedImageWidth;
    int expectedImageHeight;

    SDL_Color backgroundColor;
    SDL_Color defaultButtonColor;
    SDL_Color selectedButtonColor;
    SDL_Color buttonTextColor;

    SDL_Renderer* menuRenderer;
    SDL_Rect mainTextRectangle;
    SDL_Rect backgroundRectangle;
    TTF_Font* menuFont;
    std::string fontPath;
    int menuFontSize;

    void updateFontSize(int mainTextWidth, int mainTextHeight);
    void incrimentFontSize();
    void decrementFontSize();
    friend class Button;
};





#endif