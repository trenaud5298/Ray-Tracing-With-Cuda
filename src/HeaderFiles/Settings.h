#ifndef SETTINGS_H
#define SETTINGS_H

#include "../../lib/Other/SimpleIni.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>

class Settings {

public:
    Settings(std::string iniFile,std::string settingsFolderPath);
    ~Settings();
    std::vector< std::pair<std::string, int> > mainSettings = {
        {"HorizontalCameraSensitivity",20},
        {"VerticalCameraSensitivity",20},
        {"CameraStartPosX",0},
        {"CameraStartPosY",-2},
        {"CameraStartPosZ",-5},
        {"MainSettingsBlank4",0},
        {"MainSettingsBlank5",0}
    };


    std::vector< std::pair<std::string, int> > graphicsSettings = {
        {"ScreenResolutionX", 1200},
        {"ScreenResolutionY", 800},
        {"FieldOfView", 90},
        {"FramesPerRender", 1000},
        {"MaxRayBounces", 20},
        {"FocalDistance", 5},
        {"DepthOfFieldStrength", 0},
        {"ToneMapping",0},
        {"GraphicsSettingsBlank2",0},
        {"GraphicsSettingsBlank3",0},
        {"GraphicsSettingsBlank4",0},
        {"GraphicsSettingsBlank5",0},
        {"GraphicsSettingsBlank6",0}
    };

    std::vector< std::pair<std::string, int> > mainMenuSettings = {
        {"BackgroundColorR", 3},
        {"BackgroundColorG", 168},
        {"BackgroundColorB", 158},
        {"BackgroundOpacity", 255},
        {"DefaultButtonColorR",3},
        {"DefaultButtonColorG",168},
        {"DefaultButtonColorB",158},
        {"DefaultButtonOpacity",255},
        {"SelectedButtonColorR",0},
        {"SelectedButtonColorG",104},
        {"SelectedButtonColorB",139},
        {"SelectedButtonOpacity",255},
        {"ButtonTextColorR",255},
        {"ButtonTextColorG",255},
        {"ButtonTextColorB",255},
        {"ButtonTextOpacity",255}
    };

    std::vector< std::pair<std::string, int> > pauseMenuSettings = {
        {"BackgroundColorR", 255},
        {"BackgroundColorG", 255},
        {"BackgroundColorB", 255},
        {"BackgroundOpacity", 100},
        {"DefaultButtonColorR",3},
        {"DefaultButtonColorG",168},
        {"DefaultButtonColorB",158},
        {"DefaultButtonOpacity",150},
        {"SelectedButtonColorR",0},
        {"SelectedButtonColorG",104},
        {"SelectedButtonColorB",139},
        {"SelectedButtonOpacity",230},
        {"ButtonTextColorR",255},
        {"ButtonTextColorG",255},
        {"ButtonTextColorB",255},
        {"ButtonTextOpacity",200}
    };

    std::vector< std::pair<std::string, int> > settingsMenuSettings = {
        {"BackgroundColorR", 255},
        {"BackgroundColorG", 255},
        {"BackgroundColorB", 255},
        {"BackgroundOpacity", 255},
        {"SettingsMenuSettingsBlank1",0},
        {"SettingsMenuSettingsBlank2",0},
        {"SettingsMenuSettingsBlank3",0},
        {"SettingsMenuSettingsBlank4",0},
        {"SettingsMenuSettingsBlank5",0},
        {"SettingsMenuSettingsBlank6",0}
    };

    

    std::map<std::string, int> settingsMenusIndeces = {
        {"Main Settings",0},
        {"Graphics Settings",1},
        {"Main Menu Settings",2},
        {"Pause Menu Settings",3},
        {"Settings Menu Settings",4}
    };

    void printAllSettings();
    void saveIniFile();
    int getMainSettingsValue(std::string settingName);
    int getGraphicsSettingsValue(std::string settingName);
    int getMainMenuSettingsValue(std::string settingName);
    int getPauseMenuSettingsValue(std::string settingName);
    int getSettingsMenuSettingsValue(std::string settingName);


protected:




private:
    /*FOR SETTINGS MENU:
     * 0=MAIN SETTINGS
     * 1=GRAPHICS SETTINGS
     * 2=?
     * 3=?
     * 4=?
     * 5=?
    */
    void updateSetting(int settingsMenu, std::string setting, int newValue);
    std::string iniFile;
    std::string settingsFolderPath;


};


#endif