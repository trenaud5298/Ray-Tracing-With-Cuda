#include "../HeaderFiles/Settings.h"

#define MAIN_SETTINGS 0
#define GRAPHICS_SETTINGS 1
#define MAIN_MENU_SETTINGS 2
#define PAUSE_MENU_SETTINGS 3
#define SETTINGS_MENU_SETTINGS 4

Settings::Settings(std::string iniFile, std::string settingsFolderPath):iniFile(iniFile),settingsFolderPath(settingsFolderPath) {
    CSimpleIniA ini;
    if(iniFile != "") {
        /*INI FILE WAS FOUND*/
        ini.LoadFile(iniFile.c_str());
        std::cout<<"INI FILE FOUND"<<std::endl;
        /*UPDATE ALL SETTINGS TO MATCH INI FILE*/
        for(std::pair<std::string,int> settingsMenu : this->settingsMenusIndeces) {
            std::cout<<"Updating Section: "<<settingsMenu.first<<std::endl;
            std::list<CSimpleIniA::Entry> sectionKeys;
            if(ini.GetAllKeys(settingsMenu.first.c_str(), sectionKeys)) {
                for(const auto& key : sectionKeys) {
                    this->updateSetting(settingsMenu.second,key.pItem,ini.GetLongValue(settingsMenu.first.c_str(),key.pItem));
                }
            } else {
                std::cerr<<"FAILED TO RETRIEVE SECTION "<<settingsMenu.first<<"'s keys"<<std::endl;
            }
        }
        this->saveIniFile();

        
    } else {
        /*INI FILE WAS NOT FOUND*/
        std::cout<<"ERROR: INI FILE NOT FOUND; ATTEMPTING TO CREATE INI FILE"<<std::endl;
        this->saveIniFile();
    }

    std::cout<<"VAL: "<<this->graphicsSettings[6].second<<std::endl;
}

Settings::~Settings() {
    std::cout<<"Settings Saving And Deconstructing"<<std::endl;

    this->saveIniFile();
}


void Settings::saveIniFile() {
    CSimpleIniA ini;
    for(const auto& mainSetting : this->mainSettings) {
        ini.SetLongValue("Main Settings",mainSetting.first.c_str(),mainSetting.second);
    }
    for(const auto& graphicsSetting : this->graphicsSettings) {
        ini.SetLongValue("Graphics Settings",graphicsSetting.first.c_str(),graphicsSetting.second);
    }
    for(const auto& mainMenuSetting : this->mainMenuSettings) {
        ini.SetLongValue("Main Menu Settings",mainMenuSetting.first.c_str(),mainMenuSetting.second);
    }
    for(const auto& pauseMenuSetting : this->pauseMenuSettings) {
        ini.SetLongValue("Pause Menu Settings",pauseMenuSetting.first.c_str(),pauseMenuSetting.second);
    }
    for(const auto& settingsMenuSetting : this->settingsMenuSettings) {
        ini.SetLongValue("Settings Menu Settings",settingsMenuSetting.first.c_str(),settingsMenuSetting.second);
    }

    if(this->settingsFolderPath!="") {
        std::string iniFileSavePath;
        if(this->settingsFolderPath[0] != '.') {
            iniFileSavePath = this->settingsFolderPath + "\\Settings.ini";
        } else {
            iniFileSavePath = this->settingsFolderPath + "/Settings.ini";
        }
        std::cout<<"TRYING TO SAVE INI FILE AT PATH: "<<iniFileSavePath<<std::endl;
        ini.SaveFile(iniFileSavePath.c_str());
    } else {
        std::cout<<"ERROR: SETTINGS FOLDER NOT FOUND; SAVING INI FILE IN CURRENT WORKING DIRECTORY"<<std::endl;
        ini.SaveFile("Settings.ini");
    }
}

void Settings::updateSetting(int settingsMenu, std::string setting, int newValue) {
    switch(settingsMenu) {
        case MAIN_SETTINGS:
            for(std::pair<std::string,int>& mainSetting : this->mainSettings) {
                if(mainSetting.first == setting) {
                    mainSetting.second = newValue;
                }
            }
            break;
        
        case GRAPHICS_SETTINGS:
            for(std::pair<std::string,int>& graphicsSetting : this->graphicsSettings) {
                if(graphicsSetting.first == setting) {
                    graphicsSetting.second = newValue;
                }
            }
            break;

        case MAIN_MENU_SETTINGS:
            for(std::pair<std::string,int>& mainMenuSetting : this->mainMenuSettings) {
                if(mainMenuSetting.first == setting) {
                    mainMenuSetting.second = newValue;
                }
            }
            break;

        case PAUSE_MENU_SETTINGS:
            for(std::pair<std::string,int>& pauseMenuSetting : this->pauseMenuSettings) {
                if(pauseMenuSetting.first == setting) {
                    pauseMenuSetting.second = newValue;
                }
            }
            break;
        
        case SETTINGS_MENU_SETTINGS:
            for(std::pair<std::string,int>& settingsMenuSetting : this->settingsMenuSettings) {
                if(settingsMenuSetting.first == setting) {
                    settingsMenuSetting.second = newValue;
                }
            }
            break;
        
        default:
            std::cerr<<"ERROR: UNKNOW SETTINGS MENU"<<std::endl;
            break;

    }
}

void Settings::printAllSettings() {
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"Main Settings: "<<std::endl;
    for(std::pair<std::string,int> mainSetting : this->mainSettings) {
        std::cout<<"--"<<mainSetting.first<<": "<<mainSetting.second<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;


    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"Graphics Settings: "<<std::endl;
    for(std::pair<std::string,int> graphicsSetting : this->graphicsSettings) {
        std::cout<<"--"<<graphicsSetting.first<<": "<<graphicsSetting.second<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;


}


int Settings::getMainSettingsValue(std::string settingName) {
    /*CHECKS THROUGH MAIN SETTINGS*/
    for(std::pair<std::string,int> mainSetting : this->mainSettings) {
        if(mainSetting.first == settingName) {
            return mainSetting.second;
        }
    }

    /*IN CASE SETTING WAS NOT FOUND*/
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    std::cerr<<"ERROR SETTING NOT FOUND; RETURN A VALUE OF 0"<<std::endl;
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    return 0;
}

int Settings::getGraphicsSettingsValue(std::string settingName) {
    /*CHECKS THROUGH GRAPHICS SETTINGS*/
    for(std::pair<std::string,int> graphicsSetting : this->graphicsSettings) {
        if(graphicsSetting.first == settingName) {
            return graphicsSetting.second;
        }
    }

    /*IN CASE SETTING WAS NOT FOUND*/
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    std::cerr<<"ERROR SETTING NOT FOUND; RETURN A VALUE OF 0"<<std::endl;
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    return 0;
}

int Settings::getMainMenuSettingsValue(std::string settingName) {
    /*CHECKS THROUGH MAIN MENU SETTINGS*/
    for(std::pair<std::string,int> mainMenuSetting : this->mainMenuSettings) {
        if(mainMenuSetting.first == settingName) {
            return mainMenuSetting.second;
        }
    }

    /*IN CASE SETTING WAS NOT FOUND*/
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    std::cerr<<"ERROR SETTING NOT FOUND; RETURN A VALUE OF 0"<<std::endl;
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    return 0;
}


int Settings::getPauseMenuSettingsValue(std::string settingName) {
    /*CHECKS THROUGH PAUSE MENU SETTINGS*/
    for(std::pair<std::string,int> pauseMenuSetting : this->pauseMenuSettings) {
        if(pauseMenuSetting.first == settingName) {
            return pauseMenuSetting.second;
        }
    }

    /*IN CASE SETTING WAS NOT FOUND*/
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    std::cerr<<"ERROR SETTING NOT FOUND; RETURN A VALUE OF 0"<<std::endl;
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    return 0;
}

int Settings::getSettingsMenuSettingsValue(std::string settingName) {

    /*CHECKS THROUGH SETTINGS MENU SETTINGS*/
    for(std::pair<std::string,int> settingsMenuSetting : this->settingsMenuSettings) {
        if(settingsMenuSetting.first == settingName) {
            return settingsMenuSetting.second;
        }
    }
    /*IN CASE SETTING WAS NOT FOUND*/
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    std::cerr<<"ERROR SETTING NOT FOUND; RETURN A VALUE OF 0"<<std::endl;
    std::cerr<<"------------------\n------------------\n------------------"<<std::endl;
    return 0;
    
}