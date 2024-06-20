#ifndef SPHERE_H
#define SPHERE_H

#include "./Material.h"
#include <iostream>

class Sphere {

public:
    Sphere(float posX, float posY, float posZ, float radius, Material* sphereMaterial);
    void setMaterial(Material* newMaterial);
    float* sphereData;
    void printData();
    
protected:
    
private:

};

#endif