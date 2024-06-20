#include "../../HeaderFiles/WorldBuilding/Sphere.h"


Sphere::Sphere(float posX, float posY, float posZ, float radius, Material* sphereMaterial) {
    this->sphereData = new float[30];
    this->sphereData[0] = 0;
    this->setMaterial(sphereMaterial);
    this->sphereData[6] = posX;
    this->sphereData[7] = posY;
    this->sphereData[8] = posZ;
    this->sphereData[9] = radius;

}

void Sphere::setMaterial(Material* sphereMaterial) {
    this->sphereData[1] = sphereMaterial->colorR;
    this->sphereData[2] = sphereMaterial->colorG;
    this->sphereData[3] = sphereMaterial->colorB;
    this->sphereData[4] = sphereMaterial->smoothness;
    if(sphereMaterial->lightEmitting) {
        this->sphereData[5] = 1.0;
    } else {
        this->sphereData[5] = 0.0;
    }
}

void Sphere::printData() {
    std::cout<<"---------------------------"<<std::endl<<"Sphere"<<std::endl;
    for(int i =0; i < 9; i++) {
        std::cout<<"Index: "<<i<<"   Value: "<<this->sphereData[i]<<std::endl;
    }
}
