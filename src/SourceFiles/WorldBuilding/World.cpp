#include "../../HeaderFiles/WorldBuilding/World.h"

World::World(std::string worldName) : worldName(worldName),numOfObjects(0) {

}

void World::addSphere(Sphere* newSphere) {
    this->worldSpheres.push_back(newSphere);
    this->worldDataVector.push_back(newSphere->sphereData);
    this->numOfObjects += 1;
}

void World::addTriangle(Triangle* newTriangle) {
    this->worldTriangles.push_back(newTriangle);
    this->worldDataVector.push_back(newTriangle->triangleData);
    this->numOfObjects+=1;
}


float** World::getWorldObjectDataAsArray() {

    float** worldDataAsArray = new float*[this->numOfObjects];
    for(size_t i = 0; i < this->worldDataVector.size(); i++) {
        worldDataAsArray[i] = this->worldDataVector[i];
    }
    return worldDataAsArray;
}

void World::printWorldData() {
    for(size_t i = 0; i < this->worldDataVector.size(); i++) {
        std::cout<<"-----------------------"<<std::endl<<"WORLD OBJECT "<<i<<std::endl;
        for(size_t j = 0; j < 30; j ++) {
            std::cout<<"Index: "<<j<<"   Value: "<<this->worldDataVector[i][j]<<std::endl;
        }
    }
}