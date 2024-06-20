#ifndef WORLD_H
#define WORLD_H

#include <string>
#include "./Sphere.h"
#include "./Triangle.h"
#include <vector>
#include <iostream>
#include <vector>


class World {

public:
    World(std::string worldName);
    std::string worldName;
    int numOfObjects;
    std::vector<float*> worldDataVector;

    std::vector<Sphere*> worldSpheres;
    void addSphere(Sphere* newSphere);

    std::vector<Triangle*> worldTriangles;
    void addTriangle(Triangle* newTriangle);

    float** getWorldObjectDataAsArray();

    void printWorldData();
    
protected:

private:



};

#endif