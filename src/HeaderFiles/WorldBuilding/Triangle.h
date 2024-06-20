#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "./Vertex.h"
#include "./Material.h"
#include <math.h>
#include <iostream>

class Triangle {

public:
    Triangle(float pos1X, float pos1Y, float pos1Z, float pos2X, float pos2Y, float pos2Z, float pos3X, float pos3Y, float pos3Z, Material* triangleMaterial);
    Triangle(Vertex* vertex1, Vertex* vertex2, Vertex* vertex3, Material* triangleMaterial);
    float* triangleData;
    void setMaterial(Material* triangleMaterial);
    void updateArrayData(Vertex* vertexA, Vertex* vertexB, Vertex* vertexC);
    

protected:

private:
};


#endif