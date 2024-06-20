#ifndef VERTEX_H
#define VERTEX_H

class Vertex {
public:
    Vertex(float posX, float posY, float posZ);
    void shiftVertex(float xShift, float yShift, float zShift);
    void rotateVertexAboutOrigin(float xRotation, float yRotation, float zRotation);
    void rotateVertexAboutVertex(Vertex* rotationPoint, float xRotation, float yRotation, float zRotation);
    float posX;
    float posY;
    float posZ;
protected:


private:

};

#endif