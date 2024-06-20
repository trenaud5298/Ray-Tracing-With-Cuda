#ifndef MATERIAL_H
#define MATERIAL_H

class Material {

public:
    Material(float colorR, float colorG, float colorB, float smoothness, bool lightEmitting);
    float colorR;
    float colorG;
    float colorB;
    float smoothness;
    bool lightEmitting;
    

protected:

private:


};

#endif