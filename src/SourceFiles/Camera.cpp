#include "../HeaderFiles/Camera.h"
#include <iostream>
#include <cmath>

#define M_PI 3.141592653589793238462643383279502



Camera::Camera(Settings* settings) : fieldOfView(90),horizontalRotation(0),verticalRotation(0),focusOnPoint(false) {
    this->cameraData[0] = settings->getMainSettingsValue("CameraStartPosX");
    this->cameraData[1] = settings->getMainSettingsValue("CameraStartPosY");
    this->cameraData[2] = settings->getMainSettingsValue("CameraStartPosZ");
    std::cout<<"Camera Pos: ( "<<this->cameraData[0]<<", "<<this->cameraData[1]<<", "<<this->cameraData[2]<<" )"<<std::endl;
    this->updateFovScalar();
    this->updateRotationData();
    this->cameraData[8] = settings->getGraphicsSettingsValue("FocalDistance");
    this->cameraData[9] = static_cast<float>(settings->getGraphicsSettingsValue("DepthOfFieldStrength"))/100.0f;
    this->cameraData[10] = this->cameraData[3] * this->cameraData[8];
}

Camera::~Camera() {
    std::cout<<"Camera Deconstructing"<<std::endl;
}



void Camera::updateCameraParameters(Settings* settings) {
    int fov = settings->getGraphicsSettingsValue("FieldOfView");
    this->fieldOfView = fov < 150 ? fov : 150;
    this->horizontalCameraSensitivity = settings->getMainSettingsValue("HorizontalCameraSensitivity");
    this->verticalCameraSensitivity = settings->getMainSettingsValue("VerticalCameraSensitivity");
    this->updateFovScalar();
}


void Camera::updateFovScalar() {
    double radianAngle = static_cast<double>(M_PI * this->fieldOfView / 360);
    this->cameraData[3] = static_cast<float>(std::tan(radianAngle));
}

void Camera::updateRotationData() {
    this->cameraData[4] = std::cos(M_PI*this->horizontalRotation/180);
    this->cameraData[5] = std::cos(M_PI*this->verticalRotation/180);
    this->cameraData[6] = std::sin(M_PI*this->horizontalRotation/180);
    this->cameraData[7] = std::sin(M_PI*this->verticalRotation/180);
}

void Camera::displayInfo() {
    std::cout<<"Position: ( "<<this->cameraData[0]<<", "<<this->cameraData[1]<<", "<<this->cameraData[2]<<" )"<<std::endl;
    std::cout<<"Fov: "<<this->fieldOfView<<std::endl;
    std::cout<<"Fov Scalar: "<<this->cameraData[3]<<std::endl;
}

void Camera::rotateCamera(float horizontalRotation, float verticalRotation) {
    this->horizontalRotation += horizontalRotation;
    float newRotation = this->verticalRotation + verticalRotation;
    this->verticalRotation = (newRotation < -89) ? -89 : (newRotation > 89) ? 89 : newRotation;
    this->updateRotationData();
    this->updateFocalDistance();

}


void Camera::moveForward(float distance) {
    /*X-coord*/this->cameraData[0] += distance*this->cameraData[6];
    /*Z-coord*/this->cameraData[2] += distance*this->cameraData[4];
    this->updateFocalDistance();

}
void Camera::moveBackward(float distance) {
    /*X-coord*/this->cameraData[0] -= distance*this->cameraData[6];
    /*Z-coord*/this->cameraData[2] -= distance*this->cameraData[4];
    this->updateFocalDistance();

}
void Camera::moveRight(float distance) {
    /*X-coord*/this->cameraData[0] += distance*this->cameraData[4];
    /*Z-coord*/this->cameraData[2] -= distance*this->cameraData[6];
    this->updateFocalDistance();

}
void Camera::moveLeft(float distance) {
    /*X-coord*/this->cameraData[0] -= distance*this->cameraData[4];
    /*Z-coord*/this->cameraData[2] += distance*this->cameraData[6];
    this->updateFocalDistance();

}
void Camera::moveUp(float distance) {
    /*Y-coord*/this->cameraData[1] += distance;
    this->updateFocalDistance();

}
void Camera::moveDown(float distance) {
    /*Y-coord*/this->cameraData[1] -= distance;
    this->updateFocalDistance();
}



void Camera::increaseDepthOfFieldStrength(float amount) {
    this->cameraData[9] = this->cameraData[9] + amount;
}

void Camera::decreaseDepthOfFieldStrength(float amount) {
    this->cameraData[9] = 0 > this->cameraData[9]-amount ? 0 : this->cameraData[9]-amount;
}

void Camera::incrimentFocalDistance(float amount) {
    this->cameraData[8] = this->cameraData[8] + amount;
    this->cameraData[10] = this->cameraData[3] * this->cameraData[8];
    std::cout<<"Focal Distance: "<<this->cameraData[8]<<std::endl;
}

void Camera::decreaseFocalDistance(float amount) {
    this->cameraData[8] = 0 > this->cameraData[8]-amount ? 0 : this->cameraData[8] - amount;
    this->cameraData[10] = this->cameraData[3] * this->cameraData[8];
    std::cout<<"Focal Distance: "<<this->cameraData[8]<<std::endl;
}


void Camera::setFocalPoint(float x, float y, float z) {
    this->focusOnPoint = true;
    this->focalPointX = x;
    this->focalPointY = y;
    this->focalPointZ = z;
    this->updateFocalDistance();
}

void Camera::updateFocalDistance() {
    if(!this->focusOnPoint) {
        return;
    }
    float objectVecX,objectVecY,objectVecZ;
    float cameraDirectionVecX,cameraDirectionVecY,cameraDirectionVecZ;
    float tempVecX, tempVecY, tempVecZ;
    objectVecX = this->focalPointX - this->cameraData[0];
    objectVecY = this->focalPointY - this->cameraData[1];
    objectVecZ = this->focalPointZ - this->cameraData[2];
    tempVecX = 0;
    tempVecY = 0;
    tempVecZ = 1;
    float verticalZ = tempVecY*this->cameraData[7] + tempVecZ*this->cameraData[5];
    cameraDirectionVecX = tempVecX*this->cameraData[4]+verticalZ*this->cameraData[6];
    cameraDirectionVecY = tempVecY*this->cameraData[5]-tempVecZ*this->cameraData[7];
    cameraDirectionVecZ =  verticalZ*this->cameraData[4]-tempVecX*this->cameraData[6];
    float newFocalDistance = cameraDirectionVecX*objectVecX + cameraDirectionVecY*objectVecY + cameraDirectionVecZ*objectVecZ;
    this->cameraData[8] = 0 >= newFocalDistance ? 0.01 : newFocalDistance;
    this->cameraData[10] = this->cameraData[3] * this->cameraData[8];
}