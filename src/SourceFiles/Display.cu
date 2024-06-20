#include "../HeaderFiles/Display.cuh"
#include <iostream>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <chrono>
#include <thread>
#include "../../lib/SDL/include/SDL.h"
#include "../../lib/SDL/include/SDL_ttf.h"

#include "../HeaderFiles/Camera.h"
#include "../HeaderFiles/Button.h"


#define MAIN_MENU 0
#define RUNNING 1
#define PAUSED 2


__device__ float3 getVector(int x,int y,int imageWidth,int imageHeight,float* CameraData, float3* rayOrigin, curandState* randomGenerator) {
    float heightToWidthRatio = static_cast<float>(imageHeight)/imageWidth;
    float3 tempRayVector = make_float3(CameraData[10]*(static_cast<float>(2*x-imageWidth)/imageWidth),heightToWidthRatio*(CameraData[10]*static_cast<float>(imageHeight-2*y)/imageHeight),static_cast<float>(CameraData[8]));
    float offsetX = CameraData[9]*curand_normal(randomGenerator);
    float offsetY = CameraData[9]*curand_normal(randomGenerator);
    tempRayVector.x -= offsetX;
    tempRayVector.y -= offsetY;
    (*rayOrigin).x += offsetX;
    (*rayOrigin).y += offsetY;
    float inverseLength = rsqrtf(tempRayVector.x*tempRayVector.x+tempRayVector.y*tempRayVector.y+tempRayVector.z*tempRayVector.z);
    tempRayVector.x *= inverseLength;
    tempRayVector.y *= inverseLength;
    tempRayVector.z *= inverseLength;

    float verticalZ = tempRayVector.y*CameraData[7] + tempRayVector.z*CameraData[5];
    return make_float3(tempRayVector.x*CameraData[4]+verticalZ*CameraData[6],tempRayVector.y*CameraData[5]-tempRayVector.z*CameraData[7],verticalZ*CameraData[4]-tempRayVector.x*CameraData[6]);
}

__device__ float3 getSkyBoxColor(float3 rayVector, float angleOfSun) {
    /* FUTURE IDEA:
     * Use an Image(s) of a SkyBox and In Python
     * Convert all of the RGB Values for each
     * pixels to floats between 0 and 1 and then
     * store that data in a file to be read in
     * for this program and then store that data in
     * an array and be prepared to convert from
     * a vector to a pixel from the skybox
    */
    float3 skyColor;
    float3 sunColor = make_float3(1.0f,1.0f,1.0f);
    if(angleOfSun > 0) {
        float angleBetweenRayAndSun = acos(rayVector.x*cosf(angleOfSun) + rayVector.y*sinf(angleOfSun));
        if(angleBetweenRayAndSun < 0.05) {
            skyColor.x = sunColor.x;
            skyColor.y = sunColor.y;
            skyColor.z = sunColor.z;
            return skyColor;
        }
        float t = 0.5f*(1.0f+(sinf(acos(rsqrtf(rayVector.x*rayVector.x+rayVector.z*rayVector.z)*(rayVector.x*rayVector.x+rayVector.z*rayVector.z)))));
        float3 topSkyColor = make_float3(64.0f/255.0f,156.0f/255.0f,255.0f/255.0f);
        float3 bottomSkyColor = make_float3(1.0f,1.0f,1.0f);
        skyColor.x = topSkyColor.x*t + bottomSkyColor.x*(1-t);
        skyColor.y = topSkyColor.y*t + bottomSkyColor.y*(1-t);
        skyColor.z = topSkyColor.z*t + bottomSkyColor.z*(1-t);
        return skyColor;

    } else {
        skyColor.x = 0;
        skyColor.y = 0;
        skyColor.z = 0;
        return skyColor;
    }
}

__device__ float computeCrossProductLength(float3 vec1, float3 vec2) {
    return sqrtf(
        powf((vec1.y*vec2.z-vec2.y*vec1.z),2)+
        powf((vec1.z*vec2.x-vec2.z*vec1.x),2)+
        powf((vec1.x*vec2.y-vec2.x*vec1.y),2)
    );
}

__device__ float3 getRandomBounceVector(float3 normalVector, curandState* randomGenerator) {
    float3 randomBounceVector = make_float3(1.0f,0,0);
    float x;
    float y;
    float z;
    float squaredLength;
    for(int i = 0; i < 100; i ++) {
        x = curand_normal(randomGenerator);
        y = curand_normal(randomGenerator);
        z = curand_normal(randomGenerator);
        squaredLength = sqrtf(x*x+y*y+z*z); 
        if(squaredLength < 1) {
            randomBounceVector.x = x/squaredLength;
            randomBounceVector.y = y/squaredLength;
            randomBounceVector.z = z/squaredLength;
            break;
        }
    }
    if((normalVector.x*randomBounceVector.x + normalVector.y*randomBounceVector.y + normalVector.z*randomBounceVector.z) < 0) {
        randomBounceVector.x *= -1;
        randomBounceVector.y *= -1;
        randomBounceVector.z *= -1;

    }

    return randomBounceVector;
}

__device__ float3 getCorrectBounceVector(float3 normalVector, float3 rayVector, float intersectionObjectSmoothness,curandState* randomGenerator) {
    float3 diffuseBounceRay = make_float3(1.0f,0,0);
    float x;
    float y;
    float z;
    float squaredLength;
    for(int i = 0; i < 100; i ++) {
        x = curand_normal(randomGenerator);
        y = curand_normal(randomGenerator);
        z = curand_normal(randomGenerator);
        squaredLength = sqrtf(x*x+y*y+z*z); 
        if(squaredLength < 1) {
            x += normalVector.x;
            y += normalVector.y;
            z += normalVector.z;
            squaredLength = rsqrtf(x*x+y*y+z*z);
            diffuseBounceRay.x = x*squaredLength;
            diffuseBounceRay.y = y*squaredLength;
            diffuseBounceRay.z = z*squaredLength;
            break;
        }
    }
    
    float3 specularBounceRay;
    float doubleDotProduct = 2.0 * (rayVector.x*normalVector.x + rayVector.y*normalVector.y + rayVector.z*normalVector.z);
    specularBounceRay = make_float3 (
        rayVector.x - (doubleDotProduct*normalVector.x),
        rayVector.y - (doubleDotProduct*normalVector.y),
        rayVector.z - (doubleDotProduct*normalVector.z)
    );
    float inverseLength = 1/sqrtf(specularBounceRay.x*specularBounceRay.x + specularBounceRay.y*specularBounceRay.y + specularBounceRay.z*specularBounceRay.z);
    specularBounceRay.x *= inverseLength;
    specularBounceRay.y *= inverseLength;
    specularBounceRay.z *= inverseLength;



    float3 actualBounceRay = make_float3(
        specularBounceRay.x*intersectionObjectSmoothness + diffuseBounceRay.x*(1-intersectionObjectSmoothness),
        specularBounceRay.y*intersectionObjectSmoothness + diffuseBounceRay.y*(1-intersectionObjectSmoothness),
        specularBounceRay.z*intersectionObjectSmoothness + diffuseBounceRay.z*(1-intersectionObjectSmoothness)
    );

    return actualBounceRay;
}

__device__ float3 getCorrectBounceVector2(float3 normalVector, float3 rayVector, float intersectionObjectSmoothness,curandState* randomGenerator) {
    float3 diffuseBounceRay = make_float3(1.0f,0,0);
    float x;
    float y;
    float z;
    float squaredLength;
    float randomVal = curand_uniform(randomGenerator);
    float3 actualBounceRay;
    if(randomVal < intersectionObjectSmoothness) {
        x/*Double Dot Product*/ = 2.0 * (rayVector.x*normalVector.x + rayVector.y*normalVector.y + rayVector.z*normalVector.z);
        actualBounceRay = make_float3(
            rayVector.x - (x*normalVector.x),
            rayVector.y - (x*normalVector.y),
            rayVector.z - (x*normalVector.z)
        );
    } else {
        for(int i = 0; i < 100; i ++) {
            x = curand_normal(randomGenerator);
            y = curand_normal(randomGenerator);
            z = curand_normal(randomGenerator);
            squaredLength = sqrtf(x*x+y*y+z*z); 
            if(squaredLength < 1) {
                x += normalVector.x;
                y += normalVector.y;
                z += normalVector.z;
                squaredLength = rsqrtf(x*x+y*y+z*z);
                actualBounceRay.x = x*squaredLength;
                actualBounceRay.y = y*squaredLength;
                actualBounceRay.z = z*squaredLength;
                break;
            }
        }
    }
    return actualBounceRay;
}

__device__ float* handleRayIntersection(float3 rayOrigin, float3 rayVector, float**worldObjectData, int numWorldObjects, bool* intersected, float* shortestIntersectionDistance) {
    float* intersectionObject;
    float deltaX,deltaY,deltaZ,B,C,discriminant,distance;
    float3 intersectionPoint;
    for(int i = 0; i < numWorldObjects; i++) {
        switch(static_cast<int>(worldObjectData[i][0])) {
            case 0:
                deltaX = rayOrigin.x - worldObjectData[i][6];
                deltaY = rayOrigin.y - worldObjectData[i][7];
                deltaZ = rayOrigin.z - worldObjectData[i][8];
                B = 2*deltaX*rayVector.x + 2*deltaY*rayVector.y+2*deltaZ*rayVector.z;
                C = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ - (worldObjectData[i][9]*worldObjectData[i][9]);
                discriminant = (B*B)-(4*C);
                if(discriminant < 0){break;}
                distance = (-B-sqrtf(discriminant))*(0.5);
                
                if(distance > 0 ){
                    if(!*intersected) {
                        *intersected = true;
                        *shortestIntersectionDistance = distance;
                        intersectionObject = worldObjectData[i];
                    } else if(distance < *shortestIntersectionDistance){
                        *intersected = true;
                        *shortestIntersectionDistance = distance;
                        intersectionObject = worldObjectData[i];
                    }
                }
                break;

            
            case 1:
                distance = 
                (worldObjectData[i][9]-worldObjectData[i][6]*rayOrigin.x-worldObjectData[i][7]*rayOrigin.y-worldObjectData[i][8]*rayOrigin.z)
                /(worldObjectData[i][6]*rayVector.x + worldObjectData[i][7]*rayVector.y + worldObjectData[i][8]*rayVector.z);
                if(distance < 0 || isnan(distance) || distance > *shortestIntersectionDistance) {break;}
                //Determine if point is in triangle
                intersectionPoint.x = rayOrigin.x + rayVector.x*distance;
                intersectionPoint.y = rayOrigin.y + rayVector.y*distance;
                intersectionPoint.z = rayOrigin.z + rayVector.z*distance;
                //AB Inside Test
                if(((intersectionPoint.x-worldObjectData[i][10])*worldObjectData[i][19]  +  
                (intersectionPoint.y-worldObjectData[i][11])*worldObjectData[i][20]  +  
                (intersectionPoint.z-worldObjectData[i][12])*worldObjectData[i][21]) < 0) {break;}

                //BC Inside Test
                if(((intersectionPoint.x-worldObjectData[i][13])*worldObjectData[i][22]  +  
                (intersectionPoint.y-worldObjectData[i][14])*worldObjectData[i][23]  +  
                (intersectionPoint.z-worldObjectData[i][15])*worldObjectData[i][24]) < 0) {break;}

                //AC Inside Test
                if(((intersectionPoint.x-worldObjectData[i][16])*worldObjectData[i][25]  +  
                (intersectionPoint.y-worldObjectData[i][17])*worldObjectData[i][26]  +  
                (intersectionPoint.z-worldObjectData[i][18])*worldObjectData[i][27]) < 0) {break;}

                
                *shortestIntersectionDistance = distance;
                intersectionObject = worldObjectData[i];
                *intersected = true;
                
                break;
        }
    }
    return intersectionObject;
}





__device__ float3 getNormalVector(float* intersectionObject, float3 intersectionPoint, float3 rayVector) {
    float3 normalVector;
    switch(static_cast<int>(intersectionObject[0])) {
        case 0:
            normalVector.x = intersectionPoint.x-intersectionObject[6];
            normalVector.y = intersectionPoint.y-intersectionObject[7];
            normalVector.z = intersectionPoint.z-intersectionObject[8];
            break;

        case 1:
            normalVector.x = intersectionObject[6];
            normalVector.y = intersectionObject[7];
            normalVector.z = intersectionObject[8];
            break;
    }
    if((normalVector.x*rayVector.x + normalVector.y*rayVector.y + normalVector.z*rayVector.z)>0) {
        normalVector.x*=-1;
        normalVector.y*=-1;
        normalVector.z*=-1;
    }
    //NORMALIZE VECTOR:
    float inverseLength = rsqrtf(normalVector.x*normalVector.x + normalVector.y*normalVector.y + normalVector.z*normalVector.z);
    normalVector.x *= inverseLength;
    normalVector.y *= inverseLength;
    normalVector.z *= inverseLength;
    return normalVector;

}


__global__ void DiffuseRender(unsigned int seed, int frame,int maxRayBounces, int imageWidth, int imageHeight,int numWorldObjects, float* cameraData, float* imageDataArray, float** worldObjectData,float* troubleShootData) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < imageWidth && y < imageHeight) {
        int pixelIndex = (y * imageWidth + x) * 3;
        
        curandState randomGenerator;
        curand_init(seed+frame*79231,x+y*imageWidth*blockDim.z,0,&randomGenerator);
        //Rotate Vector:
        float3 rayOrigin = make_float3(cameraData[0],cameraData[1],cameraData[2]);
        float3 rayVector = getVector(x,y,imageWidth,imageHeight,cameraData,&rayOrigin,&randomGenerator);
        float3 tempColor;
        float3 color = make_float3(1.0,1.0,1.0);
        float3 lightColor = make_float3(0.0f,0.0f,0.0f);
        float colorScalar;
        float3 lightRecieved = {0.0f,0.0f,0.0f};
        float3 normalVector;
        bool intersected;
        float shortestIntersectionDistance;
        float* intersectionObject;
        for(int i =0; i < maxRayBounces; i ++) {
            intersected = false;
            shortestIntersectionDistance = 999999999999;
            intersectionObject = handleRayIntersection(rayOrigin,rayVector,worldObjectData,numWorldObjects,&intersected,&shortestIntersectionDistance);
            if(intersected){
                tempColor.x = intersectionObject[1];
                tempColor.y = intersectionObject[2];
                tempColor.z = intersectionObject[3];
                if(intersectionObject[5] > 0.5) {
                    lightRecieved.x += tempColor.x * color.x;
                    lightRecieved.y += tempColor.y * color.y;
                    lightRecieved.z += tempColor.z * color.z;
                }
                colorScalar = 1.0f/fmaxf(1.0f,fmaxf(tempColor.x,fmaxf(tempColor.y,tempColor.z)));
                color.x *= (tempColor.x*colorScalar);
                color.y *= (tempColor.y*colorScalar);
                color.z *= (tempColor.z*colorScalar);

               
                rayOrigin.x += rayVector.x*shortestIntersectionDistance;
                rayOrigin.y += rayVector.y*shortestIntersectionDistance;
                rayOrigin.z += rayVector.z*shortestIntersectionDistance;
                normalVector = getNormalVector(intersectionObject,rayOrigin,rayVector);
                rayOrigin.x+=normalVector.x*0.000001;
                rayOrigin.y+=normalVector.y*0.000001;
                rayOrigin.z+=normalVector.z*0.000001;
                rayVector = getCorrectBounceVector(normalVector,rayVector,intersectionObject[4],&randomGenerator);
            } else{
                // lightColor = getSkyBoxColor(rayVector,0.3);
                break;
            }
        }
        
        imageDataArray[pixelIndex] = imageDataArray[pixelIndex] + lightRecieved.x; // Red
        imageDataArray[pixelIndex + 1] = imageDataArray[pixelIndex+1] + lightRecieved.y; // Green
        imageDataArray[pixelIndex + 2] = imageDataArray[pixelIndex+2] + lightRecieved.z;
    }
};


__global__ void updateImage(int frame, int imageWidth, int imageHeight,float* imageDataArray,unsigned char* imageData,float* troubleShootData) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < imageWidth && y < imageHeight) {
        int pixelIndex = (y * imageWidth + x) * 3;
        float3 pixelColor = make_float3(
            static_cast<float>(imageDataArray[pixelIndex]/static_cast<float>(frame)),
            static_cast<float>(imageDataArray[pixelIndex+1]/static_cast<float>(frame)),
            static_cast<float>(imageDataArray[pixelIndex+2]/static_cast<float>(frame))
        );

        float maxColorVal = fmaxf(pixelColor.x,fmaxf(pixelColor.y,pixelColor.z));
        float colorScale = 1;
        if(maxColorVal > 1) {
            colorScale = 1/maxColorVal;
        }
        
        imageData[pixelIndex] = static_cast<unsigned char>(255.99f*colorScale*pixelColor.x); // Red
        imageData[pixelIndex + 1] = static_cast<unsigned char>(255.99f*colorScale*pixelColor.y); // Green
        imageData[pixelIndex + 2] = static_cast<unsigned char>(255.99f*colorScale*pixelColor.z); // Blue
    }
};

__global__ void updateImageHDR(int frame, int imageWidth, int imageHeight,float* imageDataArray,unsigned char* imageData,float* troubleShootData) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < imageWidth && y < imageHeight) {
        int pixelIndex = (y * imageWidth + x) * 3;
        float3 pixelColor = make_float3(
            static_cast<float>(imageDataArray[pixelIndex]/static_cast<float>(frame)),
            static_cast<float>(imageDataArray[pixelIndex+1]/static_cast<float>(frame)),
            static_cast<float>(imageDataArray[pixelIndex+2]/static_cast<float>(frame))
        );

        
        //Tone Map
        float luminance = 0.2126 * pixelColor.x + 0.7152 * pixelColor.y + 0.0722 * pixelColor.z;
        float mapped_luminance = luminance / (1.0f + luminance);
        float scaleFactor;
        if(luminance > 0) {
            scaleFactor = mapped_luminance / luminance;
            pixelColor.x *= scaleFactor;
            pixelColor.y *= scaleFactor;
            pixelColor.z *= scaleFactor;
        }
       
        // //Limits Max Color to 1
        // float maxColorVal = fmaxf(pixelColor.x,fmaxf(pixelColor.y,pixelColor.z));
        // float colorScale = 1;
        // if(maxColorVal > 1) {
        //     colorScale = 1/maxColorVal;
        // }
        // pixelColor.x*= colorScale;
        // pixelColor.y*= colorScale;
        // pixelColor.z*= colorScale;

        
        
        //Apply Gamma Correction
        float inverseGamma = 1.0f/2.2f;
        pixelColor.x = powf(pixelColor.x,inverseGamma);
        pixelColor.y = powf(pixelColor.y,inverseGamma);
        pixelColor.z = powf(pixelColor.z,inverseGamma);
        
        

        imageData[pixelIndex] = static_cast<unsigned char>(255.99f*fminf(1.0,pixelColor.x)); // Red
        imageData[pixelIndex + 1] = static_cast<unsigned char>(255.99f*fminf(1.0,pixelColor.y)); // Green
        imageData[pixelIndex + 2] = static_cast<unsigned char>(255.99f*fminf(1.0,pixelColor.z)); // Blue
        
    }
};


__global__ void newRender(unsigned int seed, int frame, int maxRayBounces, int imageWidth, int imageHeight,int numWorldObjects, float* cameraData, float* rawImageDataArray, unsigned char* actualImageDataArray, float** worldObjectData,float* troubleShootData) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < imageWidth && y < imageHeight) {
        int pixelIndex = (y * imageWidth + x) * 3;
        
        curandState randomGenerator;
        curand_init(seed+frame*79231,x+y*imageWidth*blockDim.z,0,&randomGenerator);
        //Rotate Vector:
        float3 rayOrigin = make_float3(cameraData[0],cameraData[1],cameraData[2]);
        float3 rayVector = getVector(x,y,imageWidth,imageHeight,cameraData,&rayOrigin,&randomGenerator);
        float3 tempColor;
        float3 color = make_float3(1.0,1.0,1.0);
        float3 lightColor = make_float3(0.0f,0.0f,0.0f);
        float colorScalar;
        float3 lightRecieved = {0.0f,0.0f,0.0f};
        float3 normalVector;
        bool intersected;
        float shortestIntersectionDistance;
        float* intersectionObject;
        for(int i =0; i < maxRayBounces; i ++) {
            intersected = false;
            shortestIntersectionDistance = 999999999999;
            intersectionObject = handleRayIntersection(rayOrigin,rayVector,worldObjectData,numWorldObjects,&intersected,&shortestIntersectionDistance);
            if(intersected){
                tempColor.x = intersectionObject[1];
                tempColor.y = intersectionObject[2];
                tempColor.z = intersectionObject[3];
                if(intersectionObject[5] > 0.5) {
                    lightRecieved.x += tempColor.x * color.x;
                    lightRecieved.y += tempColor.y * color.y;
                    lightRecieved.z += tempColor.z * color.z;
                }
                colorScalar = 1.0f/fmaxf(1.0f,fmaxf(tempColor.x,fmaxf(tempColor.y,tempColor.z)));
                color.x *= (tempColor.x*colorScalar);
                color.y *= (tempColor.y*colorScalar);
                color.z *= (tempColor.z*colorScalar);

               
                rayOrigin.x += rayVector.x*shortestIntersectionDistance;
                rayOrigin.y += rayVector.y*shortestIntersectionDistance;
                rayOrigin.z += rayVector.z*shortestIntersectionDistance;
                normalVector = getNormalVector(intersectionObject,rayOrigin,rayVector);
                rayOrigin.x+=normalVector.x*0.000001;
                rayOrigin.y+=normalVector.y*0.000001;
                rayOrigin.z+=normalVector.z*0.000001;
                rayVector = getCorrectBounceVector2(normalVector,rayVector,intersectionObject[4],&randomGenerator);
            } else{
                // lightColor = getSkyBoxColor(rayVector,0.3);
                break;
            }
        }


        float proportionalConstant = 1.0/static_cast<float>(frame);

        float3 updatedPixelColor = make_float3(
            rawImageDataArray[pixelIndex]   * (1-proportionalConstant) + lightRecieved.x * proportionalConstant,
            rawImageDataArray[pixelIndex+1] * (1-proportionalConstant) + lightRecieved.y * proportionalConstant,
            rawImageDataArray[pixelIndex+2] * (1-proportionalConstant) + lightRecieved.z * proportionalConstant
        );

        rawImageDataArray[pixelIndex]   = updatedPixelColor.x;
        rawImageDataArray[pixelIndex+1] = updatedPixelColor.y;
        rawImageDataArray[pixelIndex+2] = updatedPixelColor.z;

        // float maxColorVal = fmaxf(updatedPixelColor.x,fmaxf(updatedPixelColor.y,updatedPixelColor.z));
        // float colorScale = 1;
        // if(maxColorVal > 1) {
        //     colorScale = 1/maxColorVal;
        // }
        
        // actualImageDataArray[pixelIndex] = static_cast<unsigned char>(255.99f*colorScale*updatedPixelColor.x); // Red
        // actualImageDataArray[pixelIndex + 1] = static_cast<unsigned char>(255.99f*colorScale*updatedPixelColor.y); // Green
        // actualImageDataArray[pixelIndex + 2] = static_cast<unsigned char>(255.99f*colorScale*updatedPixelColor.z); // Blue


        //HDR TONE MAP
        // //Tone Map Colors
        float luminance = 0.2126 * updatedPixelColor.x + 0.7152 * updatedPixelColor.y + 0.0722 * updatedPixelColor.z;
        float mapped_luminance = luminance / (1.0f + luminance);
        float scaleFactor;
        if(luminance > 0) {
            scaleFactor = mapped_luminance / luminance;
            updatedPixelColor.x *= scaleFactor;
            updatedPixelColor.y *= scaleFactor;
            updatedPixelColor.z *= scaleFactor;
        }

        //Apply Gamma Correction
        float inverseGamma = 1.0f/2.2f;
        updatedPixelColor.x = powf(updatedPixelColor.x,inverseGamma);
        updatedPixelColor.y = powf(updatedPixelColor.y,inverseGamma);
        updatedPixelColor.z = powf(updatedPixelColor.z,inverseGamma);

        
        actualImageDataArray[pixelIndex] = static_cast<unsigned char>(255.99f*fminf(1.0,updatedPixelColor.x)); // Red
        actualImageDataArray[pixelIndex + 1] = static_cast<unsigned char>(255.99f*fminf(1.0,updatedPixelColor.y)); // Green
        actualImageDataArray[pixelIndex + 2] = static_cast<unsigned char>(255.99f*fminf(1.0,updatedPixelColor.z)); // Blue
        
    }
};







//Need to Redesign how function works (consider using structs)
__global__ void DiffuseRender2(unsigned int seed, int frame, int maxRayBounces, int imageWidth, int imageHeight,int numWorldObjects, float* cameraData, float* rawImageDataArray, float** worldObjectData,float* troubleShootData) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < imageWidth && y < imageHeight) {
        int pixelIndex = (y * imageWidth + x) * 3;
        
        curandState randomGenerator;
        curand_init(seed+frame*79231,x+y*imageWidth*blockDim.z,0,&randomGenerator);
        //Rotate Vector:
        float3 initialRayOrigin = make_float3(cameraData[0],cameraData[1],cameraData[2]);
        float3 initialRayVector = getVector(x,y,imageWidth,imageHeight,cameraData,&initialRayOrigin,&randomGenerator);
        float3 initialColor = {0.0f,0.0f,0.0f};
        float3 tempColor = {0.0f,0.0f,0.0f};
        float3 lightRecieved = {0.0f,0.0f,0.0f};
        float3 initialNormalVector;
        bool intersected;
        float shortestIntersectionDistance;
        float* intersectionObject;
        float colorScalar;
        //Determine the initial intersection Point:
        intersected = false;
        shortestIntersectionDistance = 999999999999;
        intersectionObject = handleRayIntersection(initialRayOrigin,initialRayVector,worldObjectData,numWorldObjects,&intersected,&shortestIntersectionDistance);
        if(intersected) {
        
            //Handle Starting Color/Light
            tempColor.x = intersectionObject[1];
            tempColor.y = intersectionObject[2];
            tempColor.z = intersectionObject[3];
            if(intersectionObject[5] > 0.5) {
                lightRecieved.x += tempColor.x;
                lightRecieved.y += tempColor.y;
                lightRecieved.z += tempColor.z;
            }
            colorScalar = 1.0f/fmaxf(1.0f,fmaxf(tempColor.x,fmaxf(tempColor.y,tempColor.z)));
            initialColor.x = (tempColor.x*colorScalar);
            initialColor.y = (tempColor.y*colorScalar);
            initialColor.z = (tempColor.z*colorScalar);
            //Shifts Origin To First Intersection Point
            initialRayOrigin = {
                initialRayOrigin.x + initialRayVector.x*shortestIntersectionDistance,
                initialRayOrigin.y + initialRayVector.y*shortestIntersectionDistance,
                initialRayOrigin.z + initialRayVector.z*shortestIntersectionDistance
            };
            initialNormalVector = getNormalVector(intersectionObject,initialRayOrigin,initialRayVector);
            initialRayOrigin.x += initialNormalVector.x*1e-6;
            initialRayOrigin.y += initialNormalVector.y*1e-6;
            initialRayOrigin.z += initialNormalVector.z*1e-6;

            float smoothness = intersectionObject[4];
            //Handles 10 Unique Bounces From Initial Point
            float3 rayVector;
            float3 rayOrigin;
            float3 normalVector;
            float3 bounceColor;
            for(int i = 0; i < 10; i ++) {
                rayVector = getCorrectBounceVector2(initialNormalVector,initialRayVector,smoothness,&randomGenerator);
                rayOrigin = initialRayOrigin;
                bounceColor = initialColor;
                for(int j = 0; j < maxRayBounces-1; j++) {
                    intersected = false;
                    shortestIntersectionDistance = 999999999999;
                    intersectionObject = handleRayIntersection(rayOrigin,rayVector,worldObjectData,numWorldObjects,&intersected,&shortestIntersectionDistance);
                    if(intersected){
                        tempColor.x = intersectionObject[1];
                        tempColor.y = intersectionObject[2];
                        tempColor.z = intersectionObject[3];
                        if(intersectionObject[5] > 0.5) {
                            lightRecieved.x += 0.1f * tempColor.x * bounceColor.x;
                            lightRecieved.y += 0.1f * tempColor.y * bounceColor.y;
                            lightRecieved.z += 0.1f * tempColor.z * bounceColor.z;
                        }
                        colorScalar = 1.0f/fmaxf(1.0f,fmaxf(tempColor.x,fmaxf(tempColor.y,tempColor.z)));
                        bounceColor.x *= (tempColor.x*colorScalar);
                        bounceColor.y *= (tempColor.y*colorScalar);
                        bounceColor.z *= (tempColor.z*colorScalar);

                    
                        rayOrigin.x += rayVector.x*shortestIntersectionDistance;
                        rayOrigin.y += rayVector.y*shortestIntersectionDistance;
                        rayOrigin.z += rayVector.z*shortestIntersectionDistance;
                        normalVector = getNormalVector(intersectionObject,rayOrigin,rayVector);
                        rayOrigin.x+=normalVector.x*0.000001;
                        rayOrigin.y+=normalVector.y*0.000001;
                        rayOrigin.z+=normalVector.z*0.000001;
                        rayVector = getCorrectBounceVector2(normalVector,rayVector,intersectionObject[4],&randomGenerator);
                    } else{
                        // lightColor = getSkyBoxColor(rayVector,0.3);
                        break;
                    }

                }
            }
        }

        rawImageDataArray[pixelIndex]    =  rawImageDataArray[pixelIndex]    +  lightRecieved.x; // Red
        rawImageDataArray[pixelIndex+1]  =  rawImageDataArray[pixelIndex+1]  +  lightRecieved.y; // Green
        rawImageDataArray[pixelIndex+2]  =  rawImageDataArray[pixelIndex+2]  +  lightRecieved.z;
        
    }
};







Display::Display(int displayWidth,int displayHeight,Camera* viewCamera): displayWidth(displayWidth),displayHeight(displayHeight),viewCamera(viewCamera),displayState(MAIN_MENU),currentFrame(0) {
    this->imageData = new unsigned char[this->displayWidth * this->displayHeight * 3];
    std::cout<<"Display Created"<<std::endl;

    //Initializes SDL VIDEO
    if (SDL_Init(SDL_INIT_VIDEO) != 0 || TTF_Init() == -1) {
        std::cout<<"MAJOR ERROR: FAILURE TO LOAD SDL"<<std::endl;
    }

    //Make Display Window
    this->displayWindow = SDL_CreateWindow("RayTracing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, this->displayWidth, this->displayHeight, SDL_WINDOW_RESIZABLE);
    if (this->displayWindow==NULL) {
        std::cout<<"MAJOR ERROR: FAILURE TO CREATE SDL WINDOW"<<std::endl;
    }
    this->resetCursor();
    
    //Initialize and Create SDL Renderer to display image
    this->imageRenderer = SDL_CreateRenderer(this->displayWindow, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetRenderDrawBlendMode(this->imageRenderer, SDL_BLENDMODE_BLEND);
    
    this->imageTexture = SDL_CreateTexture(this->imageRenderer, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, this->displayWidth, this->displayHeight);

    this->allocateGPUMemory();
    this->recieveGpuTroubleShootData = new float[30];
    std::random_device rd;
    this->randomSeedGenerator.seed(rd());
    this->seedDistribution = std::uniform_real_distribution<float>(0,1);
};



Display::~Display() {
    std::cout<<"Display Deconstructing"<<std::endl;
    if(this->rawGpuImageData)
        cudaFree(this->rawGpuImageData);
    if(this->gpuImage)
        cudaFree(this->gpuImage);
    if(this->gpuCameraData)
        cudaFree(this->gpuCameraData);
    if(this->displayWindow)
        SDL_DestroyWindow(this->displayWindow);
    if(this->imageTexture)
        SDL_DestroyTexture(this->imageTexture);
    
}




void Display::updateDisplayParameters(Settings* settings) {
    this->displayWidth = settings->getGraphicsSettingsValue("ScreenResolutionX");
    this->displayHeight = settings->getGraphicsSettingsValue("ScreenResolutionY");

}



#pragma region GPU MEMORY
void Display::allocateWorldData(World* world) {
    float** worldDataAsArray = world->getWorldObjectDataAsArray();
    this->numOfWorldObjects = world->numOfObjects;
    cudaMalloc(&this->gpuWorldObjectData, this->numOfWorldObjects * sizeof(float*));
    for(size_t i = 0; i < this->numOfWorldObjects; i++) {
        float* objectData;
        cudaMalloc(&objectData, 30 * sizeof(float));
        cudaMemcpy(objectData, worldDataAsArray[i],30*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(&this->gpuWorldObjectData[i],&objectData,sizeof(float*),cudaMemcpyHostToDevice);
    }
}



void Display::allocateGPUMemory() {
    cudaMalloc(&this->gpuTroubleShootData, 30*sizeof(float));
    cudaMalloc(&this->gpuCameraData, 30 * sizeof(float));
    cudaMalloc(&this->rawGpuImageData, this->displayWidth * this->displayHeight * 3 * sizeof(float));
    cudaMalloc(&this->gpuImage, this->displayWidth * this->displayHeight * 3 * sizeof(unsigned char));
}

void Display::reAllocateGPUMemory() {
    cudaFree(this->gpuImage);
    cudaFree(this->rawGpuImageData);
    cudaFree(this->gpuCameraData);
    cudaMalloc(&this->gpuCameraData, 30 * sizeof(float));
    cudaMalloc(&this->rawGpuImageData, this->displayWidth * this->displayHeight * 3 * sizeof(float));
    cudaMalloc(&this->gpuImage, this->displayWidth * this->displayHeight * 3 * sizeof(unsigned char));
}

void Display::updateGPUData() {
    cudaMemcpy(this->gpuCameraData, this->viewCamera->cameraData, 30 * sizeof(float), cudaMemcpyHostToDevice);
}

void Display::copyImageData() {
    cudaMemcpy(this->gpuImage, this->imageData, this->displayWidth * this->displayHeight * 3 * sizeof(unsigned char),cudaMemcpyHostToDevice);
}
#pragma endregion




void Display::getPixelInfo(int x, int y) {
    int index = (y * this->displayWidth + x) * 3;
    std::cout << "Pixel (" << x << ", " << y << "): ( " << static_cast<int>(this->imageData[index]) << ", " << static_cast<int>(this->imageData[index+1]) << ", " << static_cast<int>(this->imageData[index+2]) << " )" << std::endl;
}

void Display::resetCursor() {
    SDL_WarpMouseInWindow(this->displayWindow, this->displayWidth / 2, this->displayHeight / 2);
    this->previousMousePosX = this->displayWidth/2;
    this->previousMousePosY = this->displayHeight/2;
    this->mousePosX = this->displayWidth/2;
    this->mousePosY = this->displayHeight/2;
}

/*
NEXT STEP IS TO REORGANIZE EVENTS INTO THE EVENT MANAGER CLASS
NEEDED QUITE BADLY

*/




void Display::renderScene() {
    this->renderImage();
    SDL_RenderClear(this->imageRenderer);
    SDL_UpdateTexture(this->imageTexture, NULL, this->imageData, this->displayWidth*3);
    SDL_RenderCopy(this->imageRenderer, this->imageTexture, NULL, NULL);
    SDL_RenderPresent(this->imageRenderer);
}


void Display::handleDisplaySizeChange(int newDisplayWidth, int newDisplayHeight) {
    this->displayWidth = newDisplayWidth;
    this->displayHeight = newDisplayHeight;
    delete[] this->imageData;
    this->imageData = new unsigned char[this->displayWidth * this->displayHeight * 3];
    SDL_DestroyTexture(this->imageTexture);
    this->imageTexture = SDL_CreateTexture(this->imageRenderer, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, this->displayWidth, this->displayHeight);
    cudaFree(this->gpuImage);
    cudaMalloc(&this->gpuImage, this->displayWidth * this->displayHeight * 3 * sizeof(unsigned char));

}




void Display::renderImage() {
    dim3 threads(16, 8);
    dim3 blocksPerGrid(ceil(displayWidth / static_cast<float>(threads.x)),
                        ceil(displayHeight / static_cast<float>(threads.y)));

    //Update GPU Data
    if(this->currentFrame == 0){
        this->updateGPUData();
    }
    this->currentFrame += 1;

    int randomGenSeed = static_cast<int>(99999.0*this->seedDistribution(this->randomSeedGenerator));
    

    // Launch the CUDA kernel defined inside the class
    newRender<<<blocksPerGrid, threads>>>(randomGenSeed, this->currentFrame, 10,displayWidth, displayHeight, this->numOfWorldObjects,this->gpuCameraData,this->rawGpuImageData,this->gpuImage, this->gpuWorldObjectData,this->gpuTroubleShootData);

    cudaError_t kernelError = cudaGetLastError();
    if (kernelError != cudaSuccess) {
        std::cout << "CUDA kernel launch error: " << cudaGetErrorString(kernelError) << std::endl;
        delete this; // Free allocated memory
        return;
    }

    // Wait for kernel to finish
    cudaDeviceSynchronize();

    // Copy the processed image back from GPU to CPU
    cudaMemcpy(this->imageData, this->gpuImage, displayWidth * displayHeight * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    
    // cudaMemcpy(this->recieveGpuTroubleShootData,this->gpuTroubleShootData,30*sizeof(float),cudaMemcpyDeviceToHost);
    // std::cout<<std::endl<<"------------"<<std::endl;
    // for(int i = 0; i < 7; i ++) {
    //     std::cout<<"Value "<<i<<": "<<this->recieveGpuTroubleShootData[i]<<std::endl;
    // }
}


void Display::seriousRenderOfImage(int numOfFrames,int maxRayBounces) {
    dim3 threads(16, 8);
    dim3 blocksPerGrid(ceil(displayWidth / static_cast<float>(threads.x)),
                        ceil(displayHeight / static_cast<float>(threads.y)));

    float* imageDataArray = new float[this->displayWidth*this->displayHeight*3];
    float* gpuImageDataArray;
    cudaMalloc(&gpuImageDataArray,this->displayWidth*this->displayHeight*3*sizeof(float));

    for(int frame = 1; frame<=numOfFrames; frame++) {
        std::cout<<"Rendering Frame: "<<frame<<std::endl;
        //Update GPU Data
        this->updateGPUData();
        this->copyImageData();

        int randomGenSeed = static_cast<int>(99999.0*this->seedDistribution(this->randomSeedGenerator));
        
        
        // Launch the CUDA kernel defined inside the class
        DiffuseRender<<<blocksPerGrid, threads>>>(randomGenSeed,frame,maxRayBounces,displayWidth, displayHeight, this->numOfWorldObjects,this->gpuCameraData,gpuImageDataArray, this->gpuWorldObjectData,this->gpuTroubleShootData);
        
        cudaError_t kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            std::cout << "CUDA kernel launch error: " << cudaGetErrorString(kernelError) << std::endl;
            cudaFree(this->gpuImage); // Free allocated memory
            return;
        }

        // Wait for kernel to finish
        cudaDeviceSynchronize();

        updateImage<<<blocksPerGrid, threads>>>(frame,displayWidth, displayHeight,gpuImageDataArray,this->gpuImage,this->gpuTroubleShootData);
        kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            std::cout << "CUDA kernel launch error: " << cudaGetErrorString(kernelError) << std::endl;
            cudaFree(this->gpuImage); // Free allocated memory
            return;
        }

        // Wait for kernel to finish
        cudaDeviceSynchronize();

        // Copy the processed image back from GPU to CPU
        cudaMemcpy(this->imageData, this->gpuImage, displayWidth * displayHeight * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
        //Updates DISPLAYED Image
        SDL_RenderClear(this->imageRenderer);
        SDL_UpdateTexture(this->imageTexture, NULL, this->imageData, this->displayWidth*3);
        SDL_RenderCopy(this->imageRenderer, this->imageTexture, NULL, NULL);
        SDL_RenderPresent(this->imageRenderer);
    }
}



void Display::customResolutionRender(std::string fileName, int numOfFrames,int maxRayBounces, int resolutionX, int resolutionY) {
    dim3 threads(16, 8);
    dim3 blocksPerGrid(ceil(resolutionX / static_cast<float>(threads.x)),
                        ceil(resolutionY / static_cast<float>(threads.y)));

    float* imageDataArray = new float[resolutionX*resolutionY*3];
    float* gpuImageDataArray;
    cudaMalloc(&gpuImageDataArray,resolutionX*resolutionY*3*sizeof(float));

    unsigned char* customResImageData = new unsigned char[resolutionX*resolutionY*3];
    unsigned char* customResGPUImageData;
    cudaMalloc(&customResGPUImageData,resolutionX*resolutionY*3*sizeof(unsigned char));

    cudaError_t kernelError;
    auto startTime = std::chrono::high_resolution_clock::now();
    for(int frame = 1; frame<=numOfFrames; frame++) {
        std::cout<<"Rendering Frame: "<<frame<<std::endl;
        //Update GPU Data
        this->updateGPUData();
        this->copyImageData();

        int randomGenSeed = static_cast<int>(99999.0*this->seedDistribution(this->randomSeedGenerator));
        
        
        // Launch the CUDA kernel defined inside the class
        DiffuseRender2<<<blocksPerGrid, threads>>>(randomGenSeed,frame,maxRayBounces,resolutionX, resolutionY, this->numOfWorldObjects,this->gpuCameraData,gpuImageDataArray, this->gpuWorldObjectData,this->gpuTroubleShootData);
        
        kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            std::cout << "CUDA kernel launch error: " << cudaGetErrorString(kernelError) << std::endl;
            cudaFree(this->gpuImage); // Free allocated memory
            return;
        }

        // Wait for kernel to finish
        cudaDeviceSynchronize();
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    std::cout<<"Time taken: "<<duration.count()<<" seconds"<<std::endl;

    updateImageHDR<<<blocksPerGrid, threads>>>(numOfFrames,resolutionX, resolutionY,gpuImageDataArray,customResGPUImageData,this->gpuTroubleShootData);
        kernelError = cudaGetLastError();
        if (kernelError != cudaSuccess) {
            std::cout << "CUDA kernel launch error: " << cudaGetErrorString(kernelError) << std::endl;
            cudaFree(this->gpuImage); // Free allocated memory
            return;
        }

    // Wait for kernel to finish
    cudaDeviceSynchronize();

    // Copy the processed image back from GPU to CPU
    cudaMemcpy(customResImageData, customResGPUImageData, resolutionX*resolutionY * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);
    std::cout<<"Swapping Red and Blue Values (SDL EXPECATION OF ORDER)"<<std::endl;
    for (int y = 0; y < resolutionY; ++y) {
        for (int x = 0; x < resolutionX; ++x) {
            unsigned char temp = customResImageData[(y * resolutionX + x) * 3];
            customResImageData[(y * resolutionX + x) * 3] = customResImageData[(y * resolutionX + x) * 3 + 2];
            customResImageData[(y * resolutionX + x) * 3 + 2] = temp;
        }
    }
    std::cout<<"saving image..."<<std::endl;
    // Write image data to BMP file
    SDL_Surface* surface = SDL_CreateRGBSurfaceFrom(customResImageData, resolutionX, resolutionY, 24, resolutionX * 3, 0xFF0000, 0x00FF00, 0x0000FF, 0);
    std::string filePath = "./Saved Images/" +fileName + ".bmp";
    std::cout<<"FilePath: "<<filePath<<std::endl;
    SDL_SaveBMP(surface,filePath.c_str());
}


void Display::saveImage(const std::string& fileName) {
    std::cout<<"saving image..."<<std::endl;
    // Write image data to BMP file
    SDL_Surface* surface = SDL_CreateRGBSurface(0, displayWidth, displayHeight, 32, 0, 0, 0, 0);
    SDL_RenderReadPixels(this->imageRenderer,NULL,SDL_PIXELFORMAT_ARGB8888, surface->pixels, surface->pitch);

    std::string filePath = "./Saved Images/" +fileName + ".bmp";
    std::cout<<"FilePath: "<<filePath<<std::endl;
    SDL_SaveBMP(surface,filePath.c_str());
    
}

