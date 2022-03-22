#include "Utils.hpp"
#include "glm/ext.hpp"
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/quaternion.hpp>

//from https://learnopengl.com/code_viewer.php?code=advanced/cubemaps_skybox_data
    float skyboxVertices[] = {
        // positions          
        -1.0f,  1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f, -1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,

        -1.0f, -1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f, -1.0f,  1.0f,
        -1.0f, -1.0f,  1.0f,

        -1.0f,  1.0f, -1.0f,
        1.0f,  1.0f, -1.0f,
        1.0f,  1.0f,  1.0f,
        1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f,  1.0f,
        -1.0f,  1.0f, -1.0f,

        -1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
        1.0f, -1.0f, -1.0f,
        1.0f, -1.0f, -1.0f,
        -1.0f, -1.0f,  1.0f,
        1.0f, -1.0f,  1.0f
    };
    
float lastMouseValueX = 0;
float lastMouseValueY = 0;
bool wasLeftButtonPressed = false;


void addLineBox(Levek::LineRenderer* lineRenderer, const glm::mat4& model, const glm::vec3& half, glm::vec4 color) {

    glm::vec3 position = glm::vec4(0.0); //glm::vec3(model[3]);
    glm::vec3 v1 = position;
    v1.x -= half.x;
    v1.y += half.y;
    v1.z -= half.z;

    glm::vec3 v2 = position;
    v2.x -= half.x;
    v2.y += half.y;
    v2.z += half.z;

    glm::vec3 v3 = position;
    v3.x -= half.x;
    v3.y -= half.y;
    v3.z += half.z;

    glm::vec3 v4 = position;
    v4.x -= half.x;
    v4.y -= half.y;
    v4.z -= half.z;

    glm::vec3 v5 = position;
    v5.x += half.x;
    v5.y -= half.y;
    v5.z -= half.z;

    glm::vec3 v6 = position;
    v6.x += half.x;
    v6.y += half.y;
    v6.z -= half.z;

    glm::vec3 v7 = position;
    v7.x += half.x;
    v7.y += half.y;
    v7.z += half.z;

    glm::vec3 v8 = position;
    v8.x += half.x;
    v8.y -= half.y;
    v8.z += half.z;

 // your transformation matrix.
    glm::vec3 scale;
    glm::quat rotation;
    glm::vec3 translation;
    glm::vec3 skew;
    glm::vec4 perspective;
    glm::decompose(model, scale, rotation, translation, skew, perspective);
    //rotation = glm::conjugate(rotation);

    /*glm::vec4 v = model * v1;
    v2 = model * v2;
    v3 = model * v3;
    v4 = model * v4;
    v5 = model * v5;
    v6 = model * v6;
    v7 = model * v7;
    v8 = model * v8;*/
    
    glm::mat3 rot = glm::toMat3(rotation);//glm::mat3_cast(rotation);
    
    v1 = rot * v1;
    v2 = rot * v2;
    v3 = rot * v3;
    v4 = rot * v4;
    v5 = rot * v5;
    v6 = rot * v6;
    v7 = rot * v7;
    v8 = rot * v8;
    /*
    v1 *= rot;
    v2 *= rot;
    v3 *= rot;
    v4 *= rot;
    v5 *= rot;
    v6 *= rot;
    v7  *= rot;
*/
    v1 += translation;
    v2 += translation;
    v3 += translation;
    v4 += translation;
    v5 += translation;
    v6 += translation;
    v7 += translation;
    v8 += translation;
    //glm::vec3 v11 = glm::vec3(v[0], v[1], v[2]);
    //Levek::printVec3(v11);

    lineRenderer->AddLine(v1, v2, color);
    lineRenderer->AddLine(v2, v3, color);
    lineRenderer->AddLine(v3, v4, color);
    lineRenderer->AddLine(v4, v1, color);

    lineRenderer->AddLine(v1, v6, color);
    lineRenderer->AddLine(v4, v5, color);
    lineRenderer->AddLine(v3, v8, color);
    lineRenderer->AddLine(v2, v7, color);

    lineRenderer->AddLine(v6, v7, color);
    lineRenderer->AddLine(v7, v8, color);
    lineRenderer->AddLine(v8, v5, color);
    lineRenderer->AddLine(v5, v6, color);

};