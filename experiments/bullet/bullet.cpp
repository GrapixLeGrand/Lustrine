#include <iostream>
#include <cstdlib>
#include <ctime>

#include <string>


#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"

#define GLM_ENABLE_EXPERIMENTAL
#include "Lustrine.hpp"
#include "../Utils.hpp"
#include "LevekGL.hpp"
#include "glm/ext.hpp"
#include <glm/gtx/matrix_decompose.hpp>
#include <glm/gtx/quaternion.hpp>

int resolutionX = 1280;
int resolutionY = 720;


static void addLineBox(Levek::LineRenderer* lineRenderer, const glm::mat4& model, const glm::vec3& half, glm::vec4 color) {

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

int main(void) {

    Levek::RenderingEngine* engine = new Levek::RenderingEngine(resolutionX, resolutionY);

    Levek::Renderer* renderer = engine->getRenderer();
    Levek::LineRenderer* lineRenderer = engine->getLineRenderer();
    Levek::WindowController* windowController = engine->getWindowController();
    Levek::InputController* inputController = engine->getInputController();

    windowController->setWindowTitle("Fluid");
    #ifdef PLATFORM_UNIX
    windowController->initImGui();
    #endif
    Levek::ModelLoader* meshLoader = engine->getModelLoader();
    Levek::Model* model = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/models/billboard.obj");
    const Levek::Mesh* sphere = model->getMesh(0);
    //{20, 20, 45}
    Levek::PerspectiveCamera camera({20, 5, 20}, {5.0, 0.2, 0.2}, {0, 1, 0}, resolutionX, resolutionY);
    glm::mat4 projection = camera.getProjection();

    Levek::Shader shaderInstances = Levek::ShaderFactory::makeFromFile(
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/sphere_inst.vert",
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/sphere_inst.frag"
    );
    
    float particleScale = 1.0f;

    //////////////////////////////////////////////////////////////////////////////////////////
    Lustrine::Simulation simulation;
    Lustrine::SimulationParameters parameters;
    parameters.X = 30.0f;
    parameters.Y = 25.0f;
    parameters.Z = 30.0f;

    bool particles_shown = false;
    bool keysMovingBody = true;
    bool playerLeftGround = true;

    std::vector<Lustrine::Grid> grids (3);
    Lustrine::Grid player_grid;
    std::vector<glm::vec3> grids_positions (3);
    grids_positions[0] = {0, 0, 0};
    grids_positions[1] = {0, 0, 0};
    grids_positions[2] = {15, 0, 15};

    Lustrine::init_grid_from_magika_voxel(&grids[0], LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/models/chr_knight.vox", Lustrine::MaterialType::SOLID_STATIC);
    
    Lustrine::init_grid_box(&parameters, &player_grid, 8, 8, 8, Lustrine::MaterialType::SOLID_STATIC, glm::vec4(1.0, 0.0, 0.0, 1.0));
    Lustrine::init_grid_box(&parameters, &grids[1], 10, 20, 10, Lustrine::MaterialType::FLUID_DYNAMIC, glm::vec4(0.0, 0.2, 1.0, 1.0));
    Lustrine::init_grid_box(&parameters, &grids[2], 10, 20, 10, Lustrine::MaterialType::FLUID_DYNAMIC, glm::vec4(1.0, 0.2, 1.0, 1.0));

    Lustrine::init_simulation(&parameters, &simulation, grids, grids_positions, &player_grid, {15.0f, 17.0f, 15.0f});
    
    Lustrine::BulletPhyicsSimulation* bulletPhysics = &simulation.bullet_physics_simulation;

    //int num_objects = 10;
    
    //btTransform

    
    int box_index = Lustrine::add_box(bulletPhysics, {15, 15, 15}, true, bulletPhysics->collision_group_0, INT32_MAX);
    int box_index_2 = Lustrine::add_box(bulletPhysics, {16, 15, 16}, true);
    int box_index_3 = Lustrine::add_box(bulletPhysics, {14, 15, 14}, true);

    Lustrine::set_body_no_rotation(bulletPhysics, box_index);

    glm::vec3 half_dims_box_4 = {3.0, 1.0, 1.0};
    int box_index_4 = Lustrine::add_box(bulletPhysics, {10, 2, 10}, true, half_dims_box_4, bulletPhysics->collision_group_1, INT32_MAX);

    int ground_index = Lustrine::add_box(bulletPhysics, {parameters.X / 2, -0.5, parameters.Z / 2}, false, {parameters.X / 2, 1, parameters.Z / 2}, bulletPhysics->collision_group_0 | bulletPhysics->collision_group_1, INT32_MAX);

    glm::vec3 ground_dims = {parameters.X / 2, 1, parameters.Z / 2};

    btTransform transformBox1;
    btTransform transformBox2;
    btTransform transformBox3;
    btTransform transformBox4;

    btTransform transformGround;

    glm::mat4 boxModel(0.0);
    glm::mat4 box2Model(0.0);
    glm::mat4 box3Model(0.0);
    glm::mat4 box4Model(0.0);
    glm::mat4 groundModel(0.0);

    //////////////////////////////////////////////////////////////////////////////////////////

    Levek::VertexBuffer particlesPositionsVBO = Levek::VertexBuffer(simulation.positions.data(), simulation.positions.size() * 3 * 4);
    Levek::VertexBuffer particlesColorsVBO = Levek::VertexBuffer(simulation.colors.data(), simulation.colors.size() * 4 * 4);

    Levek::VertexBuffer sphereVBO = Levek::VertexBuffer(sphere);
    Levek::IndexBuffer sphereIBO = Levek::IndexBuffer(sphere);
    Levek::VertexBufferLayout sphereLayout = Levek::VertexBufferLayout();
    Levek::VertexBufferLayout instanceLayout = Levek::VertexBufferLayout(); 
    Levek::VertexBufferLayout colorLayout = Levek::VertexBufferLayout();
    sphereLayout.push<glm::vec3>(1); //sphere position
    sphereLayout.push<glm::vec2>(1); //sphere textures
    sphereLayout.push<glm::vec3>(1); //sphere normal 
    instanceLayout.push<glm::vec3>(1, 1); //instance offset (per instance)
    colorLayout.push<glm::vec4>(1, 1);

    Levek::VertexArray particlesVA;
    particlesVA.addBuffer(sphereVBO, sphereLayout);
    particlesVA.addBuffer(particlesPositionsVBO, instanceLayout);
    particlesVA.addBuffer(particlesColorsVBO, colorLayout);

    SkyBoxPipelineState skybox (getSkyBoxPaths());

    //plan state

    Levek::Shader planeShader = Levek::ShaderFactory::makeFromFile(
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/ground.vert",
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/ground.frag"
    );

    Levek::Model* planeModel = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/plane.obj");
    const Levek::Mesh* planeMesh = planeModel->getMesh(0);
    //Levek::Mesh planeMesh = Levek::makePlane(1.0f);

    Levek::VertexBuffer planeVBO = Levek::VertexBuffer(planeMesh);
    Levek::IndexBuffer planeIBO = Levek::IndexBuffer(planeMesh);
    Levek::VertexBufferLayout planeLayout = Levek::VertexBufferLayout();
    planeLayout.push<glm::vec3>(1);
    planeLayout.push<glm::vec2>(1);
    planeLayout.push<glm::vec3>(1);
    Levek::VertexArray planeVA;
    planeVA.addBuffer(planeVBO, planeLayout);
    
    Levek::Texture unitTexture = Levek::Texture(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/unit.png");
    unitTexture.set(Levek::TextureParameters::TextureWrapMode::REPEAT);
    unitTexture.set(Levek::TextureParameters::TextureLODFunction::LINEAR, Levek::TextureParameters::TextureLODFunction::LINEAR);
    
    //glm::mat4 planeModel = glm::mat4(1.0f);
    float factor = 1.0f;

    ////////////////////////////////////////////////////////////////////////////////
    /*
    Levek::Model* cubeModel = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/bullet/models/cube.dae");
    const Levek::Mesh* cubeMesh = model->getMesh(0);

    Levek::VertexBuffer voxelsPositionsVBO = Levek::VertexBuffer(simulation.bullet_physics_simulation.positions.data(), simulation.bullet_physics_simulation.positions.size() * 3 * 4);
    Levek::VertexBuffer voxelsColorsVBO = Levek::VertexBuffer(simulation.bullet_physics_simulation.colors.data(), simulation.bullet_physics_simulation.colors.size() * 4 * 4);

    Levek::VertexBuffer cubeVBO = Levek::VertexBuffer(cubeMesh);
    Levek::IndexBuffer cubeIBO = Levek::IndexBuffer(cubeMesh);
    Levek::VertexBufferLayout cubeLayout = Levek::VertexBufferLayout();
    Levek::VertexBufferLayout instanceCubeLayout = Levek::VertexBufferLayout(); 
    Levek::VertexBufferLayout colorCubeLayout = Levek::VertexBufferLayout();
    cubeLayout.push<glm::vec3>(1); //sphere position
    cubeLayout.push<glm::vec2>(1); //sphere textures
    cubeLayout.push<glm::vec3>(1); //sphere normal
    instanceCubeLayout.push<glm::vec3>(1, 1); //instance offset (per instance)
    colorCubeLayout.push<glm::vec4>(1, 1);

    Levek::VertexArray voxelsVA;
    voxelsVA.addBuffer(cubeVBO, cubeLayout);
    voxelsVA.addBuffer(voxelsPositionsVBO, instanceCubeLayout);
    voxelsVA.addBuffer(voxelsColorsVBO, colorCubeLayout);*/

    /////////////////////////////////////////////////////////////////////////////////

    while (!windowController->exit() && !inputController->isKeyPressed(Levek::LEVEK_KEY_Q)) {            

        //simulation.time_step = windowController->getDeltaTime();
        //sim here
        //Lustrine::simulate(&simulation, windowController->getDeltaTime());

        Lustrine::simulate_bullet(bulletPhysics, windowController->getDeltaTime());

        particlesPositionsVBO.Update(simulation.positions.data(), simulation.positions.size() * 3 * 4);
        renderer->clear();

        if (keysMovingBody == false) {
            UpdateCameraPositionWASD(inputController, camera, windowController->getDeltaTime(), 10.f);
        } else {
            
            //btVector3 velocity (0.0, 0.0, 0.0);
            float speed = 5.0f;
            glm::vec3 velocity (0.0);

            if (Lustrine::check_collision(bulletPhysics, box_index, ground_index) == true) {
                playerLeftGround = false;
            } else {
                playerLeftGround = true;
            }

            if (playerLeftGround == false) {
                if (inputController->isKeyPressed(Levek::LEVEK_KEY_W) == true) {
                    velocity.x -= speed;
                    Lustrine::set_body_velocity(bulletPhysics, box_index, velocity);
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_S) == true) {
                    velocity.x += speed;;
                    Lustrine::set_body_velocity(bulletPhysics, box_index, velocity);
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_A) == true) {
                    velocity.z += speed;
                    Lustrine::set_body_velocity(bulletPhysics, box_index, velocity);
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_D) == true) {
                    velocity.z -= speed;
                    Lustrine::set_body_velocity(bulletPhysics, box_index, velocity);
                }
            }

            if (playerLeftGround == false && inputController->isKeyPressed(Levek::LEVEK_KEY_X) == true) {
                Lustrine::apply_impulse(bulletPhysics, box_index, {0, 0.2, 0}, {0, 0, 0});
            }

        }

        UpdateCameraWithMouseOnDrag(inputController, camera, 0.2f);

        glm::mat4 view = camera.getView();
        glm::mat4 vp = camera.getProjection() * camera.getView();
        glm::mat3 view_inv = glm::inverse(glm::mat3(view)); //for billboard facing: see https://stackoverflow.com/questions/61559585/how-to-remove-rotation-from-model-view-matrix-so-that-object-always-faces-camera
        glm::vec3 lightDirection = glm::vec3(0, -1, 0); //glm::vec3(glm::normalize(view * glm::vec4(0, -1, 0, 0)));

        //render instances
        shaderInstances.bind();
        shaderInstances.setUniformMat4f("vp", vp);
        shaderInstances.setUniformMat4f("p", camera.getProjection());
        shaderInstances.setUniformMat4f("view", camera.getView());
        shaderInstances.setUniformMat3f("view_inv", view_inv);
        shaderInstances.setUniform3f("light_direction", lightDirection);
        //shaderInstances.setUniform4f("palette", palette.data(), 256);
        shaderInstances.setUniform1f("scale", particleScale);

        vp = camera.getProjection() * glm::mat4(glm::mat3(camera.getView()));
        skybox.draw(renderer, vp);
        
        if (particles_shown) {
            renderer->drawInstances(&particlesVA, &sphereIBO, &shaderInstances, simulation.num_particles);
        }
        //renderer->drawInstances(&voxelsVA, &sphereIBO, &shaderInstances, simulation.num_particles);

        //render plane
        planeShader.bind();
        unitTexture.activateAndBind(0);
        vp = camera.getProjection() * camera.getView();
        planeShader.setUniformMat4f("mvp", vp);
        planeShader.setUniform1i("tex", 0);
        renderer->draw(&planeVA, &planeIBO, &planeShader);
        
        #ifdef PLATFORM_UNIX
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGuiTabBarFlags flags = ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_Reorderable;
        ImGui::Begin("Scene");
        ImGui::BeginTabBar("Scene parameters");
        if (ImGui::BeginTabItem("Simulation")) {
            
            ImGui::Text("%d fps", (int) (1.0f / windowController->getDeltaTime()));
            ImGui::Text("%d fps sim", (int) (1.0f / simulation.time_step));

            ImGui::Text("particle radius %.3f", simulation.particleRadius);
            ImGui::Text("particle diameter %.3f", simulation.particleDiameter);
            ImGui::Text("kernel radius %.3f", simulation.kernelRadius);
            ImGui::Text("#particles: %d", simulation.num_particles);

            //ImGui::InputFloat("kernel factor", &factor, 0.001f, 100.0f, "%.3f");
            ImGui::InputFloat("cubic k", &simulation.cubic_kernel_k, 0.01f, 100.0f, "%.3f");
            ImGui::InputFloat("cubic l", &simulation.cubic_kernel_l, 0.01f, 100.0f, "%.3f");
            ImGui::InputFloat("rest density", &simulation.rest_density, 0.01f, 3000.0f, "%.3f");
            ImGui::InputFloat("mass", &simulation.mass, 0.01f, 1000.0f, "%.3f");
            ImGui::InputFloat("relaxation", &simulation.relaxation_epsilon, 0.01f, 1000.0f, "%.3f");
            ImGui::InputFloat("kernel radius", &simulation.kernelRadius, 0.01f, 1000.0f, "%.3f");
            ImGui::InputFloat("Kernel Factor", &simulation.kernelFactor, 0.001f, 1000.0f, "%.1f");
            ImGui::SliderFloat("Kernel factor", &simulation.kernelFactor, 0.001f, 1.0f, "%.3f");

            ImGui::Text("Negative pressure correction");
            ImGui::InputFloat("s_corr dq", &simulation.s_corr_dq, 0.001f, 100.0f, "%.1f");
            ImGui::InputFloat("s_corr k", &simulation.s_corr_k, 0.001f, 100.0f, "%.1f");
            ImGui::InputFloat("s_corr n", &simulation.s_corr_n, 0.001f, 100.0f, "%.1f");

            ImGui::Text("Grid");
            //float old_cell_size = simulation.cell_size;
            //ImGui::InputFloat("cell size", &simulation.cell_size, 0.001f, 100.0f, "%.1f");

            ImGui::Text("# grid cells: %d", simulation.num_grid_cells);



            if (ImGui::Button("reset")) {
                //Lustrine::init_grid_box(&simulation, &grids[0], 20, 30, 20);
                Lustrine::init_simulation(&parameters, &simulation, grids, grids_positions);
            }
            ImGui::EndTabItem();
            simulation.cubic_kernel_k *= factor;
            simulation.cubic_kernel_l *= factor;
        }
        if (ImGui::BeginTabItem("Particles")) {
            ImGui::InputFloat("particle size", &particleScale, 0.01f, 5.0f, "%.3f");
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Camera")) {
            addImGuiVec3(camera.getEye());
            addImGuiVec3(camera.getFront());
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        #endif

        lineRenderer->SetViewProjection(projection * camera.getView());
        //lineRenderer->AddLine({0, 0, 0}, {1, 0, 0}, {1.0, 0.0, 0.0, 1.0}); 
        //lineRenderer->AddLine({0, 0, 0}, {0, 1, 0}, {0.0, 1.0, 0.0, 1.0});
        //lineRenderer->AddLine({0, 0, 0}, {0, 0, 1}, {0.0, 0.0, 1.0, 1.0});
        glm::vec3 p = glm::vec3{15, 15, 15};
        glm::vec3 dims = glm::vec3{0.5, 0.5, 0.5};

        transformBox1 = bulletPhysics->rigidbodies[box_index]->getWorldTransform();
        bulletPhysics->rigidbodies[box_index]->getMotionState()->getWorldTransform(transformBox1);
        transformBox1.getOpenGLMatrix(glm::value_ptr(boxModel));

        transformBox2 = bulletPhysics->rigidbodies[box_index_2]->getWorldTransform();
        bulletPhysics->rigidbodies[box_index_2]->getMotionState()->getWorldTransform(transformBox2);
        transformBox2.getOpenGLMatrix(glm::value_ptr(box2Model));

        transformBox3 = bulletPhysics->rigidbodies[box_index_3]->getWorldTransform();
        bulletPhysics->rigidbodies[box_index_3]->getMotionState()->getWorldTransform(transformBox3);
        transformBox3.getOpenGLMatrix(glm::value_ptr(box3Model));

        transformBox4 = bulletPhysics->rigidbodies[box_index_4]->getWorldTransform();
        bulletPhysics->rigidbodies[box_index_4]->getMotionState()->getWorldTransform(transformBox4);
        transformBox4.getOpenGLMatrix(glm::value_ptr(box4Model));

        transformGround = bulletPhysics->rigidbodies[ground_index]->getWorldTransform();
        bulletPhysics->rigidbodies[ground_index]->getMotionState()->getWorldTransform(transformGround);
        transformGround.getOpenGLMatrix(glm::value_ptr(groundModel));

        //Levek::printMat4(boxModel);

        addLineBox(lineRenderer, boxModel, dims, {1.0, 0.0, 0.0, 1.0});
        addLineBox(lineRenderer, box2Model, dims, {0.0, 1.0, 0.0, 1.0});
        addLineBox(lineRenderer, box3Model, dims, {0.0, 0.0, 1.0, 1.0});
        addLineBox(lineRenderer, box4Model, half_dims_box_4, {0.0, 1.0, 1.0, 1.0});

        addLineBox(lineRenderer, groundModel, ground_dims, {1.0, 1.0, 1.0, 1.0});
    
        lineRenderer->Draw();

        inputController->poll();
        windowController->swapBuffers();
    }


    Lustrine::clean_simulation(&simulation);
    delete engine;

    return 0;
}