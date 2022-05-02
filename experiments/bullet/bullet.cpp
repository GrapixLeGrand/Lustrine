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



/**
 * @brief PipelineState that works with the LineRenderer in order to draw lines.
 * 
 */
struct LineBoxPipelineState {

    std::vector<int> body_indices;
	std::vector<glm::vec3> box_half_dims;//half dimensions of the box
	std::vector<glm::vec4> color;//color of the box
	std::vector<glm::mat4> model;//model matrix
    
    void add_body(Lustrine::Simulation* simulation, Lustrine::Grid grid) {}

	void update() {}

	void draw(Levek::LineRenderer* renderer) { /*addLineBox(renderer, model, box_half_dims, color);*/ }

};

int main(void) {

    Levek::RenderingEngine* engine = new Levek::RenderingEngine(resolutionX, resolutionY);

    Levek::Renderer* renderer = engine->getRenderer();
    Levek::LineRenderer* lineRenderer = engine->getLineRenderer();
    Levek::WindowController* windowController = engine->getWindowController();
    Levek::InputController* inputController = engine->getInputController();

    windowController->setWindowTitle("Fluid");
    //#ifdef PLATFORM_UNIX
    windowController->initImGui();
    //#endif
    Levek::ModelLoader* meshLoader = engine->getModelLoader();
    Levek::Model* model = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/models/billboard.obj");
    const Levek::Mesh* sphere = model->getMesh(0);
    //{20, 20, 45}
    Levek::PerspectiveCamera camera({20, 30, 45}, {0.2, 0.2, 0.2}, {0, 1, 0}, resolutionX, resolutionY);
    glm::mat4 projection = camera.getProjection();
    
    float particleScale = 1.0f;

    bool particles_shown = false;
    bool keysMovingBody = true;
    bool playerLeftGround = true;
    bool areParticlesDisabled = false;

    //////////////////////////////////////////////////////////////////////////////////////////
    
    Lustrine::Simulation simulation;
    Lustrine::SimulationParameters parameters;
    parameters.X = 30.0f;
    parameters.Y = 25.0f;
    parameters.Z = 30.0f;

    int subdivision = 1;
    parameters.particleRadius = 0.5f / subdivision;
    parameters.particleDiameter = 2.0f * parameters.particleRadius;

    std::vector<Lustrine::Grid> sand_grids (1);
    std::vector<glm::vec3> sand_grids_positions (1);
    sand_grids_positions[0] = {20.0f, 5.0f, 0.0f};
    //sand_grids_positions[1] = {15, 0, 15};

    std::vector<Lustrine::Grid> solid_grids (1);
    std::vector<glm::vec3> solid_grids_positions (1);
    solid_grids_positions[0] = {0, 0, 0};

    Lustrine::init_grid_from_magika_voxel(&solid_grids[0], LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_10x10.vox", { 0, 0, 0 }, Lustrine::MaterialType::SOLID);
    
    Lustrine::init_grid_box(&parameters, &sand_grids[0], 5, 10, 30, { 20.0f, 5.0f, 0.0f }, glm::vec4(0.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    //Lustrine::init_grid_box(&parameters, &sand_grids[1], 10, 20, 10, Lustrine::MaterialType::SAND, glm::vec4(1.0, 0.2, 1.0, 1.0));
    //Anemoiapolis: I need to test this game once it comes out
    Lustrine::init_simulation(
        &parameters,
        &simulation,
        sand_grids,
        solid_grids,
        subdivision
    );

    Lustrine::Bullet::Simulation* bulletPhysics = &simulation.bullet_physics_simulation;
    
    glm::vec3 half_dims_box_4 = {3.0, 1.0, 1.0};
    glm::vec3 ground_dims = {parameters.X / 2, 1, parameters.Z / 2};
    float speed = 800.0f;
    float basic_impulse = 100.0f;
    glm::vec3 gravity = glm::vec3(0.0f, -25.0f, 0.0f);
    float gravity_y = -25.0f;

    
    
    int box_index = Lustrine::Bullet::add_capsule(bulletPhysics, {15, 30, 15}, 2.0f, 3.0f); //Lustrine::Bullet::add_box(bulletPhysics, {15, 15, 15}, true);

    int box_index_2 = Lustrine::Bullet::add_box(bulletPhysics, {16, 15, 16}, true);
    int box_index_3 = Lustrine::Bullet::add_box(bulletPhysics, {14, 15, 14}, true);
    int box_index_4 = Lustrine::Bullet::add_box(bulletPhysics, {10, 2, 10}, true, half_dims_box_4); //, bulletPhysics->collision_group_1, INT32_MAX);
    int ground_index = Lustrine::Bullet::add_box(bulletPhysics, {parameters.X / 2, -0.5, parameters.Z / 2}, false, {parameters.X / 2, 1, parameters.Z / 2}); //, bulletPhysics->collision_group_0 | bulletPhysics->collision_group_1, INT32_MAX);

    int detector_block = Lustrine::Bullet::add_detector_block(bulletPhysics, {10, 10, 10}, {2, 2, 2});

    Lustrine::Bullet::set_body_no_rotation(bulletPhysics, box_index);
    

    //bulletPhysics->rigidbodies[box_index]->setActivationState(DISABLE_DEACTIVATION);

    btTransform transformBox1;
    btTransform transformBox2;
    btTransform transformBox3;
    btTransform transformBox4;

    btTransform transformDetector;

    btTransform transformGround;

    glm::mat4 boxModel(0.0);
    glm::mat4 box2Model(0.0);
    glm::mat4 box3Model(0.0);
    glm::mat4 box4Model(0.0);
    glm::mat4 groundModel(0.0);
    glm::mat4 detectorModel(0.0);

    //////////////////////////////////////////////////////////////////////////////////////////
    ParticlesPipelineSate sandParticlesPipeline(engine, simulation.positions, simulation.colors, simulation.num_sand_particles);
    ParticlesPipelineSate solidParticlesPipeline(engine, simulation.positions_solid, simulation.colors_solid, simulation.num_solid_particles);
    SkyBoxPipelineState skybox (getSkyBoxPaths());
    //WARINING BUG IN TEXTURE INIT
    GroundPipelineState groundPipelineState(engine);

    float factor = 1.0f;

    while (!windowController->exit() && !inputController->isKeyPressed(Levek::LEVEK_KEY_Q)) {            

        //simulation.time_step = windowController->getDeltaTime();
        //sim here

        simulation.bullet_physics_simulation.player_position = Lustrine::Bullet::get_body_position(&simulation.bullet_physics_simulation, box_index);
        Lustrine::simulate(&simulation, windowController->getDeltaTime());
        //if (areParticlesDisabled == false) 
        //Lustrine::Bullet::set_particles_box_colliders_positions(bulletPhysics, simulation.positions);

        /*if (Lustrine::Bullet::do_collide(bulletPhysics, box_index)) {
            std::cout << "colliding" << std::endl;
        }
        else {
            std::cout << "not colliding" << std::endl;
        }*/

        if (Lustrine::Bullet::check_collision(bulletPhysics, box_index, detector_block) == true) {
            std::cout << "detection!!!!" << std::endl;
        }

        sandParticlesPipeline.updatePositions(simulation.positions, simulation.num_sand_particles);
        renderer->clear();

        if (keysMovingBody == false) {
            UpdateCameraPositionWASD(inputController, camera, windowController->getDeltaTime(), 10.f);
        } else {
            
            glm::vec3 velocity (0.0);
            bool playerLeftGround = true;
            for (int i = 0; i < bulletPhysics->num_bodies; i++) {
                if (Lustrine::Bullet::check_collision(bulletPhysics, box_index, i) == true) {
                    playerLeftGround = false;
                    break;
                } else {
                    playerLeftGround = true;
                }
            }

            /*if (Lustrine::Bullet::check_collision(bulletPhysics, box_index, ground_index) == true) {
                playerLeftGround = false;
            } else {
                playerLeftGround = true;
            }*/

            if (playerLeftGround == false) {
                if (inputController->isKeyPressed(Levek::LEVEK_KEY_W) == true) {
                    velocity.x -= speed;
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_S) == true) {
                    velocity.x += speed;;
                    //Lustrine::Bullet::add_body_velocity(bulletPhysics, box_index, velocity);
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_A) == true) {
                    velocity.z += speed;
                    //Lustrine::Bullet::add_body_velocity(bulletPhysics, box_index, velocity);
                }

                if (inputController->isKeyPressed(Levek::LEVEK_KEY_D) == true) {
                    velocity.z -= speed;
                    //Lustrine::Bullet::add_body_velocity(bulletPhysics, box_index, velocity);
                }
                float l = glm::length(velocity); 
                if (l > 0.0f) {
                    velocity /= l;
                    velocity *= (speed * windowController->getDeltaTime());
                }
                Lustrine::Bullet::add_body_velocity(bulletPhysics, box_index, velocity);
            }
            //playerLeftGround == false &&
            if (playerLeftGround == false && inputController->isKeyPressed(Levek::LEVEK_KEY_X) == true) {
                float impulse_mag = windowController->getDeltaTime() * basic_impulse;
                Lustrine::Bullet::apply_impulse(bulletPhysics, box_index, {0, impulse_mag, 0}, {0, 0, 0});
            }

        }

        UpdateCameraWithMouseOnDrag(inputController, camera, 0.2f);

        glm::mat4 view = camera.getView();
        glm::mat4 vp = camera.getProjection() * camera.getView();
        glm::mat3 view_inv = glm::inverse(glm::mat3(view)); //for billboard facing: see https://stackoverflow.com/questions/61559585/how-to-remove-rotation-from-model-view-matrix-so-that-object-always-faces-camera
        glm::vec3 lightDirection = glm::vec3(0, -1, 0); //glm::vec3(glm::normalize(view * glm::vec4(0, -1, 0, 0)));

        sandParticlesPipeline.setUniforms(
            vp,
            camera.getProjection(),
            camera.getView(),
            view_inv,
            lightDirection,
            particleScale
        );

        vp = camera.getProjection() * glm::mat4(glm::mat3(camera.getView()));
        skybox.draw(renderer, vp);
        
        vp = camera.getProjection() * glm::mat4(glm::mat3(camera.getView()));
        skybox.draw(renderer, vp);
        sandParticlesPipeline.draw(renderer);
        solidParticlesPipeline.draw(renderer);

        vp = camera.getProjection() * camera.getView();
        groundPipelineState.setUniforms(vp);
        groundPipelineState.draw(renderer);

        
        //#ifdef PLATFORM_UNIX
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGuiTabBarFlags flags = ImGuiTabBarFlags_AutoSelectNewTabs | ImGuiTabBarFlags_Reorderable;
        ImGui::Begin("Scene");
        ImGui::BeginTabBar("Scene parameters");
        if (ImGui::BeginTabItem("Simulation")) {

            if (ImGui::Button("particles bounding box")) {
                if (bulletPhysics->particles_bounding_box_requested_state == true) {
                    std::cout << "particles enabled" << std::endl;
                    bulletPhysics->particles_bounding_box_requested_state = false;
                    //Lustrine::Bullet::enable_particles_bounding_boxes(&simulation.bullet_physics_simulation);
                } else {
                    std::cout << "particles disabled" << std::endl;
                    bulletPhysics->particles_bounding_box_requested_state = true;
                    //Lustrine::Bullet::disable_particles_bounding_boxes(&simulation.bullet_physics_simulation);
                }            
            }

            ImGui::InputFloat("speed:", &speed, 1.0f, 10000.0f, "%.3f");
            ImGui::InputFloat("gravity_y:", &gravity_y, 1.0f, 100.0f, "%.3f");

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
                Lustrine::Bullet::set_body_position(&simulation.bullet_physics_simulation, box_index, {15, 15, 15});
                Lustrine::Bullet::set_body_velocity(&simulation.bullet_physics_simulation, box_index, {0, 0, 0});
                //Lustrine::init_grid_box(&simulation, &grids[0], 20, 30, 20);
                //Lustrine::init_simulation(&parameters, &simulation, grids, grids_positions);
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
        //#endif

        gravity.y = gravity_y;
        Lustrine::Bullet::set_gravity(&simulation.bullet_physics_simulation, gravity);

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


        transformDetector = bulletPhysics->rigidbodies[detector_block]->getWorldTransform();
        bulletPhysics->rigidbodies[detector_block]->getMotionState()->getWorldTransform(transformDetector);
        transformDetector.getOpenGLMatrix(glm::value_ptr(detectorModel));

        //Levek::printMat4(boxModel);
        glm::vec3 capsule_dims = glm::vec3(1, 1.5, 1);
        addLineBox(lineRenderer, boxModel, capsule_dims, {1.0, 0.0, 0.0, 1.0});
        addLineBox(lineRenderer, box2Model, dims, {0.0, 1.0, 0.0, 1.0});
        addLineBox(lineRenderer, box3Model, dims, {0.0, 0.0, 1.0, 1.0});
        addLineBox(lineRenderer, box4Model, half_dims_box_4, {0.0, 1.0, 1.0, 1.0});

        addLineBox(lineRenderer, groundModel, ground_dims, {1.0, 1.0, 1.0, 1.0});

        addLineBox(lineRenderer, detectorModel, {2, 2, 2}, {1.0, 1.0, 1.0, 1.0});
    
        lineRenderer->Draw();

        inputController->poll();
        windowController->swapBuffers();
    }


    Lustrine::clean_simulation(&simulation);
    delete engine;

    return 0;
}