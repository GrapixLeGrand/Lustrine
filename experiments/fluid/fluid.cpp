#include <iostream>
#include <cstdlib>
#include <ctime>

#include <string>

//#ifdef PLATFORM_UNIX
//#include "GL/glew.h"

//#endif

#include "Lustrine.hpp"
#include "../Utils.hpp"
#include "LevekGL.hpp"

#include "GL/glew.h"
#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "imgui_impl_glfw.h"


int resolutionX = 1280;
int resolutionY = 720;

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
    Levek::PerspectiveCamera camera({20, 20, 45}, {0.2, 0.2, 0.2}, {0, 1, 0}, resolutionX, resolutionY);
    glm::mat4 projection = camera.getProjection();
    float particleScale = 1.0f;

    //////////////////////////////////////////////////////////////////////////////////////////
    Lustrine::Simulation simulation;
    Lustrine::SimulationParameters parameters;
    parameters.X = 40.0f;
    parameters.Y = 30.0f;
    parameters.Z = 40.0f;
    parameters.particleRadius = 0.5f;
    parameters.particleDiameter = 1.0f;

    std::vector<Lustrine::Grid> sand_grids (0);
    std::vector<glm::vec3> sand_grids_positions (4);
    sand_grids_positions[0] = {0, 0, 0};
    sand_grids_positions[1] = {20, 0, 20};
    sand_grids_positions[2] = {0, 0, 20};
    sand_grids_positions[3] = {20, 0, 0};
    
    std::vector<Lustrine::Grid> solid_grids_tmp (0);
    std::vector<Lustrine::Grid> solid_grids (1);
    std::vector<glm::vec3> solid_grids_positions (1);
    solid_grids_positions[0] = {10, 1, 10};

    Lustrine::init_grid_from_magika_voxel(&solid_grids[0], LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_physical.vox", { 10, 0, 10 }, Lustrine::MaterialType::SOLID);
    
    //Lustrine::init_grid_box(&parameters, &sand_grids[0], 20, 20, 20, glm::vec3(0, 0, 0), glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    //Lustrine::init_grid_box(&parameters, &sand_grids[1], 20, 20, 20, glm::vec3(20, 0, 20), glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    //Lustrine::init_grid_box(&parameters, &sand_grids[2], 20, 20, 20, glm::vec3(20, 0, 0), glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    //Lustrine::init_grid_box(&parameters, &sand_grids[3], 20, 20, 20, glm::vec3(0, 0, 20), glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    //Lustrine::init_grid_box(&parameters, &sand_grids[4], 20, 20, 20, glm::vec3(20, 0, 20), glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);


    Lustrine::init_simulation(
        &parameters,
        &simulation,
        sand_grids,
        solid_grids
    );

    ParticlesPipelineSate sandParticlesPipeline(engine, simulation.positions, simulation.colors, simulation.num_sand_particles);
    ParticlesPipelineSate solidParticlesPipeline(engine, simulation.positions_solid, simulation.colors_solid, simulation.num_solid_particles);
    SkyBoxPipelineState skybox (getSkyBoxPaths());
    //WARINING BUG IN TEXTURE INIT
    GroundPipelineState groundPipelineState(engine);

    //glm::mat4 planeModel = glm::mat4(1.0f);
    float factor = 1.0f;

    Lustrine::Grid test_grid;
    //init_grid_box_random(&parameters, &test_grid, 30, 1, 30, {15, 35, 15}, glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND, 0.2f);
    Lustrine::init_grid_box(&parameters, &test_grid, 20, 1, 20, {25, 30, 25}, glm::vec4(1.0, 0.2, 1.0, 1.0), Lustrine::MaterialType::SAND);
    Lustrine::add_particle_source(&simulation, &test_grid, {0, -1, 0}, 1.0f / 16.f, -1);
    Lustrine::add_particle_sink(&simulation, {0, 0, 0}, {60, 5, 60});

    while (!windowController->exit() && !inputController->isKeyPressed(Levek::LEVEK_KEY_Q)) {

        //sim here
        Lustrine::simulate(&simulation, windowController->getDeltaTime());
        //std::cout << Lustrine::query_cell_num_particles(&simulation, glm::vec3(0.0), glm::vec3(10), false) << " particles detected :)\n";
        //std::cout << simulation.ptr_sand_end << std::endl;
        sandParticlesPipeline.updatePositions(simulation.positions, simulation.num_sand_particles);
        sandParticlesPipeline.updateColors(simulation.colors, simulation.num_sand_particles);
        solidParticlesPipeline.updatePositions(simulation.positions_solid, simulation.num_solid_particles);

        renderer->clear();

        UpdateCameraPositionWASD(inputController, camera, windowController->getDeltaTime(), 10.0f);
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

        solidParticlesPipeline.setUniforms(
            vp,
            camera.getProjection(),
            camera.getView(),
            view_inv,
            lightDirection,
            particleScale
        );

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
                //Lustrine::init_simulation(&parameters, &simulation, grids, grids_positions);
                std::cout << "no reset for now" << std::endl;
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

        lineRenderer->SetViewProjection(projection * camera.getView());
        //lineRenderer->AddLine({0, 0, 0}, {1, 0, 0}, {1.0, 0.0, 0.0, 1.0}); 
        //lineRenderer->AddLine({0, 0, 0}, {0, 1, 0}, {0.0, 1.0, 0.0, 1.0});
        //lineRenderer->AddLine({0, 0, 0}, {0, 0, 1}, {0.0, 0.0, 1.0, 1.0});

        //lineRenderer->Draw();
        inputController->poll();
        windowController->swapBuffers();
    }

    delete engine;
    return 0;
}