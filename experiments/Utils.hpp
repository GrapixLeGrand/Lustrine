
#include "LevekGL.hpp"
#include <vector>
#include <string>

struct ParticlesPipelineSate {

	Levek::Model* billboardModel = nullptr;
	const Levek::Mesh* billboardMesh = nullptr;

	Levek::VertexBuffer* particlesPositionsVBO = nullptr;// = Levek::VertexBuffer((void*) simulation.positions, (size_t) simulation.num_sand_particles * 3 * 4);
    Levek::VertexBuffer* particlesColorsVBO = nullptr;// = Levek::VertexBuffer((void*) simulation.colors, (size_t) simulation.num_sand_particles * 4 * 4);

    Levek::VertexBuffer* sphereVBO = nullptr; // = Levek::VertexBuffer(sphere);
    Levek::IndexBuffer* sphereIBO = nullptr; // = Levek::IndexBuffer(sphere);
    Levek::VertexBufferLayout sphereLayout = Levek::VertexBufferLayout();
    Levek::VertexBufferLayout instanceLayout = Levek::VertexBufferLayout(); 
    Levek::VertexBufferLayout colorLayout = Levek::VertexBufferLayout();
	Levek::VertexArray* particlesVA;
	Levek::Shader shaderInstances = Levek::ShaderFactory::makeFromFile(
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/sphere_inst.vert",
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/sphere_inst.frag"
    );

	size_t size = 0;

	/**
	 * @brief Construct a new Particles Pipeline Sate object
	 * 
	 * @param positions 
	 * @param colors 
	 * @param size in num positions (NOT BYTES!)
	 */
	ParticlesPipelineSate(Levek::RenderingEngine* engine, const glm::vec3* positions, const glm::vec4* colors, size_t size_arg) {

		size = size_arg;
		particlesPositionsVBO = new Levek::VertexBuffer((const void*) positions, size * 3 * 4);
		particlesColorsVBO = new Levek::VertexBuffer((const void*) colors, (size_t) size * 4 * 4);
		Levek::ModelLoader* meshLoader = engine->getModelLoader();
    	billboardModel = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/models/billboard.obj");
    	assert(billboardModel != nullptr);
		billboardMesh = billboardModel->getMesh(0);
		assert(billboardMesh != nullptr);
		sphereVBO = new Levek::VertexBuffer(billboardMesh);
		sphereIBO = new Levek::IndexBuffer(billboardMesh);

		sphereLayout.push<glm::vec3>(1); //sphere position
    	sphereLayout.push<glm::vec2>(1); //sphere textures
    	sphereLayout.push<glm::vec3>(1); //sphere normal 
    	instanceLayout.push<glm::vec3>(1, 1); //instance offset (per instance)
    	colorLayout.push<glm::vec4>(1, 1);

		particlesVA = new Levek::VertexArray();
    	particlesVA->addBuffer(sphereVBO, &sphereLayout);
    	particlesVA->addBuffer(particlesPositionsVBO, &instanceLayout);
    	particlesVA->addBuffer(particlesColorsVBO, &colorLayout);

	}

	~ParticlesPipelineSate() {
		delete particlesPositionsVBO;
		delete particlesColorsVBO;
		delete billboardModel;
		delete billboardMesh;
		delete sphereVBO;
		delete sphereIBO;
		delete particlesVA;
	}
	
	/**
	 * @brief 
	 * 
	 * @param positions 
	 * @param new_size in num of particles (NOT bytes)
	 */
	void updatePositions(const glm::vec3* positions, size_t new_size) {
		particlesPositionsVBO->Update((const void*) positions, (size_t) new_size * 3 * 4);
	}

	void draw(Levek::Renderer* renderer) {
		renderer->drawInstances(particlesVA, sphereIBO, &shaderInstances, size);
	}

	void setUniforms(
		const glm::mat4& vp,
		const glm::mat4& p,
		const glm::mat4& v,
		const glm::mat4& v_inv,
		const glm::vec3& light_direction,
		float particle_scale 
	) {
		shaderInstances.bind();
        shaderInstances.setUniformMat4f("vp", vp);
        shaderInstances.setUniformMat4f("p", p);
        shaderInstances.setUniformMat4f("view", v);
        shaderInstances.setUniformMat3f("view_inv", v_inv);
        shaderInstances.setUniform3f("light_direction", light_direction);
        shaderInstances.setUniform1f("scale", particle_scale);
	}

};

struct GroundPipelineState {

	Levek::Shader planeShader = Levek::ShaderFactory::makeFromFile(
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/ground.vert",
        LUSTRINE_EXPERIMENTS_DIRECTORY"/fluid/shaders/ground.frag"
    );

    Levek::Model* planeModel;
    const Levek::Mesh* planeMesh;
    Levek::VertexBuffer* planeVBO;
    Levek::IndexBuffer* planeIBO;
    Levek::VertexBufferLayout planeLayout;
    Levek::VertexArray* planeVA;
    Levek::Texture* unitTexture;

	GroundPipelineState(Levek::RenderingEngine* engine): unitTexture(new Levek::Texture(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/unit.png")) {

		Levek::ModelLoader* meshLoader = engine->getModelLoader();
		planeModel = meshLoader->loadFromFile(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/plane.obj");
		planeMesh = planeModel->getMesh(0);

		planeVBO = new Levek::VertexBuffer(planeMesh);
    	planeIBO = new Levek::IndexBuffer(planeMesh);
		planeIBO->bind();
		planeIBO->unbind();
		planeLayout.push<glm::vec3>(1);
    	planeLayout.push<glm::vec2>(1);
    	planeLayout.push<glm::vec3>(1);

		planeVA = new Levek::VertexArray();
		planeVA->addBuffer(planeVBO, &planeLayout);

		//Levek::Texture* unitTexture = new Levek::Texture(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/unit.png");
    	unitTexture->set(Levek::TextureParameters::TextureWrapMode::REPEAT);
    	unitTexture->set(Levek::TextureParameters::TextureLODFunction::LINEAR, Levek::TextureParameters::TextureLODFunction::LINEAR);
	}

	void setUniforms(const glm::mat4& vp) {
		planeShader.bind();
		unitTexture->activateAndBind(0);
		planeShader.setUniformMat4f("mvp", vp);
        planeShader.setUniform1i("tex", 0);
	}

	void draw(Levek::Renderer* renderer) {
		renderer->draw(planeVA, planeIBO, &planeShader);
	}

};

//from https://learnopengl.com/code_viewer.php?code=advanced/cubemaps_skybox_data
extern float skyboxVertices[108];

struct SkyBoxPipelineState {

	Levek::CubeMap* cubeMap;
    Levek::VertexBuffer* cubeMapVbo; //= Levek::VertexBuffer(skyboxVertices, sizeof(skyboxVertices));
    Levek::VertexBufferLayout* cubeMapLayout;// = Levek::VertexBufferLayout(); cubeMapLayout.push<glm::vec3>(1);
    Levek::VertexArray* cubeMapVa;

    Levek::Shader cubeMapShader = Levek::ShaderFactory::makeFromFile(
		LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/shaders/skybox.vert",
		LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/shaders/skybox.frag"
	);

	SkyBoxPipelineState(const std::vector<std::string>& imagePaths) {
		cubeMap = new Levek::CubeMap(imagePaths, 2048, 2048);
		cubeMapVbo = new Levek::VertexBuffer(skyboxVertices, sizeof(skyboxVertices));
    	cubeMapLayout = new Levek::VertexBufferLayout(); cubeMapLayout->push<glm::vec3>(1);
		cubeMapVa = new Levek::VertexArray();
		cubeMapVa->addBuffer(cubeMapVbo, cubeMapLayout);
	}

	void draw(Levek::Renderer *renderer, glm::mat4& vp) {
		
		renderer->setDepthMask(false);
        cubeMapShader.bind();
        cubeMap->bind();
        cubeMapShader.setUniform1i("skybox", 0);
        cubeMapShader.setUniformMat4f("vp", vp);
        cubeMapShader.unbind();

        renderer->draw(cubeMapVa, &cubeMapShader);

        renderer->setDepthMask(true);
	}

	~SkyBoxPipelineState() {
		delete cubeMap;
		delete cubeMapVbo;
		delete cubeMapLayout;
		delete cubeMapVa;
	}
};

const std::vector<std::string> static getSkyBoxPaths() {

	std::vector<std::string> skyBoxImagesPaths { "right.jpg", "left.jpg", "top.jpg", "bottom.jpg", "front.jpg", "back.jpg" };
    for (auto it = skyBoxImagesPaths.begin(); it != skyBoxImagesPaths.end(); it++) {
        (*it) = LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/skybox/" + (*it);
    }

	return skyBoxImagesPaths;
}

void static UpdateCameraPositionWASD(Levek::InputController* inputController, Levek::Camera& camera, float dt, float speed) {
    
    glm::vec3 positionOffset = glm::vec3(0.0);
	const float cameraSpeed = speed * dt;

	if (inputController->isKeyPressed(Levek::LEVEK_KEY_W)) {
		positionOffset -= cameraSpeed * camera.mFront;
	}
	if (inputController->isKeyPressed(Levek::LEVEK_KEY_S)) {
		positionOffset += cameraSpeed * camera.mFront;
    }
	if (inputController->isKeyPressed(Levek::LEVEK_KEY_A)) {
		positionOffset += glm::normalize(glm::cross(camera.mFront, camera.mUp)) * cameraSpeed;
    }
	if (inputController->isKeyPressed(Levek::LEVEK_KEY_D)) {
		positionOffset -= glm::normalize(glm::cross(camera.mFront, camera.mUp)) * cameraSpeed;
    }
	//Levek::printVec3(camera.getEye());
	camera.addEye(positionOffset);
}

extern float lastMouseValueX;
extern float lastMouseValueY;
extern bool wasLeftButtonPressed;

static void UpdateCameraWithMouseOnDrag(Levek::InputController* inputController, Levek::Camera& camera, float sensivity) {

	float mouseX = inputController->getMouseX();
	float mouseY = inputController->getMouseY();
    //&& !ImGui::IsItemHovered() && !ImGui::IsWindowHovered()
	if (inputController->isLeftMouseButtonPressed()) {

		if (!wasLeftButtonPressed) {
			lastMouseValueX = mouseX;
			lastMouseValueY = mouseY;
			wasLeftButtonPressed = true;
		}

		float offsetX = (mouseX - lastMouseValueX) * sensivity;
		float offsetY = (lastMouseValueY - mouseY) * sensivity;
		lastMouseValueX = mouseX;
		lastMouseValueY = mouseY;

		camera.mYaw += offsetX;
		camera.mPitch -= offsetY;

		if (camera.mPitch > 89.0f) {
			camera.mPitch = 89.0f;
		}
		if (camera.mPitch < -89.0f) {
			camera.mPitch = -89.0f;
		}

		glm::vec3 direction;
		direction.z = cos(glm::radians(camera.mYaw)) * cos(glm::radians(camera.mPitch));
		direction.y = sin(glm::radians(camera.mPitch));
		direction.x = sin(glm::radians(-camera.mYaw)) * cos(glm::radians(camera.mPitch));
		camera.setFront(glm::normalize(direction));

	} else {
		wasLeftButtonPressed = false;
	}

	/*
	if (this->isMouseDragging || this->mouseCursorDisabled) {
		if (this->firstMouseAction) {
			this->lastX = (float)xpos;
			this->lastY = (float)ypos;
			this->firstMouseAction = false;
		}
		float xoffset = (float)xpos - this->lastX;
		float yoffset = this->lastY - (float)ypos;
		//std::cout << xoffset << " " << yoffset << std::endl;
		this->lastX = (float)xpos;
		this->lastY = (float)ypos;
		this->lastOffsetX = xoffset;
		this->lastOffsetY = yoffset;
		float sensitivity = 0.2f;
		xoffset *= sensitivity;
		yoffset *= sensitivity;
		yaw += xoffset;
		pitch += yoffset;
		if (pitch > 89.0f)
			pitch = 89.0f;
		if (pitch < -89.0f)
			pitch = -89.0f;
		glm::vec3 direction;
		direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
		direction.y = sin(glm::radians(pitch));
		direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
		target = glm::normalize(direction);
	//}
	*/
}

static void addImGuiVec3(const glm::vec3& v) {
	//ImGui::Text("[%.3f, %.3f, %.3f]", v[0], v[1], v[2]);
} 

