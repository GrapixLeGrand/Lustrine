#include "Simulation.hpp"


namespace Lustrine {


/**
 * @brief cubic kernel. from SplishSplash repository https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/SPlisHSPlasH/SPHKernels.h 
 * 
 * @param simulation 
 * @param r 
 * @return float 
 */
extern float cubic_kernel(const Simulation* simulation, float r);


/**
 * @brief cubic kernel. from SplishSplash repository https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/SPlisHSPlasH/SPHKernels.h 
 * 
 * @param simulation 
 * @param r 
 * @return float 
 */
extern float cubic_kernel(const Simulation* simulation, glm::vec3& r);

/**
 * @brief gradient of cubic kernel. from SplishSplash repository https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/SPlisHSPlasH/SPHKernels.h 
 * 
 * @param simulation 
 * @param r 
 * @return float 
 */
extern glm::vec3 cubic_kernel_grad(const Simulation* simulation, const glm::vec3& r);

extern float poly6_kernel(const Simulation* simulation, float r);
extern float poly6_kernel(const Simulation* simulation, glm::vec3& r);
extern glm::vec3 spiky_kernel(const Simulation* simulation, glm::vec3& r);

} // namespace Lustrine