
#include "Kernels.hpp"

namespace Lustrine {

float cubic_kernel(const Simulation* simulation, float r) {
    float q = (r * simulation->kernelFactor) / simulation->kernelRadius;
    float result = 0.0;
    if (q <= 1.0) {
        if (q <= 0.5) {
            float q2 = q * q;
            float q3 = q2 * q;
            result = simulation->cubic_kernel_k * (6.0f * q3 - 6.0f * q2 + 1.0f);
        } else {
            result = simulation->cubic_kernel_k * (2.0f * std::pow(1.0f - q, 3.0f));
        }
    }
    return result;
}

float cubic_kernel(const Simulation* simulation, glm::vec3& r) {
    float tmp = glm::length(r);
    return cubic_kernel(simulation, tmp);
}

glm::vec3 cubic_kernel_grad(const Simulation* simulation, const glm::vec3& r) {
    glm::vec3 result = glm::vec3(0.0);
    float rl = glm::length(r) * simulation->kernelFactor;
    float q = rl / simulation->kernelRadius;

    if (rl > 1.0e-5 && q <= 1.0) {
        const glm::vec3 grad_q = (1.0f / (rl * simulation->kernelRadius)) * r;
        if (q <= 0.5) {
            result = simulation->cubic_kernel_l * q * (3.0f * q - 2.0f) * grad_q;
        } else {
            const float f = 1.0f - q;
            result = simulation->cubic_kernel_l * (- f * f) * grad_q;
        }
    }
    return result;
}

float poly6_kernel(const Simulation* simulation, float r) {
    float result = 0.0;
    float hf = simulation->kernelRadius * simulation->kernelFactor;
    if (r <= simulation->kernelRadius) {
        result = (315.0f / (64.0f * 3.14f * std::pow(hf, 9))) * 
            std::pow(std::pow(hf, 2) - std::pow(simulation->kernelFactor * r, 2), 3);
    }
    return result;
}

float poly6_kernel(const Simulation* simulation, glm::vec3& r) {
    return poly6_kernel(simulation, r.length());
}

glm::vec3 spiky_kernel(const Simulation* simulation, glm::vec3& r) {
    glm::vec3 result = glm::vec3(0.0);
    float rl = glm::length(r);
    if (rl > 0.0 && rl <= simulation->kernelRadius) {
        float hf = simulation->kernelRadius * simulation->kernelFactor;
        float temp = ((15.0f / (3.14f * std::pow(hf, 6))) * 
            std::pow(hf - (rl * simulation->kernelFactor), 2));
        result = (r / (rl * simulation->kernelFactor)) * temp;
    }
    return result;
}

} // namespace Lustrine