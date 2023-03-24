#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        *pdf = 1;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        double numerator = exp(-pow(tan(acos(h.z)) / alpha, 2));
        double denominator = PI * pow(alpha, 2) * pow(h.z, 4);

        return numerator / denominator;
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.
        double cosThei = abs_cos_theta(wi);
        Vector3D Rs = ((eta * eta + k * k) - (2. * eta * cosThei) + (cosThei * cosThei)) / ((eta * eta + k * k) + (2. * eta * cosThei) + (cosThei * cosThei));
        Vector3D Rp = ((eta * eta + k * k) * (cosThei * cosThei) - (2. * eta * cosThei) + 1.) / ((eta * eta + k * k) * (cosThei * cosThei) + (2. * eta * cosThei) + 1.);

    return (Rs + Rp) / 2.;
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.
        if (wo.z <= 0 || wi.z <= 0) return Vector3D(0.);
        Vector3D h = (wo + wi) / (wo + wi).norm();
        return F(wi) * G(wo, wi) * D(h) / (4. * wo.z * wi.z);
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        // *wi = cosineHemisphereSampler.get_sample(pdf);
        // return MicrofacetBSDF::f(wo, *wi);

        // sample r1 and r2, two random numbers uniformly distributed within [0, 1)
        Vector2D r = sampler.get_sample();
        double the_h = atan(sqrt(-pow(alpha, 2) * log(1. - r.x)));
        double phi_h = 2. * PI * r.y;
        Vector3D h(sin(the_h) * cos(phi_h), sin(the_h) * sin(phi_h), cos(the_h));

        double p_theta = (2. * sin(the_h) * exp(-1. * pow(tan(the_h)/alpha, 2)) / (pow(alpha, 2) * pow(cos(the_h), 3)));
        double p_phi = 1. / (2. * PI);
        double pw_h = p_theta * p_phi / sin(the_h);
            
        // Assign wi based on the sampled half vector
        *wi = (2. * dot(h, wo)) * h - wo;
        wi->normalize();

        // Check sampled wi is valid
        if ((*wi).z < 0) {
            *pdf = 0.;
            return Vector3D(0.);
        } else {
            *pdf = pw_h / (4. * dot(*wi, h));
            return MicrofacetBSDF::f(wo, *wi);
        }
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        bool total_internal_reflection = !refract(wo, wi, ior);
        if (total_internal_reflection) {
            return Vector3D();
        }
        else {
            double cos_theta = wo.z;
            float eta = (cos_theta <= 0)? ior : 1 / ior;
            *pdf = 1;
            return transmittance / abs_cos_theta(*wi) / (eta * eta);
        }
        // return  Vector3D(); 
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305
        bool total_internal_reflection = !refract(wo, wi, ior);
        double cos_theta = wo.z;
        double refraction_index = (cos_theta <= 0) ? ior : 1 / ior;

        if (total_internal_reflection) {
            reflect(wo, wi);
            *pdf = 1;
            return reflectance / abs_cos_theta(*wi);
        }
        else {
            // Schlick's approximation for the Fresnel Factor
            double r0 = pow((1 - ior) / (ior + 1), 2);
            double r = r0 + (1 - r0) * pow((1 - abs(cos_theta)), 5); 
            if (coin_flip(r)) {
                reflect(wo, wi);
                *pdf = r;
                return r * reflectance / abs_cos_theta(*wi);
            }
            else {
                *pdf = 1 - r;
                return (1 - r) * transmittance / abs_cos_theta(*wi) / pow(refraction_index,2);
            }
        }
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = Vector3D(-wo.x, -wo.y, wo.z);

    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        
        double cos_theta = wo.z;
        float refraction_index = (cos_theta < 0)? ior : 1 / ior;
        float cos_theta_prime_squared = 1 - refraction_index * refraction_index * (1 - cos_theta * cos_theta);

        // total internal reflection
        if (cos_theta_prime_squared < 0) {
            return false;
        } 
        else {
            *wi = Vector3D(-refraction_index * wo.x, -refraction_index * wo.y, -wo.z/abs(wo.z) * sqrt(cos_theta_prime_squared));
            return true;
        }
    }

} // namespace CGL
