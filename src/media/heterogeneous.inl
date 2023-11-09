#include "../volume.h"

Spectrum get_majorant_op::operator()(const HeterogeneousMedium &m) {
    if (intersect(m.density, ray)) 
    {
        return get_max_value(m.density);
    } 
    else 
    {
        // std::cout<<m.density.
        // std::cout<<ray.org<<std::endl;
        return make_zero_spectrum();
    }
}

Spectrum get_sigma_s_op::operator()(const HeterogeneousMedium &m) {
    Spectrum density = lookup(m.density, p);
    Spectrum albedo = lookup(m.albedo, p);
    // std::cout<<"density: "<<density<<std::endl;
    // std::cout<<"albedo: "<<albedo<<std::endl;
    // std::cout<<density * (albedo)<<std::endl;
    return density * albedo;
}

Spectrum get_sigma_a_op::operator()(const HeterogeneousMedium &m) {
    Spectrum density = lookup(m.density, p);
    Spectrum albedo = lookup(m.albedo, p);
    // std::cout<<"density: "<<density<<std::endl;
    // std::cout<<"albedo: "<<albedo<<std::endl;
    // std::cout<<density * (Real(1) - albedo)<<std::endl;
    return density * (Real(1) - albedo);
}

Spectrum get_sigma_t_op::operator()(const HeterogeneousMedium &m) {
    Spectrum density = lookup(m.density, p);
    return density;
}