#pragma once
Spectrum nee(Vector3 p, Vector3 omega, int current_medium, int bounces, Light light, const Scene& scene, pcg32_state& rng, int light_id);

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    if(vertex_)
    {
        PathVertex vertex = *vertex_;
        Spectrum sigma_a=get_sigma_a(scene.media[vertex.exterior_medium_id],Vector3(0.f,0.f,0.f));
        Real t=distance(vertex.position, ray.org);

        Spectrum transmittance=exp(-sigma_a*t);
        Spectrum Le=make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id]))
            Le=emission(vertex, -ray.dir, scene);
        return transmittance*Le;
    }

    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    Medium medium=scene.media[scene.camera.medium_id];
    // std::cout<<"camera.medium_id: "<<scene.camera.medium_id;
    Spectrum sigma_s=get_sigma_s(medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_a=get_sigma_a(medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_t=sigma_s+sigma_a;
    Real u=next_pcg32_real<Real>(rng);
    Real t=-log(1-u)/sigma_t[0];
    Vector3 p=ray.org+t*ray.dir;

    if(vertex_)
    {
        PathVertex vertex = *vertex_;
        Real hit_t=distance(vertex.position, ray.org);
        //scatter
        if(t<hit_t)
        {
            Real trans_pdf=exp(-sigma_t[0]*t)*sigma_t[0];
            Spectrum transmittance=exp(-sigma_t*t);

            //calculate Ls_1
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            IsotropicPhase phase;
            std::optional<Vector3>dir_phase_ = sample_phase_function(phase,-ray.dir,phase_rnd_param_uv);
            if(!dir_phase_)return make_zero_spectrum();//phase sampling failed
            Vector3 dir_phase=*dir_phase_;
            Ray ray_i{p, dir_phase, get_intersection_epsilon(scene), infinity<Real>()};
            std::optional<PathVertex> vertex_i_ = intersect(scene, ray_i, ray_diff);
            if(vertex_i_)
            {
                Spectrum Ls_1=make_zero_spectrum();
                PathVertex vertex_i = *vertex_i_;
                Real t_i=distance(vertex_i.position, ray_i.org);
                Spectrum transmittance_i=exp(-sigma_t*t_i);
                if (is_light(scene.shapes[vertex_i.shape_id]))
                {
                    Ls_1=transmittance_i*emission(vertex_i, -ray_i.dir, scene)*eval(phase,-ray.dir,ray_i.dir);
                    Real Ls_1_pdf=pdf_sample_phase(phase,-ray.dir,ray_i.dir);
                    return Ls_1*transmittance*sigma_s/trans_pdf/Ls_1_pdf/(1-exp(-sigma_t[0]*hit_t));
                }
            }
            return make_zero_spectrum();
        }
        //emit
        else
        {
            Real trans_pdf=exp(-sigma_t[0]*hit_t);
            Spectrum transmittance=exp(-sigma_t*hit_t);
            if (is_light(scene.shapes[vertex.shape_id]))
            {
                return (transmittance/trans_pdf)*emission(vertex, -ray.dir, scene);
            }
            return make_zero_spectrum();
        }
    }
    else
    {
        //calculate Ls_1
        //only scatter
        Real trans_pdf=exp(-sigma_t[0]*t)*sigma_t[0];
        Spectrum transmittance=exp(-sigma_t*t);
        Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        IsotropicPhase phase;
        std::optional<Vector3>dir_phase_ = sample_phase_function(phase,-ray.dir,phase_rnd_param_uv);
        if(!dir_phase_)return make_zero_spectrum();//phase sampling failed
        Vector3 dir_phase=*dir_phase_;
        Ray ray_i{p, dir_phase, get_intersection_epsilon(scene), infinity<Real>()};
        std::optional<PathVertex> vertex_i_ = intersect(scene, ray_i, ray_diff);
        if(vertex_i_)
        {
            Spectrum Ls_1=make_zero_spectrum();
            PathVertex vertex_i = *vertex_i_;
            Real t_i=distance(vertex_i.position, ray_i.org);
            Spectrum transmittance_i=exp(-sigma_t*t_i);
            if (is_light(scene.shapes[vertex_i.shape_id]))
            {
                Ls_1=transmittance_i*emission(vertex_i, -ray_i.dir, scene)*eval(phase,-ray.dir,ray_i.dir);
                Real Ls_1_pdf=pdf_sample_phase(phase,-ray.dir,ray_i.dir);
                return Ls_1*transmittance*sigma_s/trans_pdf/Ls_1_pdf;
            }
        }
        return make_zero_spectrum();
    }
    return make_zero_spectrum();
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) 
{
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real hit_t;

    Medium current_medium=scene.media[scene.camera.medium_id];
    // std::cout<<"camera.medium_id: "<<scene.camera.medium_id;
    Spectrum sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_t=sigma_s+sigma_a;

    Spectrum current_path_throughput(1,1,1);
    Spectrum radiance(0,0,0);

    bool scatter_happen;
    int bounces=0;
    //Actually, every time in while(true) corresponds to a calculation of L_3
    while(true)
    {
        if(bounces>=scene.options.max_depth) break;
        //Russian roulette
        if(bounces>=scene.options.rr_depth)
        {
            Real rr_prob=std::min(current_path_throughput[0],0.95);
            if(next_pcg32_real<Real>(rng)>rr_prob)
                break;
            else
                current_path_throughput=current_path_throughput/rr_prob;
        } 
        //intersect
        vertex_ = intersect(scene, ray, ray_diff);
        bounces=bounces+1;//bounces+1, when shot a ray passing through scene
        //sample the transmittance
        Real u=next_pcg32_real<Real>(rng);
        Real t=-log(1-u)/sigma_t[0];
        Vector3 termination_point=ray.org+t*ray.dir;
        //the part(function) sampled corresponds to its pdf
        Real trans_pdf=exp(-sigma_t[0]*t)*sigma_t[0];
        Spectrum transmittance=exp(-sigma_t*t);
        current_path_throughput=current_path_throughput*(transmittance/trans_pdf);
        //only when intersect happends and t>hit_t, the emit part will be calculated
        scatter_happen=true;
        if(vertex_)
        {
            vertex = *vertex_;
            hit_t=distance(vertex.position, ray.org);
            if(t>hit_t) //hit the surface
            {
                scatter_happen=false;
                //the surface has no material,so pass through 
                if(vertex.material_id==-1)
                {
                    current_path_throughput=current_path_throughput*sigma_t;
                    if(dot(vertex.geometric_normal,ray.dir)<0)
                        current_medium=scene.media[vertex.interior_medium_id];
                    else
                        current_medium=scene.media[vertex.exterior_medium_id];
                    ray=Ray{vertex.position, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                    sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
                    sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
                    sigma_t=sigma_s+sigma_a;
                    continue;
                }
            }
        }
        //calculate scatter part
        if(scatter_happen)
        {
            //sample the phase function
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            PhaseFunction phase=get_phase_function(current_medium);
            std::optional<Vector3>dir_phase_ = sample_phase_function(phase,-ray.dir,phase_rnd_param_uv);
            if(!dir_phase_)return make_zero_spectrum();//phase sampling failed
            Vector3 dir_phase=*dir_phase_;
            Spectrum phase_value=eval(phase,-ray.dir,dir_phase);
            Real phase_pdf=pdf_sample_phase(phase,-ray.dir,dir_phase);

            current_path_throughput=current_path_throughput*sigma_s*(phase_value/phase_pdf);

            ray=Ray{termination_point, dir_phase, get_intersection_epsilon(scene), infinity<Real>()};
        }
        //calculate emit part
        else
        {
            // std::cout<<"emit!!!"<<std::endl;
            if (is_light(scene.shapes[vertex.shape_id]))
            {
                Spectrum l_emit=emission(vertex, -ray.dir, scene);
                radiance=radiance+current_path_throughput*sigma_t*l_emit;
            }            
            break;
        }
    }

    return radiance;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) 
{
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> vertex_;
    PathVertex vertex;
    Real hit_t;

    Medium current_medium=scene.media[scene.camera.medium_id];
    int current_medium_id=scene.camera.medium_id;
    // std::cout<<"camera.medium_id: "<<scene.camera.medium_id;
    Spectrum sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
    Spectrum sigma_t=sigma_s+sigma_a;

    Spectrum current_path_throughput(1,1,1);
    Spectrum radiance(0,0,0);

    bool scatter_happen;
    bool never_scatter_happen=1;
    int bounces=0;
    Real last_pdf_phase=0;
    Real multi_trans_pdf=1;
    //Actually, every time in while(true) corresponds to a calculation of L_3
    while(true)
    {
        if(bounces>=scene.options.max_depth) break;
        //Russian roulette
        if(bounces>=scene.options.rr_depth)
        {
            Real rr_prob=std::min(current_path_throughput[0],0.95);
            if(next_pcg32_real<Real>(rng)>rr_prob)
                break;
            else
                current_path_throughput=current_path_throughput/rr_prob;
        } 
        //intersect
        vertex_ = intersect(scene, ray, ray_diff);
        bounces=bounces+1;//bounces+1, when shot a ray passing through scene
        //sample the transmittance
        Real u=next_pcg32_real<Real>(rng);
        Real t=-log(1-u)/sigma_t[0];
        Vector3 termination_point=ray.org+t*ray.dir;
        //the part(function) sampled corresponds to its pdf
        Real trans_pdf=exp(-sigma_t[0]*t)*sigma_t[0];
        Spectrum transmittance=exp(-sigma_t*t);
        current_path_throughput=current_path_throughput*(transmittance/trans_pdf);
        //only when intersect happends and t>hit_t, the emit part will be calculated
        scatter_happen=true;
        if(vertex_)
        {            
            vertex = *vertex_;
            hit_t=distance(vertex.position, ray.org);
            if(t>hit_t) //hit the surface
            {
                scatter_happen=false;
                multi_trans_pdf=multi_trans_pdf*exp(-sigma_t[0]*hit_t);
                //the surface has no material,so pass through 
                if(vertex.material_id==-1)
                {
                    current_path_throughput*=sigma_t;

                    //update current_medium
                    if(dot(vertex.geometric_normal,ray.dir)<0)
                    {
                        current_medium=scene.media[vertex.interior_medium_id];
                        current_medium_id=vertex.interior_medium_id;
                    }
                    else
                    {
                        current_medium=scene.media[vertex.exterior_medium_id];
                        current_medium_id=vertex.exterior_medium_id;
                    }
                    ray=Ray{vertex.position+ray.dir*get_intersection_epsilon(scene), ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                    sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
                    sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
                    sigma_t=sigma_s+sigma_a;
                    continue;
                }
            }
        }
        //calculate scatter part
        if(scatter_happen)
        {
            never_scatter_happen=0;
            multi_trans_pdf=1;
            PhaseFunction phase=get_phase_function(current_medium);

            //sample the direct light
            {
                //get point from the light
                Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real light_w = next_pcg32_real<Real>(rng);
                Real shape_w = next_pcg32_real<Real>(rng);
                int light_id = sample_light(scene, light_w);
                // const Light &light = scene.lights[light_id];
                // PointAndNormal point_on_light =sample_point_on_light(light, termination_point, light_uv, shape_w, scene);
                // //some preparations
                // Vector3 dir_light = normalize(point_on_light.position - termination_point);
                // Ray shadow_ray{termination_point, dir_light, 
                //     get_shadow_epsilon(scene),
                //     (1 - get_shadow_epsilon(scene)) *
                //         distance(point_on_light.position, termination_point)};
                // std::optional<PathVertex> shadow_vertex_;
                // PathVertex shadow_vertex;
                // Vector3 current_vertex=termination_point;
                // Spectrum T(1,1,1);
                // Medium current_medium_directL=current_medium;
                // Spectrum current_sigma_t=sigma_t;
                // int shadow_bounces=0;
                // //calculate T(cumulative transmittance)，because some object can be passed through
                // while(true)
                // {
                //     if(shadow_bounces+bounces>=scene.options.max_depth) break;
                //     shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
                //     shadow_bounces++;
                //     //no intersect
                //     if(!shadow_vertex_) 
                //     {
                //         T=T*exp(-current_sigma_t*distance(current_vertex,point_on_light.position));
                //         // std::cout<<"no object occlude the direct light"<<std::endl;
                //         break;
                //     }

                //     shadow_vertex=*shadow_vertex_;
                //     //pass through the object
                //     if(shadow_vertex.material_id==-1) 
                //     {
                //         // std::cout<<"pass through the object!!!"<<std::endl;
                //         T=T*exp(-current_sigma_t*distance(current_vertex,shadow_vertex.position));
                //         shadow_ray=Ray{shadow_vertex.position, dir_light, 
                //                     get_shadow_epsilon(scene),
                //                     (1 - get_shadow_epsilon(scene)) *
                //                         distance(shadow_vertex.position,point_on_light.position)};
                //         current_vertex=shadow_vertex.position;
                //         // update current_medium_directL
                //         if(dot(shadow_vertex.geometric_normal,shadow_ray.dir)<0)
                //             current_medium_directL=scene.media[shadow_vertex.interior_medium_id];
                //         else
                //             current_medium_directL=scene.media[shadow_vertex.exterior_medium_id];

                //         current_sigma_t=get_sigma_s(current_medium_directL,Vector3(0.f,0.f,0.f))
                //                 +get_sigma_a(current_medium_directL,Vector3(0.f,0.f,0.f));
                //     }
                //     else  //light is occluded
                //     {
                //         T=make_zero_spectrum();
                //         // std::cout<<"light is occluded"<<std::endl;
                //         break;
                //     }
                // }

                // // calculate the radiance from direct light
                // Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                //     distance_squared(point_on_light.position, termination_point);
                // Spectrum phase_value=eval(phase,-ray.dir,dir_light);
                // Spectrum L_direct=emission(light, -dir_light, Real(0), point_on_light, scene);
                // Real pdf_light = light_pmf(scene, light_id) *
                //     pdf_point_on_light(light, point_on_light, termination_point, scene);
                // // if(pdf_light<=0) std::cout<<"pdf_light <=0"<<std::endl;
                // Real pdf_phase=pdf_sample_phase(phase,-ray.dir,dir_light)*G*T.x;
                // // std::cout<<T<<std::endl;
                // Real w=(pdf_light*pdf_light)/(pdf_light*pdf_light+pdf_phase*pdf_phase);
                radiance+=current_path_throughput*sigma_s*nee(termination_point, -ray.dir, current_medium_id, bounces, scene.lights[light_id], scene, rng, light_id);
            }


            //sample the phase function
            {
                Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                std::optional<Vector3>dir_phase_ = sample_phase_function(phase,-ray.dir,phase_rnd_param_uv);
                if(!dir_phase_)return make_zero_spectrum();//phase sampling failed
                Vector3 dir_phase=*dir_phase_;
                // std::cout<<dot(dir_phase,dir_phase)<<std::endl;
                Spectrum phase_value=eval(phase,-ray.dir,dir_phase);
                Real pdf_phase=pdf_sample_phase(phase,-ray.dir,dir_phase);
                last_pdf_phase=pdf_phase;
                // Real G = fabs(dot(dir_phase, bsdf_vertex->geometric_normal)) /
                //     distance_squared(bsdf_vertex->position, vertex.position);

                ray=Ray{termination_point+dir_phase*get_intersection_epsilon(scene), dir_phase, get_intersection_epsilon(scene), infinity<Real>()};

                current_path_throughput*=sigma_s*(phase_value/pdf_phase);                
            }
        }
        //calculate emit part
        else
        {
            // std::cout<<"emit!!!"<<std::endl;
            if (is_light(scene.shapes[vertex.shape_id]))
            {
                if(never_scatter_happen)
                {
                    Spectrum l_emit=emission(vertex, -ray.dir, scene);
                    radiance+=current_path_throughput*sigma_t*l_emit;
                }
            }    
            
            break;
        }
    }

    return radiance;    





    // return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
// Spectrum vol_path_tracing_5(const Scene &scene,
//                             int x, int y, /* pixel coordinates */
//                             pcg32_state &rng) {
//     // Homework 2: implememt this!
//     return make_zero_spectrum();
// }

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}


int update_medium(PathVertex isect, Ray ray){

    int medium = 0;
    if (isect.interior_medium_id != isect.exterior_medium_id) {

        if (dot(ray.dir, isect.geometric_normal) > 0) {
            medium = isect.exterior_medium_id;
        }
        else {
            medium = isect.interior_medium_id;
        }
    }
    return medium;
}
Spectrum nee(Vector3 p, Vector3 omega, int current_medium, int bounces, Light light, const Scene& scene, pcg32_state& rng, int light_id) {

    PointAndNormal p_prime = sample_point_on_light(light, p, 
                                Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                next_pcg32_real<Real>(rng), scene);

    Real T_light = 1;
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1;

    Vector3 orig_p = p;

    while (1) {

        Ray shadow_ray = Ray{p, normalize(p_prime.position - p),  get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime.position, p)};
        RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);

        if(isect) {
            PathVertex vertex = *isect;
            next_t = distance(p, vertex.position);
        }
        if(shadow_medium != -1) {
            Medium medium = scene.media[shadow_medium];
            Real sigma_a = get_sigma_a(medium, Vector3(1, 2, 3)).x;
            Real sigma_s = get_sigma_s(medium, Vector3(1, 2, 3)).x;
            Real sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }
            
        if(!isect) {
            break;
        } 
        else {
            PathVertex vertex = *isect;
            if (vertex.material_id >= 0) {
                return make_zero_spectrum();
            }
            
            shadow_bounces ++;
            int max_bounce = scene.options.max_depth;
            
            if (max_bounce != -1 && bounces + shadow_bounces + 1 >= max_bounce) {
                return make_zero_spectrum();
            }
            shadow_medium = update_medium(vertex, shadow_ray);
            p = p + next_t * (shadow_ray.dir);
        }
    }

    if (T_light > 0) {
        Vector3 omega_prime = normalize(orig_p - p_prime.position);
        Real denom = distance_squared(orig_p, p_prime.position);
        Real top = abs(dot(omega_prime, p_prime.normal));
        Real G = top / denom;
        PhaseFunction pf = get_phase_function(scene.media[current_medium]);

        Spectrum f = eval(pf, omega, -omega_prime);
        Spectrum Le = emission(light, omega_prime, Real(0), p_prime, scene);
        Real pdf_nee = light_pmf(scene, light_id) *
            pdf_point_on_light(light, p_prime, orig_p, scene);
        
        Spectrum contrib = T_light * G * f * Le / pdf_nee;
        Real pdf_phase = pdf_sample_phase(pf, omega, -omega_prime) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

        return w * contrib;
    }

    return make_zero_spectrum();

}
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };

    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = make_const_spectrum(1);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    int max_depth = scene.options.max_depth;

    bool never_scatter = true;
    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1;

    while (true) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real transmittance = 1;
        Real trans_pdf = 1;
        Real t = 0;

        if (current_medium == -1) {
            if (vertex_) {
                PathVertex vertex = *vertex_;
                ray.org = ray.org + distance(ray.org, vertex.position) * ray.dir;
            } else {
                break;
            }
        }

        if (current_medium != -1) {
            Real u = next_pcg32_real<Real>(rng);
            Medium medium = scene.media[current_medium];
            Real sigma_a = get_sigma_a(medium, Vector3(1, 2, 3)).x;
            Real sigma_s = get_sigma_s(medium, Vector3(1, 2, 3)).x;
            Real sigma_t = sigma_a + sigma_s;
            t = -log(1.0 - u) / sigma_t;

            //ray hit
            if (vertex_) {
                PathVertex vertex = *vertex_;
                Real t_hit = distance(ray.org, vertex.position);
                // if t not on surface
                if (t < t_hit) {
                    scatter = true;
                    trans_pdf = exp(-sigma_t * t) * sigma_t;
                    transmittance = exp(-sigma_t * t);
                    ray.org = ray.org + t * ray.dir;
                }
                else {
                    scatter = false;
                    trans_pdf = exp(-sigma_t * t_hit);
                    transmittance = exp(-sigma_t * t_hit);
                    ray.org = ray.org + t_hit * ray.dir;
                }
            }
            else {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                ray.org = ray.org + t * ray.dir;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (never_scatter) {
                if (is_light(scene.shapes[vertex.shape_id])) {
                    radiance += current_path_throughput * emission(vertex, -ray.dir, scene);
                }
            } 
            // else {
            //     if (is_light(scene.shapes[vertex.shape_id])) {
            //         PointAndNormal light_point = {vertex.position, vertex.geometric_normal};
            //         int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
            //         Light currlight = scene.lights[light_id];
            //         Real pdf_nee = pdf_point_on_light(currlight, light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
            //         Vector3 omega_prime = normalize(vertex.position - nee_p_cache);
            //         Real top = abs(dot(omega_prime, vertex.geometric_normal));
            //         Real bottom = length_squared(vertex.position - nee_p_cache);
                    
            //         Real G = top/bottom;
            //         Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    
            //         Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    
            //         radiance += current_path_throughput * emission(vertex, -ray.dir, scene) * w;
            //     }
            // }
        }

        if (bounces == max_depth - 1 && max_depth != -1) {
            break;
        }

        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            // if index matching surface
            if (vertex.material_id == -1) {
                current_medium = update_medium(vertex, ray);
                bounces += 1;
                ray.org = vertex.position + ray.dir * get_intersection_epsilon(scene);
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        if (scatter) {
            never_scatter = false;
            Vector3 p = ray.org;
            PhaseFunction pf = get_phase_function(scene.media[current_medium]);
            std::optional<Vector3> next_dir_ = sample_phase_function(pf, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                Real sigma_s = get_sigma_s(scene.media[current_medium], p).x;

                int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                Spectrum nee_out = nee(p, -ray.dir, current_medium, bounces, scene.lights[light_id], scene, rng, light_id);

                radiance += current_path_throughput * nee_out * sigma_s;

                dir_pdf = pdf_sample_phase(pf, -ray.dir, next_dir);
                current_path_throughput *= (eval(pf, -ray.dir, next_dir) / dir_pdf) * sigma_s;

                ray = Ray{ ray.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>() };
                nee_p_cache = p;
                multi_trans_pdf = Real(1);
            }
        }
        else {
            break;
        }

        Real rr_prob = 1;

        if (bounces >= max_depth) {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;
    }

    return radiance;
}
