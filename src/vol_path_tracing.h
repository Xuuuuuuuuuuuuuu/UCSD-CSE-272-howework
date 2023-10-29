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

    // Medium current_medium=scene.media[scene.camera.medium_id];
    int current_medium_id=scene.camera.medium_id;
    Medium current_medium;
    // std::cout<<"camera.medium_id: "<<scene.camera.medium_id;

    // Spectrum sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
    // Spectrum sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
    // Spectrum sigma_t=sigma_s+sigma_a;
    Spectrum sigma_s,sigma_a,sigma_t;

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
        //no volume ，no need to volume render
        if(current_medium_id==-1)
        {
            if(vertex_)        
            {
                vertex=*vertex_;
                if(is_light(scene.shapes[vertex.shape_id]))
                {
                    if(never_scatter_happen)
                    {
                        radiance+=current_path_throughput*emission(vertex, -ray.dir, scene);
                    }
                    else
                    {
                        Spectrum l_emit=emission(vertex, -ray.dir, scene);
                        // calculate pdf_light
                        int light_id=get_area_light_id(scene.shapes[vertex.shape_id]);
                        const Light &light = scene.lights[light_id];
                        PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                        Real pdf_light= light_pmf(scene, light_id) *
                                pdf_point_on_light(light, light_point, ray.org, scene);

                        //calculate pdf_phase_trans
                        Real G= fabs(dot(ray.dir, vertex.geometric_normal)) /
                        distance_squared(ray.org, vertex.position);
                        Real pdf_phase_trans=last_pdf_phase*G*multi_trans_pdf;

                        // Multiple importance sample：pdf_light and pdf_phase_trans
                        Real w=(pdf_phase_trans*pdf_phase_trans)/(pdf_light*pdf_light+pdf_phase_trans*pdf_phase_trans);
                        radiance+=current_path_throughput*l_emit*w;
                    }
                }

                if(vertex.material_id==-1)
                {
                    //update current_medium
                    if(dot(vertex.geometric_normal,ray.dir)<0)
                    {
                        // current_medium=scene.media[vertex.interior_medium_id];
                        current_medium_id=vertex.interior_medium_id;
                    }
                    else
                    {
                        // current_medium=scene.media[vertex.exterior_medium_id];
                        current_medium_id=vertex.exterior_medium_id;
                    }
                    ray=Ray{vertex.position+ray.dir*get_intersection_epsilon(scene), ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                    continue;
                }
                break;
            }    
            else
            {
                break;
            }
        }
        //volume render
        else
        {
            current_medium=scene.media[current_medium_id];
            sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
            sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
            sigma_t=sigma_s+sigma_a;

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
                            // current_medium=scene.media[vertex.interior_medium_id];
                            current_medium_id=vertex.interior_medium_id;
                        }
                        else
                        {
                            // current_medium=scene.media[vertex.exterior_medium_id];
                            current_medium_id=vertex.exterior_medium_id;
                        }
                        ray=Ray{vertex.position+ray.dir*get_intersection_epsilon(scene), ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                        // sigma_s=get_sigma_s(current_medium,Vector3(0.f,0.f,0.f));
                        // sigma_a=get_sigma_a(current_medium,Vector3(0.f,0.f,0.f));
                        // sigma_t=sigma_s+sigma_a;
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
                    const Light &light = scene.lights[light_id];
                    PointAndNormal point_on_light =sample_point_on_light(light, termination_point, light_uv, shape_w, scene);
                    //some preparations
                    Vector3 dir_light = normalize(point_on_light.position - termination_point);
                    Ray shadow_ray{termination_point, dir_light, 
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            distance(point_on_light.position, termination_point)};
                    std::optional<PathVertex> shadow_vertex_;
                    PathVertex shadow_vertex;
                    Vector3 current_vertex=termination_point;
                    Spectrum T(1,1,1);
                    int current_medium_id_directL=current_medium_id;
                    Medium current_medium_directL;
                    Spectrum current_sigma_t;
                    int shadow_bounces=0;
                    //calculate T(cumulative transmittance)，because some object can be passed through
                    while(true)
                    {
                        if(shadow_bounces+bounces>=scene.options.max_depth) 
                        {
                            T=make_zero_spectrum();
                            break;
                        }
                        shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
                        shadow_bounces++;
                        //no volume
                        if(current_medium_id_directL==-1)
                        {
                            current_sigma_t=make_zero_spectrum();;
                        }
                        else
                        {
                            current_medium_directL=scene.media[current_medium_id_directL];
                            current_sigma_t=get_sigma_s(current_medium_directL,Vector3(0.f,0.f,0.f))
                                +get_sigma_a(current_medium_directL,Vector3(0.f,0.f,0.f));
                        }
                        //no intersect
                        if(!shadow_vertex_) 
                        {
                            T=T*exp(-current_sigma_t*distance(current_vertex,point_on_light.position));
                            // std::cout<<"no object occlude the direct light"<<std::endl;
                            break;
                        }

                        //pass through the object
                        if(shadow_vertex_&&shadow_vertex_->material_id==-1) 
                        {
                            shadow_vertex=*shadow_vertex_;
                            // std::cout<<"pass through the object!!!"<<std::endl;
                            T=T*exp(-current_sigma_t*distance(current_vertex,shadow_vertex.position));
                            shadow_ray=Ray{shadow_vertex.position, dir_light, 
                                        get_shadow_epsilon(scene),
                                        (1 - get_shadow_epsilon(scene)) *
                                            distance(shadow_vertex.position,point_on_light.position)};
                            current_vertex=shadow_vertex.position;
                            // update current_medium_directL
                            if(dot(shadow_vertex.geometric_normal,shadow_ray.dir)<0)
                                current_medium_id_directL=shadow_vertex.interior_medium_id;
                            else
                                current_medium_id_directL=shadow_vertex.exterior_medium_id;

                            // current_sigma_t=get_sigma_s(current_medium_directL,Vector3(0.f,0.f,0.f))
                            //         +get_sigma_a(current_medium_directL,Vector3(0.f,0.f,0.f));
                        }
                        else  //light is occluded
                        {
                            T=make_zero_spectrum();
                            // std::cout<<"light is occluded"<<std::endl;
                            break;
                        }
                    }

                    // calculate the radiance from direct light
                    Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, termination_point);
                    Spectrum phase_value=eval(phase,-ray.dir,dir_light);
                    Spectrum L_direct=emission(light, -dir_light, Real(0), point_on_light, scene);
                    Real pdf_light = light_pmf(scene, light_id) *
                        pdf_point_on_light(light, point_on_light, termination_point, scene);
                    // if(pdf_light<=0) std::cout<<"pdf_light <=0"<<std::endl;
                    Real pdf_phase=pdf_sample_phase(phase,-ray.dir,dir_light)*G*T.x;
                    // std::cout<<T<<std::endl;
                    Real w=(pdf_light*pdf_light)/(pdf_light*pdf_light+pdf_phase*pdf_phase);
                    radiance+=current_path_throughput*sigma_s*phase_value*T*(L_direct/pdf_light)*G*w;            
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
                    else
                    {
                        Spectrum l_emit=emission(vertex, -ray.dir, scene);
                        // calculate pdf_light
                        int light_id=get_area_light_id(scene.shapes[vertex.shape_id]);
                        const Light &light = scene.lights[light_id];
                        PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                        Real pdf_light= light_pmf(scene, light_id) *
                                pdf_point_on_light(light, light_point, ray.org, scene);

                        //calculate pdf_phase_trans
                        Real G= fabs(dot(ray.dir, vertex.geometric_normal)) /
                        distance_squared(ray.org, vertex.position);
                        Real pdf_phase_trans=last_pdf_phase*G*multi_trans_pdf;

                        // Multiple importance sample：pdf_light and pdf_phase_trans
                        Real w=(pdf_phase_trans*pdf_phase_trans)/(pdf_light*pdf_light+pdf_phase_trans*pdf_phase_trans);
                        radiance+=current_path_throughput*sigma_t*l_emit*w;
                    }
                }    
                break;    
            }
        }
    }
    return radiance;    
    // return make_zero_spectrum();
}



//update medium
int update_medium(const std::optional<PathVertex> & surface_point_,const Ray &ray)
{
    if(dot(surface_point_->geometric_normal,ray.dir)<0)
        return surface_point_->interior_medium_id;
    else
        return surface_point_->exterior_medium_id;
}
void get_medium_sigma(const Scene &scene, const int &medium_id,Spectrum & sigma_a,Spectrum & sigma_s,Spectrum & sigma_t)
{
    if(medium_id==-1) 
    {
        sigma_a=make_zero_spectrum();
        sigma_s=make_zero_spectrum();
    }
    else
    {
        Medium medium=scene.media[medium_id];
        sigma_a=get_sigma_a(medium,Vector3(0.f,0.f,0.f));
        sigma_s=get_sigma_s(medium,Vector3(0.f,0.f,0.f));
    }
    sigma_t=sigma_a+sigma_s;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) 
{
    // Homework 2: implememt this!
    // if(x==220&&y==79)
    // {
    //     std::cout<<"debug"<<std::endl;
    // }
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);

    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_;
    // PathVertex vertex;

    int current_medium_id=scene.camera.medium_id;
    Spectrum sigma_a,sigma_s,sigma_t;
    get_medium_sigma(scene,current_medium_id,sigma_a,sigma_s,sigma_t);

    Spectrum current_path_throughput(1,1,1);
    Spectrum radiance(0,0,0);
    
    Vector3 terminate_vertex_pos=ray.org;

    bool last_is_scatter_and_reflence=false;
    Real last_f_pdf=0;
    Spectrum last_T(1,1,1);
    int max_depth = scene.options.max_depth;
    //the num_vertices initially set to 3, because the direct light sample has include one scatter（light touchs the object and change the direction）
    //3 vertices indicate 2 scatter（1 + 1 ）
    for (int num_vertices = 3; max_depth == -1 || num_vertices <= max_depth + 1; num_vertices++)
    {
        vertex_ = intersect(scene, ray, ray_diff);
        bool is_scatter=true;

        //have volume, calculate pre-factor
        if(current_medium_id!=-1)
        {
            //transmittance sample, calculate t in a straight line
            Real u=next_pcg32_real<Real>(rng);
            Real t=-log(1-u)/sigma_t[0];

            Real trans_pdf=exp(-sigma_t[0]*t)*sigma_t[0];
            Spectrum transmittance=exp(-sigma_t*t);
            current_path_throughput*=(transmittance/trans_pdf);

            if(vertex_)
            {
                Real hit_t=distance(vertex_->position,terminate_vertex_pos);
                if(t>=hit_t) 
                {
                    // if(x==122&&y==217)
                    // {
                    //     std::cout<<"debug"<<std::endl;
                    // }
                    is_scatter=false;
                }
            }

            //calculate pre-factor
            {
                //surface term
                if(!is_scatter)
                {
                    current_path_throughput*=sigma_t;
                    terminate_vertex_pos=vertex_->position;

                    //for mis
                    last_T*=exp(-sigma_t*distance(terminate_vertex_pos,ray.org));
                }
                //scatter term
                else
                {
                    current_path_throughput*=sigma_s;
                    terminate_vertex_pos=ray.org+t*ray.dir;

                    //for mis
                    last_T=Spectrum(1,1,1);
                }
            }
        }
        //no volume
        else
        {
            //no object intersect
            if(!vertex_) break;
            terminate_vertex_pos=vertex_->position;
            is_scatter=false;
        }

        //sample point is the surface point, consider the situation of emit term and vacuum
        if(!is_scatter&&vertex_) 
        {
            if(is_light(scene.shapes[vertex_->shape_id]))
            {
                if(!last_is_scatter_and_reflence)
                {
                    radiance+=current_path_throughput*emission(*vertex_, -ray.dir, scene);
                }
                else
                {
                    Spectrum l_emit=emission(*vertex_, -ray.dir, scene);

                    int light_id=get_area_light_id(scene.shapes[vertex_->shape_id]);
                    const Light &light = scene.lights[light_id];
                    PointAndNormal light_point{vertex_->position, vertex_->geometric_normal};
                    Real light_pdf= light_pmf(scene, light_id) *
                            pdf_point_on_light(light, light_point, ray.org, scene);

                    Real G= fabs(dot(ray.dir, vertex_->geometric_normal)) /
                    distance_squared(ray.org, vertex_->position);
                    Real dir_pdf=last_f_pdf*G*last_T.x;

                    // Multiple importance sample：pdf_light and pdf_phase_trans
                    Real w=(dir_pdf*dir_pdf)/(light_pdf*light_pdf+dir_pdf*dir_pdf);
                    radiance+=current_path_throughput*l_emit*w;
                }
            }
            //no material,pass through it 
            if(vertex_->material_id==-1)
            {
                ray=Ray{terminate_vertex_pos, ray.dir, get_intersection_epsilon(scene), infinity<Real>()};
                //update medium
                current_medium_id=update_medium(vertex_,ray);
                get_medium_sigma(scene,current_medium_id,sigma_a,sigma_s,sigma_t);
                // if(current_medium_id==-1)std::cout<<"cnm"<<std::endl;
                continue;
            }                
            // break;
        }

        //calculate direct light(sample the light)
        {
            //get point from the light and some preparations
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal point_on_light =sample_point_on_light(light, terminate_vertex_pos, light_uv, shape_w, scene);

            Vector3 light_dir = normalize(point_on_light.position - terminate_vertex_pos);

            //calculate T(cumulative transmittance)，because some object can be passed through
            Spectrum T(1,1,1);
            {
                Ray shadow_ray{terminate_vertex_pos, light_dir, 
                            get_shadow_epsilon(scene),
                            (1 - get_shadow_epsilon(scene)) *
                                distance(point_on_light.position, terminate_vertex_pos)};
                std::optional<PathVertex> shadow_vertex_;
                Vector3 directL_termination_point_pos=terminate_vertex_pos;

                int directL_current_medium_id=current_medium_id;
                Medium directL_current_medium;
                Spectrum directL_current_sigma_a=sigma_a;
                Spectrum directL_current_sigma_s=sigma_s;
                Spectrum directL_current_sigma_t=sigma_t;

                int directL_num_vertices = num_vertices;
                for(directL_num_vertices; max_depth == -1 || directL_num_vertices <= max_depth + 1; directL_num_vertices++)
                {
                    shadow_vertex_ = intersect(scene, shadow_ray, ray_diff);
                    //no intersect
                    if(!shadow_vertex_) 
                    {
                        T=T*exp(-directL_current_sigma_t*distance(directL_termination_point_pos,point_on_light.position));
                        // std::cout<<"no object occlude the direct light"<<std::endl;
                        break;
                    }
                    //pass through the object
                    if(shadow_vertex_&&shadow_vertex_->material_id==-1) 
                    {
                        // std::cout<<"pass through the object!!!"<<std::endl;
                        T*=exp(-directL_current_sigma_t*distance(directL_termination_point_pos,shadow_vertex_->position));
                        directL_termination_point_pos=shadow_vertex_->position;                    
                        shadow_ray=Ray{shadow_vertex_->position, light_dir, 
                                    get_shadow_epsilon(scene),
                                    (1 - get_shadow_epsilon(scene)) *
                                        distance(shadow_vertex_->position,point_on_light.position)};
                        // update current_medium_directL
                        directL_current_medium_id=update_medium(shadow_vertex_,shadow_ray);
                        get_medium_sigma(scene,directL_current_medium_id,directL_current_sigma_a,directL_current_sigma_s,directL_current_sigma_t);
                    }
                    else  //light is occluded, I think object with material can be passed through, waiting for the future!!!
                    {
                        T=make_zero_spectrum();
                        // std::cout<<"light is occluded"<<std::endl;
                        break;
                    }                
                }
                if(max_depth!=-1&&directL_num_vertices>max_depth + 1)
                {
                    T=make_zero_spectrum();
                }
            }  

            // calculate the radiance from direct light
            Real G = max(-dot(light_dir, point_on_light.normal), Real(0)) /
                distance_squared(point_on_light.position, terminate_vertex_pos);

            Spectrum f;
            Real f_pdf,dir_pdf;
            if(is_scatter)
            {
                const PhaseFunction phase=get_phase_function(scene.media[current_medium_id]);
                f=eval(phase,-ray.dir,light_dir);
                f_pdf=pdf_sample_phase(phase,-ray.dir,light_dir);
            }
            else
            {
                const Material &mat = scene.materials[vertex_->material_id];
                f=eval(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool);
                f_pdf=pdf_sample_bsdf(mat, -ray.dir, light_dir, *vertex_, scene.texture_pool);
            }
            dir_pdf=f_pdf*G*T.x;

            Spectrum L_direct=emission(light, -light_dir, Real(0), point_on_light, scene);
            Real light_pdf = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, terminate_vertex_pos, scene);

            Real w=(light_pdf*light_pdf)/(light_pdf*light_pdf+dir_pdf*dir_pdf);
            radiance+=current_path_throughput*f*(L_direct*T/light_pdf)*G*w;

        }

        //calculate indirect light(sample the phase or bsdf)
        {
            Vector3 f_dir;
            Spectrum f;
            Real f_pdf;
            if(is_scatter)
            {
                const PhaseFunction phase=get_phase_function(scene.media[current_medium_id]);
                Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                std::optional<Vector3> phase_dir_=sample_phase_function(phase,-ray.dir,phase_rnd_param_uv);
                if(!phase_dir_)
                {
                    // Phase sampling failed. Abort the loop.
                    break;
                }
                f_dir=*phase_dir_;
                f=eval(phase,-ray.dir,f_dir);
                f_pdf=pdf_sample_phase(phase,-ray.dir,f_dir);
            }
            else
            {
                const Material &mat = scene.materials[vertex_->material_id];
                Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_=
                    sample_bsdf(mat,-ray.dir,*vertex_,scene.texture_pool,bsdf_rnd_param_uv,bsdf_rnd_param_w);
                if (!bsdf_sample_) 
                {
                    // BSDF sampling failed. Abort the loop.
                    break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                f_dir= bsdf_sample.dir_out;  
                f=eval(mat, -ray.dir, f_dir, *vertex_, scene.texture_pool); 
                f_pdf=pdf_sample_bsdf(mat, -ray.dir, f_dir, *vertex_, scene.texture_pool);
            }

            current_path_throughput*=(f/f_pdf);
            ray=Ray{terminate_vertex_pos, f_dir, get_intersection_epsilon(scene), infinity<Real>()};
            last_f_pdf=f_pdf;

            //update medium when the sample point in the surface
            if(!is_scatter)
            {
                current_medium_id=update_medium(vertex_,ray);
                get_medium_sigma(scene,current_medium_id,sigma_a,sigma_s,sigma_t);
            }

            last_is_scatter_and_reflence=true;

            // Ray f_ray{terminate_vertex_pos, f_dir, get_intersection_epsilon(scene), infinity<Real>()};
            // std::optional<PathVertex> f_vertex_ = intersect(scene, f_ray);
            // Real G;
            // if (f_vertex_) 
            // {
            //     G = fabs(dot(f_dir, f_vertex_->geometric_normal)) /
            //         distance_squared(f_vertex_->position, vertex_->position);
            // } 
            // else 
            // {
            //     // We hit nothing, set G to 1 to account for the environment map contribution.
            //     G = 1;
            // }
        }

        //Russian roulette
        if(num_vertices-1>=scene.options.rr_depth)
        {
            Real rr_prob=std::min(max(current_path_throughput),0.95);
            if(next_pcg32_real<Real>(rng)>rr_prob)
                break;
            else
                current_path_throughput/=rr_prob;
        } 
    }
    return radiance;    
}

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



