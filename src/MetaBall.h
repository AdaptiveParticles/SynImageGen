//
// Created by bevanc on 26.01.17.
//

#ifndef SYNIMAGEGEN_METABALL_H
#define SYNIMAGEGEN_METABALL_H

#include "MeshDataAF.h"
#include "NoiseModel.hpp"
#include "SynImageClasses.hpp"
#include "SynImagePar.hpp"

class MetaBall {


    float x0;
    float y0;
    float z0;

    float mag;

    MetaBall(float x0,float y0,float z0,float mag): y0(y0), x0(x0), z0(z0),mag(mag) {
    }

};

float MetaBall::calculate_potential(float x,float y,float z){
    //
    // Reproduced from http://www.geisswerks.com/ryan/BLOBS/blobs.html
    //
    //

    float dx = (x - x0);
    float dy = (y - y0);
    float dz = (z - z0);
    float r_squared = (dx*dx + dy*dy + dz*dz);

    //enable this line if your blobs are of varying sizes:
    r_squared /= (mag*mag);

    // since f(r) is valid when r is in the range [0-.707],
    // r_squared should be in the range [0-.500].
    if (r_squared < 0.5f)     // same as: if (sqrtf(r_squared) < 0.707f)
    {
        q += (0.25 - r_squared + r_squared*r_squared);
    }

    return q;

}

void generate_metaball_template(Object_template& obj_template,int sample_rate,float real_size,float density = 1000000,float rad_ratio = 3.0/8.0){
    //
    //
    //  Generate a simple binary sphere for testing
    //
    //
    //
    //

    //set size
    obj_template.real_size.resize(3);
    obj_template.real_size[0] = real_size;
    obj_template.real_size[1] = real_size;
    obj_template.real_size[2] = real_size;

    obj_template.real_deltas.resize(3);

    obj_template.real_deltas[0] = real_size/sample_rate;
    obj_template.real_deltas[1] = real_size/sample_rate;
    obj_template.real_deltas[2] = real_size/sample_rate;

    obj_template.set_max_sample(density);

    //init the data
    obj_template.true_object_distribution.initialize(sample_rate, sample_rate, sample_rate, 0);


    MeshDataAF<float> potential_field;
    potential_field.initialize(sample_rate, sample_rate, sample_rate, 0);

    float center_sphere = sample_rate/2;
    float radius_sphere = sample_rate*rad_ratio;

    float curr_dist;

    int num_meta_balls = 5;

    std::vector<MetaBall> balls;

    Genrand_uni grand;

    for (int l = 0; l < num_meta_balls; ++l) {

        float mag = grand.rand_num(.2*radius_sphere,.5*radius_sphere);

        float x0 = grand.rand_num(mag,basic_sphere.true_object_distribution.x_num-mag);
        float y0 = grand.rand_num(mag,basic_sphere.true_object_distribution.y_num-mag);
        float z0 = grand.rand_num(mag,basic_sphere.true_object_distribution.z_num-mag);

        MetaBall temp_ball(x0,y0,z0,mag);

        balls.push_back(temp_ball);

    }

    // loop over and create field

    for (int i = 0; i < basic_sphere.true_object_distribution.y_num; i++) {
        for (int j = 0; j < basic_sphere.true_object_distribution.x_num; j++) {
            for (int k = 0; k < basic_sphere.true_object_distribution.z_num; k++) {

                for (int l = 0; l < num_meta_balls; ++l) {
                    potential_field(i,j,k) += balls[l].calculate_potential(i,j,k);
                }

            }
        }
    }


    //need threshold now the generated field
    float threshold = 0.15;

    for (int i = 0; i < basic_sphere.true_object_distribution.y_num; i++) {
        for (int j = 0; j < basic_sphere.true_object_distribution.x_num; j++) {
            for (int k = 0; k < basic_sphere.true_object_distribution.z_num; k++) {

                if (potential_field < threshold) {
                    basic_sphere.true_object_distribution(i,j,k) = 1*basic_sphere.max_sample;
                } else {
                    basic_sphere.true_object_distribution(i,j,k) = 0;
                }
            }
        }
    }


}


#endif //SYNIMAGEGEN_METABALL_H
