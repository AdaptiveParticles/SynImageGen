//
// Created by cheesema on 25/01/17.
//

#ifndef SYNIMAGEGEN_NOISEMODEL_HPP
#define SYNIMAGEGEN_NOISEMODEL_HPP

#include "arrayfire.h"
#include <vector>
#include <random>

class Genrand_uni {
    //generating uniform random numbers

public:

    std::default_random_engine generator; //generator

    std::uniform_real_distribution<float> dis; //distribution



    unsigned int seed;

    Genrand_uni(){
        seed = std::random_device{}();
        generator.seed(seed);
    }

    void set_seed(unsigned int seed_val){
        seed = seed_val;
        generator.seed(seed_val);

    }

    float rand_num(void) {
        //generates random number from the unifrom 0,1 dist
        return dis(generator);
    }

    float rand_num(float min,float max) {
        //generates random number from the unifrom 0,1 dist mapped to [min,max]
        return (dis(generator)*max + min);
    }



};
class Noise_Model {
    //
    // Class used to add different random distributed noise models to images
    //

public:

    std::mt19937  generator; //generator
    std::normal_distribution<float> distribution{0,1};

    unsigned long rand_seed;
    std::string noise_type;
    float gauss_var;
    float poisson_factor = 1.0f;


    Noise_Model(){
        //rand_seed = std::random_device{}();
        //af::setSeed(rand_seed);
    }

    void set_seed(uint64_t seed_val){
        rand_seed = seed_val;
        af::setSeed(rand_seed);
        generator.seed(seed_val);
    }

    void randomize_seed(){
        rand_seed = std::random_device{}();
        af::setSeed(rand_seed);
    }

    void apply_gaussian_noise(af::array& input_image,float noise_var){
        //
        //  Applied simple guassian noise to the image
        //
        //  Use I(sigma) = J + x(1)*sigma
        //

        input_image += sqrt(noise_var)*af::randn(input_image.dims(0),input_image.dims(1),input_image.dims(2));

        input_image = max(input_image,0.0);

    }

    void apply_multiplicative_noise(af::array& input_image,float noise_var){
        //
        //  Applied multiplicative guassian noise to the image
        //
        //  Use I(sigma) = J + J*x(1)*sigma
        //

        input_image += sqrt(noise_var)*af::randn(input_image.dims(0),input_image.dims(1),input_image.dims(2))*input_image;

        input_image = max(input_image,0.0);

    }

    template<typename T>
    void apply_poisson_noise_cpu(std::vector<T>& input_img){


//        std::cout << distribution.mean() << std::endl;
//        std::cout << distribution.stddev() << std::endl;

        //float val = distribution(generator);

        for (int i = 0; i < input_img.size(); ++i) {
            input_img[i] = std::abs(std::round(input_img[i] + sqrt(1.0*input_img[i])*distribution(generator)));

           // input_img[i] = std::abs(input_img[i]);
        }


    }

    void apply_poisson_noise(af::array& input_image,float noise_var){
        //
        //  Applying poisson type noise as implimented in Matlab imnoise
        //
        //
        //  Can be quite memory intensive
        //

        int th = 1;

        af::array temp;

        if (th ==1) {


            //below code has an error that doens't work for large images, therefore switched it out for this now
            for(int i = 0;i < input_image.dims(2); i++){


                //temp =    input_image(af::span,af::span,i);

                input_image(af::span,af::span,i)+= poisson_factor*sqrt(input_image(af::span,af::span,i))*af::randn(input_image.dims(0),input_image.dims(1),1);

                // temp = input_image(af::span,af::span,i);


            }
        }
        else {
            float threshold = 50;

            //those intensities above threshold
            af::array index_th = af::where(input_image >= threshold);

            // From Matlab: estimate poisson noise as guassian noise
            // For large pixel intensities the Poisson pdf becomes
            //    very similar to a Gaussian pdf of mean and of variance
            //    equal to the local pixel intensities. Ref. Mathematical Methods
            //    of Physics, 2nd Edition, Mathews, Walker (Addison Wesley)

            if (index_th.dims(0) > 0) {

                input_image(index_th) = input_image(index_th) + sqrt(input_image(index_th))*af::randn(index_th.dims(0),index_th.dims(1),index_th.dims(2));

            }

            //below threshold cases
            index_th = af::where(input_image < threshold);

            //  (Monte-Carlo Rejection Method) Ref. Numerical
            //   Recipes in C, 2nd Edition, Press, Teukolsky,
            //   Vetterling, Flannery (Cambridge Press)

            if (index_th.dims(0) > 0) {
                af::array g = exp(-1*input_image(index_th));
                af::array em = af::constant(-1,index_th.dims(0));
                af::array t = af::constant(1,index_th.dims(0));
                index_th = af::seq(0,index_th.dims(0)-1);
                af::array temp;

                //monte carlo rejection loop
                while (index_th.dims(0) > 0) {
                    em(index_th)+=1;
                    //randomize_seed();
                    t(index_th) = t(index_th)*af::randu(index_th.dims(0));

                    temp = af::where(t(index_th) > g(index_th));
                    if (temp.dims(0) > 0) {
                        index_th = index_th(temp);
                    } else {
                        break;
                    }
                }
                input_image(af::where(input_image < threshold)) = em;
            }

        }

        input_image = max(input_image,0.0);

    }

    void apply_noise_model(af::array& input_image){
        //
        //
        //  Applies the chosen noise model
        //
        //

        if (noise_type == "gaussian") {
            apply_gaussian_noise(input_image, gauss_var);
        } else if (noise_type == "multi"){
            apply_multiplicative_noise(input_image, gauss_var);
        } else if (noise_type == "poisson"){
            apply_poisson_noise(input_image, gauss_var);
        } else {
            //std::cout << "WARNING: no noise model set" << std::endl;
        }



    }





};


#endif //SYNIMAGEGEN_NOISEMODEL_HPP
