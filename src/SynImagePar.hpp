//
// Created by cheesema on 25/01/17.
//

#ifndef SYNIMAGEGEN_SYNIMAGEPAR_HPP
#define SYNIMAGEGEN_SYNIMAGEPAR_HPP

#include <iostream>

class Syn_gen_par{
    //
    //  Bevan Cheeseman
    //
    //  Synthetic image generation parameter class
    //
    //

public:


    //new parameters

    float dy;
    float dx;
    float dz;

    float psfy;
    float psfx;
    float psfz;

    float ydim;
    float xdim;
    float zdim;

    float snr_ratio;
    float background;

    //template file name location
    std::string template_name;

    //name of the syn images to be produced
    std::string name;

    //number of objects in the images
    int num_objects;

    float object_size;

    //noise model currently default is poisson
    std::string noise_model;

    //psf model
    std::string psf_model;

    //pipeline parameters
    Syn_gen_par()
    {};


};


void gen_parameter_file_txt(SynImage& syn_image,std::string save_path,std::string image_name){

    std::string image_stats_filename = save_path + image_name +  "_stats.txt";

    //output summary statistics
    std::ofstream myfile;
    myfile.open (image_stats_filename);

    myfile << "name: " << image_name << std::endl;
    myfile << "dx: " << syn_image.sampling_properties.voxel_real_dims[1] << std::endl;
    myfile << "dy: " << syn_image.sampling_properties.voxel_real_dims[0] << std::endl;
    myfile << "dz: " << syn_image.sampling_properties.voxel_real_dims[2] << std::endl;

    myfile << "xdim: " << syn_image.real_domain.size[1] << std::endl;
    myfile << "ydim: " << syn_image.real_domain.size[0] << std::endl;
    myfile << "zdim: " << syn_image.real_domain.size[2] << std::endl;

    myfile << "psfx: " << syn_image.PSF_properties.real_sigmas[1] << std::endl;
    myfile << "psfy: " << syn_image.PSF_properties.real_sigmas[0] << std::endl;
    myfile << "psfz: " << syn_image.PSF_properties.real_sigmas[2] << std::endl;

    myfile << "noise_sigma: " << sqrt(syn_image.noise_properties.gauss_var) << std::endl;

    myfile << "background: " << syn_image.global_trans.const_shift << std::endl;


}
void gen_parameter_file_txt_ground_truth(SynImage& syn_image,std::string save_path,Syn_gen_par syn_pars){

    std::string image_stats_filename = save_path + syn_pars.name +  "_syn_config.txt";

    //output summary statistics
    std::ofstream myfile;
    myfile.open (image_stats_filename);

    myfile << "name: " << syn_pars.name << std::endl;
    myfile << "dx: " << syn_image.sampling_properties.voxel_real_dims[1] << std::endl;
    myfile << "dy: " << syn_image.sampling_properties.voxel_real_dims[0] << std::endl;
    myfile << "dz: " << syn_image.sampling_properties.voxel_real_dims[2] << std::endl;

    myfile << "xdim: " << syn_image.real_domain.size[1] << std::endl;
    myfile << "ydim: " << syn_image.real_domain.size[0] << std::endl;
    myfile << "zdim: " << syn_image.real_domain.size[2] << std::endl;

    myfile << "psfx: " << syn_image.PSF_properties.real_sigmas[1] << std::endl;
    myfile << "psfy: " << syn_image.PSF_properties.real_sigmas[0] << std::endl;
    myfile << "psfz: " << syn_image.PSF_properties.real_sigmas[2] << std::endl;

    myfile << "snr_ratio: " << sqrt(syn_image.noise_properties.gauss_var) << std::endl;

    myfile << "background: " << syn_image.global_trans.const_shift << std::endl;

    myfile << "num_objects: " << syn_pars.num_objects << std::endl;

    myfile << "object_size: " << syn_pars.object_size << std::endl;

    myfile << "template_name: " << syn_pars.template_name << std::endl;


}
/*void gen_parameter_pars(SynImage& syn_image,Proc_par& pars,std::string image_name){
    //
    //
    //  Takes in the SynImage model parameters and outputs them to the APR parameter class
    //
    //
    //

    pars.name = image_name;

    pars.dy = syn_image.sampling_properties.voxel_real_dims[0];
    pars.dx = syn_image.sampling_properties.voxel_real_dims[1];
    pars.dz = syn_image.sampling_properties.voxel_real_dims[2];

    pars.psfy = syn_image.PSF_properties.real_sigmas[0];
    pars.psfx = syn_image.PSF_properties.real_sigmas[1];
    pars.psfz = syn_image.PSF_properties.real_sigmas[2];

    pars.ydim = syn_image.real_domain.size[0];
    pars.xdim = syn_image.real_domain.size[1];
    pars.zdim = syn_image.real_domain.size[2];

    pars.noise_sigma = sqrt(syn_image.noise_properties.gauss_var);
    pars.background = syn_image.global_trans.const_shift;

}*/
void get_image_par(Syn_gen_par& syn_pars,std::string output_path,std::string image_name){
    //
    //  Gets the image parameters for the specific file to be run
    //
    //  Bevan Cheeseman 2016
    //
    //
    std::cout << output_path + image_name + "_syn_config.txt" << std::endl;


    //open the files
    std::ifstream path_file;
    path_file.open (output_path + image_name + "_syn_config.txt");

    std::string out_line;

    std::size_t found;

    //get the paths

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("name: ");

    if (found!=std::string::npos){

        syn_pars.name = out_line.substr(found+6);
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }



    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("dx: ");

    if (found!=std::string::npos){

        syn_pars.dx = stof(out_line.substr(found+4));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("dy: ");

    if (found!=std::string::npos){

        syn_pars.dy = stof(out_line.substr(found+4));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("dz: ");

    if (found!=std::string::npos){

        syn_pars.dz = stof(out_line.substr(found+4));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("xdim: ");

    if (found!=std::string::npos){

        syn_pars.xdim = stof(out_line.substr(found+6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("ydim: ");

    if (found!=std::string::npos){

        syn_pars.ydim = stof(out_line.substr(found + 6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("zdim: ");

    if (found!=std::string::npos){

        syn_pars.zdim = stof(out_line.substr(found + 6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("psfx: ");

    if (found!=std::string::npos){

        syn_pars.psfx = stof(out_line.substr(found + 6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("psfy: ");

    if (found!=std::string::npos){

        syn_pars.psfy = stof(out_line.substr(found+6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("psfz: ");

    if (found!=std::string::npos){

        syn_pars.psfz = stof(out_line.substr(found+6));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("snr_ratio: ");

    if (found!=std::string::npos){

        syn_pars.snr_ratio = stof(out_line.substr(found+11));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("background: ");

    if (found!=std::string::npos){

        syn_pars.background = stof(out_line.substr(found+12));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }


    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("num_objects: ");

    if (found!=std::string::npos){

        syn_pars.num_objects = stof(out_line.substr(found+13));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("object_size: ");

    if (found!=std::string::npos){

        syn_pars.object_size = stof(out_line.substr(found+12));
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }


    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("template_name: ");

    if (found!=std::string::npos){

        syn_pars.template_name = out_line.substr(found+15);
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("noise_model: ");

    if (found!=std::string::npos){

        syn_pars.noise_model = out_line.substr(found+13);
    } else {

        std::cout << "Setting file incomplete" << std::endl;

    }

    //file name (relative path)
    std::getline(path_file,out_line);

    found = out_line.find("psf_model: ");

    if (found!=std::string::npos){

        syn_pars.psf_model = out_line.substr(found+11);
    } else {

        std::cout << "Setting file incomplete" << std::endl;
    }


}
void set_syn_image_prop(SynImage& syn_image,Syn_gen_par& syn_pars){
    //
    //
    //  Bevan Cheeseman 2016
    //
    //  Sets up the synimage properties
    //
    //
    //



    std::string image_name = syn_pars.name;

    /////////////////////////////////////////
    //////////////////////////////////////////
    // SET UP THE DOMAIN SIZE



    ///////////////////////////////////////////////////////////////////
    //
    //  sampling properties


    syn_image.sampling_properties.voxel_real_dims[0] = syn_pars.dy;
    syn_image.sampling_properties.voxel_real_dims[1] = syn_pars.dx;
    syn_image.sampling_properties.voxel_real_dims[2] = syn_pars.dz;

    //sampling rate/delta
    syn_image.sampling_properties.sampling_delta[0] = syn_pars.dy;
    syn_image.sampling_properties.sampling_delta[1] = syn_pars.dx;
    syn_image.sampling_properties.sampling_delta[2] = syn_pars.dz;

    int x_num = round(syn_pars.xdim/syn_pars.dx);
    int y_num = round(syn_pars.ydim/syn_pars.dy);
    int z_num = round(syn_pars.zdim/syn_pars.dz);


    //real size of domain
    float dom_size_y = y_num*syn_pars.dy;
    float dom_size_x = x_num*syn_pars.dx;
    float dom_size_z = z_num*syn_pars.dz;
    syn_image.real_domain.set_domain_size(0, dom_size_y, 0, dom_size_x, 0, dom_size_z);

    ///////////////////////////////////////////////////////////////////
    //PSF properties
    syn_image.PSF_properties.real_sigmas[0] = syn_pars.psfy;
    syn_image.PSF_properties.real_sigmas[1] = syn_pars.psfx;
    syn_image.PSF_properties.real_sigmas[2] = syn_pars.psfz;

    syn_image.PSF_properties.I0 = 1/(pow(2*3.14159265359,1.5)*syn_image.PSF_properties.real_sigmas[0]*syn_image.PSF_properties.real_sigmas[1]*syn_image.PSF_properties.real_sigmas[2]);

    syn_image.PSF_properties.cut_th = 0.01;

    syn_image.PSF_properties.set_guassian_window_size();

    syn_image.PSF_properties.type = syn_pars.psf_model;


    ////////////////////////////////////////////////////
    // Global Transforms

    syn_image.global_trans.const_shift = syn_pars.background;

    ///////////////////////////////////////////////////
    //Noise properties

    syn_image.noise_properties.gauss_var = sqrt(syn_pars.background);
    syn_image.noise_properties.noise_type = syn_pars.noise_model;


}


#endif //SYNIMAGEGEN_SYNIMAGEPAR_HPP
