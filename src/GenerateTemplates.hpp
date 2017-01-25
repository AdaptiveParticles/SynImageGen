//
// Created by cheesema on 25/01/17.
//

#ifndef SYNIMAGEGEN_GENERATETEMPLATES_HPP
#define SYNIMAGEGEN_GENERATETEMPLATES_HPP

#include "SynImageClasses.hpp"
#include <iostream>
#include <fstream>

void generate_sphere_template(Object_template& basic_sphere,int sample_rate,float real_size,float density = 1000000,float rad_ratio = 3.0/8.0){
    //
    //
    //  Generate a simple binary sphere for testing
    //
    //
    //
    //

    //set size
    basic_sphere.real_size.resize(3);
    basic_sphere.real_size[0] = real_size;
    basic_sphere.real_size[1] = real_size;
    basic_sphere.real_size[2] = real_size;

    basic_sphere.real_deltas.resize(3);

    basic_sphere.real_deltas[0] = real_size/sample_rate;
    basic_sphere.real_deltas[1] = real_size/sample_rate;
    basic_sphere.real_deltas[2] = real_size/sample_rate;

    basic_sphere.set_max_sample(density);



    //init the data
    basic_sphere.true_object_distribution.initialize(sample_rate, sample_rate, sample_rate, 0);



    float center_sphere = sample_rate/2;
    float radius_sphere = sample_rate*rad_ratio;

    float curr_dist;


    for (int i = 0; i < basic_sphere.true_object_distribution.y_num; i++) {
        for (int j = 0; j < basic_sphere.true_object_distribution.x_num; j++) {
            for (int k = 0; k < basic_sphere.true_object_distribution.z_num; k++) {

                curr_dist = sqrt(pow(i - center_sphere,2) + pow(j - center_sphere,2) + pow(k - center_sphere,2));

                if (curr_dist < radius_sphere) {
                    basic_sphere.true_object_distribution(i,j,k) = 1*basic_sphere.max_sample;
                } else {
                    basic_sphere.true_object_distribution(i,j,k) = 0;
                }
            }
        }
    }


    //  basic_sphere.true_object_distribution.transfer_to_arrayfire();
    //   af::array temp = basic_sphere.true_object_distribution.af_mesh;

//            std::string test_data_loc = get_path("PARTGEN_OUTPUT_PATH") ;
//        //
//            MeshDataAF<uint16_t> img_out((int)temp.dims(0),(int)temp.dims(1),(int)temp.dims(2));
//        //
//            img_out.transfer_to_arrayfire();
//            img_out.af_mesh = temp;
//            img_out.transfer_from_arrayfire();
//        //
//        //    //write data to tiff
//            std::string tiff_file_name = test_data_loc + "template.tif";
//            write_image_tiff(img_out,tiff_file_name);


}
void generate_square_template(Object_template& basic_square,int sample_rate,float real_size,float density = 1000000){
    //
    //
    //  Generate a simple binary sphere for testing
    //
    //
    //
    //


    //set size
    basic_square.real_size.resize(3);
    basic_square.real_size[0] = real_size;
    basic_square.real_size[1] = real_size;
    basic_square.real_size[2] = real_size;

    basic_square.real_deltas.resize(3);

    basic_square.real_deltas[0] = real_size/sample_rate;
    basic_square.real_deltas[1] = real_size/sample_rate;
    basic_square.real_deltas[2] = real_size/sample_rate;

    basic_square.set_max_sample(density);

    //init the data
    basic_square.true_object_distribution.initialize(sample_rate, sample_rate, sample_rate, 0);

    float center_square = sample_rate/2;
    float radius_square = 2*sample_rate/8;

    float curr_dist;

    for (int i = 0; i < basic_square.true_object_distribution.y_num; i++) {
        for (int j = 0; j < basic_square.true_object_distribution.x_num; j++) {
            for (int k = 0; k < basic_square.true_object_distribution.z_num; k++) {

                curr_dist = std::max(std::abs(i - center_square),std::max(std::abs(j - center_square), std::abs(k - center_square)));

                if (curr_dist < radius_square) {
                    basic_square.true_object_distribution(i,j,k) = 1*basic_square.max_sample;
                } else {
                    basic_square.true_object_distribution(i,j,k) = 0;
                }
            }
        }
    }

}

//
// This example program reads a .binvox file and writes
// an ASCII version of the same file called "voxels.txt"
//
// 0 = empty voxel
// 1 = filled voxel
// A newline is output after every "dim" voxels (depth = height = width = dim)
//
// Note that this ASCII version is not supported by "viewvox" and "thinvox"
//
// The x-axis is the most significant axis, then the z-axis, then the y-axis.
//


typedef unsigned char byte;



int read_binvox(std::string filespec,MeshDataAF<uint8_t>& bin_voxels)
{



    static int version;
    static int depth, height, width;
    static int size;
    static byte *voxels = 0;
    static float tx, ty, tz;
    static float scale;

    std::ifstream *input = new std::ifstream(filespec.c_str(), std::ios::in | std::ios::binary);

    //
    // read header
    //
    std::string line;
    *input >> line;  // #binvox
    if (line.compare("#binvox") != 0) {
        std::cout << "Error: first line reads [" << line << "] instead of [#binvox]" << std::endl;
        delete input;
        return 0;
    }
    *input >> version;
    std::cout << "reading binvox version " << version << std::endl;

    depth = -1;
    int done = 0;
    while(input->good() && !done) {
        *input >> line;
        if (line.compare("data") == 0) done = 1;
        else if (line.compare("dim") == 0) {
            *input >> depth >> height >> width;
        }
        else if (line.compare("translate") == 0) {
            *input >> tx >> ty >> tz;
        }
        else if (line.compare("scale") == 0) {
            *input >> scale;
        }
        else {
            std::cout << "  unrecognized keyword [" << line << "], skipping" << std::endl;
            char c;
            do {  // skip until end of line
                c = input->get();
            } while(input->good() && (c != '\n'));

        }
    }
    if (!done) {
        std::cout << "  error reading header" << std::endl;
        return 0;
    }
    if (depth == -1) {
        std::cout << "  missing dimensions in header" << std::endl;
        return 0;
    }

    size = width * height * depth;
    voxels = new byte[size];
    if (!voxels) {
        std::cout << "  error allocating memory" << std::endl;
        return 0;
    }


    bin_voxels.initialize(width, height, depth, 0);

    //
    // read voxel data
    //
    byte value;
    byte count;
    int index = 0;
    int end_index = 0;
    int nr_voxels = 0;

    input->unsetf(std::ios::skipws);  // need to read every byte now (!)
    *input >> value;  // read the linefeed char

    while((end_index < size) && input->good()) {
        *input >> value >> count;

        if (input->good()) {
            end_index = index + count;
            if (end_index > size) return 0;
            for(int i=index; i < end_index; i++) {
                bin_voxels.mesh[i] = (uint8_t) value;
            }

            if (value) nr_voxels += count;
            index = end_index;
        }  // if file still ok

    }  // while

    input->close();
    std::cout << "  read " << nr_voxels << " voxels" << std::endl;

    return 1;

}
template <typename T>
void create_template_from_file(std::string& template_path,Object_template& obj_template,std::vector<float>& obj_size,float& density,float buff_size = 0.05){
//
//  Bevan Cheeseman 2016
//
//  Takes a template file from file and generates a template
//
//  Could be a tif or an hdf5 file
//
//

//some parameters
//float buff_size = 0.05; //how much room to give either side of the original template (to allow for PSF)

//figure out what type of file it is
    std::size_t found_tif = template_path.find(".tif");
    std::size_t found_h5 = template_path.find(".h5");

//get the template
    MeshDataAF<T> template_data;

    if (found_tif!=std::string::npos){
        load_image_tiff(template_data,template_path);
    } else if (found_h5!=std::string::npos)
        load_mesh_hdf5(template_data,template_path);
    else {
        //neither tiff nor h5 found
        std::cout << "WARNING: UNKNOWN TEMPLATE FILE TYPE" << std::endl;
        return;
    }


    std::cout << template_data.x_num << std::endl;

    std::cout << round((1 + 2*buff_size)*template_data.y_num) << std::endl;

    obj_template.true_object_distribution.initialize(round((1 + 2*buff_size)*template_data.y_num),round((1 + 2*buff_size)*template_data.x_num),round((1 + 2*buff_size)*template_data.z_num),0);

    obj_template.real_size[0] = obj_size[0]*(1 + 2*buff_size);
    obj_template.real_size[1] = obj_size[1]*(1 + 2*buff_size);
    obj_template.real_size[2] = obj_size[2]*(1 + 2*buff_size);

    obj_template.real_deltas[0] = obj_template.real_size[0]/obj_template.true_object_distribution.y_num;
    obj_template.real_deltas[1] = obj_template.real_size[1]/obj_template.true_object_distribution.x_num;
    obj_template.real_deltas[2] = obj_template.real_size[2]/obj_template.true_object_distribution.z_num;

//calcualate then the max_sample
    obj_template.set_max_sample(density);

//////////////////////////////////
//copy template into the data

    template_data.transfer_to_arrayfire();
    obj_template.true_object_distribution.transfer_to_arrayfire();

//obj_template.true_object_distribution.af_mesh = template_data.af_mesh;

//fill in the array
    obj_template.true_object_distribution.af_mesh(af::seq(round(buff_size*template_data.y_num),round(buff_size*template_data.y_num) + template_data.y_num - 1) ,af::seq(round(buff_size*template_data.x_num),
                                                                                                                                                                        round(buff_size*template_data.x_num) + template_data.x_num - 1),af::seq(round(buff_size*template_data.z_num),round(buff_size*template_data.z_num) + template_data.z_num-1)) =
            template_data.af_mesh*obj_template.max_sample;

//transfer back to device
    obj_template.true_object_distribution.transfer_from_arrayfire();

//free up the GPU memory
    obj_template.true_object_distribution.free_arrayfire();

    template_data.free_arrayfire();



}

#endif //SYNIMAGEGEN_GENERATETEMPLATES_HPP
