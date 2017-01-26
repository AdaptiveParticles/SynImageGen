//
// Created by cheesema on 25/01/17.
//

#ifndef SYNIMAGEGEN_SYNIMAGECLASSES_HPP
#define SYNIMAGEGEN_SYNIMAGECLASSES_HPP

#include "arrayfire.h"
#include "MeshDataAF.h"
#include "NoiseModel.hpp"

//definitions
void perform_sep_3D_conv(af::array& main_array,af::array& y_filt, af::array& x_filt,af::array& z_filt);

class Real_domain{
    //
    //  This defines the real domain in which the image lives, its size in actual units
    //
    //
public:

    std::vector<std::vector<float>> dims; //min max of boundaries for objects in the domain
    std::vector<float> size; //length of each side of the domain


    Real_domain()
    {
        //initialize
        dims.resize(3);
        for(int i = 0;i < dims.size();i++){
            dims[i].resize(2,0);
        }
        size.resize(3,0);

    };

    void set_domain_size(float y_min,float y_max,float x_min,float x_max,float z_min,float z_max){
        //
        //  Set the domain size
        //
        dims[0][0] = y_min;
        dims[0][1] = y_max;

        dims[1][0] = x_min;
        dims[1][1] = x_max;

        dims[2][0] = z_min;
        dims[2][1] = z_max;

        size[0] = y_max - y_min;
        size[1] = x_max - x_min;
        size[2] = z_max - z_min;

    }



};
class Real_object{
    //
    //  This defines the real domain in which the image lives, its size in actual units
    //
    //
public:

    //general location and size information
    std::vector<float> location; //location gives the top corner i.e. lower value of each dimension in cartesian for real reference frame
    std::vector<float> dims; // size of the object (may remove)
    unsigned int template_id;

    float int_scale; //the intensity scaling for the object

    //individual parameters for generating the object

    Real_object()
    {
        location.resize(3,0);
        dims.resize(3,0);
        int_scale = 1;
    };


};
class Object_template{
    //
    //  The reference object template, from which images are generated
    //
public:

    MeshDataAF<float> true_object_distribution;
    MeshDataAF<float> img_object_distribution;
    std::vector<float> real_size;
    std::vector<float> real_deltas;
    float density; // flouresence density
    float max_sample;
    float max_sampled_int;

    Object_template()
    {
        real_size.resize(3,0);
        real_deltas.resize(3,0);
    };

    void set_max_sample(float d_){
        density = d_;
        max_sample = real_deltas[0]*real_deltas[1]*real_deltas[2]*density;
    }

    void free(){
        // free up memory for de_allocation
        real_size.resize(3,0);
        real_deltas.resize(3,0);
        true_object_distribution.free();
        img_object_distribution.free();

    }

};

class Sampling_prop{
    //
    //  Sampling properties object
    //
public:

    std::vector<float> voxel_real_dims; //what the voxel size is in each direction
    std::vector<float> sampling_delta; //how often you sample in each direction

    std::vector<float> rot_center; //center of rotation


    Sampling_prop()
    {
        voxel_real_dims.resize(3,0);
        sampling_delta.resize(3,0);
        rot_center.resize(3,0);
    };

    void rotate_around_z(std::vector<float> rot_center,float angle,std::vector<Real_object>& real_objects,std::vector<Object_template>& obj_templates){
        //
        //  Bevan Cheeseman 2016
        //
        //  Rotate the objects by angle around z axis
        //
        //

        // First rotate object templates
        float mid_x,mid_y;

        //this is done by rotating the midpoint, and then resetting the location to be the bottom most point in the new co_ordinate frame
        for (int i = 0; i < real_objects.size(); i++) {
            mid_y = real_objects[i].location[0] + obj_templates[real_objects[i].template_id].real_size[0]/2 - rot_center[0];
            mid_x = real_objects[i].location[1] + obj_templates[real_objects[i].template_id].real_size[1]/2  - rot_center[1];

            real_objects[i].location[0] = mid_x*cos(angle) - mid_y*sin(angle) - obj_templates[real_objects[i].template_id].real_size[0]/2 + rot_center[0];
            real_objects[i].location[1] = mid_x*sin(angle) + mid_y*cos(angle) - obj_templates[real_objects[i].template_id].real_size[1]/2 + rot_center[1];

        }

        //Now we rotate the templates by nearest neighbor interpolation

        af::array x_pos,y_pos;
        af::array x_pos_n,y_pos_n;
        af::array index;
        af::array temp;

        //need to push to device, and handle setting to zero after use as well..


        for (int i = 0; i < obj_templates.size(); i++) {

            //check object is on arrayfire device

            obj_templates[i].true_object_distribution.check_on_arrayfire();

            //rotate it
            obj_templates[i].true_object_distribution.af_mesh = af::rotate(obj_templates[i].true_object_distribution.af_mesh,angle);

            //transfer off the device
            obj_templates[i].true_object_distribution.transfer_from_arrayfire();
            //free the memory
            obj_templates[i].true_object_distribution.free_arrayfire();

        }




    }

    void rotate_around_x(std::vector<float> rot_center,float angle,std::vector<Real_object>& real_objects,std::vector<Object_template>& obj_templates){
        //
        //  Bevan Cheeseman 2016
        //
        //  Rotate the objects by angle around z axis
        //
        //

        // First rotate object templates
        float mid_z,mid_y;

        //this is done by rotating the midpoint, and then resetting the location to be the bottom most point in the new co_ordinate frame
        for (int i = 0; i < real_objects.size(); i++) {
            mid_y = real_objects[i].location[0] + obj_templates[real_objects[i].template_id].real_size[0]/2 - rot_center[0];
            mid_z = real_objects[i].location[2] + obj_templates[real_objects[i].template_id].real_size[2]/2  - rot_center[2];

            real_objects[i].location[0] = mid_z*cos(angle) - mid_y*sin(angle) - obj_templates[real_objects[i].template_id].real_size[0]/2 + rot_center[0];
            real_objects[i].location[2] = mid_z*sin(angle) + mid_y*cos(angle) - obj_templates[real_objects[i].template_id].real_size[2]/2 + rot_center[2];

        }

        //Now we rotate the templates by nearest neighbor interpolation

        af::array x_pos,y_pos;
        af::array x_pos_n,y_pos_n;
        af::array index;
        af::array temp;

        //need to push to device, and handle setting to zero after use as well..


        for (int i = 0; i < obj_templates.size(); i++) {

            //check object is on arrayfire device

            obj_templates[i].true_object_distribution.check_on_arrayfire();

            obj_templates[i].true_object_distribution.af_mesh = af::reorder(obj_templates[i].true_object_distribution.af_mesh,0,2,1);

            //rotate it
            obj_templates[i].true_object_distribution.af_mesh = af::rotate(obj_templates[i].true_object_distribution.af_mesh,angle);


            obj_templates[i].true_object_distribution.af_mesh = af::reorder(obj_templates[i].true_object_distribution.af_mesh,0,2,1);

            //transfer off the device
            obj_templates[i].true_object_distribution.transfer_from_arrayfire();
            //free the memory
            obj_templates[i].true_object_distribution.free_arrayfire();

        }




    }

    void rotate_around_y(std::vector<float> rot_center,float angle,std::vector<Real_object>& real_objects,std::vector<Object_template>& obj_templates){
        //
        //  Bevan Cheeseman 2016
        //
        //  Rotate the objects by angle around z axis
        //
        //

        // First rotate object templates
        float mid_x,mid_z;

        //this is done by rotating the midpoint, and then resetting the location to be the bottom most point in the new co_ordinate frame
        for (int i = 0; i < real_objects.size(); i++) {
            mid_z = real_objects[i].location[2] + obj_templates[real_objects[i].template_id].real_size[2]/2 - rot_center[2];
            mid_x = real_objects[i].location[1] + obj_templates[real_objects[i].template_id].real_size[1]/2  - rot_center[1];

            real_objects[i].location[2] = mid_x*cos(angle) - mid_z*sin(angle) - obj_templates[real_objects[i].template_id].real_size[2]/2 + rot_center[2];
            real_objects[i].location[1] = mid_x*sin(angle) + mid_z*cos(angle) - obj_templates[real_objects[i].template_id].real_size[1]/2 + rot_center[1];

        }

        //Now we rotate the templates by nearest neighbor interpolation

        af::array x_pos,y_pos;
        af::array x_pos_n,y_pos_n;
        af::array index;
        af::array temp;

        //need to push to device, and handle setting to zero after use as well..


        for (int i = 0; i < obj_templates.size(); i++) {

            //check object is on arrayfire device

            obj_templates[i].true_object_distribution.check_on_arrayfire();

            obj_templates[i].true_object_distribution.af_mesh = af::reorder(obj_templates[i].true_object_distribution.af_mesh,0,2,1);

            //rotate it
            obj_templates[i].true_object_distribution.af_mesh = af::rotate(obj_templates[i].true_object_distribution.af_mesh,angle);


            obj_templates[i].true_object_distribution.af_mesh = af::reorder(obj_templates[i].true_object_distribution.af_mesh,0,2,1);

            //transfer off the device
            obj_templates[i].true_object_distribution.transfer_from_arrayfire();
            //free the memory
            obj_templates[i].true_object_distribution.free_arrayfire();

        }




    }


};


class PSF_prop{
    //
    //  Point Spread Function Properties
    //
    //

public:

    std::vector<float> real_sigmas;
    std::vector<float> real_window_size;
    float I0;
    float cut_th;
    float var_del;
    int normalize;


    std::string type;

    PSF_prop(){
        real_sigmas.resize(3,0);
        real_window_size.resize(3,0);
        var_del = 0.05;
        normalize = false;

    }

    void gen_guassian_template_psf(af::array& psf_filter,Object_template& obj);

    void set_guassian_window_size(){
        //
        // Sets a real window size for a guassian PSF with threshold cut_th
        //

        for (int i = 0;  i < real_sigmas.size(); i++) {
            real_window_size[i] = sqrt(log(1/cut_th)*2*pow(real_sigmas[i],2));
        }


    }

    void apply_psf_to_template(Object_template& obj);
    void gen_sep_guassian_template_psf(af::array& psf_filter_y,af::array& psf_filter_x,af::array& psf_filter_z,Object_template& obj);

    void gen_sep_guassian_template_psf_grad(af::array& psf_filter_y,af::array& psf_filter_x,af::array& psf_filter_z,Object_template& obj,int dir);
};

void PSF_prop::apply_psf_to_template(Object_template& obj){
    //
    //  Applies a particular PSF to the template
    //

    if (type == "gauss"){

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;
        //compute PSF here
        // gen_guassian_template_psf(psf_filter,obj);

        //obj.img_object_distribution.af_mesh = convolve(obj.img_object_distribution.af_mesh,psf_filter);

        gen_sep_guassian_template_psf(psf_filter_y,psf_filter_x,psf_filter_z,obj);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);

    } else if (type == "gauss_dx") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;


        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,1);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);


    } else if (type == "gauss_dy") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;


        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,0);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);



    } else if (type == "gauss_dz") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;

        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,2);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);



    } else if (type == "gauss_dsum") {


        af::array psf_filter_x,psf_filter_y,psf_filter_z;

        af::array grad_x,grad_y,grad_z;


        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,0);
        grad_y = obj.img_object_distribution.af_mesh;

        perform_sep_3D_conv(grad_y,psf_filter_y,psf_filter_x,psf_filter_z);
        grad_y = abs(grad_y);

        grad_x = obj.img_object_distribution.af_mesh;


        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,1);
        perform_sep_3D_conv(grad_x,psf_filter_y,psf_filter_x,psf_filter_z);

        grad_x = abs(grad_x);

        grad_y = grad_x + grad_y;

        grad_x = af::constant(0,1);

        grad_z = obj.img_object_distribution.af_mesh;

        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,2);
        perform_sep_3D_conv(grad_z,psf_filter_y,psf_filter_x,psf_filter_z);

        grad_z = abs(grad_z);

        grad_y = grad_z + grad_y;

        grad_z = af::constant(0,1);

        obj.img_object_distribution.af_mesh = grad_y;


    } else if (type == "var_dy") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;

        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,0);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);


        obj.img_object_distribution.af_mesh = obj.max_sample*af::constant(1.0,obj.img_object_distribution.af_mesh.dims())*((obj.img_object_distribution.af_mesh > var_del) - (obj.img_object_distribution.af_mesh < (-var_del))  ) ;


    } else if (type == "var_dx") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;

        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,1);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);


        obj.img_object_distribution.af_mesh = obj.max_sample*af::constant(1.0,obj.img_object_distribution.af_mesh.dims())*((obj.img_object_distribution.af_mesh > var_del) - (obj.img_object_distribution.af_mesh < (-var_del))  ) ;



    } else if (type == "var_dz") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;

        gen_sep_guassian_template_psf_grad(psf_filter_y,psf_filter_x,psf_filter_z,obj,2);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);


        obj.img_object_distribution.af_mesh = obj.max_sample*af::constant(1.0,obj.img_object_distribution.af_mesh.dims())*((obj.img_object_distribution.af_mesh > var_del) - (obj.img_object_distribution.af_mesh < (-var_del))  ) ;



    }else if (type == "var") {

        af::array psf_filter;
        af::array psf_filter_x,psf_filter_y,psf_filter_z;
        //compute PSF here
        // gen_guassian_template_psf(psf_filter,obj);

        //obj.img_object_distribution.af_mesh = convolve(obj.img_object_distribution.af_mesh,psf_filter);

        gen_sep_guassian_template_psf(psf_filter_y,psf_filter_x,psf_filter_z,obj);

        perform_sep_3D_conv(obj.img_object_distribution.af_mesh,psf_filter_y,psf_filter_x,psf_filter_z);

        obj.img_object_distribution.af_mesh = obj.max_sample*af::constant(1.0,obj.img_object_distribution.af_mesh.dims())*((obj.img_object_distribution.af_mesh > var_del) - (obj.img_object_distribution.af_mesh < (-var_del))  ) ;



    } else {

        std::cout << "WARNING: NO PSF SET" << std::endl;


    }
//
//            //set your test data directory here
//    std::string test_data_loc = get_path("IMAGE_GEN_PATH") + "Test_data/";
//
//    MeshDataAF<uint16_t> img_out((int)obj.img_object_distribution.af_mesh.dims(0),(int)obj.img_object_distribution.af_mesh.dims(1),(int)obj.img_object_distribution.af_mesh.dims(2));
//
//    img_out.transfer_to_arrayfire();
//    img_out.af_mesh = abs(obj.img_object_distribution.af_mesh);
//    img_out.transfer_from_arrayfire();
//
//   // write data to tiff
//    std::string tiff_file_name = test_data_loc + "psf_image.tif";
//    write_image_tiff(img_out,tiff_file_name);





}

void PSF_prop::gen_guassian_template_psf(af::array& psf_filter,Object_template& obj){
    //
    //  Bevan Cheeseman 2016
    //
    //  Generates a guassian PSF for a given template model
    //
    //  See: https://en.wikipedia.org/wiki/Airy_disk Approximation by guassian
    //

    std::vector<float> template_window_size;
    template_window_size.resize(3);

    for (int i = 0;  i < real_sigmas.size(); i++) {
        template_window_size[i] = ceil((real_window_size[i]-.5*obj.real_deltas[i])/obj.real_deltas[i]);
    }

    psf_filter = af::array(template_window_size[0]*2 +1,template_window_size[1]*2 +1,template_window_size[2]*2 +1);

    af::array test = af::seq(-15,15);


    af::array y_coords = af::seq(-template_window_size[0],template_window_size[0]);
    af::array x_coords = af::seq(-template_window_size[1],template_window_size[1]);
    af::array z_coords = af::seq(-template_window_size[2],template_window_size[2]);

    y_coords = y_coords*obj.real_deltas[0];
    x_coords = x_coords*obj.real_deltas[1];
    x_coords = moddims(x_coords,1,x_coords.dims(0),1);
    z_coords = z_coords*obj.real_deltas[2];
    z_coords = moddims(z_coords,1,1,z_coords.dims(0));

    //create the grid for the PSF
    y_coords = af::tile(y_coords,1,(unsigned int)psf_filter.dims(1),(unsigned int)psf_filter.dims(2));
    x_coords = af::tile(x_coords,(unsigned int)psf_filter.dims(0),1,(unsigned int)psf_filter.dims(2));
    z_coords = af::tile(z_coords,(unsigned int)psf_filter.dims(0),(unsigned int)psf_filter.dims(1),1);

    psf_filter = (obj.real_deltas[0]*obj.real_deltas[1]*obj.real_deltas[2])*I0*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2))-pow(x_coords,2)/(2*pow(real_sigmas[1],2))-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));

//    //set your test data directory here
//    std::string test_data_loc = "/PHD/ImageGenData/Test_data/";
//
//    MeshDataAF<float> img_out((int)psf_filter.dims(0),(int)psf_filter.dims(1),(int)psf_filter.dims(2));
//
//    img_out.transfer_to_arrayfire();
//    img_out.af_mesh = psf_filter;
//    img_out.transfer_from_arrayfire();
//
//    //write data to tiff
//    std::string tiff_file_name = test_data_loc + "psf_image.tif";
//    write_image_tiff(img_out,tiff_file_name);
//
//    af::print("psf",psf_filter);
//
}

void PSF_prop::gen_sep_guassian_template_psf(af::array& psf_filter_y,af::array& psf_filter_x,af::array& psf_filter_z,Object_template& obj){
    //
    //  Bevan Cheeseman 2016
    //
    //  Generates a guassian PSF for a given template model
    //
    //  See: https://en.wikipedia.org/wiki/Airy_disk Approximation by guassian
    //

    std::vector<float> template_window_size;
    template_window_size.resize(3);

    for (int i = 0;  i < real_sigmas.size(); i++) {
        template_window_size[i] = ceil((real_window_size[i]-.5*obj.real_deltas[i])/obj.real_deltas[i]);
    }

    //psf_filter = af::array(template_window_size[0]*2 +1,template_window_size[1]*2 +1,template_window_size[2]*2 +1);

    af::array test = af::seq(-15,15);


    af::array y_coords = af::seq(-template_window_size[0],template_window_size[0]);
    af::array x_coords = af::seq(-template_window_size[1],template_window_size[1]);
    af::array z_coords = af::seq(-template_window_size[2],template_window_size[2]);

    y_coords = y_coords*obj.real_deltas[0];
    x_coords = x_coords*obj.real_deltas[1];
    z_coords = z_coords*obj.real_deltas[2];


    float factor = pow(I0,1.0/3.00);

    psf_filter_y = (obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));

    psf_filter_x = (obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));

    psf_filter_z = (obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));

}

void PSF_prop::gen_sep_guassian_template_psf_grad(af::array& psf_filter_y,af::array& psf_filter_x,af::array& psf_filter_z,Object_template& obj,int dir){
    //
    //  Bevan Cheeseman 2016
    //
    //  Generates a guassian PSF for a given template model
    //
    //  See: https://en.wikipedia.org/wiki/Airy_disk Approximation by guassian
    //

    std::vector<float> template_window_size;
    template_window_size.resize(3);

    for (int i = 0;  i < real_sigmas.size(); i++) {
        template_window_size[i] = ceil((real_window_size[i]-.5*obj.real_deltas[i])/obj.real_deltas[i]);
    }


    af::array y_coords = af::seq(-template_window_size[0],template_window_size[0]);
    af::array x_coords = af::seq(-template_window_size[1],template_window_size[1]);
    af::array z_coords = af::seq(-template_window_size[2],template_window_size[2]);



    y_coords = y_coords*obj.real_deltas[0];
    x_coords = x_coords*obj.real_deltas[1];
    z_coords = z_coords*obj.real_deltas[2];


    float factor = pow(I0,1.0/3.00);

    if (dir == 0){

        psf_filter_y = (-(y_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2))));
        psf_filter_x = (obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z = (obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));
    }
    else if (dir == 1){
        psf_filter_y = (obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));
        psf_filter_x = -(x_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z = (obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));

    } else if (dir == 2) {
        psf_filter_y = (obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));
        psf_filter_x = (obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z = -(z_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));
    } else if (dir == 3) {
        psf_filter_y = (obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));
        psf_filter_x = (obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z = -(z_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));

        psf_filter_y += (obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));
        psf_filter_x += -(x_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z += (obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));

        psf_filter_y += -(y_coords/(pow(real_sigmas[0],2)))*(obj.real_deltas[0])*factor*exp(-pow(y_coords,2)/(2*pow(real_sigmas[0],2)));
        psf_filter_x += (obj.real_deltas[1])*factor*exp(-pow(x_coords,2)/(2*pow(real_sigmas[1],2)));
        psf_filter_z += (obj.real_deltas[2])*factor*exp(-pow(z_coords,2)/(2*pow(real_sigmas[2],2)));
    }

}



class Global_Trans{

public:

    float const_shift;

    float grad_x;
    float grad_y;
    float grad_z;

    bool gt_ind; //indicataes if global transforms are on or off for the image generation

    //constructor
    Global_Trans()
    {
        const_shift = 0;
        grad_x = 0;
        grad_y = 0;
        grad_z = 0;
        gt_ind = true;
    };


    void apply_global_transforms(af::array& curr_image){
        af::array temp;

        curr_image = curr_image + const_shift;

        //for(int i = 0;i < curr_image.dims(2);i++){
        //temp = curr_image(af::span,af::span,i);
        //  curr_image(af::span,af::span,i) = curr_image(af::span,af::span,i) +  const_shift;
        //}
    }

    void apply_global_linear_shift(af::array& curr_image,Sampling_prop samp_prop){


        if ((grad_x + grad_y + grad_z) != 0){
            //create x grid and then use transpose for y, x_gri
            af::array grid_x = af::iota(af::dim4(1,curr_image.dims(1),1),af::dim4(curr_image.dims(0),1,1));
            grid_x = grid_x*samp_prop.voxel_real_dims[1];
            //
            //
            af::array grid_y = af::iota(af::dim4(curr_image.dims(0),1,1),af::dim4(1,curr_image.dims(1),1));
            grid_y = grid_y*samp_prop.voxel_real_dims[0];

            af::array temp;

            for(int i = 0;i < curr_image.dims(2);i++){
                curr_image(af::span,af::span,i) += grid_x*grad_x +  grid_y*grad_y + i*grad_z*samp_prop.voxel_real_dims[2];
                //curr_image(af::span,af::span,i) = temp + grid_x*grad_x +  grid_y*grad_y + i*grad_z*samp_prop.voxel_real_dims[2];
            }
        }

    }



};


class SynImage{


public:

    Real_domain real_domain;
    std::vector<Real_object> real_objects;
    std::vector<Object_template> object_templates;

    Sampling_prop sampling_properties;
    PSF_prop PSF_properties;
    Noise_Model noise_properties;
    Global_Trans global_trans;
    int ground_t;
    //1 scaling intensity response
    //To add:
    //Time date
    //Unique ID

    //constructor
    SynImage()
    {
        ground_t = 0;
    };

    void calc_voxel_convolution_rect(Object_template& obj_temp,af::array& voxel_template);

    //member function declerations
    template<typename S>
    void generate_syn_image(MeshDataAF<S>& gen_image);

    void sample_imgobj_template(af::array& sampled_imgobj,std::vector<std::vector<float>>& obj_bounds,Real_object& real_obj);

    template<typename S>
    bool apply_boundary_conditions(MeshDataAF<S>& gen_image,std::vector<std::vector<float>>& obj_bounds);

    void free_template_af(){
        for(int i = 0;i < object_templates.size();i++){
            object_templates[i].true_object_distribution.free_arrayfire();
            object_templates[i].img_object_distribution.free_arrayfire();
        }
    }

    void free_template_all(){
        for(int i = 0;i < object_templates.size();i++){
            object_templates[i].true_object_distribution.free_arrayfire();
            object_templates[i].img_object_distribution.free_arrayfire();
        }
        object_templates.clear();
        real_objects.clear();
    }

};

template<typename S>
void SynImage::generate_syn_image(MeshDataAF<S>& gen_image){
    //generates an image from the Syn template

    af::array psf_filter;

    //first load the templates to the gpu and compute the SAT
    for(int i = 0; i < object_templates.size(); i++){
        //create new img_object_distribution that starts with the original object
        object_templates[i].img_object_distribution = object_templates[i].true_object_distribution;
        object_templates[i].img_object_distribution.transfer_to_arrayfire();

        PSF_properties.apply_psf_to_template(object_templates[i]);

        calc_voxel_convolution_rect(object_templates[i], object_templates[i].img_object_distribution.af_mesh);

        object_templates[i].img_object_distribution.transfer_and_free();

    }




    //create the image data structure
    unsigned int genimage_y_num = ceil(real_domain.size[0]/sampling_properties.sampling_delta[0]);
    unsigned int genimage_x_num = ceil(real_domain.size[1]/sampling_properties.sampling_delta[1]);
    unsigned int genimage_z_num = ceil(real_domain.size[2]/sampling_properties.sampling_delta[2]);

    //init, then push to the gpu
    gen_image.initialize(genimage_y_num,genimage_x_num,genimage_z_num,0);
    gen_image.transfer_to_arrayfire();

    gen_image.af_mesh = constant(0,gen_image.af_mesh.dims(),f32);

    //pixel begin and end for template img in each direction
    std::vector<std::vector<float>> curr_obj_bounds(3);
    curr_obj_bounds[0].resize(2);
    curr_obj_bounds[1].resize(2);
    curr_obj_bounds[2].resize(2);

    af::array sampled_imgobj;


    //now sample and place the objects in the image
    for (int i = 0; i < real_objects.size(); i++) {
        //begining index in the syn_image
        curr_obj_bounds[0][0] = floor((real_objects[i].location[0] - real_domain.dims[0][0])/sampling_properties.sampling_delta[0]);
        curr_obj_bounds[1][0] = floor((real_objects[i].location[1] - real_domain.dims[1][0])/sampling_properties.sampling_delta[1]);
        curr_obj_bounds[2][0]  = floor((real_objects[i].location[2] - real_domain.dims[2][0])/sampling_properties.sampling_delta[2]);


        curr_obj_bounds[0][1] = floor((real_objects[i].location[0] + object_templates[real_objects[i].template_id].real_size[0] - real_domain.dims[0][0])/sampling_properties.sampling_delta[0]);
        curr_obj_bounds[1][1] = floor((real_objects[i].location[1] + object_templates[real_objects[i].template_id].real_size[1] - real_domain.dims[1][0])/sampling_properties.sampling_delta[1]);
        curr_obj_bounds[2][1] = floor((real_objects[i].location[2] + object_templates[real_objects[i].template_id].real_size[2] - real_domain.dims[2][0])/sampling_properties.sampling_delta[2]);

        if (apply_boundary_conditions(gen_image,curr_obj_bounds)){

            //gen_image.free_arrayfire();

            sample_imgobj_template(sampled_imgobj,curr_obj_bounds,real_objects[i]);

            //gen_image.check_on_arrayfire();

            // af::array af_temp = gen_image.af_mesh(af::seq(curr_obj_bounds[0][0],curr_obj_bounds[0][1]),af::seq(curr_obj_bounds[1][0],curr_obj_bounds[1][1]),af::seq(curr_obj_bounds[2][0],curr_obj_bounds[2][1]));
            //af::array temp = gen_image.af_mesh(af::seq(curr_obj_bounds[0][0],curr_obj_bounds[0][1]),af::seq(curr_obj_bounds[1][0],curr_obj_bounds[1][1]),af::seq(curr_obj_bounds[2][0],curr_obj_bounds[2][1]));

            //temp += sampled_imgobj;

            af::array temp = gen_image.af_mesh(af::seq(curr_obj_bounds[0][0],curr_obj_bounds[0][1]),af::seq(curr_obj_bounds[1][0],curr_obj_bounds[1][1]),af::seq(curr_obj_bounds[2][0],curr_obj_bounds[2][1]));

            gen_image.af_mesh(af::seq(curr_obj_bounds[0][0],curr_obj_bounds[0][1]),af::seq(curr_obj_bounds[1][0],curr_obj_bounds[1][1]),af::seq(curr_obj_bounds[2][0],curr_obj_bounds[2][1])) = temp +  sampled_imgobj;

            //for(int k = curr_obj_bounds[2][0];k <= curr_obj_bounds[2][1];k++){
            //    gen_image.af_mesh(af::seq(curr_obj_bounds[0][0],curr_obj_bounds[0][1]),af::seq(curr_obj_bounds[1][0],curr_obj_bounds[1][1]),k) += sampled_imgobj(af::span,af::span,k);
            //}
        }

    }
    //free up the templates from gpu memory if already hasn't hapepended
    //free_template_af();

    if (global_trans.gt_ind == true){

        global_trans.apply_global_transforms(gen_image.af_mesh);

        global_trans.apply_global_linear_shift(gen_image.af_mesh,sampling_properties);

    }

    noise_properties.apply_noise_model(gen_image.af_mesh);

    int num_byts = sizeof(S);

    //if not float or double it is uint, so get abs
    if (num_byts < 4){
        gen_image.af_mesh = abs(gen_image.af_mesh);
    }


    gen_image.transfer_from_arrayfire();


}
void SynImage::calc_voxel_convolution_rect(Object_template& obj_temp,af::array& voxel_template){
    //
    //  Bevan Cheeseman 2016
    //
    //  Calculates the voxel integral for square voxels
    //

    //size of the voxel integral
    float y_dim,x_dim,z_dim;

    y_dim = std::max((float)round(sampling_properties.voxel_real_dims[0]/obj_temp.real_deltas[0]),(float)1);
    x_dim = std::max((float)round(sampling_properties.voxel_real_dims[1]/obj_temp.real_deltas[1]),(float)1);
    z_dim = std::max((float)round(sampling_properties.voxel_real_dims[2]/obj_temp.real_deltas[2]),(float)1);

    af::array voxel_filter(y_dim,x_dim,z_dim);

    af::array voxel_filter_x(x_dim);
    af::array voxel_filter_y(y_dim);
    af::array voxel_filter_z(z_dim);

    //voxel_filter
    voxel_filter(af::span) = 1;

    voxel_filter_x(af::span) = 1;
    voxel_filter_y(af::span) = 1;
    voxel_filter_z(af::span) = 1;


    //convolve by the filter
    if ((y_dim + x_dim + z_dim) > 3) {
        //y filter
        perform_sep_3D_conv(voxel_template,voxel_filter_y,voxel_filter_x,voxel_filter_z);

    }

    obj_temp.max_sampled_int = obj_temp.max_sample*y_dim*x_dim*z_dim;

}
void SynImage::sample_imgobj_template(af::array& sampled_imgobj,std::vector<std::vector<float>>& obj_bounds,Real_object& real_obj){
    //
    //  Bevan Cheeseman 2016
    //
    //  Samples the template at the voxel locations
    //
    //


    object_templates[real_obj.template_id].img_object_distribution.check_on_arrayfire();
    //sequence of location then in local reference frame

    af::array y_points = af::seq(obj_bounds[0][0],obj_bounds[0][1]);

    y_points = y_points*sampling_properties.sampling_delta[0] + real_domain.dims[0][0] - real_obj.location[0];

    //then find the nearest index in the template
    y_points = round(y_points/object_templates[real_obj.template_id].real_deltas[0]);

    af::array x_points = af::seq(obj_bounds[1][0],obj_bounds[1][1]);

    x_points = x_points*sampling_properties.sampling_delta[1] + real_domain.dims[1][0] - real_obj.location[1];

    //then find the nearest index in the template
    x_points = round(x_points/object_templates[real_obj.template_id].real_deltas[1]);

    af::array z_points = af::seq(obj_bounds[2][0],obj_bounds[2][1]);

    z_points = z_points*sampling_properties.sampling_delta[2] + real_domain.dims[2][0] - real_obj.location[2];

    //then find the nearest index in the template
    z_points = round(z_points/object_templates[real_obj.template_id].real_deltas[2]);

    if (PSF_properties.normalize == true) {
        sampled_imgobj = (1000.0/object_templates[real_obj.template_id].max_sampled_int)*object_templates[real_obj.template_id].img_object_distribution.af_mesh(y_points,x_points,z_points);
    } else if(ground_t == 1){
        //sampled_imgobj = real_obj.template_id;
        af::array temp = object_templates[real_obj.template_id].img_object_distribution.af_mesh(y_points,x_points,z_points)/object_templates[real_obj.template_id].img_object_distribution.af_mesh(y_points,x_points,z_points);
        sampled_imgobj = real_obj.template_id*1.0*temp;
    }
    else {

        sampled_imgobj = (real_obj.int_scale)*object_templates[real_obj.template_id].img_object_distribution.af_mesh(y_points,x_points,z_points);

    }
    object_templates[real_obj.template_id].img_object_distribution.transfer_and_free();


}
template<typename S>
bool SynImage::apply_boundary_conditions(MeshDataAF<S>& gen_image,std::vector<std::vector<float>>& obj_bounds){
    //
    //  Just applies basic boundary conditions
    //

    if ((obj_bounds[0][0] >= gen_image.y_num) | (obj_bounds[1][0] >= gen_image.x_num) | (obj_bounds[2][0] >= gen_image.z_num) | (obj_bounds[0][1] < 0) | (obj_bounds[1][1] < 0)| (obj_bounds[2][1] < 0)) {
        return false;
    } else {

        //y boundaries
        if(obj_bounds[0][0] < 0){
            obj_bounds[0][0] = 0;
        }
        if(obj_bounds[0][1] >= gen_image.y_num){
            obj_bounds[0][1] = gen_image.y_num - 1;
        }

        //x boundaries
        if(obj_bounds[1][0] < 0){
            obj_bounds[1][0] = 0;
        }
        if(obj_bounds[1][1] >= gen_image.x_num){
            obj_bounds[1][1] = gen_image.x_num - 1;
        }

        //z boundaries
        if(obj_bounds[2][0] < 0){
            obj_bounds[2][0] = 0;
        }
        if(obj_bounds[2][1] >= gen_image.z_num){
            obj_bounds[2][1] = gen_image.z_num - 1;
        }
        return true;

    }

}

void perform_sep_3D_conv(af::array& main_array,af::array& y_filt, af::array& x_filt,af::array& z_filt){
//
//  Bevan Cheeseman 2016
//
//  Performs 3D filter convolution on 3D images, using seperable convolution
//

    main_array = convolve(main_array,y_filt);
//z
    main_array = af::reorder(main_array,2,1,0);
    main_array = convolve(main_array,z_filt);
    main_array = af::reorder(main_array,2,1,0);

//along x
    main_array = main_array.T();
    main_array = convolve(main_array,x_filt);
    main_array = main_array.T();


}


#endif //SYNIMAGEGEN_SYNIMAGECLASSES_HPP
