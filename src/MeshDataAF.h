//
// Created by bevanc on 25.01.17.
//

#ifndef SYNIMAGEGEN_MESHDATAAF_H
#define SYNIMAGEGEN_MESHDATAAF_H

#include <vector>
#include <algorithm>
#include "arrayfire.h"

template<class T>
class MeshDataAF {
    //Defines what a particle is and what characteristics it has
public :

    int y_num;
    int x_num;
    int z_num;

    //data on local
    std::vector<T> mesh;

    T type;

    std::vector<float> mesh_float;
    //data on arrayfire device
    af::array af_mesh;

    af::dtype data_t;


    MeshDataAF()
            :y_num(0),x_num(0),z_num(0)
    {}

    MeshDataAF(int y_num,int x_num,int z_num)
            :y_num(y_num),x_num(x_num),z_num(z_num)
    {
        mesh.resize(y_num*x_num*z_num);
        //mesh.resize(y_num,std::vector<std::vector<T> >(x_num,std::vector<T>(z_num)));
    }

    T& operator ()(int i, int j,int k){
        return mesh[y_num*(j) + i + (k)*x_num*y_num];
        //return mesh[i][j][k];
    }

    void transfer_to_arrayfire(){
        //function that transfers the current array on the pc to gpu

        //need to change it to a float datatype to work with arrayfire
        //std::copy(mesh.begin(),mesh.end(), back_inserter(mesh_float));

        af_mesh = af::array(y_num,x_num,z_num,&mesh[0]);

        data_t = af_mesh.type();


    }

    void set_size(int y_num_,int x_num_,int z_num_){

        y_num = y_num_;
        x_num = x_num_;
        z_num = z_num_;
    }

    void transfer_from_arrayfire(){
        //function that transfers the current array from the gpu to cpu


        af_mesh = af_mesh.as(data_t);

        T *h_A;
        h_A = af_mesh.host<T>();

        mesh.resize(af_mesh.elements());

        std::copy(h_A, h_A + af_mesh.elements(),mesh.begin());

        delete[] h_A;

    }

    void transfer_from_diff_array(af::array& d_af_mesh,af::dtype dt){
        //function that transfers the current array from the gpu to cpu


        d_af_mesh = d_af_mesh.as(dt);

        T *h_A;
        h_A = d_af_mesh.host<T>();



        mesh.resize(d_af_mesh.elements());


        std::copy(h_A, h_A + d_af_mesh.elements(),mesh.begin());

        delete[] h_A;


    }

    void transfer_from_diff_array_add(af::array& d_af_mesh,af::dtype dt){
        //function that transfers the current array from the gpu to cpu


        d_af_mesh = d_af_mesh.as(dt);

        T *h_A;
        h_A = d_af_mesh.host<T>();

        mesh.resize(d_af_mesh.elements());

        std::transform(h_A, h_A + d_af_mesh.elements(), mesh.begin(),
                       mesh.begin(), std::plus<T>());

        delete[] h_A;



    }

    void free(){
        //free memory
        mesh.resize(0);
        //free arrayfire
        af_mesh = af::constant(0,1);

    }

    void free_arrayfire(){
        //sets it to a zero size arraoy


        af_mesh = af::constant(0,1);

    }

    void check_on_arrayfire(){
        if ((af_mesh.dims(0) == y_num) & (af_mesh.dims(1) == x_num) & (af_mesh.dims(2) == z_num)) {

        } else {
            transfer_to_arrayfire();
        }

    }

    void transfer_and_free(){
        if ((af_mesh.dims(0) + af_mesh.dims(1) + af_mesh.dims(2)) > 0) {
            transfer_from_arrayfire();
            free_arrayfire();
        }

    }

    void initialize(T val)
    {
        mesh.resize(y_num*x_num*z_num,val);
        //mesh.insert(mesh.begin(),y_num*x_num*z_num,val);
        //mesh.resize(y_num,std::vector<std::vector<T> >(x_num,std::vector<T>(z_num)));
    }

    void initialize(int y_num_,int x_num_,int z_num_,T val)
    {
        y_num = y_num_;
        x_num = x_num_;
        z_num = z_num_;

        mesh.resize(y_num*x_num*z_num,val);
        //mesh.insert(mesh.begin(),y_num*x_num*z_num,val);
        //mesh.resize(y_num,std::vector<std::vector<T> >(x_num,std::vector<T>(z_num)));
    }


    void zero()
    {

        std::vector<T>().swap(mesh);
    }



    void setzero()
    {

        std::fill(mesh.begin(), mesh.end(), 0);
    }

    void setones()
    {

        std::fill(mesh.begin(), mesh.end(), 1.0);
    }

    void transpose(){

        std::vector<T> v2;
        std::swap(mesh, v2);

        for( unsigned int k = 0; k < z_num;k++){
            for (unsigned int i = 0; i < y_num; i++) {
                for (unsigned int j = 0; j < x_num; j++) {
                    mesh.push_back(v2[k*x_num*y_num + j * y_num + i]);
                }
            }
        }

        y_num = x_num;
        x_num = y_num;

    }


    void reorder_af(int dim0,int dim1,int dim2){


        transfer_to_arrayfire();

        af_mesh = reorder(af_mesh,dim0,dim1,dim2);

        transfer_from_arrayfire();
    }

};


#endif //SYNIMAGEGEN_MESHDATAAF_H
