#pragma once

#include <vector>

// #include <mat.h> // matlab implementation
#include <matio.h> // open source implementation
#include <string>

class MATLAB_MATIO_output {
public:

    MATLAB_MATIO_output(const std::string &name);

    ~MATLAB_MATIO_output();

    void output_vector(const std::vector<int> &output_vec,
                       const std::string &name);

    void output_vector(const std::vector<double> &output_vec,
                       const std::string &name);

    void output_scalar(double scalar,
                       const std::string &name);

    void output_scalar(long scalar,
                       const std::string &name);

    void output_scalar(int scalar,
                       const std::string &name);

    void output_string(const std::string &text,
                       const std::string &field_name_str);

    // TEMPLATES FOR OUTPUT OF arrays OF 2D OR MORE

    //

    template<size_t nb_l, size_t nb_c>
    void output_2D_matrix(const int array2d_in[nb_l][nb_c], const std::string &name) {
        int array2d[nb_l][nb_c] = {0};

        // fill 2d array
        for (int i = 0; i < nb_l; i++)
            for (int j = 0; j < nb_c; j++) array2d[i][j] = array2d_in[i][j];

        // write
        char *fieldname2d = (char *) name.c_str();
        size_t dim2d[2] = {nb_l, nb_c};
        matvar_t *variable2d = Mat_VarCreate(fieldname2d, MAT_C_INT32, MAT_T_INT32, 2, dim2d, &array2d, 0); // rank 2
        Mat_VarWrite(matfp, variable2d,
                     MAT_COMPRESSION_ZLIB);                                              // or MAT_COMPRESSION_NONE
        Mat_VarFree(variable2d);
    }

    //

    template<size_t nb_l, size_t nb_c>
    void output_2D_matrix(const double array2d_in[nb_l][nb_c], const std::string &name) {
        double array2d[nb_l][nb_c] = {0};

        // fill 2d array
        for (int i = 0; i < nb_l; i++)
            for (int j = 0; j < nb_c; j++) array2d[i][j] = array2d_in[i][j];

        // write
        char *fieldname2d = (char *) name.c_str();
        size_t dim2d[2] = {nb_l, nb_c};
        matvar_t *variable2d = Mat_VarCreate(fieldname2d, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, &array2d, 0); // rank 2
        Mat_VarWrite(matfp, variable2d,
                     MAT_COMPRESSION_ZLIB);                                                // or MAT_COMPRESSION_NONE
        Mat_VarFree(variable2d);
    }

    //

    template<size_t nb_l, size_t nb_c, size_t nb_d>
    void output_3D_matrix(const int array3d_in[nb_l][nb_c][nb_d], const std::string &name) {
        int array3d[nb_l][nb_c][nb_d] = {0};

        // fill 3d array
        for (int i = 0; i < nb_l; i++) {
            for (int j = 0; j < nb_c; j++) {
                for (int k = 0; k < nb_d; k++) {
                    array3d[i][j][k] = array3d_in[i][j][k];
                }
            }
        }

        // write
        char *fieldname3d = (char *) name.c_str();
        size_t dim3d[3] = {nb_l, nb_c, nb_d};
        matvar_t *variable3d = Mat_VarCreate(fieldname3d, MAT_C_INT32, MAT_T_INT32, 3, dim3d, &array3d, 0); // rank 3
        Mat_VarWrite(matfp, variable3d, MAT_COMPRESSION_ZLIB);
        Mat_VarFree(variable3d);
    }

    //

    template<size_t nb_l, size_t nb_c, size_t nb_d>
    void output_3D_matrix(const double array3d_in[nb_l][nb_c][nb_d], const std::string &name) {
        double array3d[nb_l][nb_c][nb_d] = {0};

        // fill 3d array
        for (int i = 0; i < nb_l; i++) {
            for (int j = 0; j < nb_c; j++) {
                for (int k = 0; k < nb_d; k++) {
                    array3d[i][j][k] = array3d_in[i][j][k];
                }
            }
        }

        // write
        char *fieldname3d = (char *) name.c_str();
        size_t dim3d[3] = {nb_l, nb_c, nb_d};
        matvar_t *variable3d = Mat_VarCreate(fieldname3d, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dim3d, &array3d, 0); // rank 3
        Mat_VarWrite(matfp, variable3d, MAT_COMPRESSION_ZLIB);
        Mat_VarFree(variable3d);
    }

private:

    mat_t *matfp = NULL; // matfp contains pointer to MAT file or NULL on failure
};
