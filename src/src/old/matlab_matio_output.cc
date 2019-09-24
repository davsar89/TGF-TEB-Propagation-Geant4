#include "matlab_matio_output.hh"

MATLAB_MATIO_output::MATLAB_MATIO_output(const std::string &name) {
    const char *filename = (char *) name.c_str();
    matfp = Mat_CreateVer(filename, NULL, MAT_FT_MAT5); // or MAT_FT_MAT4 / MAT_FT_MAT73
}

MATLAB_MATIO_output::~MATLAB_MATIO_output() {
    Mat_Close(matfp);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_vector(const std::vector<int> &output_vec, const std::string &name) {
    //
    const size_t first = output_vec.size(); // rows
    int *array1d = new int[first];

    for (uint i = 0; i < first; i++) array1d[i] = output_vec[i];

    // write
    const char *fieldname1d = (char *) name.c_str();
    size_t dim1d[1] = {first};
    matvar_t *variable1d = Mat_VarCreate(fieldname1d, MAT_C_INT32, MAT_T_INT32, 1, dim1d, array1d, 0); // rank 1
    Mat_VarWrite(matfp, variable1d, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(variable1d);
    delete[] array1d;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_vector(const std::vector<double> &output_vec, const std::string &name) {
    //
    const size_t first = output_vec.size(); // rows
    double *array1d = new double[first];

    for (uint i = 0; i < first; i++) array1d[i] = output_vec[i];

    // write
    const char *fieldname1d = (char *) name.c_str();
    size_t dim1d[1] = {first};
    matvar_t *variable1d = Mat_VarCreate(fieldname1d, MAT_C_DOUBLE, MAT_T_DOUBLE, 1, dim1d, array1d, 0); // rank 1
    Mat_VarWrite(matfp, variable1d, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(variable1d);
    delete[] array1d;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_scalar(double scalar, const std::string &name) {
    // write
    const char *fieldname = (char *) name.c_str();
    size_t dim[2] = {1, 1};
    matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim, &scalar, 0);
    Mat_VarWrite(matfp, variable, MAT_COMPRESSION_ZLIB); // or MAT_COMPRESSION_NONE
    Mat_VarFree(variable);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_scalar(long scalar, const std::string &name) {
    // write
    const char *fieldname = (char *) name.c_str();
    size_t dim[2] = {1, 1};
    matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_INT64, MAT_T_INT64, 2, dim, &scalar, 0);
    Mat_VarWrite(matfp, variable, MAT_COMPRESSION_ZLIB); // or MAT_COMPRESSION_NONE
    Mat_VarFree(variable);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_scalar(int scalar, const std::string &name) {
    // write
    const char *fieldname = (char *) name.c_str();
    size_t dim[2] = {1, 1};
    matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_INT32, MAT_T_INT32, 2, dim, &scalar, 0);
    Mat_VarWrite(matfp, variable, MAT_COMPRESSION_ZLIB); // or MAT_COMPRESSION_NONE
    Mat_VarFree(variable);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MATLAB_MATIO_output::output_string(const std::string &text_out, const std::string &field_name_str) {
    // write
    char *fieldname = (char *) field_name_str.c_str();
    char *mystring = (char *) text_out.c_str();
    size_t dim[2] = {1, sizeof(mystring) / sizeof(mystring[0])};
    matvar_t *variable = Mat_VarCreate(fieldname, MAT_C_CHAR, MAT_T_UTF8, 2, dim, mystring, 0);
    Mat_VarWrite(matfp, variable, MAT_COMPRESSION_ZLIB); // or MAT_COMPRESSION_NONE
    Mat_VarFree(variable);
}
