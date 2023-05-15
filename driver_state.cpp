#include "driver_state.h"
#include <cstring>
#include <limits>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    state.image_color = new pixel[height * width];
    state.image_depth = new float[height * width];

    for(int i = 0; i < height * width; i++){
        state.image_color[i] = make_pixel(0, 0, 0);
        state.image_depth[i] = std::numeric_limits<float>::max();

    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;

    int triangles;
    int r;
    const data_geometry* output[3];
    data_vertex input[3];
    data_geometry tempArr[3];

    switch(type)
    {
        case render_type::triangle:{
            triangles = state.num_vertices / 3;
            r = 0;

            for(int i = 0; i < triangles; i++) {
                for(int j = 0; j < 3; ++j, r += state.floats_per_vertex){
                    input[j].data = &state.vertex_data[r];
                    tempArr[j].data = input[j].data;
                    state.vertex_shader(input[j], tempArr[j], state.uniform_data);
                    output[j] = &tempArr[j];
                }
                
                clip_triangle(state, *output[0], *output[1], *output[2], 0);
            }

            break;
        }
        case render_type::indexed:{
            triangles = state.num_triangles * 3;
			r = 0;
            for(int i = 0; i < triangles; i += 3){
                for(int j = 0; j < 3; j++){
                    input[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
                    tempArr[j].data = input[j].data;
                    state.vertex_shader(input[j], tempArr[j], state.uniform_data);
                    output[j] = &tempArr[j];
                }

                clip_triangle(state, *output[0], *output[1], *output[2], 0);
            }

            break;
        }
        case render_type::fan:{
            for(int i = 0; i < state.num_vertices; i++){
                for(int j = 0; j < 3; j++){
                    if(j == 0){
                        input[j].data = state.vertex_data + (j * state.floats_per_vertex);
                    }
                    else{
                        input[j].data = state.vertex_data + ((state.floats_per_vertex) * (i + j));
                    }

                    tempArr[j].data = input[j].data;
                    state.vertex_shader(input[j], tempArr[j], state.uniform_data);
                    output[j] = &tempArr[j];
                }

                clip_triangle(state, *output[0], *output[1], *output[2], 0);
            }
            break;
        }
        case render_type::strip:{
            for(int i = 0; i < state.num_vertices - 2; i++){
                for(int j = 0; j < 3; j++){
                    input[j].data = &state.vertex_data[(i+j) * state.floats_per_vertex];
                    tempArr[j].data = input[j].data;
                    state.vertex_shader(input[j], tempArr[j], state.uniform_data);
                    output[j] = &tempArr[j];
                }

                clip_triangle(state, *output[0], *output[1], *output[2], 0);
            }

            break;
        }
        default:
            break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
   /* if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);*/

	const data_geometry* in[3];
    in[0] = &v0;
    in[1] = &v1;
    in[2] = &v2;

	if (face == 1) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    else {
        const data_geometry* temp[3] = {in[0], in[1], in[2]};
        data_geometry arr1[3];
        data_geometry arr2[3];
        vec4 P0, P1;
        vec4 a = (*in)[0].gl_Position;
        vec4 b = (*in)[1].gl_Position;
        vec4 c = (*in)[2].gl_Position;
        float tempA, tempB1, tempB2;

        // It just returns if all the vertices are inside the screen. Does nothing.
        if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3])
            return;
        else {
            if (a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3]) {
                tempB1 = (-b[3] - b[2]) / (a[2] + a[3] - b[3] - b[2]);
                tempB2 = (-a[3] - a[2]) / (c[2] + c[3] - a[3] - a[2]);

                P0 = tempB1 * a + (1 - tempB1) * b;
                P1 = tempB2 * c + (1 - tempB2) * a;

                arr1[0].data = new float[state.floats_per_vertex];
                arr1[1] = *in[1];
                arr1[2] = *in[2];

                for(int i = 0; i < state.floats_per_vertex; i++) {
                    switch (state.interp_rules[i]) {
                        case interp_type::flat:
                            arr1[0].data[i] = (*in)[0].data[i];
                            break;
                        case interp_type::smooth:
                            arr1[0].data[i] = tempB2 * (*in)[2].data[i] + (1 - tempB2) * (*in)[0].data[i];
                            break;
                        case interp_type::noperspective:
                            tempA = tempB2 * (*in)[2].gl_Position[3] / (tempB2 * (*in)[2].gl_Position[3] + (1 - tempB2) * (*in)[0].gl_Position[3]);
                            arr1[0].data[i] = tempA * (*in)[2].data[i] + (1 - tempA) * (*in)[0].data[i];
                            break;
                        default:
                            break;
                    }
                }

                arr1[0].gl_Position = P1;
                temp[0] = &arr1[0];
                temp[1] = &arr1[1];
                temp[2] = &arr1[2];

                clip_triangle(state, *temp[0], *temp[1], *temp[2], face + 1);


                arr2[0].data = new float[state.floats_per_vertex];
                arr2[1] = *in[1];
                arr2[2] = arr1[0];

                for (int i = 0; i < state.floats_per_vertex; i++) {
                    switch (state.interp_rules[i]) {
                        case interp_type::flat:
                            arr2[0].data[i] = (*in)[0].data[i];
                            break;
                        case interp_type::smooth:
                            arr2[0].data[i] = tempB1 * (*in)[0].data[i] + (1 - tempB1) * (*in)[1].data[i];
                            break;
                        case interp_type::noperspective:
                            tempA = tempB1 * (*in)[0].gl_Position[3] / (tempB1 * (*in)[0].gl_Position[3] + (1 - tempB1) * (*in)[1].gl_Position[3]);
                            arr2[0].data[i] = tempA * (*in)[0].data[i] + (1 - tempA) * (*in)[1].data[i];
                            break;
                        default:
                            break;
                    }
                }

                arr2[0].gl_Position = P0;
                temp[0] = &arr2[0];
                temp[1] = &arr2[1];
                temp[2] = &arr2[0];
            }
            clip_triangle(state, *temp[0], *temp[1], *temp[2], face + 1);
        }
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    const data_geometry* in[3];
    in[0] = &v0;
    in[1] = &v1;
    in[2] = &v2;

	float x[3];
    float y[3];
    //Implementing the z-buffering 
    float z[3];

    for(int l = 0; l < 3; l++) {
        float i = ((state.image_width / 2.0f) * ((*in)[l].gl_Position[0]/(*in)[l].gl_Position[3]) + ((state.image_width / 2.0f) - 0.5f));
        float j = ((state.image_height / 2.0f) * ((*in)[l].gl_Position[1]/(*in)[l].gl_Position[3]) + ((state.image_height / 2.0f) - 0.5f));
        float k = ((state.image_width / 2.0f) * ((*in)[l].gl_Position[2]/(*in)[l].gl_Position[3]) + ((state.image_width / 2.0f) - 0.5f));
        x[l] = i;
        y[l] = j;
        z[l] = k;
    }
     
    //Creating the boundaries
    //Minimum and maximum of triangles are calculated
    float xMin = std::min(std::min(x[0], x[1]), x[2]);
    float xMax = std::max(std::max(x[0], x[1]), x[2]);
    float yMin = std::min(std::min(y[0], y[1]), y[2]);
    float yMax = std::max(std::max(y[0], y[1]), y[2]);

    //Ensuring that the measurements are valid
    if(xMin < 0)
        xMin = 0;
    if(xMax > state.image_width)
        xMax = state.image_width;
    if(yMin < 0)
        yMin = 0;
    if(yMax > state.image_height)
        yMax = state.image_height;

    float triangleArea = (0.5f * ((x[1]*y[2] - x[2]*y[1]) - (x[0]*y[2] - x[2]*y[0]) + (x[0]*y[1] - x[1]*y[0])));

    auto *D = new float[state.floats_per_vertex];
    data_fragment fragData{D};
    data_output outData;

    //Calculate the barycentric weight with respect to vertices
    //If the pixel is in the triangle, then color it. Do it for every pixel in the bounding box of the triangle,
    for(int i = yMin; i < yMax; i++) {
        for(int j = xMin; j < xMax; j++) {
            float aP = (0.5f * ((x[1] * y[2] - x[2] * y[1]) + (y[1] - y[2])*j + (x[2] - x[1])*i)) / triangleArea;
            float bP =  (0.5f * ((x[2] * y[0] - x[0] * y[2]) + (y[2] - y[0])*j + (x[0] - x[2])*i)) / triangleArea;
            float gP = (0.5f * ((x[0] * y[1] - x[1] * y[0]) + (y[0] - y[1])*j + (x[1] - x[0])*i)) / triangleArea;

            //Ensure that alpha, beta, and gamma are not negative
            if (aP >= 0 && bP >= 0 && gP >= 0) {
                float alpha = aP;
                float beta = bP;
                float gamma = gP;
                float zB = alpha * z[0] + beta * z[1] + gamma * z[2];
    
                if(zB < state.image_depth[j + i * state.image_width]) {
                    // z-buffer value updated here
                    state.image_depth[j + i * state.image_width] = zB;
                    for (int f = 0; f < state.floats_per_vertex; f++) {
                        float k;
                        switch (state.interp_rules[f]) {
                            //Set the data = to the data of the first vertex
                            case interp_type::flat:
                                fragData.data[f] = (*in)[0].data[f];
                                break;
                            case interp_type::smooth:
                                k = (alpha / (*in)[0].gl_Position[3])
                                         + (beta / (*in)[1].gl_Position[3])
                                         + (gamma / (*in)[2].gl_Position[3]);

                                aP = alpha / (k * (*in)[0].gl_Position[3]);
                                bP = beta / (k * (*in)[1].gl_Position[3]);
                                gP = gamma / (k * (*in)[2].gl_Position[3]);
                            //interpolation with no persepctive, which means no conversion
                            case interp_type::noperspective:
                                fragData.data[f] = aP * (*in)[0].data[f]
                                                    + bP * (*in)[1].data[f]
                                                    + gP * (*in)[2].data[f];
                                break;
                            default:
                                break;
                        }
                    }

                    state.fragment_shader(fragData, outData, state.uniform_data);
                    // multiplying by 255 for every index of RGB
                    state.image_color[j + i * state.image_width] = make_pixel(static_cast<int>(outData.output_color[0] * 255),static_cast<int>(outData.output_color[1] * 255),static_cast<int>(outData.output_color[2] * 255));
                }
            }
        }
    }

    delete [] D;
    
}

