// Local Headers
#include "glitter.hpp"

// System Headers
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// Standard Headers
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <math.h>
#include <algorithm>

#include "SPHSimulator.hpp"

//application
GLFWwindow* mWindow;

//mouse input
static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
static void onRelease();
static void onPress();
static void duringPress();
static bool is_mouse_pressed = false;
static double cursor_pre_x;
static double cursor_pre_y;
static double cursor_cur_x;
static double cursor_cur_y;
static void arcball_update();
static const double CURSOR_POS_INVALID = -2;
static const double speed_factor = 2.;

//transformation matrices
static mat4 m2w;
static mat4 w2v;
static mat4 pers_proj;

//SPH Simulation
static SPHSim::BoxDef box;
static float delta = 1./60;
static SPHSim::SPHSimulator* simulator;

//forward declarations
static void init();
static void render();
static void update();
static void onExit();

static void renderBox();
static void renderParticles();
static void setTransform();

int main(int argc, char * argv[]) {

    // Load GLFW and Create a Window
    glfwInit();
    mWindow = glfwCreateWindow(mWidth, mHeight, "Parallel Fluid Simulator", nullptr, nullptr);
 
    // Check for Valid Context
    if (mWindow == nullptr) {
        fprintf(stderr, "Failed to Create OpenGL Context");
        glfwTerminate();
        return EXIT_FAILURE;
    }

    //set mouse callback
    glfwSetMouseButtonCallback(mWindow, mouse_button_callback);

    // Create Context and Load OpenGL Functions
    glfwMakeContextCurrent(mWindow);
    gladLoadGL();
    fprintf(stderr, "OpenGL %s\n", glGetString(GL_VERSION));

    // Enable Depth Test
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_GREATER);

    // Initialization
    init();

    // Rendering Loop
    while (glfwWindowShouldClose(mWindow) == false) {
        
        //exit when press esc
        if (glfwGetKey(mWindow, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(mWindow, true);

        // Background Fill Color
        glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
        glClearDepth(-1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        update();
        render();

        // Flip Buffers and Draw
        glfwSwapBuffers(mWindow);
        glfwPollEvents();
    }   glfwTerminate();

    onExit();

    return EXIT_SUCCESS;
}

static void init()
{
    //set up model to view
    //m2w = rotate(glm::radians(30.f), glm::vec3(1.f,1.f,0.f));
    //m2w = rotate(glm::radians(10.f), glm::vec3(0.f,1.f,0.f)) * m2w;

    //set up camera
    vec3 eye(0,0,50);
    vec3 center(0,0,0);
    vec3 up(0,1,0);
    w2v = lookAt(eye, center, up);

    //set up perspective projection
    pers_proj = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);    

    //setup cursor position
    cursor_cur_x = CURSOR_POS_INVALID;
    cursor_cur_y = CURSOR_POS_INVALID;
    cursor_pre_x = CURSOR_POS_INVALID;
    cursor_pre_y = CURSOR_POS_INVALID;

    //init SPH
    delta = 1./60;

    box.spanX = 25;
    box.spanY = 25;
    box.spanZ = 25;

    SPHSim::SPHConfig conf;
    conf.kernel_h = KERNEL_H;
    conf.k = K;
    conf.density0 = RESET_DENSITY;
    conf.g = G;
    conf.miu = MIU;
    conf.boundarydamping = BOUNDARYDAMPING;
    conf.damping = DAMPING;
    conf.tensionCoe = TENSION_COEF;
    conf.box = box;
    simulator = new SPHSim::SPHSimulator(conf);
    simulator->setup();

    return;
}

//
// RENDERING
//
static void render()
{
    renderBox();
    renderParticles();
}

static void renderBox()
{
    setTransform();

    double boxX = box.spanX;
    double boxY = box.spanY;
    double boxZ = box.spanZ;

  glBegin(GL_LINES);
      glColor3f(1,0,0);
      glVertex3f( -boxX/2, -boxY/2,  boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);

      glVertex3f(  boxX/2,  boxY/2,   boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,   boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2, -boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,   boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,   boxZ/2);
      glVertex3f( -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f( -boxX/2, -boxY/2,  -boxZ/2);

      glVertex3f(  boxX/2,  boxY/2,   boxZ/2);
      glVertex3f(  boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  boxX/2,  -boxY/2,  boxZ/2);
      glVertex3f(  boxX/2,  -boxY/2,  -boxZ/2);
      glVertex3f(  -boxX/2,  boxY/2,  boxZ/2);
      glVertex3f(  -boxX/2,  boxY/2,  -boxZ/2);
      glVertex3f(  -boxX/2,  -boxY/2, boxZ/2);
      glVertex3f(  -boxX/2,  -boxY/2,  -boxZ/2);

  glEnd();
}

static void renderParticles()
{
    float* data = NULL;
    int num;
    simulator->getData(&data, num);

    glPointSize(2);
    glBegin(GL_POINTS);
    glColor3f(0,1,0);
    for(int i=0; i<num; i++)
    {
        glVertex3f(data[3*i + 0],
                   data[3*i + 1],
                   data[3*i + 2]);
    }
    glEnd();
    free(data);

}

static void setTransform()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMultMatrixf(value_ptr(w2v));
    glMultMatrixf(value_ptr(m2w));

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf(value_ptr(pers_proj));

    glMatrixMode(GL_MODELVIEW);
}

//
// GAME LOGIC UPDATES
//
static void update()
{
    if(is_mouse_pressed)
        duringPress();

    vec4 g(0.f,-1.f,0.f,0.f);
    g = inverse(m2w) * g;
    simulator->setGravityDirection(g.x, g.y, g.z);
    simulator->update(delta);
}

//
// MOUSE INPUT HANDLING
//

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if(action == GLFW_PRESS)
            onPress();
        else if(action == GLFW_RELEASE)
            onRelease();
    }
}

static void onPress()
{
    is_mouse_pressed = true;
}

static void onRelease()
{
    is_mouse_pressed = false;
    cursor_cur_x = CURSOR_POS_INVALID;
    cursor_cur_y = CURSOR_POS_INVALID;
    cursor_pre_x = CURSOR_POS_INVALID;
    cursor_pre_y = CURSOR_POS_INVALID;
}

static void duringPress()
{
    //update cursor position
    cursor_pre_x = cursor_cur_x; cursor_pre_y = cursor_cur_y;
    glfwGetCursorPos(mWindow, &cursor_cur_x, &cursor_cur_y);
    if(cursor_pre_x == CURSOR_POS_INVALID)
    {
        cursor_pre_x = cursor_cur_x;
        cursor_pre_y = cursor_cur_y;
    }

    //arcball update
    arcball_update();
}


//
// ARCBALL MECHANISM
//
// ref: http://courses.cms.caltech.edu/cs171/assignments/hw3/hw3-notes/notes-hw3.html
//
static void arcball_update()
{
    //not need to update if mouse is not moved
    if(cursor_pre_x == cursor_cur_x && cursor_pre_y == cursor_cur_y)
        return;

    //compute the position in NDC space
    double ndc_pre_x = (cursor_pre_x - mWidth/2) / (mWidth/2);
    double ndc_pre_y = (mHeight/2 -cursor_pre_y) / (mHeight/2);
    double ndc_cur_x = (cursor_cur_x - mWidth/2) / (mWidth/2);
    double ndc_cur_y = (mHeight/2 -cursor_cur_y) / (mHeight/2);

    //estimate the z coordinates in NDC space using the sphere surface constraints
    double tmp, ndc_pre_z, ndc_cur_z;

    tmp = ndc_pre_x * ndc_pre_x + ndc_pre_y * ndc_pre_y;
    if(tmp > 1)
        ndc_pre_z = 0;
    else
        ndc_pre_z = std::sqrt(1 - tmp);

    tmp = ndc_cur_x * ndc_cur_x + ndc_cur_y * ndc_cur_y;
    if(tmp > 1)
        ndc_cur_z = 0;
    else
        ndc_cur_z = std::sqrt(1 - tmp);

    //compute vector u and rotate angle
    vec3 pre(ndc_pre_x, ndc_pre_y, ndc_pre_z);
    vec3 cur(ndc_cur_x, ndc_cur_y, ndc_cur_z);
    vec3 u = normalize(cross(pre, cur));
    float r = std::acos(std::min(dot(pre,cur),1.f));

    //multiplies speed factor
    r *= speed_factor;

    //update m2w matrix
    m2w = rotate(r, u) * m2w;

    return;
}

static void onExit()
{
    delete simulator;
}