#include "FIELD_2D.h"
#include <cmath>
#include <iostream>

#if _WIN32
#include <gl/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#elif __linux__
#include <GL/glut.h>
#endif

using namespace std;

// forward declarations
void apply_boundary_conditions();
void runOnce();
void runEverytime();
void recompute_timestep();
void fill_AC_field(FIELD_2D &f, FIELD_2D &pField, int dim, float delta);
float space_derivative(const FIELD_2D &f, int dim, int x, int y, float delta);

// resolution
int xRes = 250;
int yRes = 250;

// fields: phi (phase), d (diffusive mass), b (boundary mass)
FIELD_2D field(xRes, yRes);       // display field
FIELD_2D phiField(xRes, yRes);
FIELD_2D dField(xRes, yRes);
FIELD_2D bMassField(xRes, yRes);
FIELD_2D newPhi(xRes, yRes), newD(xRes, yRes), newB(xRes, yRes);
FIELD_2D AC_dpdy(xRes, yRes), AC_dpdx(xRes, yRes);  // Allen-Cahn helpers

// Screen resolution
int xScreenRes = 650;
int yScreenRes = 650;

string windowLabel("Frost Growth");

int xMouse = -1, yMouse = -1;
int mouseButton = -1, mouseState = -1, mouseModifiers = -1;
int xField = -1, yField = -1;

// UI state
bool animate = true;
bool drawingGrid = false;
bool drawingValues = true;
int displayMode = 0;  // 0=phi, 1=d, 2=b

// camera
float eyeCenter[] = {0.5, 0.5, 1};
float zoom = 1.0;

// model parameters
float rho              = 0.35f;     // supersaturation, turning up starves side branches
float kappa_eff        = 0.003f;    // freezing retention
float mu_eff           = 0.005f;    // melting
float beta_convex      = 0.05f;     // controls tip attachment, turning up speeds tip growth
float beta_concave     = 10.0f;     // concave attachment, kills side branches
float beta_scale       = 0.5f;      
float delta_beta_hex   = 0.35f;     // anisotropy
float hex_penalty      = 1.0f;      // additive angular penalty IMPORTANT
float delta_beta_sharp = 0.005f;    // binary on/off attachment
float tau_f            = 0.005f;    // d->b conversion
float tau_m            = 0.1f;      // melting timescale
float tau_a            = 0.02f;     // attachment timescale
float D_diff           = 2.0f;      

// helpers
float eps_safe      = 1e-6f;     
float seed_radius   = 5.0f;     
int   substeps      = 15;
 
// Allen-Cahn interface smoothing 
float eps0          = 0.01f;   // Interface thickness parameter
float delta_hex_ac  = 0.05f;   // Hex anisotropy strength for AC term
float ac_weight     = 0.0f;    // Disabled for now

// drift/advection parameters
float drift_strength  = 0.0f;    // vapor drift speed (0 = off, try 3-15)
float drift_angle     = M_PI/2;  // direction vapor flows (crystal grows AGAINST this)
float drift_curvature = 0.0f;    // rotate drift with distance from seed (radians/unit)
float drift_vx = 0.0f, drift_vy = 0.0f;  // computed from strength+angle

float feather_mode = 0.0f;  // 0=pure 6-fold snowflake, 1=pure 2-fold feather
float feather_axis = 0.0f;  // auto-set from drift_angle in runOnce()

// seed position
int cx_seed = 0, cy_seed = 0;

float length    = 10.0f;
float time_step = 0.001f;  // recomputed by recompute_timestep()

// step counter
int stepCount = 0;

void recompute_timestep() {
  float h = length / xRes;
  float dt_diff = h * h / (4.0f * D_diff);
  float dt_adv = (drift_strength > 0.0f) ? h / drift_strength : 1e10f;
  time_step = 0.8f * fminf(dt_diff, dt_adv);
}

// helpers

inline float sigmoid_fn(float x) {
  return 0.5f * (1.0f + tanhf(x));
}

inline float p_smooth(float phi) {
  // smooth heaviside: 0 at phi=0, 1 at phi=1
  return phi * phi * phi * (10.0f - 15.0f * phi + 6.0f * phi * phi);
}

inline float delta_interface(float phi) {
  // peaks at phi=0.5, zero at phi=0 and phi=1
  return phi * (1.0f - phi);
}

float space_derivative(const FIELD_2D &f, int dim, int x, int y, float delta) {
  if (!dim)
    return (f(x+1, y) - f(x-1, y)) / (2.0f * delta);
  else
    return (f(x, y+1) - f(x, y-1)) / (2.0f * delta);
}

// hex grid helpers (odd-r offset: odd rows shifted right)
// neighbor order: E(0), NE(1), NW(2), W(3), SW(4), SE(5)
// angles:         0,    60,    120,   180,  240,   300 

inline void hex_get_neighbors(const FIELD_2D &f, int x, int y, float n[6]) {
  if (y % 2 == 0) {
    n[0]=f(x+1,y); n[1]=f(x,  y+1); n[2]=f(x-1,y+1);
    n[3]=f(x-1,y); n[4]=f(x-1,y-1); n[5]=f(x,  y-1);
  } else {
    n[0]=f(x+1,y); n[1]=f(x+1,y+1); n[2]=f(x,  y+1);
    n[3]=f(x-1,y); n[4]=f(x,  y-1); n[5]=f(x+1,y-1);
  }
}

// convert the heat update back to the square stencil, that way connected to well-established version
// might have to backsolve for the timestep
// vary the theta0 term inside of the phase equation, to create a spiral look

inline float hex_laplacian(const float n[6], float c, float h) {
  return (2.0f/(3.0f*h*h)) * (n[0]+n[1]+n[2]+n[3]+n[4]+n[5] - 6.0f*c);
}

inline void hex_gradient(const float n[6], float h, float &gx, float &gy) {
  gx = (1.0f/(3.0f*h)) * (n[0]-n[3] + 0.5f*(n[1]-n[2]+n[5]-n[4]));
  gy = (sqrtf(3.0f)/(6.0f*h)) * (n[1]+n[2]-n[4]-n[5]);
}

inline void hex_second_derivs(const float n[6], float c, float h,
                               float &uxx, float &uyy, float &uxy) {
  float h2 = h*h;
  uxx = (n[0]+n[3]-2.0f*c) / h2;
  uyy = (-n[0]-n[3]+2.0f*(n[1]+n[2]+n[4]+n[5])-6.0f*c) / (3.0f*h2);
  uxy = (n[1]-n[2]+n[4]-n[5]) / (sqrtf(3.0f)*h2);
}

///////////////////////////////////////////////////////////////////////
// GL helper functions
///////////////////////////////////////////////////////////////////////
void refreshMouseFieldIndex(int x, int y)
{
  y = yScreenRes - y;
  float xNorm = (float)x / xScreenRes;
  float yNorm = (float)y / yScreenRes;

  float halfZoom = 0.5 * zoom;
  float xWorldMin = eyeCenter[0] - halfZoom;
  float xWorldMax = eyeCenter[0] + halfZoom;
  float xMin = (0.0 - xWorldMin) / (xWorldMax - xWorldMin);
  float xMax = (1.0 - xWorldMin) / (xWorldMax - xWorldMin);

  float yWorldMin = eyeCenter[1] - halfZoom;
  float yWorldMax = eyeCenter[1] + halfZoom;
  float yMin = (0.0 - yWorldMin) / (yWorldMax - yWorldMin);
  float yMax = (1.0 - yWorldMin) / (yWorldMax - yWorldMin);

  float xScale = 1.0, yScale = 1.0;
  if (xRes < yRes) xScale = (float)yRes / xRes;
  if (xRes > yRes) yScale = (float)xRes / yRes;

  xField = xScale * xRes * ((xNorm - xMin) / (xMax - xMin));
  yField = yScale * yRes * ((yNorm - yMin) / (yMax - yMin));

  xField = (xField < 0) ? 0 : xField;
  xField = (xField >= xRes) ? xRes - 1 : xField;
  yField = (yField < 0) ? 0 : yField;
  yField = (yField >= yRes) ? yRes - 1 : yField;
}

void printGlString(string output)
{
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

void updateTexture(FIELD_2D &texture)
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, texture.xRes(), texture.yRes(), 0, GL_LUMINANCE, GL_FLOAT, texture.data());
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glEnable(GL_TEXTURE_2D);
}

void drawGridLines()
{
  glColor4f(0.1, 0.1, 0.1, 1.0);
  float dx = 1.0 / xRes;
  float dy = 1.0 / yRes;
  if (xRes < yRes) dx *= (float)xRes / yRes;
  if (xRes > yRes) dy *= (float)yRes / xRes;

  glBegin(GL_LINES);
  for (int x = 0; x < field.xRes() + 1; x++) {
    glVertex3f(x * dx, 0, 1);
    glVertex3f(x * dx, 1, 1);
  }
  for (int y = 0; y < field.yRes() + 1; y++) {
    glVertex3f(0, y * dy, 1);
    glVertex3f(1, y * dy, 1);
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// GL callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  float halfZoom = zoom * 0.5;
  glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eyeCenter[0], eyeCenter[1], 1,
            eyeCenter[0], eyeCenter[1], 0,
            0, 1, 0);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // Hex grid: y-spacing is h*sqrt(3)/2, account for non-square grids
  float xLength = 1.0;
  float yLength = (float)yRes / xRes * sqrtf(3.0f) / 2.0f;

  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(0.0, 0.0, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(0.0, yLength, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(xLength, yLength, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(xLength, 0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  if (drawingGrid) drawGridLines();

  if (xField >= 0 && yField >= 0 && xField < field.xRes() && yField < field.yRes()) {
    glLoadIdentity();
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    float hz = 0.5 * zoom;
    glRasterPos3f(-hz * 0.95, -hz * 0.95, 0);

    char buffer[256];
    const char* modeStr[] = {"phi", "d", "b"};
    sprintf(buffer, "(%d, %d) %s=%.4f  phi=%.3f d=%.4f b=%.4f  step=%d",
            xField, yField, modeStr[displayMode], field(xField, yField),
            phiField(xField, yField), dField(xField, yField), bMassField(xField, yField),
            stepCount);
    if (drawingValues) printGlString(string(buffer));
  }

  glutSwapBuffers();
}

void printCommands()
{
  cout << "=============================================" << endl;
  cout << " GG Snowflake Growth Simulator" << endl;
  cout << "=============================================" << endl;
  cout << " a         - toggle animation" << endl;
  cout << " 1/2/3     - display phi / d / b" << endl;
  cout << " r         - reset simulation" << endl;
  cout << " i         - print parameters" << endl;
  cout << " g         - toggle grid overlay" << endl;
  cout << " v         - toggle value display" << endl;
  cout << " w         - write PPM screenshot" << endl;
  cout << " d         - toggle drift on/off" << endl;
  cout << " +/-       - adjust drift strength" << endl;
  cout << " f         - toggle feather mode (6-fold <-> 2-fold)" << endl;
  cout << " q         - quit" << endl;
  cout << " left mouse  - pan" << endl;
  cout << " right mouse - zoom" << endl;
  cout << " shift+left  - erase (set phi=0)" << endl;
}

void glutSpecial(int key, int x, int y) {}

void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key) {
  case 'a': animate = !animate; break;
  case 'g': drawingGrid = !drawingGrid; break;
  case 'v': drawingValues = !drawingValues; break;
  case '1': displayMode = 0; cout << "Display: phi" << endl; break;
  case '2': displayMode = 1; cout << "Display: d (normalized)" << endl; break;
  case '3': displayMode = 2; cout << "Display: b (normalized)" << endl; break;
  case 'r':
    stepCount = 0;
    phiField.clear(); dField.clear(); bMassField.clear();
    newPhi.clear(); newD.clear(); newB.clear();
    runOnce();
    cout << "Simulation reset." << endl;
    break;
  case 'd':
    drift_strength = (drift_strength > 0.0f) ? 0.0f : 5.0f;
    drift_vx = drift_strength * cosf(drift_angle);
    drift_vy = drift_strength * sinf(drift_angle);
    recompute_timestep();
    cout << "Drift strength: " << drift_strength << endl;
    break;
  case '+': case '=':
    drift_strength += 1.0f;
    drift_vx = drift_strength * cosf(drift_angle);
    drift_vy = drift_strength * sinf(drift_angle);
    recompute_timestep();
    cout << "Drift strength: " << drift_strength << endl;
    break;
  case '-':
    drift_strength = fmaxf(0.0f, drift_strength - 1.0f);
    drift_vx = drift_strength * cosf(drift_angle);
    drift_vy = drift_strength * sinf(drift_angle);
    recompute_timestep();
    cout << "Drift strength: " << drift_strength << endl;
    break;
  case 'f': {
    // Cycle: 0 (6-fold) -> 2 (1-fold) -> 1 (2-fold) -> 0
    if (feather_mode < 0.5f) feather_mode = 2.0f;
    else if (feather_mode > 1.5f) feather_mode = 1.0f;
    else feather_mode = 0.0f;
    const char* names[] = {"6-fold needles", "2-fold bilateral", "1-fold spike"};
    cout << "Mode: " << names[(int)feather_mode] << endl;
    if (feather_mode > 0.5f && drift_strength < 1.0f)
      cout << "  Hint: press 'd' to enable drift, then 'r' to reset seed position" << endl;
  } break;
  case 'i':
    cout << "--- Parameters ---" << endl;
    cout << "rho=" << rho << " kappa_eff=" << kappa_eff << " mu_eff=" << mu_eff << endl;
    cout << "beta_convex=" << beta_convex << " beta_concave=" << beta_concave
         << " beta_scale=" << beta_scale << endl;
    cout << "delta_beta_hex=" << delta_beta_hex << " delta_beta_sharp=" << delta_beta_sharp << endl;
    cout << "tau_f=" << tau_f << " tau_m=" << tau_m << " tau_a=" << tau_a << endl;
    cout << "D_diff=" << D_diff << " time_step=" << time_step << endl;
    cout << "drift_strength=" << drift_strength << " drift_angle=" << drift_angle
         << " drift_curvature=" << drift_curvature << endl;
    cout << "feather_mode=" << feather_mode << " feather_axis=" << feather_axis << endl;
    cout << "eps0=" << eps0 << " delta_hex_ac=" << delta_hex_ac
         << " ac_weight=" << ac_weight << endl;
    cout << "Step: " << stepCount << endl;
    break;
  case 'w': {
    static int count = 0;
    char buffer[256];
    sprintf(buffer, "output_%i.ppm", count);
    field.writePPM(buffer);
    cout << "Wrote " << buffer << endl;
    count++;
  } break;
  case '?': printCommands(); break;
  case 'q': exit(0); break;
  default: break;
  }
}

void glutMouseClick(int button, int state, int x, int y)
{
  int modifiers = glutGetModifiers();
  mouseButton = button;
  mouseState = state;
  mouseModifiers = modifiers;

  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && modifiers & GLUT_ACTIVE_SHIFT) {
    refreshMouseFieldIndex(x, y);
    for (int yy = yField - 4; yy < yField + 5; yy++)
      for (int xx = xField - 4; xx < xField + 5; xx++)
        if (xx >= 0 && xx < xRes && yy >= 0 && yy < yRes)
          phiField(xx, yy) = 0;
    return;
  }

  xMouse = x;
  yMouse = y;
}

void glutMouseMotion(int x, int y)
{
  if (mouseButton == GLUT_LEFT_BUTTON && mouseState == GLUT_DOWN && mouseModifiers & GLUT_ACTIVE_SHIFT) {
    refreshMouseFieldIndex(x, y);
    for (int yy = yField - 4; yy < yField + 5; yy++)
      for (int xx = xField - 4; xx < xField + 5; xx++)
        if (xx >= 0 && xx < xRes && yy >= 0 && yy < yRes)
          phiField(xx, yy) = 0;
    return;
  }

  float xDiff = x - xMouse;
  float yDiff = y - yMouse;
  float speed = 0.001;

  if (mouseButton == GLUT_LEFT_BUTTON) {
    eyeCenter[0] -= xDiff * speed;
    eyeCenter[1] += yDiff * speed;
  }
  if (mouseButton == GLUT_RIGHT_BUTTON)
    zoom -= yDiff * speed;

  xMouse = x;
  yMouse = y;
}

void glutPassiveMouseMotion(int x, int y)
{
  refreshMouseFieldIndex(x, y);
}

void glutIdle()
{
  if (animate)
    for (int s = 0; s < substeps; s++)
      runEverytime();
  updateTexture(field);
  glutPostRedisplay();
}

int glvuWindow()
{
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(xScreenRes, yScreenRes);
  glutInitWindowPosition(10, 10);
  glutCreateWindow(windowLabel.c_str());

  glViewport(0, 0, (GLsizei)xScreenRes, (GLsizei)yScreenRes);
  glClearColor(0.1, 0.1, 0.1, 0);

  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);
  glutPassiveMotionFunc(&glutPassiveMouseMotion);

  glutMainLoop();
  return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  if (xRes < yRes) eyeCenter[0] = (float)xRes / yRes * 0.5;
  if (yRes < xRes) eyeCenter[1] = (float)yRes / xRes * 0.5;

  phiField.clear();
  dField.clear();
  bMassField.clear();
  newPhi.clear();
  newD.clear();
  newB.clear();
  AC_dpdx.clear();
  AC_dpdy.clear();

  printCommands();
  cout << "Simulating GG continuous snowflake growth" << endl;

  runOnce();

  glutInit(&argc, argv);
  glvuWindow();
  return 1;
}

///////////////////////////////////////////////////////////////////////
// Boundary Conditions
///////////////////////////////////////////////////////////////////////
void apply_boundary_conditions()
{
  // phi: Neumann (dphi/dn = 0) - copy interior to boundary
  for (int i = 1; i < yRes - 1; i++) {
    phiField(0, i) = phiField(1, i);
    phiField(xRes - 1, i) = phiField(xRes - 2, i);
  }
  for (int i = 0; i < xRes; i++) {
    phiField(i, 0) = phiField(i, 1);
    phiField(i, yRes - 1) = phiField(i, yRes - 2);
  }
  phiField(0, 0)             = phiField(1, 1);
  phiField(xRes-1, 0)        = phiField(xRes-2, 1);
  phiField(0, yRes-1)        = phiField(1, yRes-2);
  phiField(xRes-1, yRes-1)   = phiField(xRes-2, yRes-2);

  // b: Neumann (db/dn = 0)
  for (int i = 1; i < yRes - 1; i++) {
    bMassField(0, i) = bMassField(1, i);
    bMassField(xRes - 1, i) = bMassField(xRes - 2, i);
  }
  for (int i = 0; i < xRes; i++) {
    bMassField(i, 0) = bMassField(i, 1);
    bMassField(i, yRes - 1) = bMassField(i, yRes - 2);
  }
  bMassField(0, 0)           = bMassField(1, 1);
  bMassField(xRes-1, 0)      = bMassField(xRes-2, 1);
  bMassField(0, yRes-1)      = bMassField(1, yRes-2);
  bMassField(xRes-1, yRes-1) = bMassField(xRes-2, yRes-2);

  // d: Dirichlet (d = rho at edges, infinite vapor reservoir)
  for (int i = 0; i < yRes; i++) {
    dField(0, i) = rho;
    dField(xRes - 1, i) = rho;
  }
  for (int i = 0; i < xRes; i++) {
    dField(i, 0) = rho;
    dField(i, yRes - 1) = rho;
  }
}

///////////////////////////////////////////////////////////////////////
// Allen-Cahn anisotropic helper field
// stores eps(theta) * eps'(theta) * dphi/d(dim) at each point
///////////////////////////////////////////////////////////////////////
void fill_AC_field(FIELD_2D &f, FIELD_2D &pField, int dim, float delta)
{
  for (int y = 1; y < yRes - 1; y++) {
    for (int x = 1; x < xRes - 1; x++) {
      float dp = space_derivative(pField, dim, x, y, delta);

      float phi_x = space_derivative(pField, 0, x, y, delta);
      float phi_y = space_derivative(pField, 1, x, y, delta);
      float theta = atan2f(-phi_x, -phi_y);

      float eps_val   = eps0 * (1.0f + delta_hex_ac * cosf(6.0f * theta));
      float de_dtheta = -6.0f * eps0 * delta_hex_ac * sinf(6.0f * theta);

      f(x, y) = eps_val * de_dtheta * dp;
    }
  }

  // Neumann BC on helper field
  for (int i = 1; i < yRes - 1; i++) {
    f(0, i) = f(1, i);
    f(xRes - 1, i) = f(xRes - 2, i);
  }
  for (int i = 0; i < xRes; i++) {
    f(i, 0) = f(i, 1);
    f(i, yRes - 1) = f(i, yRes - 2);
  }
}

///////////////////////////////////////////////////////////////////////
// initial conditions: hexagonal seed + uniform d = rho
///////////////////////////////////////////////////////////////////////
void runOnce()
{
  float h = length / xRes;

  // compute drift velocity components
  drift_vx = drift_strength * cosf(drift_angle);
  drift_vy = drift_strength * sinf(drift_angle);

  // feather spine grows opposite to drift direction
  feather_axis = drift_angle + (float)M_PI;

  // recompute timestep for current drift strength
  recompute_timestep();

  // seed position: if drift is on, place seed on upwind edge;
  // otherwise center as usual
  int cx, cy;
  if (drift_strength > 0.0f) {
    float margin = 0.15f;
    cx = xRes/2 + (int)(xRes * (0.5f - margin) * cosf(drift_angle));
    cy = yRes/2 + (int)(yRes * (0.5f - margin) * sinf(drift_angle));
    int sr = (int)seed_radius + 2;
    if (cx < sr) cx = sr;
    if (cx >= xRes - sr) cx = xRes - 1 - sr;
    if (cy < sr) cy = sr;
    if (cy >= yRes - sr) cy = yRes - 1 - sr;
  } else {
    cx = xRes / 2;
    cy = yRes / 2;
  }
  cx_seed = cx;
  cy_seed = cy;

  // physical center of the seed
  float pcx = (cx + 0.5f * (cy % 2)) * h;
  float pcy = cy * h * sqrtf(3.0f) / 2.0f;

  for (int y = 0; y < yRes; y++) {
    for (int x = 0; x < xRes; x++) {
      // physical position in hex grid
      float px = (x + 0.5f * (y % 2)) * h;
      float py = y * h * sqrtf(3.0f) / 2.0f;
      float dx = px - pcx;
      float dy = py - pcy;

      // signed distance to hex seed boundary (physical coords)
      float s3 = sqrtf(3.0f);
      float sr = seed_radius * h;  // seed radius in physical units
      float inradius = sr * s3 / 2.0f;
      float d1 = inradius - fabsf(dy);
      float d2 = inradius - fabsf(dx * s3 / 2.0f + dy / 2.0f);
      float d3 = inradius - fabsf(dx * s3 / 2.0f - dy / 2.0f);
      float hex_dist = fminf(d1, fminf(d2, d3));

      float phi_init = 0.5f * (1.0f + tanhf(hex_dist / h));

      phiField(x, y)  = phi_init;
      dField(x, y)    = rho * (1.0f - phi_init);
      bMassField(x, y) = 0.0f;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// main simulation update
///////////////////////////////////////////////////////////////////////
void runEverytime()
{
  apply_boundary_conditions();

  float h = length / xRes;  // hex cell spacing

  for (int y = 1; y < yRes - 1; y++) {
    for (int x = 1; x < xRes - 1; x++) {
      float phi   = phiField(x, y);
      float d_val = dField(x, y);
      float b_val = bMassField(x, y);

      float p = p_smooth(phi);

      // hex-grid phi neighbors, gradient, curvature
      float pn[6];
      hex_get_neighbors(phiField, x, y, pn);

      float phi_x, phi_y;
      hex_gradient(pn, h, phi_x, phi_y);
      float grad_sq = phi_x*phi_x + phi_y*phi_y;
      float grad_mag = sqrtf(grad_sq);

      // boundary indicator (gradient-based)
      float boundary = fminf(grad_mag * h * 2.0f, 1.0f);

      float theta = atan2f(-phi_x, -phi_y);

      // second derivatives on hex grid
      float phi_xx, phi_yy, phi_xy;
      hex_second_derivs(pn, phi, h, phi_xx, phi_yy, phi_xy);

      float denom = powf(grad_sq + eps_safe, 1.5f);
      float kappa_curv = -(phi_xx*phi_y*phi_y
                         - 2.0f*phi_x*phi_y*phi_xy
                         + phi_yy*phi_x*phi_x) / denom;

      // attachment 
      float beta_base = beta_concave
                      + (beta_convex - beta_concave) * sigmoid_fn(kappa_curv / beta_scale);
      // additive angular penalty: 0 at preferred directions, hex_penalty at off-angles
      float off_axis;
      if (feather_mode < 0.5f) {
        // mode 0: 6-fold — 6 needle directions
        off_axis = 0.5f * (1.0f + cosf(6.0f * theta));
      } else if (feather_mode < 1.5f) {
        // mode 1: 2-fold — bilateral needles along feather spine
        off_axis = 0.5f * (1.0f + cosf(2.0f * (theta - feather_axis)));
      } else {
        // mode 2: 1-fold — single needle in growth direction
        off_axis = 0.5f * (1.0f + cosf(theta - feather_axis));
      }
      float beta_eff  = beta_base + hex_penalty * off_axis;
      float attach_sig = sigmoid_fn((b_val - beta_eff) / delta_beta_sharp);
      float dphi_attach = (1.0f / tau_a) * (1.0f - phi) * attach_sig;

      float dphi_dt = dphi_attach;

      // diffusion on hex grid (vapor only)
      float dn[6];
      hex_get_neighbors(dField, x, y, dn);
      float laplace_d = hex_laplacian(dn, d_val, h);
      float dd_diff = D_diff * laplace_d * (1.0f - p);

      // advection/drift of vapor (paper eq 1c analog)
      float dd_advect = 0.0f;
      if (drift_strength > 0.0f) {
        float local_vx = drift_vx;
        float local_vy = drift_vy;

        // rotate drift direction with distance for curvy feathers
        if (drift_curvature != 0.0f) {
          float px = (x + 0.5f * (y % 2)) * h;
          float py = y * h * sqrtf(3.0f) / 2.0f;
          float scx = (cx_seed + 0.5f * (cy_seed % 2)) * h;
          float scy = cy_seed * h * sqrtf(3.0f) / 2.0f;
          float dist = sqrtf((px-scx)*(px-scx) + (py-scy)*(py-scy) + eps_safe);
          float angle_off = drift_curvature * dist;
          float ca = cosf(angle_off), sa = sinf(angle_off);
          float rvx = local_vx * ca - local_vy * sa;
          float rvy = local_vx * sa + local_vy * ca;
          local_vx = rvx;
          local_vy = rvy;
        }

        // gradient of d on hex grid (reuse dn neighbors)
        float grad_d_x, grad_d_y;
        hex_gradient(dn, h, grad_d_x, grad_d_y);

        // advection: -v . grad(d), vapor only
        dd_advect = -(local_vx * grad_d_x + local_vy * grad_d_y) * (1.0f - p);
      }

      // freezing (d -> b near interface, vapor side)
      float freeze_rate = ((1.0f - kappa_eff) / tau_f) * d_val * boundary * (1.0f - p);

      // melting (b -> d near interface) 
      float melt_rate = (mu_eff / tau_m) * b_val * boundary;

      // mass consumptioncapped
      float mass_factor = dphi_dt / fmaxf(1.0f - phi, 0.05f);

      float dd_dt = dd_diff + dd_advect - freeze_rate + melt_rate - d_val * mass_factor;
      float db_dt = freeze_rate - melt_rate - b_val * mass_factor;

      // forward euler
      newPhi(x, y) = phi   + time_step * dphi_dt;
      newD(x, y)   = d_val + time_step * dd_dt;
      newB(x, y)   = b_val + time_step * db_dt;

      if (newPhi(x, y) < 0.0f) newPhi(x, y) = 0.0f;
      if (newPhi(x, y) > 1.0f) newPhi(x, y) = 1.0f;
      if (newD(x, y) < 0.0f) newD(x, y) = 0.0f;
      if (newB(x, y) < 0.0f) newB(x, y) = 0.0f;
    }
  }

  // swap fields
  phiField  = newPhi;
  dField    = newD;
  bMassField = newB;

  // update display field
  switch (displayMode) {
    case 0:
      field = phiField;
      break;
    case 1:
      field = dField;
      if (field.max() > field.min()) field.normalize();
      break;
    case 2:
      field = bMassField;
      if (field.max() > field.min()) field.normalize();
      break;
  }

  stepCount++;
  if (stepCount % 500 == 0) {
    float maxB = 0, maxPhi = 0, maxD = 0;
    for (int i = 0; i < xRes * yRes; i++) {
      if (bMassField[i] > maxB) maxB = bMassField[i];
      if (phiField[i] > maxPhi) maxPhi = phiField[i];
      if (dField[i] > maxD) maxD = dField[i];
    }
    cout << "Step " << stepCount
         << "  maxPhi=" << maxPhi << "  maxB=" << maxB << "  maxD=" << maxD << endl;
  }
}
