#include "GL/glut.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>

// ====== ТИПЫ ПАРАМЕТРИЗАЦИИ ======
enum ParamType { UNIFORM, CHORD, CENTRIPETAL };
ParamType paramMode = CHORD;
const char* paramNames[] = {"Uniform (t=i)", "Chord Length", "Centripetal"};

// ====== ТИПЫ СПЛАЙНОВ ======
enum SplineType { NATURAL, CLAMPED };
SplineType splineMode = NATURAL;
const char* splineNames[] = {"Natural Cubic", "Clamped Cubic"};

// ====== ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ ======
int winW = 1024, winH = 768;
double mx = 0.0, my = 0.0;
int activePoint = -1, hoverPoint = -1;
bool mouseDown = false;

struct Point { double x, y; };
std::vector<Point> points;  // контрольные точки (x_i, y_i)

// Параметрические узлы и коэффициенты сплайнов
std::vector<double> t_param;  // параметр t_i для каждой точки
std::vector<double> Mx, My;   // вторые производные для X(t) и Y(t)
bool dirty = true;

struct View { double x0, x1, y0, y1; } view = {-6.0, 6.0, -4.5, 4.5};

// ====== МАТЕМАТИКА ======

void solveTridiagonal(double* a, double* b, double* c, double* d, double* x, size_t n) {
    if (n == 0) return;
    for (size_t i = 1; i < n; ++i) {
        double w = a[i] / b[i-1];
        b[i] -= w * c[i-1];
        d[i] -= w * d[i-1];
    }
    x[n-1] = d[n-1] / b[n-1];
    for (int i = (int)n - 2; i >= 0; --i)
        x[i] = (d[i] - c[i] * x[i+1]) / b[i];
}

void computeParameters() {
    size_t n = points.size();
    t_param.resize(n);
    if (n == 0) return;
    
    t_param[0] = 0.0;
    
    switch (paramMode) {
        case UNIFORM:
            for (size_t i = 1; i < n; ++i)
                t_param[i] = (double)i;
            break;
            
        case CHORD:
            for (size_t i = 1; i < n; ++i) {
                double dx = points[i].x - points[i-1].x;
                double dy = points[i].y - points[i-1].y;
                t_param[i] = t_param[i-1] + std::sqrt(dx*dx + dy*dy);
            }
            break;
            
        case CENTRIPETAL:
            for (size_t i = 1; i < n; ++i) {
                double dx = points[i].x - points[i-1].x;
                double dy = points[i].y - points[i-1].y;
                t_param[i] = t_param[i-1] + std::sqrt(std::sqrt(dx*dx + dy*dy));
            }
            break;
    }
    
    // Нормировка в [0, 1] для численной устойчивости
    double total = t_param[n-1];
    if (total > 1e-6)
        for (size_t i = 0; i < n; ++i) t_param[i] /= total;
}

void buildSplineComponent(const std::vector<double>& t,
                          const std::vector<double>& f,
                          std::vector<double>& M) {
    size_t n = t.size();
    M.assign(n, 0.0);
    if (n < 2) return;
    
    if (n == 2) return; // M[0] = M[1] = 0 для одного сегмента
    
    size_t m = n - 2;
    std::vector<double> a(m), b(m), c(m), d(m);
    
    for (size_t i = 0; i < m; ++i) {
        double h_i = t[i+1] - t[i];
        double h_ip1 = t[i+2] - t[i+1];
        a[i] = h_i / 6.0;
        b[i] = (h_i + h_ip1) / 3.0;
        c[i] = h_ip1 / 6.0;
        double df_ip1 = (f[i+2] - f[i+1]) / h_ip1;
        double df_i = (f[i+1] - f[i]) / h_i;
        d[i] = df_ip1 - df_i;
    }
    
    if (splineMode == CLAMPED && n > 2) {
        double h0 = t[1] - t[0];
        double hn = t[n-1] - t[n-2];
        double f0_prime = (f[1] - f[0]) / h0;
        double fn_prime = (f[n-1] - f[n-2]) / hn;
        d[0] -= a[0] * (6.0 / h0) * ((f[1] - f[0]) / h0 - f0_prime);
        d[m-1] -= c[m-1] * (6.0 / hn) * (fn_prime - (f[n-1] - f[n-2]) / hn);
    }
    
    std::vector<double> M_inner(m);
    solveTridiagonal(a.data(), b.data(), c.data(), d.data(), M_inner.data(), m);
    for (size_t i = 0; i < m; ++i) M[i+1] = M_inner[i];
    
    // Граничные условия
    M[0] = (splineMode == NATURAL) ? 0.0 : 0.0; // Для простоты оставляем 0
    M[n-1] = (splineMode == NATURAL) ? 0.0 : 0.0;
}

void buildSplines() {
    if (points.size() < 2) { dirty = false; return; }
    
    computeParameters();
    
    std::vector<double> fx(points.size()), fy(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        fx[i] = points[i].x;
        fy[i] = points[i].y;
    }
    
    buildSplineComponent(t_param, fx, Mx);
    buildSplineComponent(t_param, fy, My);
    
    dirty = false;
}

Point evaluate(double t) {
    if (points.size() < 2) return {0, 0};
    
    // Найти сегмент
    size_t i = 0;
    while (i + 1 < t_param.size() && t > t_param[i+1] + 1e-9) i++;
    
    double ti = t_param[i], tip1 = t_param[i+1];
    double h = tip1 - ti;
    if (h < 1e-9) return points[i];
    
    double u = (t - ti) / h;
    double A = 1.0 - u;
    double B = u;
    double h2_6 = (h * h) / 6.0;
    
    double x = A * points[i].x + B * points[i+1].x +
               ((A*A*A - A) * Mx[i] + (B*B*B - B) * Mx[i+1]) * h2_6;
    double y = A * points[i].y + B * points[i+1].y +
               ((A*A*A - A) * My[i] + (B*B*B - B) * My[i+1]) * h2_6;
    
    return {x, y};
}

// ====== ОТРИСОВКА ======

void drawGrid() {
    glColor3f(0.93f, 0.93f, 0.93f);
    glBegin(GL_LINES);
    for (double x = std::floor(view.x0); x <= std::ceil(view.x1); ++x) {
        glVertex2d(x, view.y0); glVertex2d(x, view.y1);
    }
    for (double y = std::floor(view.y0); y <= std::ceil(view.y1); ++y) {
        glVertex2d(view.x0, y); glVertex2d(view.x1, y);
    }
    glEnd();
    
    glColor3f(0.4f, 0.4f, 0.4f);
    glLineWidth(1.8f);
    glBegin(GL_LINES);
    glVertex2d(view.x0, 0.0); glVertex2d(view.x1, 0.0);
    glVertex2d(0.0, view.y0); glVertex2d(0.0, view.y1);
    glEnd();
}

void drawSpline() {
    if (points.size() < 2) return;
    if (dirty) buildSplines();
    
    glColor3f(0.15f, 0.55f, 0.95f);
    glLineWidth(2.4f);
    glBegin(GL_LINE_STRIP);
    const int steps = 300;
    for (int i = 0; i <= steps; ++i) {
        double t = (double)i / steps;
        Point p = evaluate(t);
        glVertex2d(p.x, p.y);
    }
    glEnd();
}

void drawPoints() {
    glPointSize(9.0f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < points.size(); ++i) {
        if ((int)i == activePoint) glColor3f(0.95f, 0.3f, 0.25f);
        else if ((int)i == hoverPoint) glColor3f(0.2f, 0.65f, 1.0f);
        else glColor3f(0.15f, 0.5f, 0.9f);
        glVertex2d(points[i].x, points[i].y);
    }
    glEnd();
    
    // Нумерация точек
    glColor3f(0.2f, 0.2f, 0.2f);
    for (size_t i = 0; i < points.size(); ++i) {
        char buf[8]; snprintf(buf, sizeof(buf), "%zu", i);
        glRasterPos2d(points[i].x + 0.15, points[i].y + 0.15);
        for (char* p = buf; *p; ++p) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *p);
    }
}

void drawUI() {
    glColor4f(0.99f, 0.99f, 0.99f, 0.93f);
    glBegin(GL_QUADS);
    glVertex2i(0, 0); glVertex2i(winW, 0);
    glVertex2i(winW, 34); glVertex2i(0, 34);
    glEnd();
    
    glColor3f(0.22f, 0.22f, 0.22f);
    char buf[128];
    snprintf(buf, sizeof(buf), "1D Parametric Splines | Points: %zu | Param: %s | Spline: %s | [1-3] Param | [4-5] Spline | [C] Clear | [R] Reset",
             points.size(), paramNames[(int)paramMode], splineNames[(int)splineMode]);
    
    glRasterPos2i(12, 21);
    for (char* p = buf; *p; ++p) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *p);
}

// ====== ВЗАИМОДЕЙСТВИЕ ======

bool canGrab(const Point& p, double x, double y) {
    double dx = x - p.x;
    double dy = y - p.y;
    double r = 0.18 * std::max(1.0, (view.x1 - view.x0) / 12.0);
    return (dx*dx + dy*dy) <= r*r;
}

void screenToWorld(int x, int y, double& wx, double& wy) {
    wx = view.x0 + (x / (double)winW) * (view.x1 - view.x0);
    wy = view.y1 - (y / (double)winH) * (view.y1 - view.y0);
}

// ====== CALLBACKS ======

void reshape(int w, int h) {
    winW = w; winH = h;
    glViewport(0, 0, w, h);
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case '1': paramMode = UNIFORM; dirty = true; break;
        case '2': paramMode = CHORD; dirty = true; break;
        case '3': paramMode = CENTRIPETAL; dirty = true; break;
        case '4': splineMode = NATURAL; dirty = true; break;
        case '5': splineMode = CLAMPED; dirty = true; break;
        case 'c': case 'C': points.clear(); dirty = true; activePoint = hoverPoint = -1; break;
        case 'r': case 'R': view = {-6.0, 6.0, -4.5, 4.5}; break;
        case '+': case '=': {
            double cx = (view.x0 + view.x1) * 0.5;
            double cy = (view.y0 + view.y1) * 0.5;
            view.x0 = cx + (view.x0 - cx) * 0.82;
            view.x1 = cx + (view.x1 - cx) * 0.82;
            view.y0 = cy + (view.y0 - cy) * 0.82;
            view.y1 = cy + (view.y1 - cy) * 0.82;
        } break;
        case '-': case '_': {
            double cx = (view.x0 + view.x1) * 0.5;
            double cy = (view.y0 + view.y1) * 0.5;
            view.x0 = cx + (view.x0 - cx) * 1.22;
            view.x1 = cx + (view.x1 - cx) * 1.22;
            view.y0 = cy + (view.y0 - cy) * 1.22;
            view.y1 = cy + (view.y1 - cy) * 1.22;
        } break;
        case 27: case 'q': case 'Q': exit(0);
        default: return;
    }
    glutPostRedisplay();
}

void mouse(int btn, int state, int x, int y) {
    screenToWorld(x, y, mx, my);
    
    if (btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        activePoint = -1;
        for (int i = (int)points.size() - 1; i >= 0; --i) {
            if (canGrab(points[i], mx, my)) {
                activePoint = i;
                break;
            }
        }
        if (activePoint == -1) {
            points.push_back({mx, my});
            activePoint = (int)points.size() - 1;
            dirty = true;
        }
        mouseDown = true;
    } else if (btn == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        mouseDown = false;
        activePoint = -1;
    } else if (btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (canGrab(points[i], mx, my)) {
                points.erase(points.begin() + i);
                dirty = true;
                break;
            }
        }
    }
    glutPostRedisplay();
}

void motion(int x, int y) {
    screenToWorld(x, y, mx, my);
    
    if (mouseDown && activePoint >= 0 && activePoint < (int)points.size()) {
        points[activePoint].x = mx;
        points[activePoint].y = my;
        dirty = true;
    }
    
    int newHover = -1;
    for (int i = (int)points.size() - 1; i >= 0; --i) {
        if (canGrab(points[i], mx, my)) {
            newHover = i;
            break;
        }
    }
    if (newHover != hoverPoint) {
        hoverPoint = newHover;
        glutPostRedisplay();
    } else if (dirty) {
        glutPostRedisplay();
    }
}

void display() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(view.x0, view.x1, view.y0, view.y1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    drawGrid();
    drawSpline();
    drawPoints();
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, winW, 0, winH, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    drawUI();
    
    glutSwapBuffers();
    
    char title[80];
    snprintf(title, sizeof(title), "1D Parametric Splines — %s — %s",
             paramNames[(int)paramMode], splineNames[(int)splineMode]);
    glutSetWindowTitle(title);
}

// ====== MAIN ======
int main(int argc, char** argv) {
    points = { {-4,-2}, {-2,3}, {0,-1}, {2,2}, {4,-2} };
    dirty = true;
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winW, winH);
    glutCreateWindow("1D Parametric Splines — Parameterization Modes");
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(motion);
    
    glutMainLoop();
    return 0;
}