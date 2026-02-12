#include "GL/glut.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>

// ====== СТЕПЕНИ СПЛАЙНОВ ======
enum SplineDegree { LINEAR, QUADRATIC, CUBIC_NATURAL };
SplineDegree degreeMode = CUBIC_NATURAL;
const char* degreeNames[] = {"Linear (C⁰)", "Quadratic (C¹)", "Cubic Natural (C²)"};

// ====== СХЕМЫ ПАРАМЕТРИЗАЦИИ ======
enum ParamType { UNIFORM, CHORD, CENTRIPETAL };
ParamType paramMode = CHORD;
const char* paramNames[] = {"Uniform (t=i)", "Chord Length", "Centripetal (√chord)"};

// ====== ГЛОБАЛЬНОЕ СОСТОЯНИЕ ======
int winW = 1050, winH = 720;
double mx = 0.0, my = 0.0;
int activePoint = -1, hoverPoint = -1;
bool mouseDown = false;
bool showTValues = true;

struct Point { double x, y; Point(double x=0,double y=0):x(x),y(y){} };
std::vector<Point> points;

// Параметрические узлы (нормированные в [0, 1])
std::vector<double> t_param;

// Коэффициенты сплайнов для каждого сегмента: S(u) = a + b*u + c*u^2 + d*u^3, u∈[0,1]
// где u — локальный параметр сегмента: u = (t - t_i) / (t_{i+1} - t_i)
struct Coeffs { double a, b, c, d; };
std::vector<Coeffs> coeffsX, coeffsY;
bool dirty = true;

struct View { double x0, x1, y0, y1; } view = {-6.0, 6.0, -4.5, 4.5};

// ====== ВЫЧИСЛЕНИЕ ПАРАМЕТРА t_i ======
void computeParameters() {
    size_t n = points.size();
    t_param.resize(n);
    if (n == 0) return;
    
    t_param[0] = 0.0;
    
    switch (paramMode) {
        case UNIFORM:
            // t_i = i (равномерная параметризация)
            for (size_t i = 1; i < n; ++i)
                t_param[i] = (double)i;
            break;
            
        case CHORD:
            // t_i = t_{i-1} + ||P_i - P_{i-1}|| (по длине хорд)
            for (size_t i = 1; i < n; ++i) {
                double dx = points[i].x - points[i-1].x;
                double dy = points[i].y - points[i-1].y;
                t_param[i] = t_param[i-1] + std::sqrt(dx*dx + dy*dy);
            }
            break;
            
        case CENTRIPETAL:
            // t_i = t_{i-1} + ||P_i - P_{i-1}||^{1/2} (центрипетальная)
            for (size_t i = 1; i < n; ++i) {
                double dx = points[i].x - points[i-1].x;
                double dy = points[i].y - points[i-1].y;
                double chord = std::sqrt(dx*dx + dy*dy);
                t_param[i] = t_param[i-1] + std::sqrt(chord); // chord^{1/2}
            }
            break;
    }
    
    // Нормировка в [0, 1] для численной устойчивости
    double total = t_param[n-1];
    if (total > 1e-6)
        for (size_t i = 0; i < n; ++i) t_param[i] /= total;
}

// ====== ЛИНЕЙНЫЙ СПЛАЙН (С⁰) ======
// Для сегмента [t_i, t_{i+1}] вводим локальный параметр u = (t - t_i) / h ∈ [0,1], h = t_{i+1} - t_i
// Интерполяция: S(u) = (1-u)*P_i + u*P_{i+1} = P_i + (P_{i+1} - P_i)*u
// Коэффициенты: a = P_i, b = P_{i+1} - P_i, c = d = 0
void buildLinear() {
    size_t n = points.size();
    coeffsX.clear(); coeffsY.clear();
    if (n < 2) return;
    
    coeffsX.resize(n-1); coeffsY.resize(n-1);
    
    for (size_t i = 0; i < n-1; ++i) {
        // X координата
        double x0 = points[i].x;
        double x1 = points[i+1].x;
        coeffsX[i] = {x0, x1 - x0, 0.0, 0.0};
        
        // Y координата
        double y0 = points[i].y;
        double y1 = points[i+1].y;
        coeffsY[i] = {y0, y1 - y0, 0.0, 0.0};
    }
}

// ====== КВАДРАТИЧНЫЙ СПЛАЙН (С¹) ======
// На каждом сегменте: S(t) = y_i + b_i*(t-t_i) + c_i*(t-t_i)^2
// где c_i = (δ_i - b_i) / h_i, δ_i = (y_{i+1} - y_i) / h_i
// 
// Граничное условие: задаем начальную скорость равной средней на первом сегменте
// b_0 = δ_0 (эквивалентно "нулевой начальной кривизне" для квадратичного сплайна)
// 
// Условие непрерывности первой производной в узле i дает рекуррентную формулу:
// b_{i+1} = 2*δ_i - b_i
//
// Преобразование к каноническому виду по локальному параметру u ∈ [0,1]:
// S(u) = y_i + (b_i*h)*u + ((δ_i - b_i)/h * h^2)*u^2 = y_i + (b_i*h)*u + ((δ_i - b_i)*h)*u^2
void buildQuadratic() {
    size_t n = points.size();
    coeffsX.clear(); coeffsY.clear();
    if (n < 2) return;
    
    // Шаги по параметру
    std::vector<double> h(n-1);
    for (size_t i = 0; i < n-1; ++i)
        h[i] = std::max(1e-9, t_param[i+1] - t_param[i]);
    
    // Первые производные в узлах (скорости по параметру t)
    std::vector<double> bX(n), bY(n);
    
    // Граничное условие на левом конце: начальная скорость равна средней на первом сегменте
    // Это обеспечивает "естественное" поведение в начале кривой
    if (n >= 2) {
        bX[0] = (points[1].x - points[0].x) / h[0];  // = delta_x[0]
        bY[0] = (points[1].y - points[0].y) / h[0];  // = delta_y[0]
    }
    
    // Рекуррентное вычисление производных для внутренних узлов
    // из условия C¹-непрерывности: b_i + 2*c_i*h_i = b_{i+1}
    for (size_t i = 0; i < n-2; ++i) {
        double delta_x = (points[i+2].x - points[i+1].x) / h[i+1];
        bX[i+1] = 2.0 * ((points[i+1].x - points[i].x) / h[i]) - bX[i];
        
        double delta_y = (points[i+2].y - points[i+1].y) / h[i+1];
        bY[i+1] = 2.0 * ((points[i+1].y - points[i].y) / h[i]) - bY[i];
    }
    
    // Коэффициенты для каждого сегмента в канонической форме по u ∈ [0,1]
    coeffsX.resize(n-1); coeffsY.resize(n-1);
    
    for (size_t i = 0; i < n-1; ++i) {
        // X координата
        double x0 = points[i].x;
        double bx = bX[i];
        double delta_x = (points[i+1].x - points[i].x) / h[i];
        double cx = (delta_x - bx) / h[i];  // коэффициент при (t-t_i)^2
        
        // Преобразование к форме по u: S(u) = a + b*u + c*u^2
        // где u = (t - t_i)/h, поэтому (t-t_i) = u*h, (t-t_i)^2 = u^2*h^2
        coeffsX[i] = {
            x0,           // a = начальное значение
            bx * h[i],    // b = производная * шаг = изменение за u от 0 до 1 по касательной
            cx * h[i] * h[i],  // c = коэффициент квадратичного члена, масштабированный
            0.0           // d = 0 (квадратичный сплайн)
        };
        
        // Y координата (аналогично)
        double y0 = points[i].y;
        double by = bY[i];
        double delta_y = (points[i+1].y - points[i].y) / h[i];
        double cy = (delta_y - by) / h[i];
        
        coeffsY[i] = {
            y0,
            by * h[i],
            cy * h[i] * h[i],
            0.0
        };
    }
}

// ====== КУБИЧЕСКИЙ НАТУРАЛЬНЫЙ СПЛАЙН (С²) ======
// Граничные условия: S''(t_0) = S''(t_{n-1}) = 0 (естественный сплайн)
// 
// Для вычисления моментов M_i = S''(t_i) решаем трёхдиагональную систему:
// h_{i-1}/6 * M_{i-1} + (h_{i-1}+h_i)/3 * M_i + h_i/6 * M_{i+1} = δ_i - δ_{i-1}
// где δ_i = (y_{i+1} - y_i) / h_i — разделённые разности
// 
// На сегменте [t_i, t_{i+1}] с локальным параметром u = (t-t_i)/h:
// S(u) = (1-u)*y_i + u*y_{i+1} + h^2/6 * [((1-u)^3-(1-u))*M_i + (u^3-u)*M_{i+1}]
//
// Раскладывая по степеням u ∈ [0,1]:
// S(u) = a + b*u + c*u^2 + d*u^3, где:
//   a = y_i
//   b = (y_{i+1}-y_i) - h^2/6 * (2*M_i + M_{i+1})
//   c = h^2/2 * M_i
//   d = h^2/6 * (M_{i+1} - M_i)
void buildCubicNatural() {
    size_t n = points.size();
    coeffsX.clear(); coeffsY.clear();
    if (n < 2) return;
    
    // Шаги по параметру
    std::vector<double> h(n-1);
    for (size_t i = 0; i < n-1; ++i)
        h[i] = std::max(1e-9, t_param[i+1] - t_param[i]);
    
    // Разделённые разности (скорости изменения по параметру)
    std::vector<double> deltaX(n-1), deltaY(n-1);
    for (size_t i = 0; i < n-1; ++i) {
        deltaX[i] = (points[i+1].x - points[i].x) / h[i];
        deltaY[i] = (points[i+1].y - points[i].y) / h[i];
    }
    
    // Вторые производные в узлах (моменты), граничные условия M[0] = M[n-1] = 0
    std::vector<double> MX(n, 0.0), MY(n, 0.0);
    
    if (n > 2) {
        size_t m = n - 2; // количество внутренних узлов
        std::vector<double> a(m), b(m), c(m), dX(m), dY(m);
        
        // Формирование трёхдиагональной системы для метода прогонки
        // Диагонали: нижняя a, главная b, верхняя c
        for (size_t i = 0; i < m; ++i) {
            a[i] = h[i] / 6.0;                    // коэффициент при M_i
            b[i] = (h[i] + h[i+1]) / 3.0;         // коэффициент при M_{i+1}
            c[i] = h[i+1] / 6.0;                  // коэффициент при M_{i+2}
            dX[i] = deltaX[i+1] - deltaX[i];      // правая часть для X
            dY[i] = deltaY[i+1] - deltaY[i];      // правая часть для Y
        }
        
        // Прямая прогонка (алгоритм Томаса) с проверкой на вырожденность
        const double EPS = 1e-12;
        for (size_t i = 1; i < m; ++i) {
            if (std::abs(b[i-1]) < EPS) {
                // Вырожденный случай — используем приближение
                b[i-1] = EPS;
            }
            double w = a[i] / b[i-1];
            b[i] -= w * c[i-1];
            dX[i] -= w * dX[i-1];
            dY[i] -= w * dY[i-1];
        }
        
        // Обратная подстановка
        if (std::abs(b[m-1]) < EPS) b[m-1] = EPS;
        
        MX[n-2] = dX[m-1] / b[m-1];
        MY[n-2] = dY[m-1] / b[m-1];
        
        for (int i = (int)m - 2; i >= 0; --i) {
            MX[i+1] = (dX[i] - c[i] * MX[i+2]) / b[i];
            MY[i+1] = (dY[i] - c[i] * MY[i+2]) / b[i];
        }
        // MX[0] = MX[n-1] = 0.0 (натуральные граничные условия уже установлены)
    }
    
    // Вычисление коэффициентов для каждого сегмента
    coeffsX.resize(n-1); coeffsY.resize(n-1);
    
    for (size_t i = 0; i < n-1; ++i) {
        double hh = h[i] * h[i];  // h_i^2 для компактности формул
        
        // X координата
        double x0 = points[i].x;
        double x1 = points[i+1].x;
        double M0 = MX[i];
        double M1 = MX[i+1];
        
        coeffsX[i] = {
            x0,                                           // a
            (x1 - x0) - hh/6.0 * (2.0*M0 + M1),          // b
            hh/2.0 * M0,                                  // c
            hh/6.0 * (M1 - M0)                            // d
        };
        
        // Y координата
        double y0 = points[i].y;
        double y1 = points[i+1].y;
        M0 = MY[i];
        M1 = MY[i+1];
        
        coeffsY[i] = {
            y0,                                           // a
            (y1 - y0) - hh/6.0 * (2.0*M0 + M1),          // b
            hh/2.0 * M0,                                  // c
            hh/6.0 * (M1 - M0)                            // d
        };
    }
}

// ====== ПОСТРОЕНИЕ СПЛАЙНА ======
void buildSplines() {
    if (points.size() < 2) { dirty = false; return; }
    
    // Примечание: точки используются в порядке их задания пользователем,
    // что соответствует логике параметрических кривых (порядок следования по кривой)
    
    computeParameters();
    
    switch (degreeMode) {
        case LINEAR:         buildLinear(); break;
        case QUADRATIC:      buildQuadratic(); break;
        case CUBIC_NATURAL:  buildCubicNatural(); break;
    }
    
    dirty = false;
}

// ====== ОЦЕНКА СПЛАЙНА ======
// Вычисление точки на сплайне для глобального параметра t ∈ [t_0, t_{n-1}]
// Используется бинарный поиск для определения сегмента [t_i, t_{i+1}]
// и схема Горнера для вычисления полинома по локальному параметру u ∈ [0,1]
Point evaluate(double t) {
    if (points.size() < 2) return Point(0, 0);
    if (dirty) buildSplines();
    
    // Ограничение параметра диапазоном узлов (экстраполяция не производится)
    if (t <= t_param[0]) return points[0];
    if (t >= t_param.back()) return points.back();
    
    // Бинарный поиск сегмента i такого, что t_param[i] <= t < t_param[i+1]
    size_t lo = 0, hi = t_param.size() - 2;
    while (lo < hi) {
        size_t mid = (lo + hi + 1) / 2;
        if (t_param[mid] <= t) lo = mid;
        else hi = mid - 1;
    }
    size_t i = lo;
    
    // Вычисление локального параметра u ∈ [0,1] для сегмента i
    double ti = t_param[i];
    double tip1 = t_param[i+1];
    double h = tip1 - ti;
    double u = (h < 1e-9) ? 0.0 : (t - ti) / h;
    u = std::max(0.0, std::min(1.0, u));  // защита от ошибок округления
    
    // Вычисление по схеме Горнера: S(u) = a + u*(b + u*(c + u*d))
    const auto& cx = coeffsX[i];
    const auto& cy = coeffsY[i];
    
    double x = cx.a + u * (cx.b + u * (cx.c + u * cx.d));
    double y = cy.a + u * (cy.b + u * (cy.c + u * cy.d));
    
    return Point(x, y);
}

// ====== ВИЗУАЛИЗАЦИЯ ======

void drawGrid() {
    // Фоновая сетка
    glColor3f(0.94f, 0.94f, 0.94f);
    glBegin(GL_LINES);
    for (double x = std::floor(view.x0); x <= std::ceil(view.x1) + 1e-9; x += 0.5) {
        glVertex2d(x, view.y0); glVertex2d(x, view.y1);
    }
    for (double y = std::floor(view.y0); y <= std::ceil(view.y1) + 1e-9; y += 0.5) {
        glVertex2d(view.x0, y); glVertex2d(view.x1, y);
    }
    glEnd();
    
    // Оси координат
    glColor3f(0.5f, 0.5f, 0.5f);
    glLineWidth(1.6f);
    glBegin(GL_LINES);
    glVertex2d(view.x0, 0.0); glVertex2d(view.x1, 0.0);
    glVertex2d(0.0, view.y0); glVertex2d(0.0, view.y1);
    glEnd();
    glLineWidth(1.0f);  // сброс толщины
}

void drawSpline() {
    if (points.size() < 2) return;
    if (dirty) buildSplines();
    
    // Цвет по степени гладкости
    switch (degreeMode) {
        case LINEAR:        glColor3f(0.65f, 0.65f, 0.65f); break; // серый — C⁰
        case QUADRATIC:     glColor3f(0.25f, 0.75f, 0.35f); break; // зелёный — C¹
        case CUBIC_NATURAL: glColor3f(0.18f, 0.52f, 0.92f); break; // синий — C²
    }
    
    glLineWidth(2.6f);
    glBegin(GL_LINE_STRIP);
    
    const int subdivisions = 250;  // плотность разбиения для гладкости
    double t0 = t_param.front();
    double t1 = t_param.back();
    
    for (int i = 0; i <= subdivisions; ++i) {
        double t = t0 + (t1 - t0) * i / subdivisions;
        Point p = evaluate(t);
        glVertex2d(p.x, p.y);
    }
    glEnd();
    
    // Цветовая шкала параметра t (градиент от фиолетового к оранжевому)
    if (showTValues) {
        glLineWidth(4.2f);
        glBegin(GL_LINE_STRIP);
        const int colorSteps = 40;
        for (int i = 0; i <= colorSteps; ++i) {
            double t = t0 + (t1 - t0) * i / (double)colorSteps;
            Point p = evaluate(t);
            double alpha = (double)i / colorSteps;  // нормированное значение t
            // Градиент: фиолетовый (0.2, 0.0, 0.9) -> оранжевый (0.9, 0.5, 0.2)
            glColor3f(
                (float)(0.2 + 0.7*alpha), 
                (float)(0.0 + 0.5*alpha), 
                (float)(0.9 - 0.7*alpha)
            );
            glVertex2d(p.x, p.y);
        }
        glEnd();
    }
    glLineWidth(1.0f);  // сброс
}

void drawPoints() {
    // Контрольные точки
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < points.size(); ++i) {
        if ((int)i == activePoint) glColor3f(0.95f, 0.25f, 0.2f);      // красная — активная
        else if ((int)i == hoverPoint) glColor3f(0.2f, 0.6f, 1.0f);    // голубая — под курсором
        else glColor3f(0.15f, 0.45f, 0.85f);                           // синяя — обычная
        glVertex2d(points[i].x, points[i].y);
    }
    glEnd();
    
    // Подписи значений параметра t_i
    if (showTValues && t_param.size() == points.size()) {
        glColor3f(0.18f, 0.18f, 0.18f);
        for (size_t i = 0; i < points.size(); ++i) {
            char buf[16];
            snprintf(buf, sizeof(buf), "t=%.2f", t_param[i]);
            glRasterPos2d(points[i].x + 0.16, points[i].y - 0.26);
            for (char* p = buf; *p; ++p) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *p);
        }
    }
    
    // Нумерация точек
    glColor3f(0.32f, 0.32f, 0.32f);
    for (size_t i = 0; i < points.size(); ++i) {
        char buf[8];
        snprintf(buf, sizeof(buf), "%zu", i);
        glRasterPos2d(points[i].x - 0.13, points[i].y + 0.22);
        for (char* p = buf; *p; ++p) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *p);
    }
    
    // Ломаная через контрольные точки (опорная линия)
    glColor4f(0.75f, 0.75f, 0.75f, 0.75f);
    glLineWidth(1.1f);
    glBegin(GL_LINE_STRIP);
    for (const auto& p : points) glVertex2d(p.x, p.y);
    glEnd();
    glLineWidth(1.0f);
}

void drawUI() {
    // Фон статус-бара
    glColor4f(0.99f, 0.99f, 0.99f, 0.95f);
    glBegin(GL_QUADS);
    glVertex2i(0, 0); glVertex2i(winW, 0);
    glVertex2i(winW, 38); glVertex2i(0, 38);
    glEnd();
    
    // Рамка
    glColor3f(0.88f, 0.88f, 0.88f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(0, 0); glVertex2i(winW, 0);
    glVertex2i(winW, 38); glVertex2i(0, 38);
    glEnd();
    
    // Текст статуса
    glColor3f(0.18f, 0.18f, 0.18f);
    char status[256];
    snprintf(status, sizeof(status),
        "1D Parametric Splines | Points: %zu | Degree: %s | Param: %s | "
        "[1-3] Degree | [4-6] Param | [V] t-values | [C] Clear | [R] Reset | "
        "LMB: add/drag | RMB: delete | [+/-] Zoom",
        points.size(), degreeNames[(int)degreeMode], paramNames[(int)paramMode]);
    
    glRasterPos2i(14, 24);
    for (char* p = status; *p; ++p) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *p);
}

// ====== ВЗАИМОДЕЙСТВИЕ ======

// Проверка попадания курсора в окрестность точки
bool canGrab(const Point& p, double x, double y) {
    double dx = x - p.x;
    double dy = y - p.y;
    // Адаптивный радиус захвата в зависимости от масштаба
    double r = 0.21 * std::max(1.0, (view.x1 - view.x0) / 13.0);
    return (dx*dx + dy*dy) <= r*r;
}

// Преобразование экранных координат в мировые
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
        case '1': degreeMode = LINEAR; dirty = true; break;
        case '2': degreeMode = QUADRATIC; dirty = true; break;
        case '3': degreeMode = CUBIC_NATURAL; dirty = true; break;
        case '4': paramMode = UNIFORM; dirty = true; break;
        case '5': paramMode = CHORD; dirty = true; break;
        case '6': paramMode = CENTRIPETAL; dirty = true; break;
        case 'v': case 'V': showTValues = !showTValues; break;
        case 'c': case 'C': 
            points.clear(); 
            dirty = true; 
            activePoint = hoverPoint = -1; 
            break;
        case 'r': case 'R': 
            view = {-6.0, 6.0, -4.5, 4.5}; 
            break;
        case '+': case '=': {  // приближение
            double cx = (view.x0 + view.x1) * 0.5;
            double cy = (view.y0 + view.y1) * 0.5;
            double scale = 0.83;  // коэффициент масштабирования
            view.x0 = cx + (view.x0 - cx) * scale;
            view.x1 = cx + (view.x1 - cx) * scale;
            view.y0 = cy + (view.y0 - cy) * scale;
            view.y1 = cy + (view.y1 - cy) * scale;
        } break;
        case '-': case '_': {  // отдаление
            double cx = (view.x0 + view.x1) * 0.5;
            double cy = (view.y0 + view.y1) * 0.5;
            double scale = 1.20;  // обратный коэффициент
            view.x0 = cx + (view.x0 - cx) * scale;
            view.x1 = cx + (view.x1 - cx) * scale;
            view.y0 = cy + (view.y0 - cy) * scale;
            view.y1 = cy + (view.y1 - cy) * scale;
        } break;
        case 27: case 'q': case 'Q': exit(0);
        default: return;
    }
    glutPostRedisplay();
}

void mouse(int btn, int state, int x, int y) {
    screenToWorld(x, y, mx, my);
    
    if (btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        // Поиск существующей точки для захвата (в обратном порядке — верхние сверху)
        activePoint = -1;
        for (int i = (int)points.size() - 1; i >= 0; --i) {
            if (canGrab(points[i], mx, my)) {
                activePoint = i;
                break;
            }
        }
        
        // Если не захвачена существующая — добавляем новую
        if (activePoint == -1) {
            // Проверка минимального расстояния до существующих точек
            const double MIN_DIST_SQ = 0.25;  // (0.5)^2
            bool tooClose = false;
            for (const auto& p : points) {
                double dx = mx - p.x;
                double dy = my - p.y;
                if (dx*dx + dy*dy < MIN_DIST_SQ) { 
                    tooClose = true; 
                    break; 
                }
            }
            if (!tooClose) {
                points.push_back(Point(mx, my));
                activePoint = (int)points.size() - 1;
                dirty = true;
            }
        }
        mouseDown = true;
    } 
    else if (btn == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        mouseDown = false;
        activePoint = -1;
    } 
    else if (btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        // Удаление точки под курсором
        for (size_t i = 0; i < points.size(); ++i) {
            if (canGrab(points[i], mx, my)) {
                points.erase(points.begin() + i);
                dirty = true;
                if (hoverPoint == (int)i) hoverPoint = -1;
                break;
            }
        }
    }
    glutPostRedisplay();
}

void motion(int x, int y) {
    screenToWorld(x, y, mx, my);
    
    // Перетаскивание активной точки
    if (mouseDown && activePoint >= 0 && activePoint < (int)points.size()) {
        const double MIN_DIST_SQ = 0.25;
        bool tooClose = false;
        for (size_t i = 0; i < points.size(); ++i) {
            if ((int)i == activePoint) continue;
            double dx = mx - points[i].x;
            double dy = my - points[i].y;
            if (dx*dx + dy*dy < MIN_DIST_SQ) { 
                tooClose = true; 
                break; 
            }
        }
        if (!tooClose) {
            points[activePoint].x = mx;
            points[activePoint].y = my;
            dirty = true;
        }
    }
    
    // Обновление состояния наведения
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
    
    // Мировые координаты для графики
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(view.x0, view.x1, view.y0, view.y1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    drawGrid();
    drawSpline();
    drawPoints();
    
    // Экранное пространство для интерфейса
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, winW, 0, winH, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    drawUI();
    
    glutSwapBuffers();
    
    // Обновление заголовка окна
    char title[128];
    snprintf(title, sizeof(title), "1D Parametric Splines | %s | %s",
             degreeNames[(int)degreeMode], paramNames[(int)paramMode]);
    glutSetWindowTitle(title);
}

// ====== MAIN ======
int main(int argc, char** argv) {
    // Начальные точки, демонстрирующие влияние параметризации
    points = {
        Point(-5.0, -1.5), Point(-3.0, 3.0), Point(-1.5, -2.8),
        Point(0.5, 2.2), Point(2.5, -2.5), Point(4.5, 1.8)
    };
    dirty = true;
    
    // Резервирование памяти для оптимизации
    points.reserve(50);
    t_param.reserve(50);
    coeffsX.reserve(50);
    coeffsY.reserve(50);
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(winW, winH);
    glutCreateWindow("1D Parametric Splines — Verified Mathematics");
    
    // Настройка сглаживания
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
    // Регистрация callback-функций
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(motion);
    
    glutMainLoop();
    return 0;
}