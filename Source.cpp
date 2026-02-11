#include "E:\git_work\Computer_Graphics_Course\rgz\GL\glut.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <array>



// ====== Константы ======
const double POINT_GRAB_RADIUS = 0.5;
const int MIN_POINTS = 4; // Для сплайна Лагранжа 3-й степени нужно минимум 4 точки
const double LABEL_BASE_HEIGHT = 0.25;
const double TEXT_DEFAULT_HEIGHT = 119.05;

// ====== Глобальные переменные ======
int windowWidth = 800, windowHeight = 600;
struct Point {
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};
std::vector<Point> points;
int moving_point = -1;
int hover_point = -1;
Point mouseWorld;
Point cameraPos = { 0.0, 0.0 };
Point cameraSize = { 10.0, 10.0 };
double step = 0.01;
double label_step = 1.0;

// ====== Полином Лагранжа 3-й степени по 4 точкам ======
double lagrange3(double x, const std::array<Point, 4>& pts) {
    double L = 0.0;
    for (int i = 0; i < 4; ++i) {
        double term = pts[i].y;
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                double denom = pts[i].x - pts[j].x;
                if (std::abs(denom) < 1e-12) return 0.0;
                term *= (x - pts[j].x) / denom;
            }
        }
        L += term;
    }
    return L;
}

double lagrange3Derivative(double x, const std::array<Point, 4>& pts) {
    double dL = 0.0;
    for (int i = 0; i < 4; ++i) {
        double sum = 0.0;
        for (int k = 0; k < 4; ++k) {
            if (k == i) continue;
            double prod = 1.0;
            for (int j = 0; j < 4; ++j) {
                if (j == i || j == k) continue;
                double denom = pts[i].x - pts[j].x;
                if (std::abs(denom) < 1e-12) return 0.0;
                prod *= (x - pts[j].x) / denom;
            }
            double denom = pts[i].x - pts[k].x;
            if (std::abs(denom) < 1e-12) return 0.0;
            sum += prod / denom;
        }
        dL += pts[i].y * sum;
    }
    return dL;
}

// ====== Вспомогательные функции ======
bool canGrab(const Point& p) {
    double dx = mouseWorld.x - p.x;
    double dy = mouseWorld.y - p.y;
    return std::sqrt(dx * dx + dy * dy) <= POINT_GRAB_RADIUS;
}

bool canAdd(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy) > 2 * POINT_GRAB_RADIUS;
}

void updateMouseWorld(int x, int y) {
    double aspect = (double)windowWidth / windowHeight;
    cameraSize.x = cameraSize.y * aspect;
    mouseWorld.x = (double)x / windowWidth * (2 * cameraSize.x) - cameraSize.x + cameraPos.x;
    mouseWorld.y = (1.0 - (double)y / windowHeight) * (2 * cameraSize.y) - cameraSize.y + cameraPos.y;
}

void calcCamera() {
    if (cameraSize.y > 100.0) cameraSize.y = 100.0;
    if (cameraSize.y < 0.1) cameraSize.y = 0.1;
    double aspect = (double)windowWidth / windowHeight;
    cameraSize.x = cameraSize.y * aspect;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(cameraPos.x - cameraSize.x, cameraPos.x + cameraSize.x,
        cameraPos.y - cameraSize.y, cameraPos.y + cameraSize.y,
        -1, 1);
    glMatrixMode(GL_MODELVIEW);

    step = cameraSize.x / windowWidth;
    label_step = 100.0 * cameraSize.x / windowWidth;
    double log_step = std::floor(std::log10(label_step));
    double mantissa = label_step / std::pow(10.0, log_step);
    if (mantissa < 1.5) mantissa = 1.0;
    else if (mantissa < 3.0) mantissa = 2.0;
    else if (mantissa < 7.0) mantissa = 5.0;
    else mantissa = 10.0;
    label_step = mantissa * std::pow(10.0, log_step);
}

void updateWindowTitle() {
    char title[256];
    bool duplicateX = false;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        if (std::abs(points[i].x - points[i + 1].x) < 1e-9) {
            duplicateX = true;
            snprintf(title, sizeof(title), "КГ Лаб 4 | Ошибка: несколько точек в одном x! (%.2f)", points[i].x);
            break;
        }
    }
    if (!duplicateX) {
        if ((int)points.size() < MIN_POINTS) {
            snprintf(title, sizeof(title), "КГ Лаб 4 | Недостаточно точек! (%d / %d)", (int)points.size(), MIN_POINTS);
        }
        else {
            snprintf(title, sizeof(title), "КГ Лаб 4 | Точек: %d; красным - сплайн, синим - производная", (int)points.size());
        }
    }
    glutSetWindowTitle(title);
}

void renderText(double x, double y, double height, const char* text) {
    double scale = height / TEXT_DEFAULT_HEIGHT;
    glPushMatrix();
    glTranslated(x, y - height, 0.0);
    glScaled(scale, scale, 1.0);
    for (const char* c = text; *c != '\0'; ++c) {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, *c);
    }
    glPopMatrix();
}

// ====== Обработчики меню ======
void menuCallback(int value) {
    switch (value) {
    case 1: // Добавить точку
        if (canAdd(Point(), mouseWorld)) {
            bool tooClose = false;
            for (const auto& p : points) {
                if (!canAdd(p, mouseWorld)) {
                    tooClose = true;
                    break;
                }
            }
            if (!tooClose) {
                points.push_back(mouseWorld);
                std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
                    return a.x < b.x;
                    });
            }
        }
        break;
    case 2: // Удалить точку под курсором
        for (size_t i = 0; i < points.size(); ++i) {
            if (canGrab(points[i])) {
                points.erase(points.begin() + i);
                break;
            }
        }
        break;
    case 3: // Очистить все
        points.clear();
        break;
    case 10: cameraPos.y += 1.0; break;
    case 11: cameraPos.y -= 1.0; break;
    case 12: cameraPos.x -= 1.0; break;
    case 13: cameraPos.x += 1.0; break;
    case 14: cameraSize.y *= 1.2; break; // Отдалить
    case 15: cameraSize.y *= 0.8; break; // Приблизить
    }
    calcCamera();
    glutPostRedisplay();
}

void createMenu() {
    int cameraMenu = glutCreateMenu(menuCallback);
    glutAddMenuEntry("Сдвиг вверх", 10);
    glutAddMenuEntry("Сдвиг вниз", 11);
    glutAddMenuEntry("Сдвиг влево", 12);
    glutAddMenuEntry("Сдвиг вправо", 13);
    glutAddMenuEntry("Отдалить", 14);
    glutAddMenuEntry("Приблизить", 15);

    int mainMenu = glutCreateMenu(menuCallback);
    glutAddSubMenu("Область отображения", cameraMenu);
    glutAddMenuEntry("Добавить точку", 1);
    glutAddMenuEntry("Удалить точку под курсором", 2);
    glutAddMenuEntry("Очистить", 3);

    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

// ====== Callbacks ======
void reshape(int w, int h) {
    windowWidth = w;
    windowHeight = h;
    glViewport(0, 0, w, h);
    calcCamera();
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
    case 'w': cameraPos.y += 1.0; break;
    case 's': cameraPos.y -= 1.0; break;
    case 'a': cameraPos.x -= 1.0; break;
    case 'd': cameraPos.x += 1.0; break;
    case 'q': cameraSize.y *= 1.2; break;
    case 'e': cameraSize.y *= 0.8; break;
    case 'z': points.clear(); break;
    default: return;
    }
    calcCamera();
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    updateMouseWorld(x, y);
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (canGrab(points[i])) {
                moving_point = (int)i;
                return;
            }
        }
        bool tooClose = false;
        for (const auto& p : points) {
            if (!canAdd(p, mouseWorld)) {
                tooClose = true;
                break;
            }
        }
        if (!tooClose) {
            points.push_back(mouseWorld);
            std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
                return a.x < b.x;
                });
        }
    }
    else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        moving_point = -1;
    }
    else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
        for (size_t i = 0; i < points.size(); ++i) {
            if (canGrab(points[i])) {
                points.erase(points.begin() + i);
                break;
            }
        }
    }
    glutPostRedisplay();
}

void motion(int x, int y) {
    updateMouseWorld(x, y);
    if (moving_point != -1 && moving_point < (int)points.size()) {
        bool conflict = false;
        for (size_t i = 0; i < points.size(); ++i) {
            if ((int)i == moving_point) continue;
            if (!canAdd(points[i], mouseWorld)) {
                conflict = true;
                break;
            }
        }
        if (!conflict) {
            points[moving_point] = mouseWorld;
            if (moving_point > 0 && points[moving_point].x < points[moving_point - 1].x) {
                std::swap(points[moving_point], points[moving_point - 1]);
                moving_point--;
            }
            else if (moving_point < (int)points.size() - 1 && points[moving_point].x > points[moving_point + 1].x) {
                std::swap(points[moving_point], points[moving_point + 1]);
                moving_point++;
            }
        }
    }
    glutPostRedisplay();
}

void passiveMotion(int x, int y) {
    updateMouseWorld(x, y);
    int new_hover = -1;
    for (size_t i = 0; i < points.size(); ++i) {
        if (canGrab(points[i])) {
            new_hover = (int)i;
            break;
        }
    }
    if (new_hover != hover_point) {
        hover_point = new_hover;
        glutPostRedisplay();
    }
}

void display() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Оси и сетка
    glColor3f(0.0f, 0.0f, 0.0f);
    glLineWidth(2.0f);
    glBegin(GL_LINES);
    glVertex2d(cameraPos.x - cameraSize.x, 0.0);
    glVertex2d(cameraPos.x + cameraSize.x, 0.0);
    glVertex2d(0.0, cameraPos.y - cameraSize.y);
    glVertex2d(0.0, cameraPos.y + cameraSize.y);
    glEnd();

    glColor3f(0.8f, 0.8f, 0.8f);
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    double from_x = std::floor((cameraPos.x - cameraSize.x) / label_step) * label_step;
    double to_x = std::ceil((cameraPos.x + cameraSize.x) / label_step) * label_step;
    for (double x = from_x; x <= to_x; x += label_step) {
        glVertex2d(x, cameraPos.y - cameraSize.y);
        glVertex2d(x, cameraPos.y + cameraSize.y);
    }
    double from_y = std::floor((cameraPos.y - cameraSize.y) / label_step) * label_step;
    double to_y = std::ceil((cameraPos.y + cameraSize.y) / label_step) * label_step;
    for (double y = from_y; y <= to_y; y += label_step) {
        glVertex2d(cameraPos.x - cameraSize.x, y);
        glVertex2d(cameraPos.x + cameraSize.x, y);
    }
    glEnd();

    // Подписи
    glColor3f(0.0f, 0.0f, 0.0f);
    from_x = std::floor((cameraPos.x - cameraSize.x) / label_step) * label_step;
    to_x = std::ceil((cameraPos.x + cameraSize.x) / label_step) * label_step;
    for (double x = from_x; x <= to_x; x += label_step) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", std::abs(x) < 0.1 ? 0.0 : x);
        renderText(x, 0.0, LABEL_BASE_HEIGHT, buf);
    }
    from_y = std::floor((cameraPos.y - cameraSize.y) / label_step) * label_step;
    to_y = std::ceil((cameraPos.y + cameraSize.y) / label_step) * label_step;
    for (double y = from_y; y <= to_y; y += label_step) {
        char buf[32];
        snprintf(buf, sizeof(buf), "%.1f", std::abs(y) < 0.1 ? 0.0 : y);
        renderText(0.0, y, LABEL_BASE_HEIGHT, buf);
    }

    // Точки
    glPointSize(8.0f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < points.size(); ++i) {
        if ((int)i == hover_point) {
            glColor3f(0.7f, 0.2f, 0.2f);
        }
        else {
            glColor3f(1.0f, 0.2f, 0.2f);
        }
        glVertex2d(points[i].x, points[i].y);
    }
    glEnd();

    // Сплайн и производная
    if ((int)points.size() >= MIN_POINTS) {
        // Сплайн (красный)
        glColor3f(1.0f, 0.2f, 0.2f);
        glBegin(GL_LINE_STRIP);
        for (int seg = 0; seg <= (int)points.size() - 4; ++seg) {
            std::array<Point, 4> quad = { points[seg], points[seg + 1], points[seg + 2], points[seg + 3] };
            double x_start = quad[1].x;
            double x_end = quad[2].x;
            int steps = std::max(1, (int)((x_end - x_start) / step));
            for (int i = 0; i <= steps; ++i) {
                double t = x_start + (x_end - x_start) * i / steps;
                double y = lagrange3(t, quad);
                glVertex2d(t, y);
            }
        }
        glEnd();

        // Производная (синий)
        glColor3f(0.0f, 0.0f, 1.0f);
        glBegin(GL_LINE_STRIP);
        for (int seg = 0; seg <= (int)points.size() - 4; ++seg) {
            std::array<Point, 4> quad = { points[seg], points[seg + 1], points[seg + 2], points[seg + 3] };
            double x_start = quad[1].x;
            double x_end = quad[2].x;
            int steps = std::max(1, (int)((x_end - x_start) / step));
            for (int i = 0; i <= steps; ++i) {
                double t = x_start + (x_end - x_start) * i / steps;
                double y = lagrange3Derivative(t, quad);
                glVertex2d(t, y);
            }
        }
        glEnd();
    }
    else {
        // Ломаная
        glColor3f(1.0f, 0.2f, 0.2f);
        glBegin(GL_LINE_STRIP);
        for (const auto& p : points) {
            glVertex2d(p.x, p.y);
        }
        glEnd();
    }

    glFinish();
    glutSwapBuffers();
    updateWindowTitle();
}

// ====== Main ======
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("КГ Лаб 4");

    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(passiveMotion);

    createMenu();

    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

    calcCamera();
    glutMainLoop();
    return 0;
}